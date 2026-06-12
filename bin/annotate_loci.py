#!/usr/bin/env python3
"""Annotate locus summary with NEAREST_GENE (from credsets) and optionally
FLAMES gene prioritization.  Produces locus_annotated.txt with three new columns:
  NEAREST_GENE         – modal nearest gene across credset variants at the locus
  PRIORITIZED_GENE     – FLAMES top gene when available, else NEAREST_GENE
  PRIORITIZATION_METHOD – "FLAMES" or "NEAREST_GENE"
"""

import argparse
import glob
import logging
import os
import sys

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description="Annotate locus summary")
    p.add_argument("--locus_summary", required=True)
    p.add_argument("--credsets", nargs="+", required=True)
    p.add_argument("--flames_scores", nargs="*", default=[],
                   help="FLAMES scores directories (one per population)")
    p.add_argument("--output", required=True)
    return p.parse_args()


def add_nearest_gene(locus_summary, credset_paths):
    """Range-join credset variants onto locus_summary to get NEAREST_GENE.

    credset LOCUS format: "CHR:LEAD_POS" (e.g. "1:12345678").
    locus_summary LOCUS: integer cluster ID from collect_lead_variants.py.
    Join key: credset LEAD_POS must fall within locus_summary [LOCUS_START, LOCUS_END]
    on the same chromosome.
    """
    frames = []
    for f in credset_paths:
        try:
            df = pd.read_csv(f, sep="\t", usecols=["LOCUS", "NEAREST_GENE"])
            frames.append(df)
        except Exception as e:
            log.warning(f"Could not read credset {f}: {e}")

    if not frames:
        locus_summary = locus_summary.copy()
        locus_summary["NEAREST_GENE"] = pd.NA
        return locus_summary

    credset_loci = (
        pd.concat(frames, ignore_index=True)
        .drop_duplicates(subset=["LOCUS", "NEAREST_GENE"])
    )

    # Parse "CHR:POS" into integer columns
    split = credset_loci["LOCUS"].str.split(":", expand=True)
    credset_loci = credset_loci.copy()
    credset_loci["LOCUS_CHR"] = split[0].astype(int)
    credset_loci["LOCUS_POS"] = split[1].astype(int)

    loci = locus_summary[["LOCUS", "CHR", "LOCUS_START", "LOCUS_END"]].copy()
    loci["CHR"] = loci["CHR"].astype(int)
    loci["LOCUS_START"] = loci["LOCUS_START"].astype(int)
    loci["LOCUS_END"] = loci["LOCUS_END"].astype(int)

    # Cross-join on CHR, then filter by position window
    merged = credset_loci.merge(loci, left_on="LOCUS_CHR", right_on="CHR", suffixes=("_cred", "_locus"))
    merged = merged[
        (merged["LOCUS_POS"] >= merged["LOCUS_START"]) &
        (merged["LOCUS_POS"] <= merged["LOCUS_END"])
    ]

    nearest = (
        merged.groupby("LOCUS_locus")["NEAREST_GENE"]
        .agg(lambda x: x.mode().iloc[0] if len(x) > 0 else pd.NA)
        .rename("NEAREST_GENE")
        .reset_index()
        .rename(columns={"LOCUS_locus": "LOCUS"})
    )

    return locus_summary.merge(nearest, on="LOCUS", how="left")


def load_flames_top_genes(flames_score_dirs, locus_summary, credset_paths):
    """Return {locus_int: gene_name} for the top FLAMES-prioritised gene per locus.

    FLAMES GenomicLocus integers are assigned by prepare_flames_inputs.R as
    row_number() of credset loci sorted by (CHR, POS).  We reconstruct that
    mapping here to translate FLAMES integers back to our locus_summary LOCUS ids.
    """
    # Collect all files in the score directories
    all_files = []
    for d in flames_score_dirs:
        if not os.path.isdir(d):
            log.warning(f"FLAMES scores directory not found: {d}")
            continue
        for pattern in ("*.txt", "*.tsv"):
            all_files.extend(glob.glob(os.path.join(d, pattern)))

    if not all_files:
        log.warning("No FLAMES score files found – skipping FLAMES gene prioritisation")
        return {}

    frames = []
    for f in all_files:
        try:
            df = pd.read_csv(f, sep="\t")
            frames.append(df)
        except Exception as e:
            log.warning(f"Could not read FLAMES file {f}: {e}")

    if not frames:
        log.warning("FLAMES files could not be parsed – falling back to NEAREST_GENE")
        return {}

    flames_df = pd.concat(frames, ignore_index=True)

    # Identify columns by name convention (case-insensitive)
    cols_lower = {c.lower(): c for c in flames_df.columns}
    gene_col = next((cols_lower[k] for k in cols_lower if "gene" in k), None)
    locus_col = next(
        (cols_lower[k] for k in cols_lower if "genomic" in k or ("locus" in k and "gene" not in k)),
        None
    )
    if gene_col is None or locus_col is None:
        log.warning(f"Cannot identify gene/locus columns in FLAMES output. Columns: {flames_df.columns.tolist()}")
        return {}

    # Score column: first numeric column that isn't the locus column
    numeric_cols = [c for c in flames_df.columns
                    if pd.api.types.is_numeric_dtype(flames_df[c]) and c != locus_col]
    score_col = numeric_cols[0] if numeric_cols else None

    if score_col:
        top_genes = (
            flames_df.sort_values(score_col, ascending=False)
            .drop_duplicates(subset=[locus_col])
            [[locus_col, gene_col]]
        )
    else:
        top_genes = flames_df.drop_duplicates(subset=[locus_col])[[locus_col, gene_col]]

    log.info(f"FLAMES: identified {len(top_genes)} locus–gene pairs (locus={locus_col}, gene={gene_col})")

    # Reconstruct GenomicLocus → locus_summary LOCUS mapping
    # prepare_flames_inputs.R assigns GenomicLocus = row_number() of credset
    # loci ordered by (chr, pos) where locus = "chr:pos"
    cred_frames = []
    for f in credset_paths:
        try:
            df = pd.read_csv(f, sep="\t", usecols=["LOCUS"])
            cred_frames.append(df)
        except Exception:
            pass

    if not cred_frames:
        log.warning("Could not read credset loci for FLAMES mapping")
        return {}

    cred_loci_raw = pd.concat(cred_frames).drop_duplicates(subset=["LOCUS"])
    split = cred_loci_raw["LOCUS"].str.split(":", expand=True)
    cred_loci_raw = cred_loci_raw.copy()
    cred_loci_raw["LOCUS_CHR"] = split[0].astype(int)
    cred_loci_raw["LOCUS_POS"] = split[1].astype(int)
    cred_loci_raw = cred_loci_raw.sort_values(["LOCUS_CHR", "LOCUS_POS"]).reset_index(drop=True)
    cred_loci_raw["GenomicLocus"] = range(1, len(cred_loci_raw) + 1)

    # Range-join credset loci → locus_summary integer LOCUS
    loci = locus_summary[["LOCUS", "CHR", "LOCUS_START", "LOCUS_END"]].copy()
    loci["CHR"] = loci["CHR"].astype(int)
    loci["LOCUS_START"] = loci["LOCUS_START"].astype(int)
    loci["LOCUS_END"] = loci["LOCUS_END"].astype(int)

    mapped = cred_loci_raw.merge(loci, left_on="LOCUS_CHR", right_on="CHR")
    mapped = mapped[
        (mapped["LOCUS_POS"] >= mapped["LOCUS_START"]) &
        (mapped["LOCUS_POS"] <= mapped["LOCUS_END"])
    ][["GenomicLocus", "LOCUS"]].drop_duplicates()

    result = top_genes.merge(mapped, left_on=locus_col, right_on="GenomicLocus", how="inner")
    return dict(zip(result["LOCUS"], result[gene_col]))


def main():
    args = parse_args()

    locus_summary = pd.read_csv(args.locus_summary, sep="\t")

    locus_summary = add_nearest_gene(locus_summary, args.credsets)
    locus_summary["PRIORITIZED_GENE"] = locus_summary["NEAREST_GENE"]
    locus_summary["PRIORITIZATION_METHOD"] = "NEAREST_GENE"

    if args.flames_scores:
        flames_genes = load_flames_top_genes(args.flames_scores, locus_summary, args.credsets)
        if flames_genes:
            mask = locus_summary["LOCUS"].isin(flames_genes)
            locus_summary.loc[mask, "PRIORITIZED_GENE"] = (
                locus_summary.loc[mask, "LOCUS"].map(flames_genes)
            )
            locus_summary.loc[mask, "PRIORITIZATION_METHOD"] = "FLAMES"
            log.info(f"FLAMES prioritised {mask.sum()} of {len(locus_summary)} loci")

    locus_summary.to_csv(args.output, sep="\t", index=False)
    log.info(f"Written: {args.output}")


if __name__ == "__main__":
    main()
