#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gwaslab as gl
import argparse
import json
import gzip
import re
import os
import numpy as np
import pandas as pd

DPI = 300
POP_ORDER = ["all", "eur", "afr", "amr", "eas", "csa", "mid"]

GWASLAB_KEEP = {
    "SNPID", "rsID", "CHR", "POS", "EA", "NEA", "EAF",
    "BETA", "SE", "Z", "P", "MLOG10P",
    "N", "N_CASE", "N_CONTROL", "I2", "N_CONTRIBUTIONS", "MAF",
}


def load_sumstats(path):
    df = pd.read_csv(
        path, sep="\t", compression="gzip",
        usecols=lambda c: c in GWASLAB_KEEP,
        low_memory=False,
    )
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df.dropna(subset=["CHR", "POS", "BETA", "SE"])
    df["CHR"] = df["CHR"].astype(int)
    df["POS"] = df["POS"].astype(int)
    return df


def load_gene_index(gtf_path):
    """Build {chrom_str: [(midpoint, gene_name)]} from Ensembl GTF."""
    index = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 9)
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0].lstrip("chr")
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            attrs = parts[8]
            m = re.search(r'gene_name "([^"]+)"', attrs) or \
                re.search(r'gene_id "([^"]+)"', attrs)
            gene_name = m.group(1) if m else f"{chrom}:{(start + end) // 2}"
            index.setdefault(chrom, []).append(((start + end) // 2, gene_name))
    return index


def nearest_gene(gene_index, chrom, pos):
    entries = gene_index.get(str(chrom), [])
    if not entries:
        return f"{chrom}:{pos}"
    mids, names = zip(*entries)
    idx = int(np.argmin(np.abs(np.array(mids) - pos)))
    return names[idx]


def register_gtf_with_gwaslab(gtf_path):
    """Register local GTF in gwaslab's config so it's used instead of downloading."""
    config_path = gl.options.paths["config"]
    os.makedirs(os.path.dirname(config_path), exist_ok=True)

    config = {}
    if os.path.exists(config_path):
        try:
            with open(config_path) as f:
                config = json.load(f)
        except Exception:
            pass

    config.setdefault("downloaded", {})
    config["downloaded"]["ensembl_hg38_gtf"] = {"local_path": gtf_path}

    with open(config_path, "w") as f:
        json.dump(config, f)


def main():
    parser = argparse.ArgumentParser(description="Create gwaslab locus plots for GWS loci")
    parser.add_argument("--locus_summary", required=True)
    parser.add_argument("--sumstats", nargs="+", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--population_vcf_map", required=True,
                        help="JSON mapping population → VCF path")
    args = parser.parse_args()

    locus_df = pd.read_csv(args.locus_summary, sep="\t")
    if locus_df.empty:
        print("No loci to plot.")
        return

    vcf_map = json.loads(args.population_vcf_map)
    eur_vcf = vcf_map.get("eur")

    print("Registering local GTF with gwaslab...")
    register_gtf_with_gwaslab(args.gtf)

    sumstats_map = {}
    for f in args.sumstats:
        basename = os.path.basename(f)
        pop = basename.replace(f"{args.phenotype}-", "").replace(".txt.gz", "")
        sumstats_map[pop] = f

    print(f"Loading {len(sumstats_map)} population sumstats...")
    all_data = {pop: load_sumstats(path) for pop, path in sumstats_map.items()}

    print("Loading gene index from GTF...")
    gene_index = load_gene_index(args.gtf)

    ordered_pops = sorted(
        sumstats_map.keys(),
        key=lambda x: POP_ORDER.index(x.lower()) if x.lower() in POP_ORDER else len(POP_ORDER),
    )

    for _, locus in locus_df.iterrows():
        chrom = int(locus["CHR"])
        start = int(locus["LOCUS_START"])
        end = int(locus["LOCUS_END"])
        lead_pos = int(locus["POS"])

        gene_name = nearest_gene(gene_index, chrom, lead_pos)
        gene_clean = re.sub(r"[^A-Za-z0-9_]", "_", gene_name)
        out_prefix = f"{chrom}-{start}-{end}-{gene_clean}"

        panels, vcfs, titles = [], [], []

        for pop in ordered_pops:
            region_df = all_data[pop]
            region_df = region_df[
                (region_df["CHR"] == chrom)
                & (region_df["POS"] >= start)
                & (region_df["POS"] <= end)
            ].copy()

            if region_df.empty:
                print(f"  No data in {chrom}:{start}-{end} for {pop}, skipping panel.")
                continue

            # Fallback to EUR VCF for any population without a dedicated entry
            vcf_path = vcf_map.get(pop.lower(), eur_vcf)

            panels.append(region_df)
            vcfs.append(vcf_path)
            titles.append(pop.upper())

        if not panels:
            print(f"No data for locus {chrom}:{start}-{end}, skipping.")
            continue

        print(f"Plotting {out_prefix} ({len(panels)} populations)...")
        try:
            fig, _ = gl.plot_stacked_mqq(
                panels,
                vcfs=vcfs,
                mode="r",
                region=(chrom, start, end),
                titles=titles,
                build="38",
                save=None,
                verbose=False,
            )
            fig.savefig(f"{out_prefix}.locus.png", dpi=DPI, bbox_inches="tight")
            fig.savefig(f"{out_prefix}.locus.pdf", bbox_inches="tight")
            print(f"  Saved {out_prefix}.locus.png/pdf")
        except Exception as e:
            print(f"  WARNING: Failed to plot {out_prefix}: {e}")
        finally:
            plt.close("all")


if __name__ == "__main__":
    main()
