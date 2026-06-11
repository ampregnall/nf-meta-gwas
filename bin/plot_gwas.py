#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")

import argparse
import gzip
import re
import numpy as np
import pandas as pd
from scipy import stats
import plotnine as p9

GWS_THRESHOLD = 5e-8
COLOR_ODD = "#045ea7"
COLOR_EVEN = "#82afd3"
COLOR_SIG = "#990000"
FIG_WIDTH = 8
FIG_HEIGHT = 5
DPI = 300
CHR_GAP = 5_000_000


# ---------------------------------------------------------------------------
# Core statistics
# ---------------------------------------------------------------------------

def mlog10p_from_z(beta, se):
    """
    Numerically stable -log10(p) computed from Z = beta/se.
    Uses the log-scale survival function to avoid underflow for large Z.
    p = 2 * P(N(0,1) > |Z|)  →  log(p) = log(2) + logsf(|Z|)
    """
    z = np.abs(np.asarray(beta, dtype=float) / np.asarray(se, dtype=float))
    log_p = np.log(2) + stats.norm.logsf(z)
    return -log_p / np.log(10)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_sumstats(path):
    keep = {"SNPID", "rsID", "CHR", "POS", "BETA", "SE", "P", "Z", "MLOG10P"}
    df = pd.read_csv(
        path, sep="\t", compression="gzip",
        usecols=lambda c: c in keep,
        dtype={"CHR": str},
        low_memory=False,
    )
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df.dropna(subset=["CHR", "POS", "BETA", "SE"])
    df["CHR"] = df["CHR"].astype(int)
    df["POS"] = df["POS"].astype(int)
    df["BETA"] = df["BETA"].astype(float)
    df["SE"] = df["SE"].astype(float)
    return df


def compute_mlog10p(df):
    """
    Use pre-computed MLOG10P column when present (meta-analysis output).
    Otherwise compute from BETA/SE via log-scale logsf to handle underflow
    (issue #12), with a final fallback to stored P for any non-finite results.
    """
    if "MLOG10P" in df.columns:
        vals = df["MLOG10P"].values.astype(float)
    else:
        vals = mlog10p_from_z(df["BETA"].values, df["SE"].values)
        if "P" in df.columns:
            p_stored = df["P"].values.astype(float)
            fallback = -np.log10(np.clip(p_stored, a_min=np.finfo(float).tiny, a_max=None))
            bad = ~np.isfinite(vals)
            vals[bad] = fallback[bad]
    df = df.copy()
    df["mlog10p"] = np.clip(vals, a_min=0.0, a_max=None)
    return df


# ---------------------------------------------------------------------------
# Plot geometry helpers
# ---------------------------------------------------------------------------

def add_cumulative_pos(df):
    """Lay chromosomes end-to-end with a small gap between them."""
    chrom_max = df.groupby("CHR")["POS"].max().sort_index()
    offsets = {}
    cumulative = 0
    for chrom in sorted(chrom_max.index):
        offsets[chrom] = cumulative
        cumulative += chrom_max[chrom] + CHR_GAP
    df = df.copy()
    offset_series = pd.Series(offsets)
    df["cumpos"] = df["POS"] + df["CHR"].map(offset_series)
    return df, offset_series


def assign_colors(df):
    """
    Whole chromosomes that contain at least one GWS variant → COLOR_SIG.
    Remaining odd chromosomes → COLOR_ODD, even → COLOR_EVEN.
    """
    gws_line = -np.log10(GWS_THRESHOLD)
    sig_chroms = set(df.loc[df["mlog10p"] >= gws_line, "CHR"].unique())

    def _color(chrom):
        if chrom in sig_chroms:
            return COLOR_SIG
        return COLOR_ODD if chrom % 2 == 1 else COLOR_EVEN

    df = df.copy()
    df["color"] = df["CHR"].map(_color)
    return df


def chr_tick_midpoints(df, offsets):
    """Cumulative midpoint per chromosome for x-axis tick placement."""
    ticks = {}
    for chrom, grp in df.groupby("CHR"):
        ticks[chrom] = grp["POS"].median() + offsets.get(chrom, 0)
    return ticks


# ---------------------------------------------------------------------------
# Gene annotation from GTF
# ---------------------------------------------------------------------------

def nearest_gene_from_gtf(positions, gtf_path):
    """
    For each (chr, pos) in *positions*, return the name of the nearest gene
    midpoint found in the Ensembl GTF.  Reads only 'gene' feature lines for
    the relevant chromosomes to keep memory usage low.
    """
    target_chroms = {str(c) for c, _ in positions}
    genes = []
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 9)
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0].lstrip("chr")
            if chrom not in target_chroms:
                continue
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            attrs = parts[8]
            m = re.search(r'gene_name "([^"]+)"', attrs) or \
                re.search(r'gene_id "([^"]+)"', attrs)
            gene_name = m.group(1) if m else f"{chrom}:{(start + end) // 2}"
            genes.append((chrom, (start + end) // 2, gene_name))

    if not genes:
        return {pos: f"{pos[0]}:{pos[1]}" for pos in positions}

    gene_df = pd.DataFrame(genes, columns=["chr", "mid", "gene"])
    result = {}
    for chrom, pos in positions:
        sub = gene_df[gene_df["chr"] == str(chrom)]
        if sub.empty:
            result[(chrom, pos)] = f"{chrom}:{pos}"
            continue
        idx = (sub["mid"] - pos).abs().idxmin()
        result[(chrom, pos)] = sub.at[idx, "gene"]
    return result


def get_gws_annotations(df):
    """
    Return the most significant GWS variant per chromosome.
    Returns an empty DataFrame if no variants pass the GWS threshold.
    """
    gws_line = -np.log10(GWS_THRESHOLD)
    gws = df[df["mlog10p"] >= gws_line]
    if gws.empty:
        return pd.DataFrame()
    return gws.loc[gws.groupby("CHR")["mlog10p"].idxmax()].copy().reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def manhattan_plot(df, offsets, annotations, gtf_path):
    ticks = chr_tick_midpoints(df, offsets)
    x_breaks = [ticks[c] for c in sorted(ticks)]
    x_labels = [str(c) for c in sorted(ticks)]
    sig_line = -np.log10(GWS_THRESHOLD)
    color_vals = sorted(df["color"].unique())

    plot = (
        p9.ggplot(df, p9.aes(x="cumpos", y="mlog10p", color="color"))
        + p9.geom_point(size=0.35, alpha=0.75, stroke=0)
        + p9.geom_hline(yintercept=sig_line, linetype="dashed", color="#444444", size=0.4)
        + p9.scale_color_manual(values={c: c for c in color_vals}, guide=False)
        + p9.scale_x_continuous(
            breaks=x_breaks,
            labels=x_labels,
            expand=(0.01, 0),
            name="Chromosome",
        )
        + p9.scale_y_continuous(
            expand=(0.01, 0),
            name=r"$-\log_{10}(p)$",
        )
        + p9.theme_classic()
        + p9.theme(
            figure_size=(FIG_WIDTH, FIG_HEIGHT),
            axis_text_x=p9.element_text(size=7),
            axis_title=p9.element_text(size=10),
        )
    )

    if not annotations.empty:
        positions = list(zip(annotations["CHR"].tolist(), annotations["POS"].tolist()))
        gene_map = nearest_gene_from_gtf(positions, gtf_path)
        annotations = annotations.copy()
        annotations["label"] = annotations.apply(
            lambda r: gene_map.get((r["CHR"], r["POS"]), ""), axis=1
        )
        plot = (
            plot
            + p9.geom_point(data=annotations, color="black", size=1.5, stroke=0)
            + p9.geom_text(
                data=annotations,
                mapping=p9.aes(label="label"),
                color="black",
                size=6,
                nudge_y=0.25,
                va="bottom",
                fontweight="bold",
            )
        )

    return plot


def qq_plot(df):
    """
    QQ plot of observed vs expected -log10(p).
    Downsamples points below -log10(p) = 1 for plotting efficiency while
    preserving the correct expected quantile for each kept point.
    """
    n = len(df)
    df_sorted = df.sort_values("mlog10p", ascending=False).reset_index(drop=True)
    df_sorted["expected"] = -np.log10((df_sorted.index + 1) / n)

    above = df_sorted[df_sorted["mlog10p"] >= 1.0]
    below = df_sorted[df_sorted["mlog10p"] < 1.0]
    if len(below) > 0:
        below = below.sample(frac=0.05, random_state=42)
    qq_df = pd.concat([above, below], ignore_index=True)

    max_val = max(qq_df["mlog10p"].max(), qq_df["expected"].max()) * 1.05

    return (
        p9.ggplot(qq_df, p9.aes(x="expected", y="mlog10p"))
        + p9.geom_point(size=0.35, alpha=0.75, color=COLOR_ODD, stroke=0)
        + p9.geom_abline(slope=1, intercept=0, color="#cc0000", linetype="dashed", size=0.5)
        + p9.scale_x_continuous(
            limits=(0, max_val),
            expand=(0.01, 0),
            name=r"Expected $-\log_{10}(p)$",
        )
        + p9.scale_y_continuous(
            limits=(0, max_val),
            expand=(0.01, 0),
            name=r"Observed $-\log_{10}(p)$",
        )
        + p9.theme_classic()
        + p9.theme(figure_size=(FIG_WIDTH, FIG_HEIGHT))
    )


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

def save_plot(plot, prefix, suffix):
    for ext in ("png", "pdf"):
        p9.ggsave(
            plot,
            filename=f"{prefix}.{suffix}.{ext}",
            dpi=DPI,
            width=FIG_WIDTH,
            height=FIG_HEIGHT,
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate Manhattan and QQ plots for GWAS summary statistics"
    )
    parser.add_argument("--input",  required=True, help="Gzipped TSV sumstats (gwaslab fmt)")
    parser.add_argument("--output", required=True, help="Output filename prefix")
    parser.add_argument("--gtf",    required=True, help="Ensembl GTF for nearest-gene annotation")
    args = parser.parse_args()

    df = load_sumstats(args.input)
    df = compute_mlog10p(df)
    df, offsets = add_cumulative_pos(df)
    df = assign_colors(df)

    annotations = get_gws_annotations(df)
    save_plot(manhattan_plot(df, offsets, annotations, args.gtf), args.output, "manhattan")
    save_plot(qq_plot(df), args.output, "qq")


if __name__ == "__main__":
    main()
