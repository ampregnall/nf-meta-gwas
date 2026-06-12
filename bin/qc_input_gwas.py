#!/usr/bin/env python3
"""Per-cohort GWAS QC plots: precision, EAF, Z-score, and effect size distributions."""

import matplotlib
matplotlib.use("Agg")

import argparse
import numpy as np
import pandas as pd
import warnings
from scipy import stats
import plotnine as p9
from plotnine import (
    ggplot, aes, geom_violin, geom_histogram, geom_density, geom_line,
    geom_hline, scale_fill_identity, scale_x_discrete, scale_x_continuous,
    scale_y_continuous, labs, theme_classic, theme, element_text,
)

FIG_WIDTH = 8
FIG_HEIGHT = 5
DPI = 300
COLOR_ODD = "#045ea7"
COLOR_EVEN = "#82afd3"
N_DENSITY_MAX = 200_000
N_VIOLIN_PER_CHR = 5_000


def parse_args():
    p = argparse.ArgumentParser(description="Per-cohort GWAS QC plots")
    p.add_argument("--input",      required=True, help="Munged sumstats (.sumstats.munged.txt.gz)")
    p.add_argument("--output",     required=True, help="Output filename prefix")
    p.add_argument("--phenotype",  required=True)
    p.add_argument("--cohort",     required=True)
    p.add_argument("--population", required=True)
    return p.parse_args()


def load_sumstats(path):
    want = {"CHR", "POS", "BETA", "SE", "Z", "EAF", "MAF", "N", "N_CASE", "N_CONTROL"}
    header = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns.tolist()
    usecols = [c for c in header if c in want]
    df = pd.read_csv(path, sep="\t", compression="gzip", usecols=usecols, low_memory=False)
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df.dropna(subset=["CHR", "BETA", "SE"]).copy()
    df = df[df["CHR"].between(1, 22)]
    df["CHR"] = df["CHR"].astype(int)
    df["BETA"] = df["BETA"].astype(float)
    df["SE"] = df["SE"].astype(float)
    if "N" not in df.columns and "N_CASE" in df.columns and "N_CONTROL" in df.columns:
        df["N"] = pd.to_numeric(df["N_CASE"], errors="coerce") + \
                  pd.to_numeric(df["N_CONTROL"], errors="coerce")
    elif "N" in df.columns:
        df["N"] = pd.to_numeric(df["N"], errors="coerce")
    if "EAF" in df.columns:
        df["EAF"] = pd.to_numeric(df["EAF"], errors="coerce")
    if "Z" in df.columns:
        df["Z"] = pd.to_numeric(df["Z"], errors="coerce")
    return df


def precision_sample_size_plot(df):
    """Violin of log10(SE * sqrt(N)) grouped by chromosome. Constant SE*sqrt(N) is expected."""
    sub = df.dropna(subset=["N", "SE"]).copy()
    sub = sub[(sub["N"] > 0) & (sub["SE"] > 0)]
    sub["metric"] = np.log10(sub["SE"] * np.sqrt(sub["N"]))
    sub = sub.groupby("CHR", group_keys=False).apply(
        lambda g: g.sample(min(len(g), N_VIOLIN_PER_CHR), random_state=42)
    ).reset_index(drop=True)
    chrs = sorted(sub["CHR"].unique())
    sub["CHR_f"] = pd.Categorical(sub["CHR"].astype(str), categories=[str(c) for c in chrs], ordered=True)
    sub["fill"] = sub["CHR"].map(lambda c: COLOR_ODD if c % 2 == 1 else COLOR_EVEN)
    return (
        ggplot(sub, aes(x="CHR_f", y="metric", fill="fill"))
        + geom_violin(color="none", alpha=0.85, scale="width")
        + scale_fill_identity(guide=False)
        + scale_x_discrete(name="Chromosome")
        + scale_y_continuous(name=r"log₁₀(SE × √N)")
        + geom_hline(
            yintercept=float(sub["metric"].median()),
            linetype="dashed", color="#444444", size=0.4
        )
        + theme_classic()
        + theme(
            figure_size=(FIG_WIDTH, FIG_HEIGHT),
            axis_text_x=element_text(size=7),
            axis_title=element_text(size=10),
        )
    )


def eaf_distribution_plot(df):
    """Histogram of effect allele frequency."""
    col = "EAF" if "EAF" in df.columns else "MAF"
    sub = df[[col]].rename(columns={col: "eaf"}).dropna()
    return (
        ggplot(sub, aes(x="eaf"))
        + geom_histogram(bins=50, fill=COLOR_ODD, color="white", alpha=0.85)
        + scale_x_continuous(name="Effect allele frequency", limits=(0, 1))
        + scale_y_continuous(name="Variant count")
        + theme_classic()
        + theme(figure_size=(FIG_WIDTH, FIG_HEIGHT), axis_title=element_text(size=10))
    )


def zscore_distribution_plot(df):
    """Density of Z-scores with standard N(0,1) overlay."""
    sub = df[["Z"]].dropna().copy()
    if len(sub) > N_DENSITY_MAX:
        sub = sub.sample(N_DENSITY_MAX, random_state=42)
    z_min = float(sub["Z"].quantile(0.001))
    z_max = float(sub["Z"].quantile(0.999))
    zrange = np.linspace(z_min, z_max, 400)
    normal_df = pd.DataFrame({"Z": zrange, "density": stats.norm.pdf(zrange)})
    return (
        ggplot(sub, aes(x="Z"))
        + geom_density(fill=COLOR_EVEN, color=COLOR_ODD, alpha=0.5, size=0.6)
        + geom_line(
            data=normal_df, mapping=aes(x="Z", y="density"),
            color="#cc0000", linetype="dashed", size=0.7
        )
        + scale_x_continuous(name="Z-score", limits=(z_min, z_max))
        + scale_y_continuous(name="Density")
        + theme_classic()
        + theme(figure_size=(FIG_WIDTH, FIG_HEIGHT), axis_title=element_text(size=10))
    )


def beta_distribution_plot(df):
    """Histogram of effect sizes (BETA)."""
    sub = df[["BETA"]].dropna()
    if len(sub) > N_DENSITY_MAX:
        sub = sub.sample(N_DENSITY_MAX, random_state=42)
    return (
        ggplot(sub, aes(x="BETA"))
        + geom_histogram(bins=100, fill=COLOR_ODD, color="white", alpha=0.85)
        + scale_x_continuous(name="Effect size (BETA)")
        + scale_y_continuous(name="Variant count")
        + theme_classic()
        + theme(figure_size=(FIG_WIDTH, FIG_HEIGHT), axis_title=element_text(size=10))
    )


def write_qc_summary(df, prefix, phenotype, cohort, population):
    beta = df["BETA"].dropna()
    q1 = float(beta.quantile(0.25))
    q3 = float(beta.quantile(0.75))
    iqr = q3 - q1
    row = {
        "phenotype": phenotype,
        "cohort": cohort,
        "population": population,
        "n_variants": len(df),
        "median_beta": float(beta.median()),
        "q1_beta": q1,
        "q3_beta": q3,
        "whisker_lo": float(max(q1 - 1.5 * iqr, float(beta.min()))),
        "whisker_hi": float(min(q3 + 1.5 * iqr, float(beta.max()))),
        "median_se": float(df["SE"].dropna().median()),
    }
    if "EAF" in df.columns:
        row["median_eaf"] = float(df["EAF"].dropna().median())
    if "N" in df.columns:
        row["median_n"] = float(df["N"].dropna().median())
    pd.DataFrame([row]).to_csv(f"{prefix}.qc-summary.tsv", sep="\t", index=False)


def save_plot(plot, prefix, suffix):
    for ext in ("png", "pdf"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p9.ggsave(
                plot, filename=f"{prefix}.{suffix}.{ext}",
                dpi=DPI, width=FIG_WIDTH, height=FIG_HEIGHT,
            )


def main():
    args = parse_args()
    df = load_sumstats(args.input)

    if "N" in df.columns and df["N"].notna().any():
        save_plot(precision_sample_size_plot(df), args.output, "precision-sample-size")

    if "EAF" in df.columns or "MAF" in df.columns:
        save_plot(eaf_distribution_plot(df), args.output, "eaf-distribution")

    if "Z" in df.columns and df["Z"].notna().any():
        save_plot(zscore_distribution_plot(df), args.output, "zscore-distribution")

    save_plot(beta_distribution_plot(df), args.output, "effect-size-distribution")

    write_qc_summary(df, args.output, args.phenotype, args.cohort, args.population)


if __name__ == "__main__":
    main()
