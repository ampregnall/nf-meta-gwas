#!/usr/bin/env python3
"""Meta-analysis QC plots: I² heterogeneity distribution and N_CONTRIBUTIONS breakdown,
plus 99% credible set size distribution violin plot."""

import matplotlib
matplotlib.use("Agg")

import argparse
import os
import warnings
import numpy as np
import pandas as pd
import plotnine as p9
from plotnine import (
    ggplot, aes, geom_histogram, geom_col, geom_violin, geom_jitter,
    geom_vline, scale_x_continuous, scale_x_discrete, scale_y_continuous,
    scale_y_log10, scale_fill_manual, theme_classic, theme, element_text,
)

FIG_WIDTH = 8
FIG_HEIGHT = 5
DPI = 300
COLOR_MAIN = "#045ea7"
COLOR_LIGHT = "#82afd3"

POP_ORDER = ["all", "eur", "afr", "amr", "eas", "csa", "mid"]


def parse_args():
    p = argparse.ArgumentParser(description="Meta-analysis QC plots")
    p.add_argument("--input",  help="Collected meta sumstats (.txt.gz)")
    p.add_argument("--output", required=True, help="Output prefix")
    p.add_argument("--credsets", nargs="+",
                   help="Credset .txt files (all populations) for CS size violin plot")
    p.add_argument("--phenotype", help="Phenotype name (used to parse population from credset filenames)")
    return p.parse_args()


def load_meta(path):
    want = {"I2", "N_CONTRIBUTIONS"}
    header = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns.tolist()
    usecols = [c for c in header if c in want]
    if not usecols:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", compression="gzip", usecols=usecols, low_memory=False)


def heterogeneity_plot(df):
    """Histogram of I² with reference lines at 25, 50, and 75%."""
    sub = df[["I2"]].copy()
    sub["I2"] = pd.to_numeric(sub["I2"], errors="coerce")
    sub = sub.dropna()
    # I2 is only meaningful for multi-cohort variants; 0 includes single-cohort
    return (
        ggplot(sub, aes(x="I2"))
        + geom_histogram(bins=50, fill=COLOR_MAIN, color="white", alpha=0.85)
        + geom_vline(xintercept=25, linetype="dashed", color="#cc0000", size=0.5)
        + geom_vline(xintercept=50, linetype="dashed", color="#cc0000", size=0.5)
        + geom_vline(xintercept=75, linetype="dashed", color="#cc0000", size=0.5)
        + scale_x_continuous(name="I² (%)", limits=(0, 100))
        + scale_y_continuous(name="Variant count")
        + theme_classic()
        + theme(figure_size=(FIG_WIDTH, FIG_HEIGHT), axis_title=element_text(size=10))
    )


def n_contributions_plot(df):
    """Bar chart of variants by number of contributing cohorts."""
    sub = df[["N_CONTRIBUTIONS"]].copy()
    sub["N_CONTRIBUTIONS"] = pd.to_numeric(sub["N_CONTRIBUTIONS"], errors="coerce")
    sub = sub.dropna()
    counts = sub["N_CONTRIBUTIONS"].astype(int).value_counts().reset_index()
    counts.columns = ["n_cohorts", "n_variants"]
    counts = counts.sort_values("n_cohorts")
    counts["n_cohorts_f"] = pd.Categorical(
        counts["n_cohorts"].astype(str),
        categories=[str(c) for c in sorted(counts["n_cohorts"].unique())],
        ordered=True,
    )
    return (
        ggplot(counts, aes(x="n_cohorts_f", y="n_variants"))
        + geom_col(fill=COLOR_MAIN, alpha=0.85)
        + scale_x_discrete(name="Number of contributing cohorts")
        + scale_y_continuous(name="Variant count")
        + theme_classic()
        + theme(figure_size=(FIG_WIDTH, FIG_HEIGHT), axis_title=element_text(size=10))
    )


def load_credsets(paths, phenotype):
    """Load credset files and return a DataFrame with per-locus CS sizes and population labels."""
    frames = []
    for f in paths:
        basename = os.path.basename(f)
        pop = basename.replace(f"{phenotype}-", "").replace(".credset.txt", "")
        df = pd.read_csv(f, sep="\t", usecols=["LOCUS"])
        if df.empty:
            continue
        sizes = df.groupby("LOCUS").size().reset_index(name="cs_size")
        sizes["population"] = pop.upper()
        frames.append(sizes)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def credset_violin_plot(df):
    """Violin plot of 99% credible set sizes per population (log y-axis)."""
    ordered = sorted(
        df["population"].unique(),
        key=lambda x: POP_ORDER.index(x.lower()) if x.lower() in POP_ORDER else len(POP_ORDER),
    )
    df = df.copy()
    df["population"] = pd.Categorical(df["population"], categories=ordered, ordered=True)

    return (
        ggplot(df, aes(x="population", y="cs_size"))
        + geom_violin(fill=COLOR_MAIN, color="white", alpha=0.80, trim=True, scale="width")
        + geom_jitter(width=0.12, size=0.9, alpha=0.55, color=COLOR_LIGHT, stroke=0)
        + scale_y_log10(
            name="99% credible set size (variants)",
            breaks=[1, 2, 5, 10, 25, 50, 100, 250, 500],
        )
        + scale_x_discrete(name="Population")
        + theme_classic()
        + theme(
            figure_size=(FIG_WIDTH, FIG_HEIGHT),
            axis_title=element_text(size=10),
            axis_text_x=element_text(size=9),
        )
    )


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

    if args.credsets:
        if not args.phenotype:
            raise ValueError("--phenotype is required when --credsets is provided")
        cs_df = load_credsets(args.credsets, args.phenotype)
        if not cs_df.empty:
            save_plot(credset_violin_plot(cs_df), args.output, "credset-sizes")
        return

    df = load_meta(args.input)
    if df.empty:
        return

    if "I2" in df.columns:
        save_plot(heterogeneity_plot(df), args.output, "heterogeneity")

    if "N_CONTRIBUTIONS" in df.columns:
        save_plot(n_contributions_plot(df), args.output, "n-contributions")


if __name__ == "__main__":
    main()
