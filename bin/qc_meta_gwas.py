#!/usr/bin/env python3
"""Meta-analysis QC plots: I² heterogeneity distribution and N_CONTRIBUTIONS breakdown."""

import matplotlib
matplotlib.use("Agg")

import argparse
import warnings
import numpy as np
import pandas as pd
import plotnine as p9
from plotnine import (
    ggplot, aes, geom_histogram, geom_col, geom_vline,
    scale_x_continuous, scale_x_discrete, scale_y_continuous,
    theme_classic, theme, element_text,
)

FIG_WIDTH = 8
FIG_HEIGHT = 5
DPI = 300
COLOR_MAIN = "#045ea7"
COLOR_LIGHT = "#82afd3"


def parse_args():
    p = argparse.ArgumentParser(description="Meta-analysis QC plots")
    p.add_argument("--input",  required=True, help="Collected meta sumstats (.txt.gz)")
    p.add_argument("--output", required=True, help="Output prefix")
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
    df = load_meta(args.input)
    if df.empty:
        return

    if "I2" in df.columns:
        save_plot(heterogeneity_plot(df), args.output, "heterogeneity")

    if "N_CONTRIBUTIONS" in df.columns:
        save_plot(n_contributions_plot(df), args.output, "n-contributions")


if __name__ == "__main__":
    main()
