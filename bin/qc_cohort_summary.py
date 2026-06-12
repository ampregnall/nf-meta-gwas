#!/usr/bin/env python3
"""Cross-cohort QC summary plots: filter cascade, LDSC summary, effect size comparison."""

import matplotlib
matplotlib.use("Agg")

import argparse
import numpy as np
import pandas as pd
import warnings
import plotnine as p9
from plotnine import (
    ggplot, aes, geom_col, geom_crossbar, geom_linerange, geom_hline,
    geom_errorbar, geom_blank, facet_wrap, scale_fill_manual, scale_x_discrete,
    scale_y_continuous, theme_classic, theme, element_text, element_blank,
)

FIG_WIDTH = 8
FIG_HEIGHT = 5
DPI = 300
PALETTE = ["#045ea7", "#82afd3", "#990000", "#f4a261", "#2a9d8f", "#e9c46a"]


def parse_args():
    p = argparse.ArgumentParser(description="Cross-cohort QC summary plots")
    p.add_argument("--filter_stats",  required=True, help="Combined filter_stats.txt (per phenotype)")
    p.add_argument("--ldsc_h2",       required=True, help="Combined ldsc_h2.txt (per phenotype)")
    p.add_argument("--qc_summaries",  required=True, nargs="+", help="Per-cohort qc-summary.tsv files")
    p.add_argument("--output",        required=True, help="Output prefix")
    return p.parse_args()


def _cohort_label(row):
    return f"{row['cohort']}\n({row['population']})"


def filter_cascade_plot(path):
    df = pd.read_csv(path, sep="\t")
    df["label"] = df.apply(_cohort_label, axis=1)
    steps = [
        ("After basic check",   "n_after_basic_check"),
        ("After palindromic",   "n_after_palindromic"),
        ("After status filter", "n_after_status_filter"),
        ("After MAC filter",    "n_after_mac_filter"),
    ]
    available = [(label, col) for label, col in steps if col in df.columns]
    if not available:
        return None
    long = []
    for step_label, col in available:
        tmp = df[["label", col]].copy()
        tmp.columns = ["label", "n"]
        tmp["step"] = step_label
        long.append(tmp)
    plot_df = pd.concat(long, ignore_index=True)
    step_order = [s for s, _ in available]
    plot_df["step"] = pd.Categorical(plot_df["step"], categories=step_order, ordered=True)
    color_map = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(step_order)}
    return (
        ggplot(plot_df, aes(x="label", y="n", fill="step"))
        + geom_col(position="dodge", alpha=0.85)
        + scale_fill_manual(values=color_map, name="Filter step")
        + scale_x_discrete(name="Cohort (population)")
        + scale_y_continuous(name="Variant count")
        + theme_classic()
        + theme(
            figure_size=(FIG_WIDTH, FIG_HEIGHT),
            axis_text_x=element_text(size=8, angle=20, ha="right"),
            legend_title=element_text(size=9),
            legend_text=element_text(size=8),
        )
    )


def ldsc_summary_plot(path):
    df = pd.read_csv(path, sep="\t")
    df["label"] = df.apply(_cohort_label, axis=1)
    panels_spec = [
        ("Lambda_GC",  "Lambda GC",     1.0,  None),
        ("Intercept",  "LDSC Intercept", 1.0, "Intercept_se"),
        ("h2",         "SNP h²",         None, "h2_se"),
        ("Mean_Chi2",  "Mean Chi²",      1.0,  None),
    ]
    long = []
    for col, panel_label, ref, se_col in panels_spec:
        if col not in df.columns:
            continue
        tmp = df[["label", col]].copy()
        tmp.columns = ["label", "value"]
        tmp["se"] = df[se_col].fillna(0.0) if se_col and se_col in df.columns else 0.0
        tmp["ref"] = ref
        tmp["panel"] = panel_label
        long.append(tmp)
    if not long:
        return None
    plot_df = pd.concat(long, ignore_index=True)
    panel_order = [l for c, l, _, _ in panels_spec if c in df.columns]
    plot_df["panel"] = pd.Categorical(plot_df["panel"], categories=panel_order, ordered=True)
    ref_lines = plot_df[plot_df["ref"].notna()][["panel", "ref"]].drop_duplicates()
    return (
        ggplot(plot_df, aes(x="label", y="value"))
        + geom_col(fill=PALETTE[0], alpha=0.8)
        + geom_errorbar(
            aes(ymin="value - se", ymax="value + se"),
            width=0.25, color="#333333"
        )
        + (geom_hline(
               data=ref_lines, mapping=aes(yintercept="ref"),
               linetype="dashed", color="#cc0000", size=0.5
           ) if not ref_lines.empty else geom_blank())
        + facet_wrap("panel", scales="free_y", ncol=2)
        + scale_x_discrete(name="Cohort (population)")
        + scale_y_continuous(name="Value")
        + theme_classic()
        + theme(
            figure_size=(FIG_WIDTH, FIG_HEIGHT),
            axis_text_x=element_text(size=7, angle=20, ha="right"),
            strip_text=element_text(size=9),
        )
    )


def effect_size_comparison_plot(qc_paths):
    dfs = [pd.read_csv(f, sep="\t") for f in qc_paths]
    df = pd.concat(dfs, ignore_index=True).sort_values(["phenotype", "cohort", "population"])
    df["label"] = df.apply(_cohort_label, axis=1)
    return (
        ggplot(df, aes(x="label"))
        + geom_linerange(
            aes(ymin="whisker_lo", ymax="whisker_hi"),
            color=PALETTE[0], size=0.6
        )
        + geom_crossbar(
            aes(y="median_beta", ymin="q1_beta", ymax="q3_beta"),
            fill=PALETTE[1], color=PALETTE[0], alpha=0.8, width=0.5
        )
        + geom_hline(yintercept=0, linetype="dashed", color="#444444", size=0.4)
        + scale_x_discrete(name="Cohort (population)")
        + scale_y_continuous(name="Effect size (BETA)")
        + theme_classic()
        + theme(
            figure_size=(FIG_WIDTH, FIG_HEIGHT),
            axis_text_x=element_text(size=8, angle=20, ha="right"),
            axis_title=element_text(size=10),
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

    cascade = filter_cascade_plot(args.filter_stats)
    if cascade:
        save_plot(cascade, args.output, "filter-cascade")

    ldsc = ldsc_summary_plot(args.ldsc_h2)
    if ldsc:
        save_plot(ldsc, args.output, "ldsc-summary")

    save_plot(effect_size_comparison_plot(args.qc_summaries), args.output, "effect-size-comparison")


if __name__ == "__main__":
    main()
