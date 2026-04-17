#!/usr/bin/env python

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(
    description="Aggregate lead variants across populations, cluster into loci, and flag population-specific hits"
)
parser.add_argument("--inputs", nargs="+", required=True, help="Lead variant .txt files for all populations")
parser.add_argument("--phenotype", type=str, required=True, help="Phenotype name (used to parse population from filename)")
parser.add_argument("--output", type=str, required=True, help="Output file prefix")
args = parser.parse_args()

OUTPUT_COLS = ["LOCUS", "CHR", "LOCUS_START", "LOCUS_END", "POS", "EA", "NEA", "BETA", "SE", "P", "POPULATIONS", "POPULATION_SPECIFIC"]


def genome_cluster(df, distance=0, cluster_col="LOCUS"):
    """
    Re-implementation of tidygenomics::genome_cluster.

    Groups by CHR, clusters overlapping intervals (within max_distance),
    then assigns globally unique 0-indexed integer cluster IDs matching
    the as.factor(paste0(chr, '-', per_chr_id)) - 1 behaviour in R.
    """
    df = df.copy().sort_values(["CHR", "START"])
    keys = []
    current_chr, current_end, chr_cluster_id = None, -1, 0

    for _, row in df.iterrows():
        if row["CHR"] != current_chr:
            current_chr = row["CHR"]
            current_end = row["END"]
            chr_cluster_id = 1
        elif row["START"] <= current_end + distance:
            current_end = max(current_end, row["END"])
        else:
            chr_cluster_id += 1
            current_end = row["END"]
        keys.append(f"{row['CHR']}-{chr_cluster_id}")

    unique_keys = {k: i for i, k in enumerate(dict.fromkeys(keys))}
    df[cluster_col] = [unique_keys[k] for k in keys]
    return df


def is_population_specific(pops_str):
    pops = [p.lower() for p in pops_str.split(";")]
    return "eur" not in pops and "all" not in pops


# Load and tag each file with its population
frames = []
for f in args.inputs:
    pop = os.path.basename(f).replace(f"{args.phenotype}-", "").replace(".lead.variants.txt", "")
    df = pd.read_csv(f, sep="\t")
    if df.empty:
        continue
    df["POP"] = pop
    frames.append(df)

out_file = f"{args.output}.locus_summary.txt"

if not frames:
    pd.DataFrame(columns=OUTPUT_COLS).to_csv(out_file, sep="\t", index=False)
else:
    df = pd.concat(frames, ignore_index=True)

    # Compute 500 kb windows
    df["START"] = (df["POS"] - 500_000).clip(lower=0)
    df["END"] = df["POS"] + 500_000

    # Cluster overlapping windows into loci
    df = genome_cluster(df, distance=0, cluster_col="LOCUS")

    # Lead variant per locus (minimum P)
    lead_per_locus = df.loc[df.groupby("LOCUS")["P"].idxmin()][
        ["LOCUS", "POS", "EA", "NEA", "BETA", "SE", "P"]
    ]

    # Locus-level summary
    summary = (
        df.groupby("LOCUS")
        .agg(
            CHR=("CHR", "first"),
            LOCUS_START=("START", "min"),
            LOCUS_END=("END", "max"),
            POPULATIONS=("POP", lambda x: ";".join(sorted(x.unique()))),
        )
        .reset_index()
    )
    summary = summary.merge(lead_per_locus, on="LOCUS")
    summary["POPULATION_SPECIFIC"] = summary["POPULATIONS"].apply(is_population_specific)
    summary[OUTPUT_COLS].to_csv(out_file, sep="\t", index=False)
