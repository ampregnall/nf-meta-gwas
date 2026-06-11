#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(
    description="Prepare MAGMA SNP-level input (SNP, P, N) from meta-analysis sumstats"
)
parser.add_argument("--input",  required=True, help="Gzipped TSV sumstats")
parser.add_argument("--output", required=True, help="Output file path")
args = parser.parse_args()

keep = {"rsID", "P", "N", "MLOG10P"}
df = pd.read_csv(
    args.input, sep="\t", compression="gzip",
    usecols=lambda c: c in keep,
    low_memory=False,
)

# Reconstruct P from MLOG10P when available — avoids underflow for very
# significant variants where the stored P column rounds to zero (issue #12)
if "MLOG10P" in df.columns:
    df["P"] = np.power(10.0, -df["MLOG10P"].values.astype(float))

df = df.rename(columns={"rsID": "SNP"})

# Remove variants with missing or blank rsID (MAGMA requires a valid identifier)
df = df.dropna(subset=["SNP", "P", "N"])
df = df[df["SNP"].astype(str).str.strip() != ""]

# Clip P at float minimum to prevent MAGMA receiving exact zero
df["P"] = df["P"].clip(lower=np.finfo(float).tiny)

df[["SNP", "P", "N"]].to_csv(args.output, sep="\t", index=False)
