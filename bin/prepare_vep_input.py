#!/usr/bin/env python3
"""Convert a credible set file to VEP's 5-column headerless input format."""

import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Format credible set for VEP input")
    p.add_argument("--credset", required=True, help="Credible set TSV from ABF fine-mapping")
    p.add_argument("--output",  required=True, help="Output VEP input file path")
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.credset, sep="\t", usecols=["SNPID", "CHR", "POS", "NEA", "EA"])
    df = df.dropna(subset=["CHR", "POS", "NEA", "EA", "SNPID"]).copy()
    df["CHR"] = df["CHR"].astype(int)
    df["POS"] = df["POS"].astype(int)
    # VEP format: chr  pos  REF/ALT  strand  id
    # NEA = reference (REF), EA = alternative (ALT)
    df["ALLELES"] = df["NEA"] + "/" + df["EA"]
    df["STRAND"]  = "+"
    df[["CHR", "POS", "ALLELES", "STRAND", "SNPID"]].to_csv(
        args.output, sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    main()
