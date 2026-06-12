#!/usr/bin/env python3
"""Join VEP annotation results back to the credible set."""

import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Annotate credible set with VEP consequences")
    p.add_argument("--credset",    required=True, help="Original credible set TSV")
    p.add_argument("--vep_output", required=True, help="VEP --tab output file")
    p.add_argument("--output",     required=True, help="Output annotated credible set path")
    return p.parse_args()


def read_vep_tab(path):
    """
    Parse VEP --tab output.
    Comment lines start with ##; the column header line starts with a single #.
    """
    header = None
    rows   = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            elif line.startswith("#"):
                header = line.lstrip("#").split("\t")
            else:
                rows.append(line.split("\t"))
    if header is None:
        return pd.DataFrame()
    return pd.DataFrame(rows, columns=header) if rows else pd.DataFrame(columns=header)


def main():
    args   = parse_args()
    credset = pd.read_csv(args.credset, sep="\t")
    vep     = read_vep_tab(args.vep_output)

    vep_cols = ["VEP_SYMBOL", "VEP_Consequence", "VEP_IMPACT", "VEP_VARIANT_CLASS"]

    if vep.empty or "Uploaded_variation" not in vep.columns:
        for col in vep_cols:
            credset[col] = pd.NA
    else:
        vep = vep.rename(columns={
            "Uploaded_variation": "SNPID",
            "SYMBOL":             "VEP_SYMBOL",
            "Consequence":        "VEP_Consequence",
            "IMPACT":             "VEP_IMPACT",
            "VARIANT_CLASS":      "VEP_VARIANT_CLASS",
        })
        keep = ["SNPID"] + [c for c in vep_cols if c in vep.columns]
        credset = credset.merge(vep[keep], on="SNPID", how="left")

    credset.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
