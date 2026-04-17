#!/usr/bin/env python

import gwaslab as gl
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Extract lead variants and create Manhattan plot"
)
parser.add_argument("--input", type=str, required=True, help="Path to input file")
parser.add_argument("--output", type=str, required=True, help="Output filename")
parser.add_argument("--phenotype", type=str, required=True, help="Phenotype name")
parser.add_argument("--cohort", type=str, required=True, help="Cohort name")
parser.add_argument("--population", type=str, required=True, help="Population label")
parser.add_argument("--gtf", type=str, required=True, help="GTF file")
args = parser.parse_args()


sumstats = gl.Sumstats(args.input, fmt="gwaslab", build="38")
sumstats.basic_check()

# Extract lead vars
lead_vars = sumstats.get_lead(anno=True, gtf_path=args.gtf)
lead_out = f"{args.output}.lead.variants.txt"

if lead_vars.empty:
    # Determine trait type from input file columns to build the correct empty header
    header_cols = pd.read_csv(args.input, sep="\t", nrows=0, compression="gzip").columns
    is_binary = "N_CASE" in header_cols
    base_cols = ["SNPID", "rsID", "CHR", "POS", "EA", "NEA", "STATUS",
                 "EAF", "BETA", "SE", "Z", "P", "I2", "N"]
    if is_binary:
        base_cols += ["N_CASE", "N_CONTROL"]
    base_cols += ["LOCATION", "GENE"]
    pd.DataFrame(columns=base_cols).to_csv(lead_out, sep="\t", index=False)
else:
    lead_vars.to_csv(lead_out, sep="\t", index=False)
