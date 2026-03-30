#!/usr/bin/env python

import argparse
import gwaslab as gl
import numpy as np
from scipy.stats import norm

parser = argparse.ArgumentParser(description="Test gwaslab installation")
parser.add_argument("--input", type=str, required=True, help="Path to input file")
parser.add_argument("--output", type=str, required=True, help="Output filenmame")
parser.add_argument(
    "--ldsc", type=str, required=True, help="Path to Pan-UKBB LD reference panels"
)
parser.add_argument("--phenotype", type=str, required=True, help="Phenotype name")
parser.add_argument("--cohort", type=str, required=True, help="Cohort name")
parser.add_argument("--population", type=str, required=True, help="Population label")
parser.add_argument("--gtf", type=str, required=True, help="GTF file")
args = parser.parse_args()


sumstats = gl.Sumstats(args.input, fmt="gwaslab", build="38")

# Perform LDSC correction
sumstats_hapmap3 = sumstats.filter_hapmap3(inplace=False)
sumstats_hapmap3.estimate_h2_by_ldsc(ref_ld=args.ldsc, w_ld=args.ldsc)

# Save LDSC h2 results
ldsc_df = sumstats_hapmap3.ldsc_h2.copy()
ldsc_df.insert(0, "population", args.population)
ldsc_df.insert(0, "cohort", args.cohort)
ldsc_df.insert(0, "phenotype", args.phenotype)
ldsc_out = f"{args.output}.ldsc_h2.txt"
ldsc_df.to_csv(ldsc_out, index=False, sep="\t")

