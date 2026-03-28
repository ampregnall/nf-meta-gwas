#!/usr/bin/env python

import argparse
import gwaslab as gl
import numpy as np
from scipy.stats import norm

parser = argparse.ArgumentParser(description="Test gwaslab installation")
parser.add_argument("--input", type=str, required=True, help="Path to input file")
parser.add_argument("--output", type=str, required=True, help="Output filenmame")
parser.add_argument("--ldsc", type=str, required=True, help = "Path to Pan-UKBB LD reference panels")
args = parser.parse_args()


sumstats = gl.Sumstats(args.input, fmt="gwaslab", build="38")

# Perform LDSC correction
sumstats_hapmap3 = sumstats.filter_hapmap3(inplace=False)
sumstats_hapmap3.estimate_h2_by_ldsc(ref_ld = args.ldsc,  w_ld = args.ldsc)


if np.float64(sumstats_hapmap3.ldsc_h2['Intercept'][0]) > 1:
    sumstats.data['SE'] = sumstats.data['SE'] * np.sqrt(np.float64(sumstats_hapmap3.ldsc_h2['Intercept'][0]))
    sumstats.data['Z'] = sumstats.data['BETA'] / sumstats.data['SE']
    sumstats.data['P'] = 2 * norm.sf(abs(sumstats.data['Z']))
    
# Save results
sumstats_out = f"{args.output}.sumstats.processed.txt.gz"
parquet_out = f"{args.output}.sumstats.parquet"
sumstats.data.to_csv(sumstats_out, index = False, compression="gzip", sep = "\t")
sumstats.data.to_parquet(parquet_out)
