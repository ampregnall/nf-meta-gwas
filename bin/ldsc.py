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


if np.float64(sumstats_hapmap3.ldsc_h2["Intercept"][0]) > 1:
    sumstats.data["SE"] = sumstats.data["SE"] * np.sqrt(
        np.float64(sumstats_hapmap3.ldsc_h2["Intercept"][0])
    )
    sumstats.data["Z"] = sumstats.data["BETA"] / sumstats.data["SE"]
    sumstats.data["P"] = 2 * norm.sf(abs(sumstats.data["Z"]))

# Save LDSC h2 results
ldsc_df = sumstats_hapmap3.ldsc_h2.copy()
ldsc_df.insert(0, "population", args.population)
ldsc_df.insert(0, "cohort", args.cohort)
ldsc_df.insert(0, "phenotype", args.phenotype)
ldsc_out = f"{args.output}.ldsc_h2.txt"
ldsc_df.to_csv(ldsc_out, index=False, sep="\t")

# Create manhattan and qq plot
fig = sumstats.plot_mqq(
    skip=3,
    anno="GENENAME",
    build="38",
    title=f"{args.cohort} {args.phenotype} {args.population}".upper(),
    anno_style="expand",
    fontsize=8,
    anno_fontsize=8,
    anno_gtf_path=args.gtf,
    font_family="DejaVu Sans",
    fig_kwargs={"figsize": (7.5, 5), "dpi": 400},
)

fig.savefig(f"{args.output}-manhattan-qq.png", dpi=400, bbox_inches="tight")
fig.savefig(f"{args.output}-manhattan-qq.pdf", dpi=400, bbox_inches="tight")

# Save sumstats results
sumstats_out = f"{args.output}.sumstats.processed.txt.gz"
parquet_out = f"{args.output}.sumstats.parquet"
sumstats.data.to_csv(sumstats_out, index=False, compression="gzip", sep="\t")
sumstats.data.to_parquet(parquet_out)
