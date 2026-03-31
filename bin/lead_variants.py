#!/usr/bin/env python

import gwaslab as gl
import argparse

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
var_list = lead_vars.sort_values(by="P").head(10)

fig = sumstats.plot_mqq(
    skip=3,
    anno="GENENAME",
    build="38",
    highlight=list(var_list.GENE),
    anno_set=list(var_list.GENE),
    title=f"{args.cohort} {args.phenotype} {args.population}".upper(),
    anno_style="expand",
    fontsize=8,
    anno_fontsize=8,
    anno_max_rows=10,
    anno_gtf_path=args.gtf,
    font_family="DejaVu Sans",
    fig_kwargs={"figsize": (7.5, 5), "dpi": 400},
)

# Save results
fig.savefig(f"{args.output}-manhattan-qq.png", dpi=400, bbox_inches="tight")
fig.savefig(f"{args.output}-manhattan-qq.pdf", dpi=400, bbox_inches="tight")
lead_out = f"{args.output}.lead.variants.txt"
lead_vars.to_csv(lead_out, sep="\t", index=False)
