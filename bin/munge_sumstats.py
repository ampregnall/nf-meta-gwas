#!/usr/bin/env python

import argparse
import gwaslab as gl

parser = argparse.ArgumentParser(description="Test gwaslab installation")
parser.add_argument("--input", type=str, required=True, help="Path to input file")
parser.add_argument("--output", type=str, required=True, help="Output filenmame")
parser.add_argument(
    "--type", type=str, required=True, help="Trait type (binary or continuous)"
)
parser.add_argument("--mac", type=int, required=True, help="Minimum minor allele count")
parser.add_argument("--chain", type=str, required=True, help="Path to chainfile")
parser.add_argument(
    "--fasta",
    type=str,
    required=True,
    help="Path to UCSD Human reference genome (GRCh38",
)
parser.add_argument(
    "--dbsnp",
    type=str,
    required=True,
    help="Path to dbSNP VCF file for rsID assignment",
)
parser.add_argument(
    "--popvcf",
    type=str,
    required=True,
    help="Path to population specific VCF file for strand inference",
)
parser.add_argument(
    "--ldsc", type=str, required=True, help="Path to Pan-UKBB LD reference panels"
)
parser.add_argument("--threads", type=int, required=True, help="Number of threads")
parser.add_argument("--phenotype", type=str, required=True, help="Phenotype name")
parser.add_argument("--cohort", type=str, required=True, help="Cohort name")
parser.add_argument("--population", type=str, required=True, help="Population label")
args = parser.parse_args()


sumstats = gl.Sumstats(args.input, fmt="gwaslab")
sumstats.infer_build(verbose=True)

# Liftover if statistics are on hg19
if sumstats.build == "19":
    sumstats.liftover(to_build="38", from_build="19", chain_path=args.chain)

sumstats.basic_check(
    remove=True,
    remove_dup=True,
    normalize=True,
    remove_dup_kwargs={"mode": "md", "keep": "first", "keep_col": "P"},
)

# Harmonize summary statistics
sumstats.harmonize(
    basic_check=False,
    ref_seq=args.fasta,
    ref_rsid_vcf=args.dbsnp,
    ref_infer=args.popvcf,  # Ancestry specific. Logic handled by Nextflow
    ref_alt_freq="AF",
    threads=args.threads,  # pass threads from nextflow process,
    sweep_mode=True,
)

status_str = sumstats.data["STATUS"].astype(str)
mask = ~(status_str.str[6].isin(["7", "8"]) | status_str.str[5].isin(["6"]))
sumstats.data = sumstats.data[mask].copy()

if args.type == "binary":
    sumstats.data["MAC"] = sumstats.data["MAF"] * sumstats.data["N_CASE"] * 2
else:
    sumstats.data["MAC"] = sumstats.data["MAF"] * sumstats.data["N"] * 2

sumstats.data = sumstats.data[sumstats.data["MAC"] > args.mac]
sumstats.data["DAF"] = sumstats.data["EAF"] - sumstats.data["RAF"]
sumstats.filter_value('abs(DAF) < 0.20', inplace=True)

# Create DAF plot
fig, _ = sumstats.plot_daf(
    title=f"{args.cohort} {args.phenotype} {args.population}".upper(),
    fontsize=8,
    font_family="DejaVu Sans",
    fig_kwargs={"figsize": (7.5, 5), "dpi": 400},
)

fig.savefig(f"{args.output}-daf.png", dpi=400, bbox_inches="tight")
fig.savefig(f"{args.output}-daf.pdf", dpi=400, bbox_inches="tight")

# Save results
sumstats_out = f"{args.output}.sumstats.munged.txt.gz"
sumstats.data.to_csv(sumstats_out, index=False, compression="gzip", sep="\t")
