#!/usr/bin/env python

import argparse
import pandas as pd
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
parser.add_argument(
    "--col_overrides",
    type=str,
    required=False,
    default=None,
    help="Comma-separated key=value column name overrides (e.g. snpid=rsid,chrom=chr,pos=bp,ea=A1,nea=A2,beta=b,se=se,eaf=af). "
         "When provided, all 7 required keys (chrom, pos, ea, nea, beta, se, eaf) must be included.",
)
args = parser.parse_args()

REQUIRED_OVERRIDE_KEYS = {"chrom", "pos", "ea", "nea", "beta", "se", "eaf"}
ALLOWED_OVERRIDE_KEYS = REQUIRED_OVERRIDE_KEYS | {
    "snpid", "rsid", "p", "n", "n_case", "n_con", "or", "z", "info", "i2", "phet"
}


def parse_col_overrides(override_str):
    if not override_str or not override_str.strip():
        return {}
    overrides = {}
    for pair in (p.strip() for p in override_str.split(",")):
        if "=" not in pair:
            raise ValueError(f"Invalid col_overrides entry (missing '='): '{pair}'")
        key, _, value = pair.partition("=")
        key = key.strip().lower()
        if key not in ALLOWED_OVERRIDE_KEYS:
            raise ValueError(
                f"Unknown override key '{key}'. Allowed keys: {sorted(ALLOWED_OVERRIDE_KEYS)}"
            )
        overrides[key] = value.strip()
    missing = REQUIRED_OVERRIDE_KEYS - overrides.keys()
    if missing:
        raise ValueError(
            f"col_overrides is missing required keys: {sorted(missing)}"
        )
    return overrides


col_overrides = parse_col_overrides(args.col_overrides)

if col_overrides:
    sumstats = gl.Sumstats(args.input, **col_overrides)
else:
    sumstats = gl.Sumstats(args.input, fmt="gwaslab")
sumstats.infer_build(verbose=True)

# Liftover if statistics are on hg19
if sumstats.build == "19":
    sumstats.liftover(to_build="38", from_build="19", chain_path=args.chain)

n_initial = len(sumstats.data)

sumstats.basic_check(
    remove=True,
    remove_dup=True,
    normalize=True,
    remove_dup_kwargs={"mode": "md", "keep": "first", "keep_col": "P"},
)
n_after_basic_check = len(sumstats.data)

# Additional filtering
sumstats.filter_palindromic(mode="out", inplace=True)
n_after_palindromic = len(sumstats.data)

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
n_after_status = len(sumstats.data)

if args.type == "binary":
    sumstats.data["MAC"] = sumstats.data["MAF"] * sumstats.data["N_CASE"] * 2
else:
    sumstats.data["MAC"] = sumstats.data["MAF"] * sumstats.data["N"] * 2

sumstats.data = sumstats.data[sumstats.data["MAC"] > args.mac]
n_after_mac = len(sumstats.data)

# Save variant count tracking table
filter_stats = pd.DataFrame([{
    "phenotype": args.phenotype,
    "cohort": args.cohort,
    "population": args.population,
    "n_initial": n_initial,
    "n_after_basic_check": n_after_basic_check,
    "n_removed_basic_check": n_initial - n_after_basic_check,
    "n_after_palindromic": n_after_palindromic,
    "n_removed_palindromic": n_after_basic_check - n_after_palindromic,
    "n_after_status_filter": n_after_status,
    "n_removed_status_filter": n_after_palindromic - n_after_status,
    "n_after_mac_filter": n_after_mac,
    "n_removed_mac_filter": n_after_status - n_after_mac,
}])
filter_stats.to_csv(f"{args.output}.filter_stats.txt", index=False, sep="\t")

# Save results
sumstats_out = f"{args.output}.sumstats.munged.txt.gz"
sumstats.data.to_csv(sumstats_out, index=False, compression="gzip", sep="\t")
