import gwaslab as gl
import argparse

parser = argparse.ArgumentParser(description="Argument parser for tissue and cell type enrichment")
parser.add_argument("--input", type=str, required=True, help="Path to input file")
parser.add_argument("--output", type=str, required=True, help="Output filenmame")
parser.add_argument("--ld", type=str, required=True, help="Path to Pan-UKBB LD reference panels")
parser.add_argument("--ldcts", type=str, required=True, help="Phenotype name")
parser.add_argument("--weights", type=str, required=True, help="Cohort name")
parser.add_argument("--output", type=str, required=True, help="Population label")
args = parser.parse_args()

sumstats = gl.Sumstats(args.input, fmt="gwaslab", build="38")
sumstats.basic_check()

sumstats.estimate_h2_cts_by_ldsc(
    ref_ld_chr=args.ld,
    ref_ld_chr_cts=args.ldcts,
    w_ld_chr=args.weights
)

sumstats.ldsc_h2_cts.to_csv(args.output, index=False, sep="\t")
