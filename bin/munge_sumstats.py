#!/usr/bin/env python

import argparse
import gwaslab as gl
import numpy as np
from scipy.stats import norm

parser = argparse.ArgumentParser(description="Test gwaslab installation")
parser.add_argument("--input", type=str, required=True, help="Path to input file")
parser.add_argument("--output", type=str, required=True, help="Output filenmame")
parser.add_argument("--type", type=str, required=True, help="Trait type (binary or continuous)")
parser.add_argument("--mac", type=int, required=True, help="Minimum minor allele count")
parser.add_argument("--chain", type=str, required=True, help="Path to chainfile")
parser.add_argument("--fasta", type=str, required=True, help="Path to UCSD Human reference genome (GRCh38")
parser.add_argument("--dbsnp", type=str, required=True, help="Path to dbSNP VCF file for rsID assignment")
parser.add_argument("--popvcf", type=str, required=True, help = "Path to population specific VCF file for strand inference")
parser.add_argument("--ldsc", type=str, required=True, help = "Path to Pan-UKBB LD reference panels")
parser.add_argument("--threads", type=int, required=True, help = "Number of threads")
args = parser.parse_args()


sumstats = gl.Sumstats(args.input, fmt="gwaslab")
sumstats.infer_build(verbose=True)

# Liftover if statistics are on hg19
if sumstats.build == '19':
    sumstats.liftover(to_build='38', from_build='19', chain_path=args.chain)
    
# Harmonize summary statistics
sumstats.harmonize(
    ref_seq=args.fasta,
    ref_rsid_vcf=args.dbsnp, 
    ref_infer=args.popvcf, # Ancestry specific. Logic handled by Nextflow
    ref_alt_freq="AF", 
    threads=args.threads, # pass threads from nextflow process, 
    sweep_mode=True 
    )

sumstats.data.to_csv("debug.txt.gz", index = False, compression="gzip")

# Perform LDSC correction
sumstats_hapmap3 = sumstats.filter_hapmap3(inplace=False)
sumstats_hapmap3.estimate_h2_by_ldsc(ref_ld = args.ldsc,  w_ld = args.ldsc)

if np.float64(sumstats_hapmap3.ldsc_h2['Intercept'][0]) > 1:
    sumstats.data['SE'] = sumstats.data['SE'] * np.sqrt(np.float64(sumstats_hapmap3.ldsc_h2['Intercept'][0]))
    sumstats.data['Z'] = sumstats.data['BETA'] / sumstats.data['SE']
    sumstats.data['P'] = 2 * norm.sf(abs(sumstats.data['Z']))
    
    
if args.type == "binary":
  sumstats.data['MAC'] = sumstats.data['MAF'] * sumstats.data['N_CASE'] * 2
else:
  sumstats.data['MAC'] = sumstats.data['MAF'] * sumstats.data['N'] * 2

sumstats = sumstats.filter_value('MAC > @args.mac', inplace=True)

# Save results
parquet_out = f"{args.output}.sumstats.processed.parquet"
sumstats_out = f"{args.output}.sumstats.txt.gz"

sumstats.data.to_parquet(parquet_out)
sumstats.data.to_csv(sumstats_out, index = False, compression="gzip")
