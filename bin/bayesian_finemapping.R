library(data.table)
library(tidyverse)

# Load custom scripts
source(snakemake@params[["creds"]])
source(snakemake@params[["locus"]])

sumstats <- fread(snakemake@input[["sumstats"]])
lead <- fread(snakemake@input[["lead"]])

# Extract all variants within 500 kb of lead variants
sumstats <- locus_extract(sumstats,
                          sumstats_chr_col = CHR, 
                          sumstats_pos_col = POS, 
                          locus_df = lead,
                          locus_chr_col = CHR,
                          locus_pos_col = POS,
                          locus_gene_col = GENE,
                          locus_size = 1e6)

credset <- calc_credset(sumstats, 
                        locus_marker_col = locus_marker, 
                        effect_col = BETA, 
                        se_col = SE, 
                        samplesize_col = N, 
                        cred_interval = 0.99)

credset <- credset %>% dplyr::select(SNPID, rsID, LOCUS = locus_marker, 
                                     CHR = chromosome, POS = position,
                                     EA, NEA, BETA, SE, NEAREST_GENE = gene,
                                     BF = bf, BF_PIP = posterior_prob)

fwrite(credset, snakemake@output[["credset"]], sep = "\t")