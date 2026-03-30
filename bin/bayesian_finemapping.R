#!/usr/bin/env Rscript

box::use(
  data.table[fread, fwrite],
  dplyr[select, rename],
  glue[glue],
  argparse[ArgumentParser]
)

box::use(./utilities[calc_credset, locus_extract])

parser <- ArgumentParser(description = "Perform Approximate Bayes Factor fine-mapping")
parser$add_argument("--input",      type = "character", required = TRUE,
                    help = "Path to meta-analysis summary statistics (gz)")
parser$add_argument("--lead",       type = "character", required = TRUE,
                    help = "Path to lead variants file")
parser$add_argument("--output",     type = "character", required = TRUE,
                    help = "Output prefix")
parser$add_argument("--phenotype",  type = "character", required = TRUE,
                    help = "Phenotype name")
parser$add_argument("--population", type = "character", required = TRUE,
                    help = "Population label")
args <- parser$parse_args()

sumstats <- fread(args$input)
lead     <- fread(args$lead)

# Extract all variants within 500 kb of each lead locus
sumstats_loci <- locus_extract(
  sumstats,
  sumstats_chr_col = CHR,
  sumstats_pos_col = POS,
  locus_df         = lead,
  locus_chr_col    = CHR,
  locus_pos_col    = POS,
  locus_gene_col   = GENENAME,
  locus_size       = 1e6
)

credset <- calc_credset(
  sumstats_loci,
  locus_marker_col = locus_marker,
  effect_col       = BETA,
  se_col           = SE,
  samplesize_col   = N,
  cred_interval    = 0.99
)

credset <- credset |>
  dplyr::select(
    SNPID, rsID,
    LOCUS        = locus_marker,
    CHR          = chromosome,
    POS          = position,
    EA, NEA, BETA, SE,
    NEAREST_GENE = gene,
    BF           = bf,
    BF_PIP       = posterior_prob
  )

fwrite(credset, glue("{args$output}.credset.txt"), sep = "\t")
