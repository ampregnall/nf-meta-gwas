#!/usr/bin/env Rscript

# Load required libraries
box::use(
  vroom[vroom, vroom_write],
  dplyr,
  MungeSumstats[format_sumstats],
  logger[log_info],
  argparse[ArgumentParser],
  BSgenome.Hsapiens.NCBI.GRCh38,
  SNPlocs.Hsapiens.dbSNP144.GRCh38,
  BSgenome.Hsapiens.1000genomes.hs37d5
)

parser <- ArgumentParser(
  description = "Argument parser for munging summary statistics"
)

parser$add_argument(
  "--input",
  type = "character",
  required = TRUE,
  help = "Input summary statistics file"
)

parser$add_argument(
  "--output",
  type = "character",
  required = TRUE,
  help = "Output file"
)

parser$add_argument(
  "--type",
  type = "character",
  required = TRUE,
  help = "Trait type: binary or quantitative"
)

parser$add_argument(
  "--dbsnp",
  type = "character",
  required = TRUE
)

parser$add_argument(
  "--cpus",
  type = "integer",
  default = 1
)

parser$add_argument(
  "--mac",
  type = "integer",
  default = 50
)

args <- parser$parse_args()

# Load summary statistics
log_info("Loading summary statistic file: {args$input}")
df <- vroom(args$input)

log_info("Summary statistic file contains {nrow(df)} rows before munging")
log_info("Running MungeSumstats")

df_formatted <- format_sumstats(
  df,
  ref_genome = "GRCh38",
  #convert_ref_genome = "GRCh38",
  dbSNP = 144,
  dbSNP_tarball = args$dbsnp,
  return_data = TRUE,
  return_format = "data.table",
  nThread = args$cpus
)

log_info(
  "Summary statistic file contains {nrow(df_formatted)} rows after munging"
)

# Apply minor allele count filter
if (args$type == "binary") {
  df_formatted <- df_formatted |>
    mutate(MAF = pmin(FRQ, 1 - FRQ), MAC = MAF * N_CAS * 2) |>
    filter(MAC > args$mac)
} else {
  df_formatted <- df_formatted |>
    mutate(MAF = pmin(FRQ, 1 - FRQ), MAC = MAF * N * 2) |>
    filter(MAC >= args$mac)
}

log_info(
  "Summary statistic file contains {nrow(df_formatted)} rows after applying minimum allele count"
)

vroom_write(df_formatted, args$output, sep = "\t")
