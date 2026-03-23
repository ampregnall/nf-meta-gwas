#!/usr/bin/env Rscript

# Load required libraries
box::use(
  vroom[vroom, vroom_write],
  dplyr,
  MungeSumstats[format_sumstats],
  logger[log_info],
  argparse[ArgumentParser]
)

parser <- ArgumentParser(
  description = "Argument parser for munging summary statistics"
)

parser$add_argument(
  "--input",
  type = "string",
  required = TRUE,
  help = "Input summary statistics file"
)

parser$add_argument(
  "--output",
  type = "string",
  required = TRUE,
  help = "Output file"
)

parser$add_argument(
  "--type",
  type = "string",
  required = TRUE,
  help = "Trait type: binary or quantitative"
)

parser$add_argument(
  "--dbsnp",
  type = "string",
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
  convert_ref_genome = "GRCh38",
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
    mutate(MAF = pmin(FRQ, 1 - FRQ), MAC = MAC * N_CASE * 2) |>
    filter(MAC > args$mac)
} else {
  df_formatted <- df_formatted |>
    mutate(MAF = pmin(FRQ, 1 - FRQ), MAC = MAC * N * 2) |>
    filter(MAC >= args$mac)
}

log_info(
  "Summary statistic file contains {nrow(df_formatted)} rows after applying minimum allele count"
)

vroom_write(df_formatted, args$output, sep = "\t")
