#!/usr/bin/env Rscript

box::use(
  vroom[vroom, vroom_write],
  arrow,
  dplyr[
    filter,
    mutate,
    collect,
    group_by,
    summarise,
    select,
    if_else,
    across,
    arrange,
    any_of,
    all_of, 
    rename,
    everything
  ],
  glue[glue],
  logger[log_info, log_warn],
  argparse[ArgumentParser]
)


# -----------------------------------------------------------------------------
# 1. Argument parsing
# -----------------------------------------------------------------------------
parser <- ArgumentParser(
  description = "Fixed-effects IVW meta-analysis for one chromosome"
)

parser$add_argument(
  "--chr",
  type = "integer",
  required = TRUE,
  metavar = "INT",
  help = "Chromosome to process (e.g. 1)"
)

parser$add_argument(
  "--type",
  dest = "trait_type",
  type = "character",
  required = TRUE,
  choices = c("binary", "continuous"),
  metavar = "STR",
  help = "Trait type: 'binary' or 'continuous' [required]"
)

parser$add_argument(
  "--output",
  dest = "out_prefix",
  type = "character",
  default = "meta",
  metavar = "STR",
  help = "Output file prefix [default: meta]"
)

args <- parser$parse_args()

# -----------------------------------------------------------------------------
# 2. Hard-coded column schema
# -----------------------------------------------------------------------------
CHR_COL <- "CHR"
POS_COL <- "POS"
RSID_COL <- "rsID"
EA_COL <- "EA"
NEA_COL <- "NEA"
EAF_COL <- "EAF"
BETA_COL <- "BETA"
SE_COL <- "SE"
N_COL <- "N"
N_CASE_COL <- "N_CASE"
N_CONTROL_COL <- "N_CONTROL"

by_cols <- c(CHR_COL, POS_COL, RSID_COL, EA_COL, NEA_COL)
eaf_n_col <- paste0("N_", EAF_COL) # temporary column for EAF weighting

if (args$trait_type == "binary") {
  n_cols_to_sum <- c(EAF_COL, N_CASE_COL, N_CONTROL_COL, N_COL)
} else {
  n_cols_to_sum <- c(EAF_COL, N_COL)
}

# -----------------------------------------------------------------------------
# 3. Discover input files
# -----------------------------------------------------------------------------
parquet_files <- list.files(".", pattern = args$pattern, full.names = TRUE)
log_info("Found {length(parquet_files)} parquet files")

# -----------------------------------------------------------------------------
# 4. Open dataset and run IVW meta-analysis
# -----------------------------------------------------------------------------
log_info("Processing chromosome {args$chr} ({args$trait_type})")

parquet_schema <- arrow::open_dataset(parquet_files[1])$schema

# Identify any dictionary-encoded fields and override them to utf8
unified_fields <- lapply(parquet_schema$fields, function(field) {
  if (field$type$id == arrow::Type$DICTIONARY) {
    arrow::field(field$name, arrow::utf8())
  } else {
    field
  }
})

unified_schema <- arrow::schema(unified_fields)
ds <- arrow::open_dataset(parquet_files, schema = unified_schema)

meta_results <- ds |>
  filter(.data[[CHR_COL]] == args$chr) |>
  filter(
    !is.na(.data[[BETA_COL]]),
    !is.na(.data[[SE_COL]]),
    .data[[SE_COL]] > 0
  ) |>
  mutate(
    W = 1 / (.data[[SE_COL]])^2,
    B = .data[[BETA_COL]] * W,
    WB2 = (B^2) / W,
    # Weight EAF by N before grouping so we can compute a N-weighted mean later
    across(
      any_of(EAF_COL),
      ~ .x * .data[[N_COL]]
    )
  ) |>
  group_by(across(all_of(by_cols))) |>
  summarise(
    N_CONTRIBUTIONS = n(),
    # Accumulate total N used for EAF weighting
    across(
      any_of(EAF_COL),
      ~ sum(if_else(is.na(.x), 0, .data[[N_COL]]), na.rm = TRUE),
      .names = "N_{.col}"
    ),
    across(any_of(POS_COL), ~ min(.x, na.rm = TRUE)),
    WB2 = sum(WB2, na.rm = TRUE),
    across(any_of(c("W", "B")), ~ sum(.x, na.rm = TRUE)),
    across(any_of(n_cols_to_sum), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  collect() |>
  mutate(
    B_w = B,
    B = B / W,
    SE = 1 / sqrt(W),
    # Recover N-weighted mean EAF
    across(
      any_of(EAF_COL),
      ~ .x / .data[[eaf_n_col]]
    )
  ) |>
  select(-any_of(eaf_n_col)) |>
  mutate(
    Z = B / SE,
    P = 2 * stats::pnorm(-abs(Z)),
    # Cochran's Q heterogeneity
    Q = WB2 - (B_w^2) / W,
    Q_DF = pmax(N_CONTRIBUTIONS - 1L, 0L),
    Q_PVAL = stats::pchisq(Q, df = Q_DF, lower.tail = FALSE),
    Q_PVAL = if_else(Q_DF == 0L, 1, Q_PVAL),
    I2 = if_else(Q > 0, pmax((Q - Q_DF) / Q, 0), 0)
  )

if (nrow(meta_results) == 0) {
  log_warn("No variants remaining after filtering for chromosome {args$chr}")
}

# -----------------------------------------------------------------------------
# 5. Final column selection and ordering
# -----------------------------------------------------------------------------
meta_results <- meta_results |> arrange(CHR, POS)

if (args$trait_type == "binary") {
  meta_results <- meta_results |>
    select(
      CHR,
      POS = all_of(POS_COL),
      rsID = all_of(RSID_COL),
      EA = all_of(EA_COL),
      NEA = all_of(NEA_COL),
      B,
      SE,
      P,
      Z,
      EAF = all_of(EAF_COL),
      N_CONTRIBUTIONS,
      N_CASE,
      N_CONTROL,
      N,
      Q,
      Q_DF,
      Q_PVAL,
      I2
    )
} else {
  meta_results <- meta_results |>
    select(
      CHR,
      POS = all_of(POS_COL),
      rsID = all_of(RSID_COL),
      EA = all_of(EA_COL),
      NEA = all_of(NEA_COL),
      B,
      SE,
      P,
      Z,
      EAF = all_of(EAF_COL),
      N_CONTRIBUTIONS,
      N = all_of(N_COL),
      Q,
      Q_DF,
      Q_PVAL,
      I2
    )
}

meta_results <- meta_results |>
  rename(BETA = B) |> 
  mutate(SNPID = glue("{CHR}:{POS}:{NEA}:{EA}")) |> 
  select(SNPID, everything())

# -----------------------------------------------------------------------------
# 6. Write output
# -----------------------------------------------------------------------------
out_file <- glue("{args$out_prefix}.chr{args$chr}.txt")
vroom_write(meta_results, out_file)
