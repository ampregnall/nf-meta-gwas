#!/usr/bin/env Rscript

box::use(
  vroom[vroom, vroom_write],
  dplyr[filter, select, mutate, distinct, pull, group_by, ungroup,
        bind_rows, between, n, everything],
  tidyr[pivot_wider, drop_na],
  tibble[column_to_rownames],
  purrr[list_rbind],
  glue[glue],
  argparse[ArgumentParser],
  arrow[open_dataset, collect],
  logger[log_info, log_warn]
)

box::use(./utilities[locus_extract])

# =============================================================================
# eQTL parquet column mapping
# UPDATE THESE NAMES to match your parquet file schema before running.
# Use pyarrow or R arrow::read_parquet() with nrow=1 to inspect column names.
# =============================================================================
EQTL_CHR  <- "chromosome"          # chromosome (string: "1", "2", ... "22")
EQTL_POS  <- "position"            # genomic position (integer)
EQTL_GENE <- "molecular_trait_id"  # gene / trait identifier
EQTL_RSID <- "rsid"               # rsID
EQTL_BETA <- "beta"               # effect size
EQTL_SE   <- "se"                 # standard error
# =============================================================================

parser <- ArgumentParser(description = "HyPrColoc eQTL colocalization")
parser$add_argument("--sumstats",     required = TRUE, help = "Meta-analysis sumstats (.txt.gz)")
parser$add_argument("--lead",         required = TRUE, help = "Lead variants file")
parser$add_argument("--eqtl_parquet", required = TRUE, help = "eQTL catalog parquet file or dir")
parser$add_argument("--output",       required = TRUE, help = "Output prefix")
parser$add_argument("--phenotype",    required = TRUE)
parser$add_argument("--population",   required = TRUE)
parser$add_argument("--eqtl_study",   required = TRUE, help = "Study label (from filename)")
args <- parser$parse_args()

# ---------------------------------------------------------------------------
# Empty output helpers
# ---------------------------------------------------------------------------
empty_results <- function() {
  data.frame(
    phenotype = character(), eqtl_study = character(), locus = character(),
    eqtl_gene = character(), traits = character(), posterior_prob = double(),
    regional_prob = double(), candidate_snp = character(),
    posterior_explained_by_snp = double(), dropped_trait = character()
  )
}

empty_failures <- function() {
  data.frame(
    phenotype = character(), eqtl_study = character(), locus = character(),
    eqtl_gene = character(), n_gwas_snps = integer(), n_eqtl_snps = integer(),
    n_overlap_snps = integer(), error = character()
  )
}

write_outputs <- function(results, failures, prefix) {
  vroom_write(results,  glue("{prefix}.hyprcoloc.tsv"),          delim = "\t")
  vroom_write(failures, glue("{prefix}.hyprcoloc.failures.tsv"), delim = "\t")
}

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
log_info("Loading lead variants: {args$lead}")
lead <- vroom(args$lead, show_col_types = FALSE) |>
  filter(!is.na(POS), !is.na(CHR))

if (nrow(lead) == 0) {
  log_info("No lead variants — writing empty outputs")
  write_outputs(empty_results(), empty_failures(), args$output)
  quit(status = 0)
}

log_info("Loading GWAS sumstats: {args$sumstats}")
sumstats_raw <- vroom(
  args$sumstats,
  col_select = c(rsID, CHR, POS, BETA, SE, EA, NEA),
  show_col_types = FALSE
) |>
  filter(!is.na(rsID), !is.na(BETA), !is.na(SE)) |>
  mutate(trait = args$phenotype)

# Extract GWAS variants within all locus windows in one pass
sumstats_loci <- locus_extract(
  sumstats_raw,
  sumstats_chr_col = CHR,
  sumstats_pos_col = POS,
  locus_df         = lead,
  locus_chr_col    = CHR,
  locus_pos_col    = POS,
  locus_gene_col   = GENE,
  locus_size       = 1e6
)

log_info("Opening eQTL dataset: {args$eqtl_parquet}")
eqtl_ds <- open_dataset(args$eqtl_parquet)

# ---------------------------------------------------------------------------
# Process each locus
# ---------------------------------------------------------------------------
results_list  <- list()
failures_list <- list()

locus_ids <- unique(sumstats_loci$locus_marker)
log_info("{length(locus_ids)} loci to process")

for (locus_id in locus_ids) {
  gwas_locus <- sumstats_loci |> filter(locus_marker == locus_id)
  locus_parts <- strsplit(locus_id, ":")[[1]]
  chr_val   <- as.integer(locus_parts[1])
  pos_val   <- as.integer(locus_parts[2])
  start_val <- pos_val - 500000L
  end_val   <- pos_val + 500000L

  valid_rsids <- gwas_locus |> distinct(rsID) |> pull(rsID)
  n_gwas      <- length(valid_rsids)

  # eQTL window — predicate pushdown via arrow row-group statistics
  eqtl_locus <- eqtl_ds |>
    filter(
      .data[[EQTL_CHR]] == as.character(chr_val),
      .data[[EQTL_POS]] >= start_val,
      .data[[EQTL_POS]] <= end_val
    ) |>
    select(
      gene = .data[[EQTL_GENE]],
      rsID = .data[[EQTL_RSID]],
      BETA = .data[[EQTL_BETA]],
      SE   = .data[[EQTL_SE]]
    ) |>
    collect() |>
    filter(!is.na(rsID), rsID %in% valid_rsids, !is.na(BETA), !is.na(SE))

  if (nrow(eqtl_locus) == 0) {
    log_warn("Locus {locus_id}: no overlapping eQTL variants — skipping")
    next
  }

  genes <- unique(eqtl_locus$gene)
  log_info("Locus {locus_id}: {length(genes)} eQTL genes with overlapping variants")

  for (gene_id in genes) {
    gene_df <- eqtl_locus |>
      filter(gene == gene_id) |>
      distinct(rsID, .keep_all = TRUE) |>
      mutate(trait = gene_id)

    n_eqtl <- nrow(gene_df)

    combined <- bind_rows(
      gwas_locus |> select(rsID, BETA, SE, trait),
      gene_df    |> select(rsID, BETA, SE, trait)
    )

    .betas <- combined |>
      distinct(rsID, trait, .keep_all = TRUE) |>
      select(rsID, BETA, trait) |>
      pivot_wider(names_from = trait, values_from = BETA) |>
      column_to_rownames("rsID") |>
      drop_na() |>
      as.matrix()

    .ses <- combined |>
      distinct(rsID, trait, .keep_all = TRUE) |>
      select(rsID, SE, trait) |>
      pivot_wider(names_from = trait, values_from = SE) |>
      column_to_rownames("rsID") |>
      drop_na() |>
      as.matrix()

    n_overlap <- nrow(.betas)

    result <- tryCatch({
      hyprcoloc::hyprcoloc(
        .betas, .ses,
        trait.names = colnames(.betas),
        snp.id      = rownames(.betas)
      )[[1]]
    }, error = function(e) {
      msg <- conditionMessage(e)
      log_warn("HyPrColoc failed — locus {locus_id} / gene {gene_id}: {msg}")
      failures_list[[length(failures_list) + 1L]] <<- data.frame(
        phenotype      = args$phenotype,
        eqtl_study     = args$eqtl_study,
        locus          = locus_id,
        eqtl_gene      = gene_id,
        n_gwas_snps    = n_gwas,
        n_eqtl_snps    = n_eqtl,
        n_overlap_snps = n_overlap,
        error          = msg,
        stringsAsFactors = FALSE
      )
      NULL
    })

    if (!is.null(result) && nrow(result) > 0) {
      result$phenotype  <- args$phenotype
      result$eqtl_study <- args$eqtl_study
      result$locus      <- locus_id
      result$eqtl_gene  <- gene_id
      results_list[[length(results_list) + 1L]] <- result
    }
  }
}

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
results_out  <- if (length(results_list)  > 0) {
  list_rbind(results_list) |> select(phenotype, eqtl_study, locus, eqtl_gene, everything())
} else {
  empty_results()
}

failures_out <- if (length(failures_list) > 0) {
  list_rbind(failures_list)
} else {
  empty_failures()
}

write_outputs(results_out, failures_out, args$output)
log_info("Done: {nrow(results_out)} colocalization results, {nrow(failures_out)} failures logged")
