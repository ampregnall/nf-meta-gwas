#!/usr/bin/env Rscript

box::use(
  vroom[vroom, vroom_write],
  dplyr[mutate, arrange, filter, distinct, row_number, select, desc, group_by, ungroup],
  glue[glue],
  purrr[walk2],
  tibble[tibble],
  argparse[ArgumentParser]
)

parser <- ArgumentParser(
  description = "Prepare FLAMES inputs from ABF fine-mapping credible sets"
)
parser$add_argument(
  "--credset",
  required = TRUE,
  help = "Credible set file from ABF_FINEMAPPING (TSV)"
)
parser$add_argument(
  "--output",
  required = TRUE,
  help = "Output file prefix"
)
args <- parser$parse_args()

credset <- vroom(args$credset, show_col_types = FALSE)

# FLAMES requires CHR:POS:EA:NEA — the reverse of our pipeline convention
credset <- credset |>
  mutate(SNPID_FLAMES = glue("{CHR}:{POS}:{EA}:{NEA}"))

# Derive ordered loci from the LOCUS column (format: "{chr}:{lead_pos}")
locus_parts  <- do.call(rbind, strsplit(unique(credset$LOCUS), ":", fixed = TRUE))
loci_ordered <- data.frame(
  LOCUS = unique(credset$LOCUS),
  chr   = as.integer(locus_parts[, 1]),
  pos   = as.integer(locus_parts[, 2]),
  stringsAsFactors = FALSE
) |>
  arrange(chr, pos) |>
  mutate(GenomicLocus = row_number())

# Per-locus credset files (index, cred1, prob1) in FLAMES format
dir.create("credsets", showWarnings = FALSE)

walk2(loci_ordered$LOCUS, loci_ordered$GenomicLocus, function(loc, locus_num) {
  locus_credset <- credset |>
    filter(LOCUS == loc) |>
    arrange(desc(BF_PIP)) |>
    mutate(index = row_number()) |>
    select(index, cred1 = SNPID_FLAMES, prob1 = BF_PIP)
  vroom_write(locus_credset, glue("credsets/{loc}.txt"), delim = "\t")
})

# Genomic locus file: GenomicLocus, chr, start, end (±500 kb window)
genomic_locus <- loci_ordered |>
  mutate(start = pos - 500000L, end = pos + 500000L) |>
  select(GenomicLocus, chr, start, end)

vroom_write(genomic_locus, glue("{args$output}.genomic-locus.txt"), delim = "\t")

# Index file: relative paths so they resolve correctly inside the Nextflow
# work directory when staged for the FLAMES process
df_index <- tibble(
  Filename    = paste0("credsets/", loci_ordered$LOCUS, ".txt"),
  GenomicLocus = loci_ordered$GenomicLocus,
  Annotfiles  = paste0("flames-annotate/FLAMES_annotated_", loci_ordered$LOCUS, ".txt")
)

vroom_write(df_index, glue("{args$output}.index.txt"), delim = "\t")
