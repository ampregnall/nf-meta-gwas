#' @export
calc_credset <- function(df, 
                         locus_marker_col = locus_marker, 
                         effect_col = effect, 
                         se_col = std_err, 
                         samplesize_col = samplesize, 
                         cred_interval = 0.99) {
  df %>%
    dplyr::group_by({{ locus_marker_col }}) %>%
    dplyr::mutate(bf = exp(0.5 * ({{ effect_col }}^2 / {{ se_col }}^2 - log({{ samplesize_col }})))) %>%
    dplyr::mutate(posterior_prob = bf / sum(bf)) %>%
    dplyr::arrange(dplyr::desc(posterior_prob)) %>%
    dplyr::mutate(cum_sum = cumsum(posterior_prob)) %>%
    dplyr::group_by({{ locus_marker_col }}) %>%
    dplyr::filter(cum_sum <= cred_interval | posterior_prob > cred_interval) %>%
    dplyr::select(-cum_sum) %>%
    dplyr::ungroup()
}

#' @export
locus_extract <- function(sumstats_df, # df containing summary statistics
                          sumstats_chr_col = chromosome,
                          sumstats_pos_col = position,
                          locus_df, # df containing gws loci
                          locus_chr_col = chromosome,
                          locus_pos_col = position,
                          locus_gene_col = gene,
                          locus_size = 1e6) {
  
  locus_df <- locus_df %>%
    dplyr::select(chromosome = {{ locus_chr_col }}, lead_pos = {{ locus_pos_col }}, gene = {{ locus_gene_col }}) %>%
    dplyr::distinct(chromosome, lead_pos, .keep_all = TRUE)
  
  sumstats_df <- sumstats_df %>%
    dplyr::rename(chromosome = {{ sumstats_chr_col }}, position = {{ sumstats_pos_col }}) %>%
    dplyr::mutate(chromosome = as.numeric(chromosome), position = as.numeric(position)) %>%
    dplyr::inner_join(locus_df, by = "chromosome", relationship = "many-to-many") %>%
    dplyr::filter(between(position, lead_pos - locus_size / 2, lead_pos + locus_size / 2)) %>%
    dplyr::collect() %>%
    dplyr::mutate(locus_marker = glue::glue("{chromosome}:{lead_pos}")) %>%
    dplyr::select(locus_marker, chromosome, lead_pos, gene, everything())
  
  return(sumstats_df)
}