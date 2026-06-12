include { samplesheetToList } from 'plugin/nf-schema'
include { MUNGE_SUMSTATS } from "${projectDir}/modules/munge_sumstats.nf"
include { LDSC_CORRECTION } from "${projectDir}/modules/ldsc_correction.nf"
include { COLLECT_LDSC_RESULTS } from "${projectDir}/modules/collect_ldsc.nf"
include { META_ANALYZE } from "${projectDir}/modules/meta_analysis.nf"
include { COLLECT_META_RESULTS } from "${projectDir}/modules/collect_meta.nf"
include { EXTRACT_LEAD_VARIANTS } from "${projectDir}/modules/lead_variants.nf"
include { HERITABILITY } from "${projectDir}/modules/heritability.nf"
include { ABF_FINEMAPPING } from "${projectDir}/modules/abf_finemapping.nf"
include { COLLECT_FILTER_STATS } from "${projectDir}/modules/collect_filter_stats.nf"
include { TISSUE_ENRICHMENT } from "${projectDir}/modules/tissue_enrichment.nf"
include { COLLECT_LEAD_VARIANTS } from "${projectDir}/modules/collect_lead_variants.nf"
include { PLOT_GWAS as PLOT_INPUT_GWAS  } from "${projectDir}/modules/plot_gwas.nf"
include { PLOT_GWAS as PLOT_META_GWAS  } from "${projectDir}/modules/plot_gwas.nf"
include { QC_INPUT_GWAS                } from "${projectDir}/modules/qc_input_gwas.nf"
include { QC_COHORT_SUMMARY            } from "${projectDir}/modules/qc_cohort_summary.nf"
include { QC_META_GWAS                 } from "${projectDir}/modules/qc_meta_gwas.nf"
include { PREPARE_MAGMA_INPUT          } from "${projectDir}/modules/prepare_magma_input.nf"
include { MAGMA_GENE                   } from "${projectDir}/modules/magma_gene.nf"
include { MAGMA_TISSUE                 } from "${projectDir}/modules/magma_tissue.nf"
include { POPS                         } from "${projectDir}/modules/pops.nf"
include { PREPARE_FLAMES_INPUTS        } from "${projectDir}/modules/prepare_flames_inputs.nf"
include { FLAMES                       } from "${projectDir}/modules/flames.nf"

workflow {
    ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    ch_sumstats = ch_input.map { phenotype, trait_type, population, cohort, col_overrides ->
            def input_file = file("data/raw/${phenotype}/${cohort}-${phenotype}-${population}.sumstats.txt.gz")
            if (!input_file.exists()) error "Input file not found: ${input_file}"
            [[phenotype: phenotype, type: trait_type, population: population, cohort: cohort, col_overrides: col_overrides ?: ""], input_file]
        }

    ch_sumstats_munged = MUNGE_SUMSTATS(ch_sumstats)
    ch_ldsc_out = LDSC_CORRECTION(ch_sumstats_munged.sumstats_munged)

    // Gather filter stats per phenotype and publish combined file
    ch_filter_stats_collected = COLLECT_FILTER_STATS(
        ch_sumstats_munged.filter_stats
            .map { meta, txt -> tuple([phenotype: meta.phenotype], txt) }
            .groupTuple(by: 0)
    )

    // Gather LDSC h2 results per phenotype and publish combined file
    ch_ldsc_collected = COLLECT_LDSC_RESULTS(
        ch_ldsc_out.ldsc_h2
            .map { meta, txt -> tuple([phenotype: meta.phenotype], txt) }
            .groupTuple(by: 0)
    )

    ch_ldsc_corrected = ch_ldsc_out.sumstats_parquet

    ch_chroms = Channel.of(1..22)

    ch_across_populations = ch_ldsc_corrected
        .map { meta, sumstats ->
            tuple([phenotype: meta.phenotype, population: "all", type: meta.type], sumstats)
        }
        .groupTuple(by: 0)
        .combine(ch_chroms)

    ch_within_populations = ch_ldsc_corrected
        .map { meta, sumstats ->
            tuple([phenotype: meta.phenotype, population: meta.population, type: meta.type], sumstats)
        }
        .groupTuple(by: 0)
        .combine(ch_chroms)

    ch_meta = META_ANALYZE(ch_across_populations.mix(ch_within_populations))
    ch_collected_meta = COLLECT_META_RESULTS(
        ch_meta
            .map { meta, txt ->
                tuple([phenotype: meta.phenotype, population: meta.population], txt)
            }
            .groupTuple(by: 0)
    )

    ch_lead = EXTRACT_LEAD_VARIANTS(ch_collected_meta.sumstats)

    // Aggregate lead variants across all populations per phenotype, cluster into loci, flag population-specific hits
    ch_lead.lead_variants
        .map { meta, txt -> tuple([phenotype: meta.phenotype], txt) }
        .groupTuple(by: 0)
        | COLLECT_LEAD_VARIANTS

    // Heritability estimation on per-population meta results only (ldsc_reference is per-population)
    HERITABILITY(
        ch_collected_meta.sumstats.filter { it[0].population != "all" }
    )

    // ABF fine-mapping: join collected meta sumstats with lead variants
    // Filter out cases with no GWS hits (header-only lead variants file = 1 line)
    ch_abf = ABF_FINEMAPPING(
        ch_collected_meta.sumstats
            .join(ch_lead.lead_variants, by: 0)
            .filter { meta, sumstats, leads -> leads.countLines() > 1 }
    )

    // Tissue and cell type enrichment on all meta-analysis results
    TISSUE_ENRICHMENT(ch_collected_meta.sumstats)

    // Manhattan and QQ plots for each input GWAS (post-munging)
    PLOT_INPUT_GWAS(ch_sumstats_munged.sumstats_munged)

    // Manhattan and QQ plots for each meta-analysis result
    PLOT_META_GWAS(ch_collected_meta.sumstats)

    // ---------------------------------------------------------------------------
    // QC figures (Issue #4)
    // ---------------------------------------------------------------------------

    // Per-cohort QC plots + summary stats for cross-cohort comparison
    ch_qc_input = QC_INPUT_GWAS(ch_sumstats_munged.sumstats_munged)

    // Join per-phenotype collected tables with grouped per-cohort summaries
    ch_qc_summaries_by_phenotype = ch_qc_input.qc_summary
        .map { meta, tsv -> [[phenotype: meta.phenotype], tsv] }
        .groupTuple(by: 0)

    QC_COHORT_SUMMARY(
        ch_filter_stats_collected.filter_stats
            .join(ch_ldsc_collected.ldsc_h2,          by: 0)
            .join(ch_qc_summaries_by_phenotype,        by: 0)
    )

    // Meta-analysis heterogeneity and N_CONTRIBUTIONS QC plots
    QC_META_GWAS(ch_collected_meta.sumstats)

    // ---------------------------------------------------------------------------
    // MAGMA / PoPS / FLAMES gene-prioritisation pipeline — opt-in
    //
    // Set params.magma_bfile (and related MAGMA/PoPS params) to enable MAGMA
    // and PoPS. Set params.flames_annotation_data to additionally enable FLAMES.
    // Leave any param null to skip that tier. The core pipeline (harmonisation,
    // meta-analysis, fine-mapping) always runs regardless of these params.
    // ---------------------------------------------------------------------------

    if (params.magma_bfile) {
        ch_magma_input = PREPARE_MAGMA_INPUT(ch_collected_meta.sumstats)

        ch_magma = MAGMA_GENE(ch_magma_input.magma_input)

        // MAGMA tissue expression and PoPS both branch from the genes.raw output
        ch_magma_tissue = MAGMA_TISSUE(ch_magma.genes_raw)
        ch_pops         = POPS(ch_magma.genes_raw)

        // FLAMES additionally requires annotation data; ch_abf is always defined
        // (gated upstream by GWS hits) so this reference is always safe
        if (params.flames_annotation_data) {
            ch_flames_prep = PREPARE_FLAMES_INPUTS(ch_abf.credset)

            ch_for_flames = ch_flames_prep.flames_inputs
                .join(ch_pops.pops_preds,      by: 0)
                .join(ch_magma.genes_out,      by: 0)
                .join(ch_magma_tissue.gsa_out, by: 0)

            FLAMES(ch_for_flames)
        }
    }
}
