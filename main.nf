include { samplesheetToList } from 'plugin/nf-schema'
include { MUNGE_SUMSTATS } from "${projectDir}/modules/munge_sumstats.nf"
include { LDSC_CORRECTION } from "${projectDir}/modules/ldsc_correction.nf"
include { COLLECT_LDSC_RESULTS } from "${projectDir}/modules/collect_ldsc.nf"
include { META_ANALYZE } from "${projectDir}/modules/meta_analysis.nf"

workflow {
    ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    ch_sumstats = ch_input.map { phenotype, trait_type, population, cohort ->
            def input_file = file("data/raw/${phenotype}/${cohort}-${phenotype}-${population}.sumstats.txt.gz")
            if (!input_file.exists()) error "Input file not found: ${input_file}"
            [[phenotype: phenotype, type: trait_type, population: population, cohort: cohort], input_file]
        }

    ch_sumstats_munged = MUNGE_SUMSTATS(ch_sumstats)
    ch_ldsc_out = LDSC_CORRECTION(ch_sumstats_munged)

    // Gather LDSC h2 results per phenotype and publish combined file
    ch_ldsc_out.ldsc_h2
        .map { meta, txt -> tuple([phenotype: meta.phenotype], txt) }
        .groupTuple(by: 0)
        | COLLECT_LDSC_RESULTS

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
}
