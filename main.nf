include { samplesheetToList } from 'plugin/nf-schema'
include { MUNGE_SUMSTATS } from "${projectDir}/modules/munge_sumstats.nf"

workflow {
    ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    ch_sumstats = ch_input.map { phenotype, trait_type, population, cohort ->
            def input_file = file("data/raw/${phenotype}/${cohort}-${phenotype}-${population}.sumstats.txt.gz")
            if (!input_file.exists()) error "Input file not found: ${input_file}"
            [[phenotype: phenotype, type: trait_type, population: population, cohort: cohort], input_file]
        }

   ch_sumstats_munged = MUNGE_SUMSTATS(ch_sumstats) 

}