include { samplesheetToList } from 'plugin/nf-schema'

workflow {
    ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    ch_input
        .map { phenotype, trait_type, population, cohort ->
            def input_file = file("data/raw/${phenotype}/${cohort}.${population}.sumstats.txt.gz")
            if (!input_file.exists()) error "Input file not found: ${input_file}"
            [[phenotype: phenotype, type: trait_type, population: population, cohort: cohort], input_file]
        }
        .show()

   //MUNGE_SUMSTATS(ch_input)     
}