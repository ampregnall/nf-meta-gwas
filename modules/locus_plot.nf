process LOCUS_PLOT {
    cpus 4
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=${32000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    publishDir { "${launchDir}/results/plots/locus/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta), path(locus_summary), path(sumstats_files)

    output:
        tuple val(meta), path("*.locus.png"), optional: true, emit: png
        tuple val(meta), path("*.locus.pdf"), optional: true, emit: pdf

    script:
    def vcf_map = params.population_vcf + [all: params.population_vcf.eur]
    def vcf_map_json = groovy.json.JsonOutput.toJson(vcf_map)
    """
    locus_plot.py \
        --locus_summary ${locus_summary} \
        --sumstats ${sumstats_files} \
        --phenotype ${meta.phenotype} \
        --gtf ${params.gtf} \
        --population_vcf_map '${vcf_map_json}'
    """
}
