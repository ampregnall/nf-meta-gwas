process MUNGE_SUMSTATS {
    cpus 4
    // conda '${projectDir}/envs/munge.yaml'
    
    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.munged.sumstats.gz"), emit: sumstats_munged
        // tuple val(meta), path("*.statistics.txt"),     emit: statistics

    script:
    """
    munge_sumstats.R \
        --input ${sumstats} \
        --output ${meta.phenotype}-${meta.cohort}-${meta.population} \
        --type ${meta.type} \
        --dbsnp ${params.dbsnp} \
        --cpus ${task.cpus} \
        --mac ${params.mac}
    """
}