process MUNGE_SUMSTATS {
    cpus 4
    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.munged.sumstats.gz"), emit: sumstats_munged
        // tuple val(meta), path("*.statistics.txt"),     emit: statistics

    script:
    """
    Rscript munge_sumstats.R \
        --input ${sumstats} \
        --output ${meta.cohort}.${meta.population} \
        --type ${meta.type} \
        --dbsnp ${params.dbsnp_tarball} \
        --cpus ${cpus}
        --mac ${params.mac}
    """
}