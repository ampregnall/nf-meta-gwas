process MUNGE_SUMSTATS {
    cpus 4
    container 'ghcr.io/ampregnall/nf-meta-gwas:0.1.0'

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.munged.sumstats.gz"), emit: sumstats_munged

    script:
    def prefix = "${meta.phenotype}-${meta.cohort}-${meta.population}"
    """
    test_gwaslab.py --output ${prefix}
    """
}