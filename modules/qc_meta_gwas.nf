process QC_CREDSET_SIZES {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = "-R \"rusage[mem=4000]\""
    publishDir { "${launchDir}/results/plots/qc/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta), path(credset_files)

    output:
        tuple val(meta), path("*.{png,pdf}")

    script:
    def prefix = "${meta.phenotype}"
    """
    qc_meta_gwas.py \
        --credsets ${credset_files} \
        --phenotype ${meta.phenotype} \
        --output ${prefix}
    """
}

process QC_META_GWAS {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=8000]\"" }

    publishDir(
        { "${launchDir}/results/plots/qc/${meta.phenotype}" },
        mode: 'copy'
    )

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.{png,pdf}")

    script:
    def prefix = "${meta.phenotype}-${meta.population}-meta"
    """
    qc_meta_gwas.py \
        --input  ${sumstats} \
        --output ${prefix}
    """
}
