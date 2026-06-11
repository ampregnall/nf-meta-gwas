process MAGMA_TISSUE {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/magma:latest'
    clusterOptions = { "-R \"rusage[mem=${8000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/magma/${meta.phenotype}" },
        mode: 'copy'
    )

    input:
        tuple val(meta), path(genes_raw)

    output:
        tuple val(meta), path("*.gsa.out"), emit: gsa_out

    script:
    def prefix = "${meta.phenotype}-${meta.population}-gtex"
    """
    magma \
        --gene-results ${genes_raw} \
        --gene-covar   ${params.magma_gtex_covar} \
        --out          ${prefix}
    """
}
