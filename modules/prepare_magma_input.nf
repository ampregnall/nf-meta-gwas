process PREPARE_MAGMA_INPUT {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/magma:latest'
    clusterOptions = { "-R \"rusage[mem=8000]\"" }

    publishDir(
        { "${launchDir}/results/magma/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.magma-input.txt"
    )

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.magma-input.txt"), emit: magma_input

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    prepare_magma_input.py \
        --input  ${sumstats} \
        --output ${prefix}.magma-input.txt
    """
}
