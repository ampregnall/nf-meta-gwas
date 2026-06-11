process PREPARE_FLAMES_INPUTS {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/meta-analysis:latest'
    clusterOptions = { "-R \"rusage[mem=8000]\"" }

    publishDir(
        { "${launchDir}/results/flames/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.{genomic-locus.txt,index.txt}"
    )

    input:
        tuple val(meta), path(credset)

    output:
        tuple val(meta),
              path("credsets/"),
              path("*.genomic-locus.txt"),
              path("*.index.txt"),
              emit: flames_inputs

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    prepare_flames_inputs.R \
        --credset ${credset} \
        --output  ${prefix}
    """
}
