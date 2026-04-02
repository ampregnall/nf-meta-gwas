process ABF_FINEMAPPING {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/meta-analysis:latest'
    publishDir(
        { "${launchDir}/results/finemapping/${meta.phenotype}" },
        mode: 'copy'
    )
    clusterOptions = { "-R \"rusage[mem=${32000 * task.attempt}]\"" }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    input:
        tuple val(meta), path(sumstats), path(lead_variants)

    output:
        tuple val(meta), path("*.credset.txt")

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    bayesian_finemapping.R \
      --input ${sumstats} \
      --lead ${lead_variants} \
      --output ${prefix} \
      --phenotype ${meta.phenotype} \
      --population ${meta.population}
    """
}
