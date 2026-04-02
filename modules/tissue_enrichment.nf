process TISSUE_ENRICHMENT {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    publishDir(
        { "${launchDir}/results/enrichment/${meta.phenotype}" },
        mode: 'copy'
    )
    clusterOptions = { "-R \"rusage[mem=${24000 * task.attempt}]\"" }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.enrichment.txt")

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    def ldbaseline = params.ldbaseline
    """
    enrichment.py \
      --input ${sumstats} \
      --output ${prefix}.enrichment.txt \
      --ldbaseline ${ldbaseline} \
      --ldcts ${params.ldcts} \
      --weights ${params.ldcts_weights}
    """
}
