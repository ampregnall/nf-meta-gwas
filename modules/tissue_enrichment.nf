process TISSUE_ENRICHMENT {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    publishDir(
        { "${launchDir}/results/enrichment/${meta.phenotype}" },
        mode: 'copy'
    )
    memory { 16.GB + (12.GB * (task.attempt - 1)) }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.enrichment.txt")

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    def ldsc = params.ldsc_reference[meta.population]
    """
    enrichment.py \
      --input ${sumstats} \
      --output ${prefix}.enrichment.txt \
      --ld ${ldsc} \
      --ldcts ${params.ldcts} \
      --weights ${params.ldcts_weights}
    """
}
