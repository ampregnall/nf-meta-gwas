process EXTRACT_LEAD_VARIANTS {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    publishDir(
        { "${launchDir}/results/lead-variants/${meta.phenotype}" },
        mode: 'copy'
    )
    memory { 16.GB + (12.GB * (task.attempt - 1)) }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.lead.variants.txt"), emit: lead_variants
        tuple val(meta), path("*.{pdf,png}"),         emit: plots

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    lead_variants.py \
      --input ${sumstats} \
      --output ${prefix} \
      --phenotype ${meta.phenotype} \
      --cohort meta \
      --population ${meta.population} \
      --gtf ${params.gtf}
    """
}
