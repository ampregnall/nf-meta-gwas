process HERITABILITY {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    publishDir(
        { "${launchDir}/results/heritability/${meta.phenotype}" },
        mode: 'copy'
    )
    memory { 16.GB + (12.GB * (task.attempt - 1)) }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.ldsc_h2.txt")

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    def ldsc = params.ldsc_reference[meta.population]
    """
    heritability.py \
      --input ${sumstats} \
      --output ${prefix} \
      --ldsc ${ldsc} \
      --phenotype ${meta.phenotype} \
      --cohort meta \
      --population ${meta.population}
    """
}
