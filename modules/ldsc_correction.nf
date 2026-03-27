process LDSC_CORRECTION {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    publishDir { "${launchDir}/data/sumstats-processed/${meta.phenotype}" }, mode: 'copy'
    memory { 16.GB + (12.GB * (task.attempt - 1)) }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    
    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.sumstats.processed.txt.gz"), emit: sumstats_munged

    script:
    def prefix = "${meta.phenotype}-${meta.cohort}-${meta.population}"
    def ldsc = params.ldsc_reference[meta.population]
    """
    ldsc.py \
      --input ${sumstats} \
      --output ${prefix} \
      --ldsc ${ldsc}
    """
}