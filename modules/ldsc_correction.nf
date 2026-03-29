process LDSC_CORRECTION {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    publishDir(
        { "${launchDir}/data/sumstats-processed/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.sumstats.processed.txt.gz"
    )
    publishDir(
        { "${launchDir}/results/manhattan-qq/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.{pdf,png}"
    )
    memory { 16.GB + (12.GB * (task.attempt - 1)) }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    input:
        tuple val(meta), path(sumstats)
    output:
        tuple val(meta), path("*.sumstats.parquet"),          emit: sumstats_parquet
        tuple val(meta), path("*.sumstats.processed.txt.gz"), emit: sumstats_txt
        tuple val(meta), path("*.ldsc_h2.txt"),               emit: ldsc_h2
        tuple val(meta), path("*.{pdf,png}"),                 emit: plots
    script:
    def prefix = "${meta.phenotype}-${meta.cohort}-${meta.population}"
    def ldsc = params.ldsc_reference[meta.population]
    """
    ldsc.py \
      --input ${sumstats} \
      --output ${prefix} \
      --ldsc ${ldsc} \
      --phenotype ${meta.phenotype} \
      --cohort ${meta.cohort} \
      --population ${meta.population}
    """
}