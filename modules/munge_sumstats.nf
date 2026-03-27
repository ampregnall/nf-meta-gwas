process MUNGE_SUMSTATS {
    cpus 4
    container 'ghcr.io/ampregnall/nf-meta-gwas:0.2.0'
    publishDir { "${launchDir}/data/sumstats-processed/${meta.phenotype}" }, mode: 'copy'
    memory { 48.GB + (48.GB * (task.attempt - 1)) }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    
    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.sumstats.processed.txt.gz"), emit: sumstats_munged

    script:
    def prefix = "${meta.phenotype}-${meta.cohort}-${meta.population}"
    def popvcf = params.population_vcf[meta.population]
    def ldsc = params.ldsc_reference[meta.population]
    """
    munge_sumstats.py \
      --input ${sumstats} \
      --output ${prefix} \
      --type ${meta.type} \
      --mac ${params.mac} \
      --chain ${params.chain} \
      --fasta ${params.fasta} \
      --dbsnp ${params.dbsnp} \
      --popvcf ${popvcf} \
      --ldsc ${ldsc} \
      --threads ${task.cpus} \
    """
}