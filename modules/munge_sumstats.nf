process MUNGE_SUMSTATS {
    cpus 4
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=${48000 * task.attempt}]\"" }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/daf/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.{pdf,png}"
    )

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.sumstats.munged.txt.gz"), emit: sumstats_munged
        tuple val(meta), path("*.filter_stats.txt"),       emit: filter_stats

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
      --phenotype ${meta.phenotype} \
      --cohort ${meta.cohort} \
      --population ${meta.population}
    """
}
