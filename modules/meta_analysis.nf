process META_ANALYZE {
    cpus 1
    publishDir { "${launchDir}/data/meta-analysis/${meta.phenotype}" }, mode: 'copy'
    clusterOptions = { "-R \"rusage[mem=${48000 * task.attempt}]\"" }
    maxRetries 4
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    
    input:
        tuple val(meta), path(sumstats), val(chrom)

    output:
        tuple val(meta), path("*.txt.gz")

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    meta_analysis.R \
      --output ${prefix} \
      --type ${meta.type} \
      --chr ${chrom}
    """
}