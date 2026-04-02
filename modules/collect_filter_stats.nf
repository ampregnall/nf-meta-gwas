process COLLECT_FILTER_STATS {
    cpus 1
    publishDir { "${launchDir}/results/filter_stats" }, mode: 'copy'

    input:
        tuple val(meta), path(filter_stats_files)

    output:
        tuple val(meta), path("*.filter_stats.txt")

    script:
    def prefix = "${meta.phenotype}"
    """
    awk 'NR==1 || FNR>1' ${filter_stats_files} > ${prefix}.filter_stats.txt
    """
}
