process COLLECT_META_RESULTS {
    cpus 1
    publishDir { "${launchDir}/data/meta-analysis/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta), path(meta_files)

    output:
        tuple val(meta), path("*.txt.gz"), emit: sumstats

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    sorted_files=\$(ls ${meta_files} | sort -V)
    awk 'NR==1 || FNR>1' \${sorted_files} > ${prefix}.txt
    gzip ${prefix}.txt
    """
}