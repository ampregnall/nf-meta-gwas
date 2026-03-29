process COLLECT_LDSC_RESULTS {
    cpus 1
    publishDir { "${launchDir}/results/ldsc/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta), path(ldsc_files)

    output:
        tuple val(meta), path("*.ldsc_h2.txt")

    script:
    def prefix = "${meta.phenotype}"
    """
    awk 'NR==1 || FNR>1' ${ldsc_files} > ${prefix}.ldsc_h2.txt
    """
}
