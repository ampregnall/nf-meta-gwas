process COLLECT_LEAD_VARIANTS {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    publishDir { "${launchDir}/results/lead-variants/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta), path(lead_files)

    output:
        tuple val(meta), path("*.locus_summary.txt")

    script:
    def prefix = "${meta.phenotype}"
    """
    collect_lead_variants.py \
      --inputs ${lead_files} \
      --phenotype ${meta.phenotype} \
      --output ${prefix}
    """
}
