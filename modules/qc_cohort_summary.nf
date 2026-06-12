process QC_COHORT_SUMMARY {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=8000]\"" }

    publishDir(
        { "${launchDir}/results/plots/qc/${meta.phenotype}" },
        mode: 'copy'
    )

    input:
        tuple val(meta), path(filter_stats), path(ldsc_h2), path(qc_summaries)

    output:
        tuple val(meta), path("*.{png,pdf}")

    script:
    def prefix = "${meta.phenotype}"
    """
    qc_cohort_summary.py \
        --filter_stats ${filter_stats} \
        --ldsc_h2      ${ldsc_h2} \
        --qc_summaries ${qc_summaries} \
        --output       ${prefix}
    """
}
