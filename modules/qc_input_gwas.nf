process QC_INPUT_GWAS {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=${8000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/plots/qc/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.{png,pdf}"
    )

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.{png,pdf}"),       emit: qc_plots,   optional: true
        tuple val(meta), path("*.qc-summary.tsv"),  emit: qc_summary

    script:
    def prefix = "${meta.phenotype}-${meta.cohort}-${meta.population}"
    """
    qc_input_gwas.py \
        --input      ${sumstats} \
        --output     ${prefix} \
        --phenotype  ${meta.phenotype} \
        --cohort     ${meta.cohort} \
        --population ${meta.population}
    """
}
