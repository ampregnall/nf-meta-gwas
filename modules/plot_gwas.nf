process PLOT_GWAS {
    cpus 2
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=${8000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/plots/manhattan/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.manhattan.*"
    )
    publishDir(
        { "${launchDir}/results/plots/qc/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.qq.*"
    )

    input:
        tuple val(meta), path(sumstats)

    output:
        tuple val(meta), path("*.manhattan.png"), emit: manhattan_png
        tuple val(meta), path("*.manhattan.pdf"), emit: manhattan_pdf
        tuple val(meta), path("*.qq.png"),        emit: qq_png
        tuple val(meta), path("*.qq.pdf"),        emit: qq_pdf

    script:
    def prefix = meta.cohort
        ? "${meta.phenotype}-${meta.cohort}-${meta.population}"
        : "${meta.phenotype}-${meta.population}-meta"
    """
    plot_gwas.py \
      --input  ${sumstats} \
      --output ${prefix} \
      --gtf    ${params.gtf}
    """
}
