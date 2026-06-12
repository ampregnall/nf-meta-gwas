process RENDER_REPORT {
    cpus 2
    container 'ghcr.io/ampregnall/nf-meta-gwas/reports:latest'
    clusterOptions = { "-R \"rusage[mem=8000]\"" }
    publishDir { "${launchDir}/results/reports/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta),
              path(locus_annotated),
              path(credsets),
              path(filter_stats),
              path(ldsc_h2),
              val(has_vep),
              val(has_hyprcoloc),
              val(has_flames)

    output:
        tuple val(meta), path("*.qc.html"),      emit: qc_report
        tuple val(meta), path("*.results.html"), emit: results_report

    script:
    def prefix      = "${meta.phenotype}"
    def results_dir = "${launchDir}/results"
    def vepFlag     = has_vep       ? "--has_vep"       : ""
    def hyprFlag    = has_hyprcoloc ? "--has_hyprcoloc" : ""
    def flamesFlag  = has_flames    ? "--has_flames"    : ""
    """
    render_reports.py \\
        --phenotype       ${prefix} \\
        --results_dir     ${results_dir} \\
        --locus_annotated ${locus_annotated} \\
        --credsets        ${credsets} \\
        --filter_stats    ${filter_stats} \\
        --ldsc_h2         ${ldsc_h2} \\
        ${vepFlag} \\
        ${hyprFlag} \\
        ${flamesFlag} \\
        --templates_dir   ${projectDir}/assets/reports \\
        --output_prefix   ${prefix}
    """
}
