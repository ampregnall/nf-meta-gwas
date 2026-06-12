process HYPRCOLOC {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/meta-analysis:latest'
    clusterOptions = { "-R \"rusage[mem=${16000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    input:
        tuple val(meta), path(sumstats), path(lead_variants), path(eqtl_parquet)

    output:
        tuple val(meta), path("*.hyprcoloc.tsv"),          emit: results
        tuple val(meta), path("*.hyprcoloc.failures.tsv"), emit: failures

    script:
    def study_name = eqtl_parquet.baseName
    def prefix     = "${meta.phenotype}-${meta.population}-${study_name}"
    """
    hyprcoloc_eqtl.R \\
        --sumstats     ${sumstats} \\
        --lead         ${lead_variants} \\
        --eqtl_parquet ${eqtl_parquet} \\
        --output       ${prefix} \\
        --phenotype    ${meta.phenotype} \\
        --population   ${meta.population} \\
        --eqtl_study   ${study_name}
    """
}
