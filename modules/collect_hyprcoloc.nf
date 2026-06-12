process COLLECT_HYPRCOLOC {
    cpus 1
    publishDir { "${launchDir}/results/colocalization/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta), path(result_files), path(failure_files)

    output:
        tuple val(meta), path("*.hyprcoloc.tsv"),          emit: results
        tuple val(meta), path("*.hyprcoloc.failures.tsv"), emit: failures

    script:
    def prefix = "${meta.phenotype}"
    """
    awk 'NR==1 || FNR>1' ${result_files}  > ${prefix}.hyprcoloc.tsv
    awk 'NR==1 || FNR>1' ${failure_files} > ${prefix}.hyprcoloc.failures.tsv
    """
}
