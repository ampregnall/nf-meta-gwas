process ANNOTATE_CREDSET_VEP {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=4000]\"" }

    publishDir(
        { "${launchDir}/results/vep/${meta.phenotype}" },
        mode: 'copy'
    )

    input:
        tuple val(meta), path(credset), path(vep_annotation)

    output:
        tuple val(meta), path("*.vep_annotated_credset.txt"), emit: annotated_credset

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    annotate_credset_vep.py \\
        --credset    ${credset} \\
        --vep_output ${vep_annotation} \\
        --output     ${prefix}.vep_annotated_credset.txt
    """
}
