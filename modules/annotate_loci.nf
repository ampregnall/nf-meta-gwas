process ANNOTATE_LOCI {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = "-R \"rusage[mem=4000]\""
    publishDir { "${launchDir}/results/lead-variants/${meta.phenotype}" }, mode: 'copy'

    input:
        tuple val(meta), path(locus_summary), path(credsets), val(flames_scores_dirs)

    output:
        tuple val(meta), path("*.locus_annotated.txt"), emit: locus_annotated

    script:
    def prefix = "${meta.phenotype}"
    def flamesArg = (flames_scores_dirs instanceof List && flames_scores_dirs.size() > 0)
        ? "--flames_scores ${flames_scores_dirs.join(' ')}"
        : ""
    """
    annotate_loci.py \
        --locus_summary ${locus_summary} \
        --credsets ${credsets} \
        ${flamesArg} \
        --output ${prefix}.locus_annotated.txt
    """
}
