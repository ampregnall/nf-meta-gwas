process POPS {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/pops:latest'
    clusterOptions = { "-R \"rusage[mem=${8000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/pops/${meta.phenotype}" },
        mode: 'copy'
    )

    input:
        tuple val(meta), path(genes_raw)

    output:
        tuple val(meta), path("*.preds"), emit: pops_preds

    script:
    // PoPS expects the MAGMA prefix (basename without .genes.raw extension)
    def magma_prefix = "${meta.phenotype}-${meta.population}"
    def out_prefix   = "${meta.phenotype}-${meta.population}-pops"
    """
    pops.py \
        --gene_annot_path      ${params.pops_gene_annot} \
        --feature_mat_prefix   ${params.pops_feature_prefix} \
        --num_feature_chunks   ${params.pops_num_feature_chunks} \
        --magma_prefix         ${magma_prefix} \
        --control_features     ${params.pops_control_features} \
        --out_prefix           ${out_prefix}
    """
}
