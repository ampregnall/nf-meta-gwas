process MAGMA_GENE {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/magma:latest'
    clusterOptions = { "-R \"rusage[mem=${16000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/magma/${meta.phenotype}" },
        mode: 'copy'
    )

    input:
        tuple val(meta), path(magma_input)

    output:
        tuple val(meta), path("*.genes.raw"), emit: genes_raw
        tuple val(meta), path("*.genes.out"), emit: genes_out

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    magma \
        --bfile        ${params.magma_bfile} \
        --gene-annot   ${params.magma_gene_annot} \
        --pval         ${magma_input} use=SNP,P ncol=N \
        --gene-model   snp-wise=mean \
        --out          ${prefix}
    """
}
