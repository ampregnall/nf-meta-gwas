process PREPARE_VEP_INPUT {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest'
    clusterOptions = { "-R \"rusage[mem=4000]\"" }

    input:
        tuple val(meta), path(credset)

    output:
        tuple val(meta), path("*.vep_input.txt"), emit: vep_input

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    prepare_vep_input.py \
        --credset ${credset} \
        --output  ${prefix}.vep_input.txt
    """
}
