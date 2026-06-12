process VEP {
    cpus 4
    container 'docker://ensemblorg/ensembl-vep'
    clusterOptions = { "-R \"rusage[mem=${8000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/vep/${meta.phenotype}" },
        mode: 'copy',
        pattern: "*.vep_annotation.tab"
    )

    input:
        tuple val(meta), path(vep_input)

    output:
        tuple val(meta), path("*.vep_annotation.tab"), emit: vep_annotation

    script:
    def prefix = "${meta.phenotype}-${meta.population}"
    """
    vep \\
        --input_file  ${vep_input} \\
        --output_file ${prefix}.vep_annotation.tab \\
        --cache \\
        --offline \\
        --pick \\
        --tab \\
        --fields "Uploaded_variation,SYMBOL,Consequence,IMPACT,VARIANT_CLASS" \\
        --dir         ${params.vep_cache} \\
        --species     homo_sapiens \\
        --assembly    GRCh38 \\
        --force_overwrite \\
        --fork        ${task.cpus}
    """
}
