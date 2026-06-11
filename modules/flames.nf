process FLAMES {
    cpus 1
    container 'ghcr.io/ampregnall/nf-meta-gwas/flames:latest'
    clusterOptions = { "-R \"rusage[mem=${8000 * task.attempt}]\"" }
    maxRetries 2
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    publishDir(
        { "${launchDir}/results/flames/${meta.phenotype}" },
        mode: 'copy'
    )

    input:
        tuple val(meta),
              path(credsets_dir),
              path(genomic_locus),
              path(index_file),
              path(pops_preds),
              path(genes_out),
              path(gsa_out)

    output:
        tuple val(meta), path("flames-annotate/"), emit: annotated
        tuple val(meta), path("flames-scores/"),   emit: scores

    script:
    """
    mkdir -p flames-annotate flames-scores

    FLAMES.py annotate \
        -o  flames-annotate \
        -a  ${params.flames_annotation_data} \
        -p  ${pops_preds} \
        -m  ${genes_out} \
        -mt ${gsa_out} \
        --GenomicRiskLoci ${genomic_locus} \
        -id ${index_file} \
        --build      GRCh38 \
        --credset_95 False

    FLAMES.py FLAMES \
        -id ${index_file} \
        -o  flames-scores
    """
}
