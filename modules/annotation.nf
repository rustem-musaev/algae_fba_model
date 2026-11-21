// Step 7: functional annotation

process EGGNOG_MAPPER {
    tag "${meta.bin}"
    label 'process_medium'
    publishDir "${params.outdir}/07_annotation/eggnog", mode: params.publish_mode

    input:
    tuple val(meta), path(proteins)
    path db

    output:
    tuple val(meta), path("${meta.bin}.emapper.annotations"), emit: annotations
    path "${meta.bin}.emapper.*"

    script:
    """
    emapper.py -i ${proteins} --itype proteins --tax_scope ${params.eggnog_tax_scope} \\
        --data_dir ${db} --cpu ${task.cpus} -o ${meta.bin} --output_dir . ${params.eggnog_args}
    """
}

process INTERPROSCAN {
    tag "${meta.bin}"
    label 'process_medium'
    publishDir "${params.outdir}/07_annotation/interproscan", mode: params.publish_mode

    input:
    tuple val(meta), path(proteins)

    output:
    tuple val(meta), path("${meta.bin}.interproscan.tsv"), emit: tsv

    script:
    """
    ${params.interproscan_cmd} -i ${proteins} -f tsv --goterms --pathways \\
        --cpu ${task.cpus} -o ${meta.bin}.interproscan.tsv ${params.interproscan_args}
    """
}
