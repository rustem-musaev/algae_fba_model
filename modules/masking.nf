// Step 5: shorten headers and soft-mask repeats

process SHORTEN_HEADERS {
    tag "${meta.bin}"
    label 'process_single'
    publishDir "${params.outdir}/05_masking/${meta.bin}", mode: params.publish_mode, pattern: '*.tsv'

    input:
    tuple val(meta), path(bin)

    output:
    tuple val(meta), path("${meta.bin}.renamed.fa"), emit: fasta
    path "${meta.bin}.header_map.tsv"

    script:
    // RepeatMasker refuses headers over 50 characters. Instead of cutting names
    // (which can make two contigs end up with the same name) every contig gets a
    // short new name (c000001, c000002, ...) and the old->new names are saved.
    """
    rename_headers.sh ${bin} ${meta.bin}.renamed.fa ${meta.bin}.header_map.tsv
    """
}

process REPEATMASKER {
    tag "${meta.bin}"
    label 'process_medium'
    publishDir "${params.outdir}/05_masking/${meta.bin}", mode: params.publish_mode

    input:
    tuple val(meta), path(fasta)
    path lib

    output:
    tuple val(meta), path("${meta.bin}.masked.fa"), emit: masked
    path "*.{out,tbl}", optional: true

    script:
    // each RepeatMasker "-pa" job uses ~4 threads
    def pa = Math.max(1, (task.cpus / 4) as int)
    def what = lib ? "-lib ${lib}" : "-species ${params.repeatmasker_species}"
    """
    RepeatMasker ${what} -xsmall -pa ${pa} -dir . ${params.repeatmasker_args} ${fasta}

    # if no repeats are found RepeatMasker writes no .masked file; then the bin is used as it is
    if [ -f ${fasta}.masked ]; then
        cp ${fasta}.masked ${meta.bin}.masked.fa
    else
        echo "No repeats found in ${meta.bin}; using it unmasked." >&2
        cp ${fasta} ${meta.bin}.masked.fa
    fi
    """
}
