// Step 3: read mapping + MetaBAT2 binning

process MAP_AND_DEPTH {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/03_binning/${meta.id}", mode: params.publish_mode, pattern: '*.depth.txt'

    input:
    tuple val(meta), path(contigs), path(reads)

    output:
    tuple val(meta), path("${meta.id}.depth.txt"), emit: depth

    script:
    def r = reads instanceof List ? reads : [reads]
    def readArgs = r.size() == 2 ? "-1 ${r[0]} -2 ${r[1]}" : "-U ${r[0]}"
    """
    bowtie2-build --threads ${task.cpus} ${contigs} contigs_index
    bowtie2 -p ${task.cpus} -x contigs_index ${readArgs} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam
    jgi_summarize_bam_contig_depths --outputDepth ${meta.id}.depth.txt ${meta.id}.bam
    """
}

process METABAT2 {
    tag "${meta.id}:${label}"
    label 'process_medium'
    publishDir "${params.outdir}/03_binning/${meta.id}/${label}", mode: params.publish_mode

    input:
    tuple val(meta), path(fasta), val(label), path(depth)

    output:
    tuple val(meta), val(label), path("bins/*.fa"), emit: bins, optional: true

    script:
    """
    mkdir -p bins
    if [ ! -s ${fasta} ]; then
        echo "No ${label} contigs for ${meta.id}; nothing to bin." >&2
        exit 0
    fi

    # MetaBAT2 needs the depth table to list exactly the same contigs, in the same order
    subset_depth.sh ${fasta} ${depth} ${label}.depth.txt

    metabat2 -i ${fasta} -a ${label}.depth.txt -o bins/${meta.id}_${label}_bin \\
        -m ${params.metabat_min_contig} -t ${task.cpus} ${params.metabat_args}
    """
}
