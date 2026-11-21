// Step 1: metagenome assembly

process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/01_assembly", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: contigs
    path "${meta.id}_megahit.log"

    script:
    def r = reads instanceof List ? reads : [reads]
    def readArgs = r.size() == 2 ? "-1 ${r[0]} -2 ${r[1]}" : "-r ${r[0]}"
    """
    megahit ${readArgs} -t ${task.cpus} -o megahit_out ${params.megahit_args}
    cp megahit_out/final.contigs.fa ${meta.id}.contigs.fa
    cp megahit_out/log ${meta.id}_megahit.log
    """
}

process METASPADES {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/01_assembly", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: contigs
    path "${meta.id}_spades.log"

    script:
    def r = reads instanceof List ? reads : [reads]
    if (r.size() != 2) error "metaSPAdes needs paired reads (sample ${meta.id} has only one file). Use --assembler megahit instead."
    def memGb = task.memory ? task.memory.toGiga() : 64
    """
    metaspades.py -1 ${r[0]} -2 ${r[1]} -t ${task.cpus} -m ${memGb} -o spades_out ${params.spades_args}
    cp spades_out/contigs.fasta ${meta.id}.contigs.fa
    cp spades_out/spades.log ${meta.id}_spades.log
    """
}
