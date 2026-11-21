// Step 6: gene prediction (GeneMark-ES, then AUGUSTUS with GeneMark hints)

process GENEMARK_ES {
    tag "${meta.bin}"
    label 'process_medium'
    label 'genemark'
    publishDir "${params.outdir}/06_gene_prediction/${meta.bin}", mode: params.publish_mode

    input:
    tuple val(meta), path(masked)

    output:
    tuple val(meta), path("${meta.bin}.genemark.gtf"), emit: gtf

    script:
    """
    gmes_petap.pl --seq ${masked} --ES --soft_mask auto --cores ${task.cpus} \\
        --min_contig ${params.genemark_min_contig} ${params.genemark_args}
    cp genemark.gtf ${meta.bin}.genemark.gtf
    """
}

process AUGUSTUS {
    tag "${meta.bin}"
    label 'process_single'
    publishDir "${params.outdir}/06_gene_prediction/${meta.bin}", mode: params.publish_mode

    input:
    tuple val(meta), path(masked), path(genemark_gtf)

    output:
    tuple val(meta), path("${meta.bin}.faa"),          emit: proteins
    tuple val(meta), path("${meta.bin}.augustus.gff"), emit: gff
    path "${meta.bin}.hints.gtf"
    path "${meta.bin}.codingseq", optional: true

    script:
    def extrinsic = params.augustus_extrinsic_cfg ? "--extrinsicCfgFile=${params.augustus_extrinsic_cfg}" : ''
    """
    # make the GeneMark GTF readable as AUGUSTUS hints (gene/mRNA -> genicpart, UTRs -> UTR, source=M)
    fix_genemark_gtf.sh ${genemark_gtf} ${meta.bin}.hints.gtf

    augustus --species=${params.augustus_species} --hintsfile=${meta.bin}.hints.gtf \\
        --softmasking=1 ${extrinsic} ${params.augustus_args} ${masked} > ${meta.bin}.augustus.gff

    # protein sequences from the AUGUSTUS genes
    getAnnoFasta.pl --seqfile=${masked} ${meta.bin}.augustus.gff
    [ -f ${meta.bin}.augustus.codingseq ] && mv ${meta.bin}.augustus.codingseq ${meta.bin}.codingseq
    [ -f ${meta.bin}.augustus.aa ] || touch ${meta.bin}.augustus.aa

    # AUGUSTUS calls genes g1, g2, ... in every bin, so put the bin name in front
    # to keep protein names unique once bins are pooled (e.g. for EUKulele)
    tag_proteins.sh ${meta.bin}.augustus.aa ${meta.bin}.faa ${meta.bin}
    """
}
