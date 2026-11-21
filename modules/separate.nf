// Step 2: split eukaryotic and prokaryotic contigs with DeepMicroClass

process DEEPMICROCLASS {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/02_euk_prok_split/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.filtered.fa"),        emit: filtered
    tuple val(meta), path("${meta.id}.eukaryote.fa"),       emit: euk
    tuple val(meta), path("${meta.id}.prokaryote.fa"),      emit: prok
    tuple val(meta), path("${meta.id}.other.fa"),           emit: other
    tuple val(meta), path("${meta.id}.deepmicroclass.tsv"), emit: tsv
    path "${meta.id}.class_summary.tsv"

    script:
    """
    # drop very short contigs first (they cannot be binned anyway and slow everything down)
    filter_contigs.sh ${contigs} ${meta.id}.filtered.fa ${params.min_contig_length}

    DeepMicroClass predict -i ${meta.id}.filtered.fa -o dmc_out -d ${params.deepmicroclass_device} ${params.deepmicroclass_args}

    # DeepMicroClass writes one table of scores; take it whatever its exact name is
    tsv=\$(find dmc_out -name '*.tsv' | head -n 1)
    if [ -z "\$tsv" ]; then echo "DeepMicroClass did not write a .tsv table" >&2; exit 1; fi
    cp "\$tsv" ${meta.id}.deepmicroclass.tsv

    # each contig goes to the class with the highest score
    split_by_class.py --fasta ${meta.id}.filtered.fa --tsv ${meta.id}.deepmicroclass.tsv --prefix ${meta.id}
    """
}
