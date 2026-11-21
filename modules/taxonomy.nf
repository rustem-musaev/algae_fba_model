// Step 8: taxonomy (EUKulele on predicted proteins, optional Kraken2 on contigs)

process EUKULELE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/08_taxonomy/eukulele", mode: params.publish_mode

    input:
    // only protein files go in the sample folder (EUKulele has a bug when the
    // genome FASTA sits in the same folder)
    tuple val(meta), path(proteins, stageAs: 'proteins/*')
    path reference_dir

    output:
    tuple val(meta), path("${meta.id}.estimated-taxonomy.tsv"), emit: taxonomy
    path "eukulele_${meta.id}"

    script:
    def ref = reference_dir ? "--reference_dir ${reference_dir}" : ''
    """
    EUKulele --sample_dir proteins -m mags --p_ext .faa -o eukulele_${meta.id} \\
        --database ${params.eukulele_database} ${ref} --CPUs ${task.cpus} ${params.eukulele_args}

    # pool the per-bin tables into one (keep a single header line)
    awk 'FNR == 1 && NR != 1 { next } { print }' eukulele_${meta.id}/taxonomy_estimation/*estimated-taxonomy.out \\
        > ${meta.id}.estimated-taxonomy.tsv
    """
}

process EXTRACT_TAXA_PROTEINS {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir}/08_taxonomy/proteins_by_taxon/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(taxonomy), path(proteins, stageAs: 'proteins/*')
    path taxa_list

    output:
    path "*_proteins.faa", optional: true

    script:
    """
    cat proteins/*.faa > all_proteins.faa
    extract_taxa_proteins.sh ${taxonomy} all_proteins.faa ${taxa_list}
    """
}

process KRAKEN2 {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/08_taxonomy/kraken2", mode: params.publish_mode

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    path "${meta.id}.kraken2.*"

    script:
    """
    kraken2 --db ${db} --threads ${task.cpus} --output ${meta.id}.kraken2.out \\
        --report ${meta.id}.kraken2.report ${params.kraken2_args} ${contigs}
    """
}
