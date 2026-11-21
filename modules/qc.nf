// Step 4: quality check of eukaryotic bins

process QUAST {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/04_quality/quast", mode: params.publish_mode

    input:
    tuple val(meta), path(bins, stageAs: 'bins/*')

    output:
    path "quast_${meta.id}"

    script:
    """
    quast.py -o quast_${meta.id} -t ${task.cpus} bins/*
    """
}

process EUKCC {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/04_quality/eukcc", mode: params.publish_mode

    input:
    tuple val(meta), path(bins, stageAs: 'bins/*')
    path db

    output:
    path "eukcc_${meta.id}"

    script:
    // EukCC only looks at files ending in .fa (MetaBAT2 bins already do)
    """
    eukcc folder --out eukcc_${meta.id} --threads ${task.cpus} --db ${db} bins
    """
}

process BUSCO {
    tag "${meta.bin}"
    label 'process_medium'
    publishDir "${params.outdir}/04_quality/busco", mode: params.publish_mode

    input:
    tuple val(meta), path(bin)

    output:
    tuple val(meta), path("short_summary.*.txt"), emit: summary
    path "${meta.bin}"

    script:
    def lineage = params.busco_lineage == 'auto'     ? '--auto-lineage'
                : params.busco_lineage == 'auto_euk' ? '--auto-lineage-euk'
                : "-l ${params.busco_lineage}"
    def offline = params.busco_download_path ? "--offline --download_path ${params.busco_download_path}" : ''
    """
    busco -i ${bin} -m genome ${lineage} -o ${meta.bin} -c ${task.cpus} ${offline} -f ${params.busco_args}

    # keep one summary per bin for the plot: the most specific lineage if there is one
    f=\$(ls ${meta.bin}/short_summary.specific.*.txt 2>/dev/null | head -n 1 || true)
    [ -z "\$f" ] && f=\$(ls ${meta.bin}/short_summary.generic.*.txt 2>/dev/null | head -n 1 || true)
    [ -z "\$f" ] && f=\$(ls ${meta.bin}/short_summary*.txt | head -n 1)
    cp "\$f" .
    """
}

process BUSCO_PLOT {
    label 'process_low'
    publishDir "${params.outdir}/04_quality/busco", mode: params.publish_mode

    input:
    path(summaries, stageAs: 'BUSCO_summaries/*')

    output:
    path "BUSCO_summaries"

    script:
    """
    generate_plot.py -wd BUSCO_summaries
    """
}
