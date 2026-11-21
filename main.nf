#!/usr/bin/env nextflow
/*
 * eukprok-nf : eukaryotic + prokaryotic metagenome pipeline
 *
 * reads -> assembly -> euk/prok split (DeepMicroClass) -> binning (MetaBAT2)
 *       -> bin quality (QUAST, BUSCO, EukCC)
 *       -> eukaryotic bins: header shortening -> RepeatMasker -> GeneMark-ES
 *          -> AUGUSTUS (GeneMark hints) -> proteins
 *          -> eggNOG-mapper / InterProScan -> EUKulele taxonomy -> per-taxon protein sets
 *       -> optional Kraken2 on the assembled contigs
 *
 * AuReMe / MetaDraft (metabolic model drafting) stay manual: see README.
 */
nextflow.enable.dsl = 2

include { MEGAHIT ; METASPADES             } from './modules/assembly'
include { DEEPMICROCLASS                   } from './modules/separate'
include { MAP_AND_DEPTH ; METABAT2         } from './modules/binning'
include { QUAST ; BUSCO ; BUSCO_PLOT ; EUKCC } from './modules/qc'
include { SHORTEN_HEADERS ; REPEATMASKER   } from './modules/masking'
include { GENEMARK_ES ; AUGUSTUS           } from './modules/genepred'
include { EGGNOG_MAPPER ; INTERPROSCAN     } from './modules/annotation'
include { EUKULELE ; EXTRACT_TAXA_PROTEINS ; KRAKEN2 } from './modules/taxonomy'

// ---------------------------------------------------------------------------
// Check settings before anything starts
// ---------------------------------------------------------------------------
def qcTools    = params.qc_tools.toString().tokenize(',')*.trim()
def annotTools = params.annotation_tools.toString().tokenize(',')*.trim().findAll { it && it != 'none' }

if (!params.input) error "Please give a samplesheet with --input samples.csv (see README)"
if (!(params.assembler in ['megahit', 'metaspades'])) error "--assembler must be 'megahit' or 'metaspades'"
if ('eukcc' in qcTools && !params.eukcc_db) error "EukCC is switched on but --eukcc_db is not set. Give the EukCC database folder, or run only BUSCO with --qc_tools busco"
if ('eggnog' in annotTools && !params.eggnog_db) error "eggNOG-mapper is switched on but --eggnog_db is not set. Give the eggNOG data folder, or use --annotation_tools none (or interproscan)"
if (params.taxa_list && !params.run_eukulele) error "--taxa_list needs EUKulele (--run_eukulele true)"

def withBin(meta, f) {
    // give every bin its own short, safe name, e.g. S1_euk_bin_3
    def name = f.baseName.replace('.', '_')
    return [ meta + [bin: name], f ]
}

workflow {

    // -----------------------------------------------------------------------
    // 0. Read the samplesheet: sample,fastq_1,fastq_2[,contigs]
    // -----------------------------------------------------------------------
    samples = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq_1) error "Samplesheet rows need at least 'sample' and 'fastq_1': ${row}"
            def meta  = [id: row.sample]
            def reads = row.fastq_2 ? [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
                                    : [file(row.fastq_1, checkIfExists: true)]
            def contigs = row.contigs ? file(row.contigs, checkIfExists: true) : []
            [meta, reads, contigs]
        }

    reads = samples.map { meta, r, c -> [meta, r] }

    by_input = samples.branch { meta, r, c ->
        given:    c
        assemble: true
    }

    // -----------------------------------------------------------------------
    // 1. Assembly (skipped for samples that already have a 'contigs' column)
    // -----------------------------------------------------------------------
    to_assemble = by_input.assemble.map { meta, r, c -> [meta, r] }
    if (params.assembler == 'megahit') {
        MEGAHIT(to_assemble)
        assembled = MEGAHIT.out.contigs
    } else {
        METASPADES(to_assemble)
        assembled = METASPADES.out.contigs
    }
    contigs = assembled.mix(by_input.given.map { meta, r, c -> [meta, c] })

    if (params.kraken2_db) {
        KRAKEN2(contigs, file(params.kraken2_db, checkIfExists: true))
    }

    // -----------------------------------------------------------------------
    // 2. Split eukaryotic / prokaryotic contigs
    // -----------------------------------------------------------------------
    DEEPMICROCLASS(contigs)

    // -----------------------------------------------------------------------
    // 3. Binning: map reads once to all kept contigs, then bin euk and prok separately
    // -----------------------------------------------------------------------
    MAP_AND_DEPTH(DEEPMICROCLASS.out.filtered.join(reads))

    to_bin = DEEPMICROCLASS.out.euk.map  { meta, f -> [meta, f, 'euk'] }
        .mix(DEEPMICROCLASS.out.prok.map { meta, f -> [meta, f, 'prok'] })
        .combine(MAP_AND_DEPTH.out.depth, by: 0)          // -> meta, fasta, label, depth

    METABAT2(to_bin)

    all_bins = METABAT2.out.bins.map { meta, label, b -> [meta, label, b instanceof List ? b : [b]] }
    euk_bins_per_sample = all_bins.filter { it[1] == 'euk' }.map { meta, label, b -> [meta, b] }
    euk_bins = euk_bins_per_sample.flatMap { meta, b -> b.collect { f -> withBin(meta, f) } }

    // -----------------------------------------------------------------------
    // 4. Quality check of eukaryotic bins
    // -----------------------------------------------------------------------
    QUAST(euk_bins_per_sample)

    if ('eukcc' in qcTools) {
        EUKCC(euk_bins_per_sample, file(params.eukcc_db, checkIfExists: true))
    }
    if ('busco' in qcTools) {
        BUSCO(euk_bins)
        BUSCO_PLOT(BUSCO.out.summary.map { meta, s -> s }.collect())
    }

    if (!params.skip_gene_prediction) {

        // -----------------------------------------------------------------------
        // 5. Masking: shorten headers (RepeatMasker limit) then soft-mask repeats
        // -----------------------------------------------------------------------
        SHORTEN_HEADERS(euk_bins)
        rm_lib = params.repeatmasker_lib ? file(params.repeatmasker_lib, checkIfExists: true) : []
        REPEATMASKER(SHORTEN_HEADERS.out.fasta, rm_lib)

        // -----------------------------------------------------------------------
        // 6. Gene prediction: GeneMark-ES -> hints -> AUGUSTUS -> proteins
        // -----------------------------------------------------------------------
        GENEMARK_ES(REPEATMASKER.out.masked)
        AUGUSTUS(REPEATMASKER.out.masked.join(GENEMARK_ES.out.gtf))

        proteins = AUGUSTUS.out.proteins.filter { meta, f -> f.size() > 0 }

        // -----------------------------------------------------------------------
        // 7. Functional annotation
        // -----------------------------------------------------------------------
        if ('eggnog' in annotTools) {
            EGGNOG_MAPPER(proteins, file(params.eggnog_db, checkIfExists: true))
        }
        if ('interproscan' in annotTools) {
            INTERPROSCAN(proteins)
        }

        // -----------------------------------------------------------------------
        // 8. Taxonomy of predicted proteins (EUKulele), then optional per-taxon protein sets
        // -----------------------------------------------------------------------
        if (params.run_eukulele) {
            proteins_per_sample = proteins
                .map { meta, f -> [[id: meta.id], f] }
                .groupTuple()
            eukulele_ref = params.eukulele_reference_dir ? file(params.eukulele_reference_dir, checkIfExists: true) : []
            EUKULELE(proteins_per_sample, eukulele_ref)

            if (params.taxa_list) {
                EXTRACT_TAXA_PROTEINS(
                    EUKULELE.out.taxonomy.join(proteins_per_sample),
                    file(params.taxa_list, checkIfExists: true)
                )
            }
        }
    }
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nDone. Results are in: ${params.outdir}\n"
        : "\nThe pipeline stopped with an error. Look at the message above and at .nextflow.log\n")
}
