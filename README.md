# eukprok-nf

This is the "Eukaryotic + Prokaryotic Metagenomes" protocol turned into a Nextflow pipeline.
You give it sequencing reads. It runs every step from assembly to annotated, taxonomically labelled eukaryotic bins. Each step's results go into its own numbered folder. I made the Nextflow pipeline after the fact to make it reproducible, the actual analyses were performed by running the scripts separately.

## What it does, step by step

| # | Step | Tool | Results folder |
|---|------|------|----------------|
| 1 | Put reads together into contigs | MEGAHIT (default) or metaSPAdes | `01_assembly` |
| 2 | Sort contigs into eukaryotic, prokaryotic and "other" | DeepMicroClass | `02_euk_prok_split` |
| 3 | Group contigs into bins (genomes), eukaryotes and prokaryotes separately | Bowtie2 + MetaBAT2 | `03_binning` |
| 4 | Check how good the eukaryotic bins are | QUAST, BUSCO (+ plot), EukCC | `04_quality` |
| 5 | Mask repeats (soft-masking) | RepeatMasker | `05_masking` |
| 6 | Find genes | GeneMark-ES, then AUGUSTUS using the GeneMark genes as hints | `06_gene_prediction` |
| 7 | Guess what the proteins do | eggNOG-mapper and/or InterProScan | `07_annotation` |
| 8 | Guess which organism each protein comes from | EUKulele (+ optional Kraken2 on contigs) | `08_taxonomy` |

Steps 4 to 8 run on the **eukaryotic** bins. The prokaryotic bins are made and saved in `03_binning/<sample>/prok/`, but nothing more is done with them, just as in the original protocol.

**Not included:** AuReMe and MetaDraft (building the metabolic model). AuReMe needs template model files that you prepare by hand, and some of its steps are interactive. MetaDraft is a point-and-click program. Use the proteins from step 6 or step 8 as input for those, as before.

## Changes from the original protocol (and why)

- **Contig names before RepeatMasker.** The old AWK script kept the last 20 characters of each name. That can give two contigs the same name. Now every contig gets a short new name (`c000001`, `c000002`, …), and a table (`*.header_map.tsv`) records which new name belongs to which original contig.
- **Protein names.** AUGUSTUS names genes `g1`, `g2`, … in *every* bin, so the same names repeat across bins. The pipeline adds the bin name in front (for example `S1_euk_bin_3__g1.t1`). This keeps names unique when proteins from many bins are pooled for EUKulele.
- **Bins with no repeats.** If RepeatMasker finds nothing to mask, it writes no output file. The pipeline then uses the bin as it is instead of stopping.
- **Bins that are too small for GeneMark-ES.** GeneMark needs enough DNA to train on. When it fails on a bin, that bin is skipped and the rest continue. To stop on failures instead, use `--genemark_skip_failures false`.
- **Per-taxon protein files** (your `all_genera.txt` loop) are built in, and no longer need seqkit. Pass the list with `--taxa_list`.
- **BUSCO lineage.** By default BUSCO searches only eukaryotic lineages (`--busco_lineage auto_euk`), because these are eukaryotic bins.

One scientific point to keep in mind: the hints file uses `source=M`, as in your protocol. With AUGUSTUS's default settings, "M" hints are treated as almost certainly correct, so AUGUSTUS will mostly copy GeneMark's gene models. If you want AUGUSTUS to weigh the hints less, give it your own settings file with `--augustus_extrinsic_cfg`.

## What you need before you start

1. **Nextflow** (version 23.04 or newer) and **Java**.
2. **Conda or mamba** (the pipeline installs most tools by itself). Your conda channels should include `conda-forge` and `bioconda`.
3. Tools you have to install yourself and put on your PATH:
   - **GeneMark-ES** (needs a free licence key from the GeneMark website)
   - **InterProScan**, only if you want it (it needs a very large data download)
4. Databases (download once, then point to the folder):
   - **eggNOG-mapper** data folder → `--eggnog_db`
   - **EukCC** database → `--eukcc_db` (or skip EukCC with `--qc_tools busco`)
   - **BUSCO** lineages: downloaded automatically if the compute nodes have internet. If they don't, download them first and use `--busco_download_path`
   - **RepeatMasker**: the Dfam library for your group (for example Chlorophyta) has to be installed in RepeatMasker
   - **EUKulele** database: downloaded automatically, or point to one with `--eukulele_reference_dir`
   - **Kraken2** database: optional; Kraken2 runs only if you give `--kraken2_db`

## How to run it

**1. Make a samplesheet.** This is a CSV file with one line per sample. See `assets/samplesheet_example.csv`:

```
sample,fastq_1,fastq_2,contigs
siteA,/path/siteA_R1.fastq.gz,/path/siteA_R2.fastq.gz,
siteB,/path/siteB_R1.fastq.gz,/path/siteB_R2.fastq.gz,/path/siteB_assembly.fa
```

- `fastq_2` can be empty if you have single-end reads.
- `contigs` is optional. If you fill it in, that sample skips assembly and uses your contigs. The reads are still needed, because binning uses read coverage.

**2. Run it.** Each of these is one line you can paste.

Test the setup on an interactive node:

```
nextflow run main.nf -profile conda,local --input samples.csv --eggnog_db /path/eggnog_data --eukcc_db /path/eukcc_db --outdir results
```

On Dardel (each step goes to the Slurm queue as its own job):

```
nextflow run main.nf -profile conda,dardel --slurm_account naiss2025-x-xx --input samples.csv --eggnog_db /path/eggnog_data --eukcc_db /path/eukcc_db --outdir results -resume
```

Run Nextflow itself inside `screen` or `tmux` (or as a small long Slurm job), so it keeps going when you log out.

**3. If something fails,** fix the problem and run the *same command again* with `-resume` added at the end. Steps that already finished are not run again.

## Settings you will most likely change

You can set each of these on the command line (for example `--repeatmasker_species alveolata`), or change them in `nextflow.config`.

| Setting | Default | What it means |
|---------|---------|---------------|
| `--assembler` | `megahit` | or `metaspades` (needs paired reads) |
| `--min_contig_length` | `1000` | shorter contigs are thrown away before sorting |
| `--qc_tools` | `busco,eukcc` | which quality tools to run |
| `--busco_lineage` | `auto_euk` | or `auto`, or a fixed lineage such as `chlorophyta_odb10` |
| `--repeatmasker_species` | `chlorophyta` | your group's name in NCBI/Dfam |
| `--repeatmasker_lib` | none | use your own repeat library instead (for example from RepeatModeler) |
| `--genemark_min_contig` | `1000` | GeneMark ignores contigs shorter than this |
| `--augustus_species` | `chlamydomonas` | the closest species AUGUSTUS has a model for |
| `--annotation_tools` | `eggnog` | `eggnog`, `interproscan`, `eggnog,interproscan`, or `none` |
| `--eggnog_tax_scope` | `3041` | NCBI taxon ID to limit eggNOG to (3041 = Chlorophyta) |
| `--taxa_list` | none | text file, one taxon name per line → one protein file per taxon |
| `--skip_gene_prediction` | `false` | set to `true` to stop after the quality check |

**Important for breviates:** the defaults (`chlorophyta`, `chlamydomonas`, `3041`) come from the green-algae project. For other organisms, change all three.

**CPUs, memory and time** for each job are set in the `process` section of `nextflow.config`. Edit the numbers there if your jobs need more, or are waiting too long in the queue.

## Folder layout

```
main.nf                  the pipeline (order of steps)
nextflow.config          all settings, resources, and profiles (conda, dardel, slurm, local)
conf/conda.config        which conda package each step uses
modules/                 one file per group of steps
bin/                     small helper scripts (sorting contigs, renaming headers, fixing the GTF, ...)
envs/deepmicroclass.yml  DeepMicroClass is installed with pip, so it has its own environment file
assets/                  example samplesheet and example taxa list
```

After a run, `results/pipeline_info/report.html` shows how long each step took and how much memory it used.
