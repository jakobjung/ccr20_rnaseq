# RNA-Seq analysis of PNA-mediated knockdown of the colibactin regulator ClbR in *E. coli* CCR20

- Project name: Antisense PNA targeting of the colibactin transcriptional activator ClbR in the
  UPEC strain CCR20

- Experiments/samples: Sarah Nentwich, Linda Popella

- Supervision: Jörg Vogel, Lars Barquist

- Data analysis/RNA-Seq analysis: Jakob J. Jung

- Start: 2026-08

## Introduction

RNA-Seq analysis of the uropathogenic *E. coli* (UPEC) strain CCR20, treated with peptide nucleic
acids (PNAs) targeting translation of ClbR, the transcriptional activator of the colibactin (*pks*)
biosynthetic gene cluster. Two backgrounds were compared: CCR20+ (JVS-12599, wild-type-like) and
CCR20 ΔclbR (JVS-13616, clean deletion of *clbR*), each in 3 independent clones. Samples were
treated for 1 h or 4 h with:

- **H2O** (untreated control)
- **PNA 992** (PPNA-clbR) — the on-target PNA, targets the translation of *clbR*
- **PNA 22** (PPNA-nt-control, scr-*acpP*) — scrambled/non-targeting control
- **PNA 2107** (PPNA-clbR-mm) — sequence-matched mismatch control for PNA 992

The ΔclbR background was not treated with PNA 992, since it has no *clbR* mRNA to target; it
instead serves as a genetic control for on-target specificity. PNA 2107 is nominally the more
rigorous control (same sequence as PNA 992 but mismatched), but in other experiments it has shown
some unspecific, mildly growth-inhibitory effects — so part of this analysis is checking whether
that holds here, in which case PNA 22 (scr-*acpP*) should be used as the main control for the
transcriptomic comparison instead.

In total 42 samples (2 backgrounds x up to 4 treatments x 2 time points x 3 clones, unbalanced —
see [data/20260709_RNAseq-Samples_Label.xlsx](data/20260709_RNAseq-Samples_Label.xlsx) for the
exact design). Sample 35 (ΔclbR, PNA 2107, 4 h) failed Novogene's QC but was sequenced anyway, so
treat its results with caution.

The main biological question: does PNA 992 cause broad downregulation of the colibactin gene
cluster while keeping off-target transcriptional effects to a minimum?

Sequencing was done by Novogene (batch `X208SC26074449-Z01-F001`), each sample split across two
flowcells/lanes (a top-up run), paired-end.

## Directory structure

The project is divided in 3 main directories:

-   [data](data) : contains all raw, intermediate and final data of the project.
-   [analysis](analysis): contains analysis files, such as figures, plots, tables, etc.
-   [scripts](scripts): contains scripts to process and analyze data from the data directory.

Some directories have their own README.md file with information on the respective files.

## Workflow

Here I describe the workflow to reproduce the results. Mapping is submitted as cluster jobs,
everything downstream (DE analysis) runs locally in R.

### 1. Prerequisites

For running the whole analysis, one needs the following packages/tools/software (make sure
they're on `$PATH`, e.g. via `module load` or a conda env — edit the placeholder near the top of
`trimm_map_BB.sh`):

- **BBMap** (v38.84) & **BBDuk**
- **samtools** (v1.12)
- **featureCounts** (Subread, v2.0.3+ for `--countReadPairs` support; e.g. 2.0.6)
- **FastQC**
- **R** (v4.1.1+) with Bioconductor/CRAN packages for the downstream DE analysis
- Access to a cluster job scheduler (the script submits jobs via `bsub`; adapt this if yours uses
  something else)

### 2. Mapping

All raw FastQ files live under
[./data/fastq/X208SC26074449-Z01-F001/01.RawData/](./data/fastq/X208SC26074449-Z01-F001/01.RawData/),
one subdirectory per sample, each containing paired-end reads split across two
flowcells/lanes. Sample-to-condition mapping is in
[./data/20260709_RNAseq-Samples_Label.xlsx](data/20260709_RNAseq-Samples_Label.xlsx).

Reads are mapped against *E. coli* 536 (CP000247.1), a UPEC reference strain carrying the
colibactin island (reference fasta/gff in
[./data/reference_sequences/](./data/reference_sequences/)). Note the annotation does not carry
*clbA-clbS* gene symbols — locate the colibactin/pks cluster by searching the gff for
`polyketide`/`non-ribosomal peptide` in the product field, or by coordinates from the literature.
Since featureCounts doesn't reliably parse GFF3's `key=value` attribute syntax, the counting step
first derives a plain SAF file (`GeneID`/`Chr`/`Start`/`End`/`Strand`) from the gff3's
CDS/sRNA/tRNA/rRNA rows rather than passing the gff3 to featureCounts directly.

[./scripts/trimm_map_BB.sh](./scripts/trimm_map_BB.sh) loops through all samples: concatenating
lanes, trimming with bbduk, mapping with bbmap, and sorting/indexing the resulting bam, then counts
reads over CDS/sRNA/tRNA/rRNA features once all samples are mapped. All its paths are resolved
relative to the script's own location, so it can be submitted as a single job from anywhere, e.g.
from the repo root:

```bash
bsub -q long -n 4 -M 10000 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" \
     -o scripts/stdout.%J -e scripts/stderr.%J bash scripts/trimm_map_BB.sh
```

(`normal`'s 12h run limit may be tight running 42 samples sequentially; `long` gives 48h.)

This produces, per sample, a sorted+indexed BAM in [./data/rna_align](./data/rna_align), and a
combined count table [./data/rna_align/counttable_ccr20.txt](./data/rna_align/counttable_ccr20.txt).
Trimmed libraries and per-sample FastQC reports are kept in [./data/libs](./data/libs).

The job's stdout log doubles as bbduk/bbmap's own log output, so
[./scripts/get_mapping_stats.sh](./scripts/get_mapping_stats.sh) can be run directly against it to
get a summary trimming/mapping table:

```bash
bash scripts/get_mapping_stats.sh scripts/stdout.<jobid>
```

### 3. Differential expression analysis

To run the differential expression analysis, run the R script
[./scripts/fuso_rnaseq.R](./scripts/fuso_rnaseq.R), adapted to import
[./data/rna_align/counttable_ccr20.txt](./data/rna_align/counttable_ccr20.txt) and the sample sheet
above, and to set up contrasts for CCR20+ vs. CCR20 ΔclbR crossed with PNA 992/22/2107/H2O at 1 h
and 4 h. This outputs figures/tables to the [./analysis](./analysis) directory, including a check
of PNA 2107 off-target effects (to decide whether PNA 22 should be used as the main control) and
the pks/colibactin cluster response to PNA 992 in the CCR20+ background.
