# GEO submission (4 h samples only)

## What's prepared

Run `scripts/build_geo_submission.R` (from the repo root) to (re)generate:

- **`processed_counttable_4h.txt`** — the featureCounts read-count matrix, restricted to the 21
  4 h samples used in the analysis. This is the "processed data" file GEO expects.
- **`GEO_metadata_4h.xlsx`** — a metadata workbook matching GEO's real `seq_template.xlsx` column
  layout, fully filled in:
  - `Series` — study title, summary, overall design, contributor list
  - `Samples` — one row per 4 h sample: library name, title, library strategy, organism,
    tissue/cell line/cell type, genotype, treatment, batch, molecule, single-or-paired-end,
    instrument model, description, processed data file, and the raw file columns (already
    matching the real filenames below)
  - `Protocols` — growth, treatment, extraction, library construction, 6 data-processing steps,
    genome build/assembly, and processed-file format/content, all with real text now (no
    placeholders left)

And run `scripts/gather_raw_fastq_4h.sh` **on the farm** (it needs `data/fastq/.../01.RawData/` or
`data/libs/`, neither of which exist locally) to produce, in `GEO_submission/raw_fastq_4h/`:

- 42 raw FASTQ files (21 samples x R1/R2), named `<unique_identifier>_R1.fastq.gz`/`_R2.fastq.gz` —
  already done, confirmed present on the farm.
- `md5sums_4h.txt` — MD5 checksums for all 42, already generated.

## What's still needed — the GEO portal side

Nothing left to fill in on this end. What's left is entirely on GEO's submission portal:

1. **Register a GEO account** (if you don't have one):
   https://www.ncbi.nlm.nih.gov/geo/info/submission.html → "GEO submission portal" (this now
   handles both the GEO metadata and the linked SRA raw-read deposit in one flow).
2. **Start a submission** and **upload** the 42 files in `GEO_submission/raw_fastq_4h/` plus
   `GEO_submission/processed_counttable_4h.txt` to the FTP/Aspera destination the portal gives you.
   Use `md5sums_4h.txt` to verify nothing got corrupted in transit.
3. **Transcribe `GEO_metadata_4h.xlsx`** into the portal's submission form (or attach it directly,
   if the portal accepts an uploaded template).
4. GEO curators review the submission and assign a GSE accession before it's released (either
   immediately or held private until a release date you specify, e.g. matching a paper's
   publication).
