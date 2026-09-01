## Title: Prepare a GEO submission package for the 4h-only CCR20 PNA/ClbR RNA-seq data
## Author: Jakob Jung
## Description: Builds the two pieces of a GEO high-throughput-sequencing
##              submission that can be assembled without farm/cluster access:
##              (1) the processed data file - the featureCounts count matrix,
##              restricted to the 21 4h samples actually used downstream, and
##              (2) the metadata workbook (Series/Samples/Protocols) matching
##              GEO's seq_template.xlsx layout. Raw FASTQ files, their
##              checksums, and the actual FTP/SRA upload happen on the farm -
##              this only prepares what can be built from files in this repo.
##              See GEO_submission/README.md for the remaining manual steps.

library(readr)
library(dplyr)
library(writexl)

setwd("~/Documents/ccr20_PNA_rnaseq")
dir.create("GEO_submission", showWarnings = FALSE)

## 1. Processed data file - featureCounts matrix, 4h samples only, with clean
## sample-label column names instead of full bam paths
count_table <- read.delim("data/counttable_ccr20.txt", skip = 1)
colnames(count_table) <- gsub(".*[./]([^./]+)\\.bam$", "\\1", colnames(count_table))

sample_info <- read_csv("data/20260709_RNAseq_sample_table.csv", show_col_types = FALSE) %>%
  filter(treatment_time_h == 4)

stopifnot(all(sample_info$label %in% colnames(count_table)))
processed <- count_table[, c("Geneid", "Chr", "Start", "End", "Strand", "Length", sample_info$label)]
write.table(processed, "GEO_submission/processed_counttable_4h.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote processed count matrix:", nrow(processed), "features x", nrow(sample_info), "4h samples\n")

## 2. Metadata workbook (Series / Samples / Protocols)
strain_display <- c(wt = "CCR20+ (JVS-12599)", delta_clbR = "CCR20 ΔclbR (JVS-13616)")
genotype_display <- c(wt = "wild type", delta_clbR = "ΔclbR (clean deletion)")
pna_display <- c(H2O = "untreated (H2O)",
                 PPNA_clbR = "PNA 992 (on-target, targets clbR translation)",
                 PPNA_nt_control = "PNA 22 (non-targeting scrambled control)",
                 PPNA_clbR_mm = "PNA 2107 (sequence-matched mismatch control for PNA 992)")

## columns match GEO's actual seq_template.xlsx "Samples" tab exactly (no
## dedicated strain/time/replicate columns there - strain+time fold into
## genotype/title/description, and replicate is only distinguishable via
## library name/title, same as GEO itself expects for bacterial submissions).
## processed_file_2/raw_file_3/raw_file_4 stay empty - we only have one
## processed file (shared across all samples) and one raw file pair (R1/R2,
## lanes already merged) per sample - and are relabelled below to match the
## template's repeated "processed data file"/"raw file" headers, since a
## tibble can't hold duplicate column names but the xlsx sheet can.
samples <- sample_info %>%
  transmute(
    library_name = unique_identifier,
    title = paste0(strain_display[strain], " clone ", sub("clone_", "", clone),
                   ", ", pna_display[PNA], ", 4 h"),
    library_strategy = "RNA-Seq",
    organism = "Escherichia coli",
    tissue = "",
    cell_line = "",
    cell_type = "bacterial culture",
    genotype = paste0(strain_display[strain], ", ", genotype_display[strain]),
    treatment = paste0(pna_display[PNA], ", 4 h"),
    batch = "Novogene X208SC26074449-Z01-F001",
    molecule = "total RNA",
    single_or_paired_end = "paired-end",
    instrument_model = "[FILL IN: sequencer model - check the Novogene report]",
    description = paste0("RNA-seq of E. coli CCR20, ", pna_display[PNA],
                         ", 4 h post-treatment, biological replicate (clone) ",
                         sub("clone_", "", clone), "."),
    processed_file_1 = "processed_counttable_4h.txt",
    processed_file_2 = "",
    raw_file_1 = paste0(unique_identifier, "_R1.fastq.gz"),
    raw_file_2 = paste0(unique_identifier, "_R2.fastq.gz"),
    raw_file_3 = "",
    raw_file_4 = ""
  )
colnames(samples) <- c("library name", "title", "library strategy", "organism", "tissue",
                       "cell line", "cell type", "genotype", "treatment", "batch", "molecule",
                       "single or paired-end", "instrument model", "description",
                       "processed data file", "processed data file",
                       "raw file", "raw file", "raw file", "raw file")

series <- data.frame(
  field = c("title", "summary", "overall design", "contributor(s)", "organism", "supplementary file"),
  value = c(
    "RNA-Seq analysis of PNA-mediated knockdown of the colibactin regulator ClbR in uropathogenic Escherichia coli CCR20",
    paste("Uropathogenic Escherichia coli (UPEC) strains carrying the colibactin (pks) biosynthetic",
          "gene cluster produce colibactin, a genotoxin implicated in colorectal cancer. We used",
          "peptide nucleic acids (PNAs) to knock down translation of ClbR, the transcriptional",
          "activator of the pks cluster, in the UPEC strain CCR20, and profiled the transcriptome by",
          "RNA-seq to assess on-target colibactin cluster repression and off-target effects. Two",
          "genetic backgrounds were compared: CCR20+ (wild-type-like) and a CCR20 ΔclbR clean",
          "deletion mutant (a genetic specificity control), each treated for 4 h with H2O",
          "(untreated), PNA 992 (on-target, targets clbR translation; CCR20+ only, since ΔclbR",
          "lacks the target mRNA), PNA 22 (a non-targeting scrambled control), or PNA 2107 (a",
          "sequence-matched mismatch control for PNA 992)."),
    paste("RNA was extracted from E. coli CCR20+ and CCR20 ΔclbR cultures, in triplicate",
          "independent-clone biological replicates, 4 h after treatment with H2O, PNA 992, PNA 22,",
          "or PNA 2107 (CCR20+: all four treatments; CCR20 ΔclbR: H2O, PNA 22 and PNA 2107 only,",
          "since it lacks the PNA 992 target). 21 libraries were sequenced paired-end and analyzed",
          "to identify PNA 992-induced changes in colibactin cluster expression and off-target",
          "transcriptional effects."),
    "Nentwich, Sarah; Popella, Linda; Ghosh, Chandradhish; Adelmann, Juliane; Vogel, Joerg; Barquist, Lars; Hoebartner, Claudia; Jung, Jakob J",
    "Escherichia coli",
    "processed_counttable_4h.txt"
  ),
  stringsAsFactors = FALSE
)

protocols <- data.frame(
  protocol = c("growth protocol", "treatment protocol", "extraction protocol",
              "library construction protocol", "data processing steps"),
  text = c(
    "[FILL IN: growth medium, temperature, growth phase at treatment - not recorded in this repo]",
    "[FILL IN: PNA concentration, delivery method, exposure duration before harvest - not recorded in this repo]",
    "[FILL IN: RNA extraction kit/method, rRNA depletion if any - not recorded in this repo]",
    "[FILL IN: library prep kit, strandedness - not recorded in this repo]",
    paste(
      "1. Adapters removed and bases with Phred quality score <10 trimmed with BBDuk (BBTools v39.01).",
      "2. Trimmed reads mapped to the E. coli 536 (UPEC, GenBank CP000247.1) genome with BBMap (v39.01).",
      "3. Reads quantified against CDS, sRNA, tRNA and rRNA features with featureCounts (Subread v2.0.6; -p --countReadPairs).",
      "4. Processed data file: featureCounts read-count matrix, restricted to the 4 h samples.",
      sep = "\n"
    )
  ),
  stringsAsFactors = FALSE
)

write_xlsx(list(Series = series, Samples = samples, Protocols = protocols),
          "GEO_submission/GEO_metadata_4h.xlsx")
cat("Wrote metadata workbook: GEO_submission/GEO_metadata_4h.xlsx\n")
