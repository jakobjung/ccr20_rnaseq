## Title: locus_tag -> gene_name mapping for E. coli 536 (ECP_ locus tags)
## Author: Jakob Jung
## Description: E. coli 536's own GFF3 carries no gene= symbols on any CDS
##              (confirmed: 0/4685). Real gene names need an ortholog-mapping
##              step against a strain that HAS them - E. coli K-12 MG1655's
##              RefSeq annotation is complete/curated and 536 is close enough
##              to K-12 for reciprocal-best-hit BLASTP to work well. Run this
##              once; fuso_rnaseq.R reads the resulting CSV (falling back to
##              the locus tag itself for anything without a reciprocal hit).

suppressMessages({
  library(Biostrings)
  library(rtracklayer)
})

setwd("~/Documents/ccr20_PNA_rnaseq")

REF_DIR <- "data/reference_sequences"
k12_url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_translated_cds.faa.gz"
k12_faa_gz <- file.path(REF_DIR, "ecoli_k12_proteins.faa.gz")
k12_faa <- file.path(REF_DIR, "ecoli_k12_proteins.faa")
ecp_faa <- file.path(REF_DIR, "ecoli536_proteins.faa")

# --- 1. extract + translate every 536 CDS from the genome fasta ------------
cds <- readGFF(file.path(REF_DIR, "ecoli536_sRNAs_modified.gff3"),
               filter = list(type = "CDS"),
               tags = c("locus_tag"),
               columns = c("type", "start", "end", "strand", "seqid"))
cds <- as.data.frame(cds)

genome <- readDNAStringSet(file.path(REF_DIR, "ecoli536.fasta"))
contig <- genome[[1]]  # single-contig genome (CP000247.1)

extract_protein <- function(i) {
  s <- subseq(contig, start = cds$start[i], end = cds$end[i])
  if (cds$strand[i] == "-") s <- reverseComplement(s)
  as.character(translate(s, genetic.code = getGeneticCode("11")))
}
proteins <- vapply(seq_len(nrow(cds)), extract_protein, character(1))
names(proteins) <- cds$locus_tag
proteins <- sub("\\*$", "", proteins)  # drop the trailing stop codon

writeXStringSet(AAStringSet(proteins), ecp_faa)
cat("Wrote", length(proteins), "536 proteins to", ecp_faa, "\n")

# --- 2. get the K-12 MG1655 reference proteome (has real gene symbols) -----
if (!file.exists(k12_faa_gz)) {
  cat("Downloading E. coli K-12 MG1655 (GCF_000005845.2) protein FASTA...\n")
  download.file(k12_url, k12_faa_gz, mode = "wb")
}
if (!file.exists(k12_faa)) system(paste("gunzip -k", shQuote(k12_faa_gz)))

# K-12 headers look like:
#   >lcl|NC_000913.3_prot_NP_414542.1_1 [gene=thrL] [locus_tag=b0001] ...
# BLAST's qseqid/sseqid columns use the FULL first whitespace-delimited token
# ("lcl|NC_000913.3_prot_NP_414542.1_1"), not just the bare accession - the
# lookup key here must match that exactly, or the final merge silently joins
# on nothing.
k12_headers <- names(readAAStringSet(k12_faa))
k12_id <- sub("\\s.*", "", k12_headers)
k12_gene <- sub(".*\\[gene=([^]]+)\\].*", "\\1", k12_headers)
k12_locus <- sub(".*\\[locus_tag=([^]]+)\\].*", "\\1", k12_headers)
k12_lookup <- data.frame(k12_protein_id = k12_id, k12_locus_tag = k12_locus, gene_name = k12_gene)

# --- 3. reciprocal best-hit BLASTP (536 <-> K-12) ---------------------------
blast_dir <- file.path(REF_DIR, "blast_db")
dir.create(blast_dir, showWarnings = FALSE)
k12_db <- file.path(blast_dir, "ecoli_k12")
ecp_db <- file.path(blast_dir, "ecoli536")

system(paste("makeblastdb -in", shQuote(k12_faa), "-dbtype prot -out", shQuote(k12_db)))
system(paste("makeblastdb -in", shQuote(ecp_faa), "-dbtype prot -out", shQuote(ecp_db)))

fwd_out <- file.path(blast_dir, "536_vs_k12.tsv")  # 536 query -> K-12 subject
rev_out <- file.path(blast_dir, "k12_vs_536.tsv")  # K-12 query -> 536 subject
blast_cols <- "qseqid sseqid pident length evalue bitscore"

system(paste("blastp -query", shQuote(ecp_faa), "-db", shQuote(k12_db),
             "-out", shQuote(fwd_out), "-outfmt", shQuote(paste("6", blast_cols)),
             "-max_target_seqs 1 -evalue 1e-10 -num_threads 4"))
system(paste("blastp -query", shQuote(k12_faa), "-db", shQuote(ecp_db),
             "-out", shQuote(rev_out), "-outfmt", shQuote(paste("6", blast_cols)),
             "-max_target_seqs 1 -evalue 1e-10 -num_threads 4"))

col_names <- strsplit(blast_cols, " ")[[1]]
fwd <- read.table(fwd_out, col.names = col_names, stringsAsFactors = FALSE)
rev <- read.table(rev_out, col.names = col_names, stringsAsFactors = FALSE)

# best hit per query in each direction (blast already sorts by bitscore desc
# and -max_target_seqs 1 keeps only the top one, but dedupe defensively)
fwd_best <- fwd[!duplicated(fwd$qseqid), ]
rev_best <- rev[!duplicated(rev$qseqid), ]

# reciprocal best hit: 536 gene X's best K-12 hit is Y, AND K-12 gene Y's own
# best 536 hit is X again
merged <- merge(fwd_best, rev_best, by.x = c("qseqid", "sseqid"), by.y = c("sseqid", "qseqid"),
                suffixes = c("_fwd", "_rev"))
names(merged)[names(merged) == "qseqid"] <- "locus_tag"       # 536 locus tag
names(merged)[names(merged) == "sseqid"] <- "k12_protein_id"  # K-12 protein accession

out <- merge(merged[, c("locus_tag", "k12_protein_id", "pident_fwd", "evalue_fwd")],
            k12_lookup, by = "k12_protein_id")
out <- out[, c("locus_tag", "k12_locus_tag", "gene_name", "pident_fwd", "evalue_fwd")]
names(out)[names(out) %in% c("pident_fwd", "evalue_fwd")] <- c("pident", "evalue")
out <- out[order(out$locus_tag), ]

out_path <- file.path(REF_DIR, "ecp_gene_names.csv")
write.csv(out, out_path, row.names = FALSE)
cat("Wrote", nrow(out), "reciprocal best-hit gene names to", out_path,
    "(", round(100 * nrow(out) / nrow(cds), 1), "% of", nrow(cds), "536 CDS genes )\n")
