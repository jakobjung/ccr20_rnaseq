## Title: RNA-seq analysis of PNA-mediated ClbR knockdown in E. coli CCR20 - 4h only
## Author: Jakob Jung
## Description: same DE analysis as scripts/fuso_rnaseq.R, but the 1h samples
##              are dropped before anything else runs (not just before PCA) -
##              filterByExpr, normalization, dispersion estimation and every
##              contrast below are all computed from the 4h subset alone.
##              Scoped to just what Sarah/Lars's email asked for: PCA, the 8
##              volcano comparisons, and the 6-comparison heatmap/supplement
##              table (no KEGG enrichment, no comprehensive multi-panel
##              exports - see fuso_rnaseq.R for the full analysis).

# Load libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(edgeR)
library(ggpubr)
library(readr)
library(MetBrewer)
library(viridis)
library(RColorBrewer)
library(BiocGenerics)
library(cowplot)
library(circlize)
library(ComplexHeatmap)
library("rtracklayer")

setwd("~/Documents/ccr20_PNA_rnaseq")

# PDF device with proper Unicode support (needed for the "Δ" in ΔclbR) - see
# fuso_rnaseq.R for why this shadows cairo_pdf() rather than using it directly.
cairo_pdf <- function(filename, width, height, ...) {
  before <- dev.cur()
  if (identical(Sys.info()[["sysname"]], "Darwin") && exists("quartz")) {
    tryCatch(quartz(type = "pdf", file = filename, width = width, height = height, ...),
             error = function(e) NULL)
    if (dev.cur() != before) return(invisible(NULL))
  }
  tryCatch(grDevices::cairo_pdf(filename, width = width, height = height, ...),
           error = function(e) NULL)
  if (dev.cur() != before) return(invisible(NULL))
  pdf(filename, width = width, height = height, ...)
}

dir.create("analysis_4h/volcano_plots", showWarnings = FALSE, recursive = TRUE)
dir.create("analysis_4h/heatmaps", showWarnings = FALSE, recursive = TRUE)

# Differential gene expression analysis
count_table <- read.delim("data/counttable_ccr20.txt", skip = 1)
gwc <- as.data.frame(count_table[,6:length(count_table)])
colnames(gwc) <- gsub(".*[./]([^./]+)\\.bam$", "\\1", colnames(gwc))
rownames(gwc) <- count_table$Geneid

# Sample metadata, exactly as in fuso_rnaseq.R
sample_info <- read_csv("data/20260709_RNAseq_sample_table.csv", col_names = TRUE) %>%
  transmute(label,
            background = recode(strain, wt = "CCR20", delta_clbR = "dClbR"),
            pna = recode(PNA, H2O = "H2O", PPNA_clbR = "PPNA_clbR",
                         PPNA_nt_control = "scr_acpP", PPNA_clbR_mm = "PPNAmm_ctrl_clbR"),
            time = paste0(treatment_time_h, "h"),
            clone = clone) %>%
  mutate(group = paste(background, pna, time, sep = "_"))

# drop the 1h samples HERE, before anything downstream (mapping stats, PCA,
# filterByExpr, normalization, dispersion) ever sees them - not just filtered
# out of individual plots afterward
sample_info <- sample_info[sample_info$time == "4h", ]
stopifnot(all(sample_info$label %in% colnames(gwc)))
gwc <- gwc[, c("Length", sample_info$label)]

# Mapping statistics (now naturally 4h-only, since gwc only has those columns)
gene_types <- readGFF("./data/reference_sequences/ecoli536_sRNAs_modified.gff3",
                      filter = list(type = c("CDS", "sRNA", "tRNA", "rRNA")),
                      tags = c("locus_tag"),
                      columns = c("type")) %>%
  as_tibble() %>%
  distinct(locus_tag, .keep_all = TRUE) %>%
  mutate(feature_group = ifelse(type %in% c("CDS", "sRNA"), "CDS & sRNA", as.character(type))) %>%
  select(locus_tag, feature_group)

mapping_stats <- gwc %>% select(-Length)
mapping_stats$locus_tag <- rownames(gwc)
mapping_stats <- mapping_stats %>%
  left_join(gene_types, by = "locus_tag") %>%
  filter(!is.na(feature_group)) %>%
  pivot_longer(-c(locus_tag, feature_group), names_to = "sample", values_to = "counts") %>%
  group_by(sample, feature_group) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  mutate(feature_group = factor(feature_group, levels = c("rRNA", "tRNA", "CDS & sRNA")))

g_mapping <- mapping_stats %>%
  ggplot(aes(x = sample, y = counts, fill = feature_group)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_pubr() +
  scale_fill_manual(values = c("rRNA" = "darkorange", "tRNA" = "steelblue", "CDS & sRNA" = "grey40")) +
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6), name = "reads (in million)") +
  theme(axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1)) +
  labs(x = "Sample", title = "CCR20 (4h) - reads assigned to CDS/sRNA/tRNA/rRNA features") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"),
        axis.title.x = element_blank(),
        legend.title = element_blank())

cairo_pdf("analysis_4h/RNA_mapping_stats.pdf", width = 12, height = 5)
g_mapping
dev.off()

# One DGEList/design/normalization/dispersion model, 4h samples only
y_all <- DGEList(counts = gwc[, sample_info$label],
                 group = factor(sample_info$group),
                 genes = gwc[, "Length", drop = FALSE])
keep_all <- filterByExpr(y_all)
print(table(keep_all))
y_all <- y_all[keep_all, , keep.lib.sizes = FALSE]

design_all <- model.matrix(~0 + group, data = y_all$samples)
colnames(design_all) <- levels(y_all$samples$group)
rownames(design_all) <- colnames(y_all)

y_all <- calcNormFactors(y_all, method = "TMM")
y_all <- estimateDisp(y_all, design_all, robust = TRUE)
fit_all <- glmQLFit(y_all, design_all, robust = TRUE)

logCPM_all <- cpm(y_all, log = TRUE, prior.count = 2)

# plot theme, treatment colours, display-name helper - same as fuso_rnaseq.R
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA, size = 1.5),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             title=element_text(colour="black", size=23),
             axis.text=element_text(colour="black", size=18),axis.ticks=element_line(colour="black"),
             axis.title=element_text(colour="black", size=21),
             plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = "top",
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 15),
             plot.title = element_text(size = 23, face = "bold", hjust=0.5),
                axis.line = element_line(size = 1)
              )

mycols <- c(H2O = "grey40", PPNA_clbR = "firebrick", scr_acpP = "steelblue", PPNAmm_ctrl_clbR = "darkorange")

pretty_name <- function(x) {
  x <- gsub("PPNAmm_ctrl_clbR", "PPNA-clbR-mm", x)
  x <- gsub("PPNA_clbR", "PPNA-clbR", x)
  x <- gsub("scr_acpP", "PPNA-nt-control", x)
  x <- gsub("dClbR", "ΔclbR", x)
  x <- gsub("_", " ", x)
  x <- gsub("vs", "vs.", x)
  x
}

## PCA - every timepoint is 4h already, so shape now marks BACKGROUND
## (CCR20 vs ΔclbR) instead of time, and there's no more hollow/filled trick
## for background - colour (treatment, mycols) plus shape (background) is
## enough on its own, so a plain solid-colour point shape works directly with
## `colour` (no separate `fill` aesthetic at all).
make_pca_plot <- function(pcx, pcy, sample_subset, suffix) {
  logCPM_sub <- logCPM_all[, sample_subset, drop = FALSE]
  info_sub <- sample_info[sample_subset, ]
  pca <- prcomp(t(logCPM_sub))
  pca_df <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0( colnames(pca_df), " (", paste( as.character(percentage), "%", ")", sep="") )
  pca_df$treatment <- factor(info_sub$pna, levels = names(mycols))
  pca_df$background <- factor(info_sub$background, levels = c("CCR20", "dClbR"))
  pca_df$clone <- gsub("clone_", "", info_sub$clone)
  present_treatments <- levels(droplevels(pca_df$treatment))
  xcol <- paste0("PC", pcx); ycol <- paste0("PC", pcy)
  y_range <- diff(range(pca_df[[ycol]]))

  p <- ggplot(pca_df, aes(x = .data[[xcol]], y = .data[[ycol]], colour = treatment, shape = background)) +
    geom_point(size = 10, stroke = 1.4) +
    geom_text_repel(aes(x = .data[[xcol]], y = .data[[ycol]], label = clone), inherit.aes = FALSE,
                    size = 3.5, colour = "black", nudge_y = 0.035 * y_range,
                    box.padding = 0.3, min.segment.length = 0, max.overlaps = Inf,
                    segment.color = "grey50", segment.size = 0.3, seed = 1) +
    theme +
    xlab(percentage[pcx]) +
    ylab(percentage[pcy]) +
    scale_color_manual(values = mycols[present_treatments], labels = pretty_name(present_treatments)) +
    scale_shape_manual(values = c(CCR20 = 17, dClbR = 16), labels = pretty_name(c("CCR20", "dClbR"))) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = 16)),
           shape = guide_legend(override.aes = list(colour = "black"))) +
    theme(plot.title = element_text(face = "bold", size = 40), legend.box = "vertical")

  file_name <- paste0("PCA_", suffix, if (!(pcx == 1 && pcy == 2)) paste0("_PC", pcx, "_vs_PC", pcy) else "")
  cairo_pdf(paste0("analysis_4h/", file_name, ".pdf"), width = 11, height = 12)
  print(p)
  dev.off()
  p
}

pca_variants <- list(all = rep(TRUE, nrow(sample_info)),
                     no_mm = sample_info$pna != "PPNAmm_ctrl_clbR")

for (suffix in names(pca_variants)) {
  subset_mask <- pca_variants[[suffix]]
  g_pca_12 <- make_pca_plot(pcx = 1, pcy = 2, subset_mask, suffix)
  g_pca_13 <- make_pca_plot(pcx = 1, pcy = 3, subset_mask, suffix)
  g_pca_23 <- make_pca_plot(pcx = 2, pcy = 3, subset_mask, suffix)

  cairo_pdf(paste0("analysis_4h/PCA_", suffix, "_PC1_till_3.pdf"), width = 11, height = 30)
  print(plot_grid(g_pca_12, g_pca_13, g_pca_23, ncol = 1,
                  scale = 0.95, labels = c("A", "B", "C"), label_size = 40))
  dev.off()
}

## Colibactin (pks) cluster genes and gene-name lookup - same curated list and
## same BLAST-vs-K12 mapping (scripts/build_gene_name_mapping.R) as
## fuso_rnaseq.R; not re-derived here.
clb_genes <- tibble::tribble(
  ~clb_gene, ~locus_tag, ~clb_product,
  "clbA", "ECP_1982", "4'-phosphopantetheinyl transferase",
  "clbR", "ECP_1981", "putative transcriptional regulator",
  "clbB", "ECP_1979", "putative polyketide synthase (NRPS-PKS hybrid, largest gene ~9.6 kb)",
  "clbC", "ECP_1978", "putative polyketide synthase",
  "clbD", "ECP_1977", "probable 3-hydroxybutyryl-CoA dehydrogenase",
  "clbE", "ECP_1976", "putative D-alanyl carrier protein",
  "clbF", "ECP_1975", "putative acyl-CoA dehydrogenase",
  "clbG", "ECP_1974", "putative malonyl-CoA transacylase",
  "clbH", "ECP_1973", "probable peptide synthetase protein",
  "clbI", "ECP_1972", "probable polyketide synthase protein",
  "clbJ", "ECP_1971", "putative non-ribosomal peptide synthase",
  "clbK", "ECP_1970", "putative non-ribosomal peptide synthase",
  "clbL", "ECP_1969", "probable amidase",
  "clbM", "ECP_1968", "hypothetical drug/sodium antiporter (MATE efflux)",
  "clbN", "ECP_1967", "non-ribosomal peptide synthetase",
  "clbO", "ECP_1966", "putative polyketide synthase",
  "clbP", "ECP_1965", "hypothetical protein",
  "clbQ", "ECP_1964", "putative thioesterase",
  "clbS", "ECP_1963", "hypothetical protein"
)
clb_gene_lookup <- setNames(clb_genes$clb_gene, clb_genes$locus_tag)

blast_gene_names_path <- "data/reference_sequences/ecp_gene_names.csv"
blast_gene_lookup <- if (file.exists(blast_gene_names_path)) {
  bgn <- read.csv(blast_gene_names_path, stringsAsFactors = FALSE)
  setNames(bgn$gene_name, bgn$locus_tag)
} else {
  c()
}

gene_name_lookup <- function(locus_tags) {
  out <- locus_tags
  in_clb <- locus_tags %in% names(clb_gene_lookup)
  out[in_clb] <- clb_gene_lookup[locus_tags[in_clb]]
  in_blast <- !in_clb & locus_tags %in% names(blast_gene_lookup)
  out[in_blast] <- blast_gene_lookup[locus_tags[in_blast]]
  out
}

## Contrasts - only the 8 comparisons from the email, flat (no more
## background/timepoint grouping needed - every sample here is already 4h)
lvl <- function(bg, pna) paste(bg, pna, "4h", sep = "_")

contrast_specs_4h <- list(
  PPNA_clbR_vs_H2O                          = c(lvl("CCR20", "PPNA_clbR"),  lvl("CCR20", "H2O")),
  PPNA_clbR_vs_scr_acpP                     = c(lvl("CCR20", "PPNA_clbR"),  lvl("CCR20", "scr_acpP")),
  PPNA_clbR_vs_PPNAmm_ctrl_clbR             = c(lvl("CCR20", "PPNA_clbR"),  lvl("CCR20", "PPNAmm_ctrl_clbR")),
  scr_acpP_vs_H2O                           = c(lvl("CCR20", "scr_acpP"),   lvl("CCR20", "H2O")),
  PPNAmm_ctrl_clbR_vs_H2O                   = c(lvl("CCR20", "PPNAmm_ctrl_clbR"), lvl("CCR20", "H2O")),
  PPNA_clbR_CCR20_vs_scr_acpP_dClbR         = c(lvl("CCR20", "PPNA_clbR"),  lvl("dClbR", "scr_acpP")),
  PPNA_clbR_CCR20_vs_PPNAmm_ctrl_clbR_dClbR = c(lvl("CCR20", "PPNA_clbR"),  lvl("dClbR", "PPNAmm_ctrl_clbR")),
  dClbR_vs_CCR20_H2O                        = c(lvl("dClbR", "H2O"),        lvl("CCR20", "H2O"))
)

volcano_contrasts <- names(contrast_specs_4h)  # all 8, in the requested order
heatmap_contrasts <- c("PPNA_clbR_vs_H2O", "PPNA_clbR_vs_scr_acpP", "PPNA_clbR_vs_PPNAmm_ctrl_clbR",
                       "scr_acpP_vs_H2O", "PPNAmm_ctrl_clbR_vs_H2O", "dClbR_vs_CCR20_H2O")

exprs <- sapply(contrast_specs_4h, function(pair) paste0(pair[1], " - ", pair[2]))
cm <- makeContrasts(contrasts = exprs, levels = design_all)
colnames(cm) <- names(contrast_specs_4h)

res_list <- lapply(colnames(cm), function(cn) {
  r <- glmQLFTest(fit_all, contrast = cm[, cn])
  r$table$FDR <- p.adjust(r$table$PValue, method = "fdr")
  r
})
names(res_list) <- colnames(cm)

# parse a group-level string (e.g. "dClbR_scr_acpP_4h") back into its
# (background, pna, time) parts - same approach as fuso_rnaseq.R
parse_level <- function(level_str) {
  for (bg in c("CCR20", "dClbR")) {
    prefix <- paste0(bg, "_"); suffix <- "_4h"
    if (startsWith(level_str, prefix) && endsWith(level_str, suffix)) {
      return(list(bg = bg, pna = substr(level_str, nchar(prefix) + 1, nchar(level_str) - nchar(suffix))))
    }
  }
  stop("parse_level: could not parse '", level_str, "'")
}

# every contrast here is 4h, so the time is never shown in the title - it'd
# just be redundant on every single panel
contrast_title <- function(pair) {
  a <- parse_level(pair[1]); b <- parse_level(pair[2])
  paste0(pretty_name(a$pna), " ", pretty_name(a$bg), " vs. ", pretty_name(b$pna), " ", pretty_name(b$bg))
}

## Volcano plots - do_volcano(), copied verbatim from fuso_rnaseq.R
do_volcano <- function(restab, targetgene = NULL, pointsize = 2, x_limit = F, y_limit = F, show_sig = F,
                       alpha = 0.05, color_sig = T, marked_gene_names = NULL, marked_gene_title = NULL,
                       minlogfc = 1, title = "Volcano", marked_genes = NULL, add_labels = T, gene_names = NULL,
                       color_threshold_lines = "black",
                       cols = c("target" = "darkorange", "marked_sig" = "red", "marked_nonsig" = "darkred",
                                "up" = "darkorange", "down" = "steelblue", "sRNA" = "darkred",
                                "other" = "darkgrey")) {
  if (!is.null(gene_names)) {
      rownames(restab) <- gene_names
  }

  g <- ggplot(restab) +
    geom_point(
      data = restab,
      aes(x = logFC, y = -log10(FDR), fill = "other"),shape = 21, color="darkgrey",
      cex = pointsize, alpha = 0.4
    ) +
    theme_bw() +
    geom_hline(yintercept = -log10(alpha),
               color = color_threshold_lines, linetype = 3) +
    geom_vline(xintercept = c(-minlogfc, minlogfc),
               color = color_threshold_lines, linetype = 3) +
    theme(axis.title.x = element_text(size = 15),
          legend.position = "none",
          axis.title.y = element_text(size = 15),
          axis.text = element_text(size = 10, colour = "black"),
          panel.background = element_rect(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
    ggtitle(title) +
    xlab(expression("log"[2] * " fold change")) +
    ylab(expression("- log"[10] * " P-value (FDR)"))

  if (x_limit) {
    g <- g + scale_x_continuous(expand = c(0, 0), breaks = c(-10, -5, -2, 0, 2, 5, 10),
                                limits = c(-x_limit, x_limit), oob = scales::squish)
  }
  if (y_limit) {
    g <- g + scale_y_continuous(expand = expansion(mult = c(0, 0.03)), breaks = seq(0, 40, 5),
                                limits = c(0, y_limit), oob = scales::squish)
  }

  if (color_sig) {
    g <- g +
      geom_point(
        data = restab[restab$FDR < alpha & restab$logFC < -minlogfc,],
        aes(x = logFC, y = -log10(FDR), fill = "down"),alpha = 0.5,
        cex = pointsize,shape=21) +
      geom_point(
        data = restab[restab$FDR < alpha & restab$logFC > minlogfc,],
        aes(x = logFC, y = -log10(FDR), fill = "up"),alpha = 0.5,
        cex = pointsize,shape=21)
  }

  if (!is.null(marked_genes)) {
    marked_data <- restab[marked_genes,]
    marked_data$sig_marked <- ifelse(marked_data$FDR < alpha & abs(marked_data$logFC) > minlogfc,
                                     "marked_sig", "marked_nonsig")
    g <- g + geom_point(
      data = marked_data,
      aes(x = logFC, y = -log10(FDR), fill = sig_marked),
      cex = pointsize * 1.4, shape=21,alpha = 0.85)
    if (!is.null(marked_gene_names)) {
      marked_data$label <- marked_gene_names[rownames(marked_data)]
      marked_up <- marked_data[marked_data$logFC > 0, ]
      marked_down <- marked_data[marked_data$logFC <= 0, ]
      g <- g + geom_text_repel(
        data = marked_up,
        aes(x = logFC, y = -log10(FDR), label = label),
        color = "darkred", fontface = "italic",
        nudge_x = x_limit - marked_up$logFC,
        min.segment.length = 0, direction = "y", hjust = 1, size = 4,
        segment.color = "darkred", segment.alpha = 0.5,
        max.overlaps = Inf, parse = F) +
      geom_text_repel(
        data = marked_down,
        aes(x = logFC, y = -log10(FDR), label = label),
        color = "darkred", fontface = "italic",
        nudge_x = -x_limit - marked_down$logFC,
        min.segment.length = 0, direction = "y", hjust = 0, size = 4,
        segment.color = "darkred", segment.alpha = 0.5,
        max.overlaps = Inf, parse = F)
    }
  }

  if (show_sig) {
    top_up <- restab[which(restab$FDR < alpha & restab$logFC > minlogfc),]
    top_down <- restab[which(restab$FDR < alpha & restab$logFC < -(minlogfc)),]
    if (length(rownames(top_up)) > 0 && (length(rownames(top_up)) > 3)) {
      top_up <- top_up[order(-top_up$logFC),][1:10,]
    }
    if (length(rownames(top_down)) > 0 && (length(rownames(top_down)) > 3)) {
      top_down <- top_down[order(top_down$logFC),][1:10,]
    }
    top_peaks <- rbind(top_up, top_down)
    top_peaks <- na.omit(top_peaks)
    g_labels <- c(targetgene, rownames(top_peaks))
    if (!is.null(marked_genes)) g_labels <- setdiff(g_labels, marked_genes)
  }

  if (add_labels == T) {
    g_labels <- unique(g_labels)
    rup <- restab[g_labels,][restab[g_labels,]$logFC > 0,]
    g <- g + geom_text_repel(
      data = rup,
      aes(x = logFC, y = -log10(FDR), label = rownames(rup)),
      fontface = "italic", nudge_x = x_limit - rup$logFC,
      min.segment.length = 0, direction = "y", hjust = 1, size = 3,
      segment.color = "black", segment.alpha = 0.5, parse = F
    )
    rdown <- restab[g_labels,][restab[g_labels,]$logFC < 0,]
    g <- g + geom_text_repel(
      data = rdown,
      aes(x = logFC, y = -log10(FDR), label = rownames(rdown)),
      fontface = "italic", nudge_x = -x_limit - rdown$logFC,
      min.segment.length = 0, direction = "y", hjust = 0, size = 3,
      segment.color = "black", segment.alpha = 0.5, parse = F)
  }
  g + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
}

# 01_ .. 08_ prefixes preserve the email's requested order in a plain file listing
volcano_file_stubs <- c(
  "01_PPNA-clbR_vs_H2O_CCR20", "02_PPNA-clbR_vs_scr-acpP_CCR20", "03_PPNA-clbR_vs_PPNAclbRmm_CCR20",
  "04_scr-acpP_vs_H2O_CCR20", "05_PPNAclbRmm_vs_H2O_CCR20", "06_PPNA-clbR_CCR20_vs_scr-acpP_dClbR",
  "07_PPNA-clbR_CCR20_vs_PPNAclbRmm_dClbR", "08_dClbR_vs_CCR20_H2O"
)
names(volcano_file_stubs) <- volcano_contrasts

for (cn in volcano_contrasts) {
  restab <- res_list[[cn]]$table
  ttle <- contrast_title(contrast_specs_4h[[cn]])
  volc <- do_volcano(restab, pointsize = 2, x_limit = 12, y_limit = 10, show_sig = F,
                     alpha = 0.001, color_sig = T, title = ttle,
                     minlogfc = 2, add_labels = F,
                     color_threshold_lines = "black",
                     marked_genes = intersect(clb_genes$locus_tag, rownames(restab)),
                     marked_gene_names = clb_gene_lookup)
  base_path <- paste0("analysis_4h/volcano_plots/", volcano_file_stubs[[cn]])
  ggsave(paste0(base_path, ".png"), volc, width = 12, height = 12, units = "cm", dpi = 300)
  ggsave(paste0(base_path, ".pdf"), volc, device = cairo_pdf, width = 12, height = 12, units = "cm")
}

## Heatmaps + supplementary table - same split/gene-name/adaptive-padding
## design as fuso_rnaseq.R's curated 4h deliverable, just against this
## script's flat (not group/cn-nested) contrast list
top_updown_genes <- function(restab, exclude, n = 10) {
  tab <- restab[!(rownames(restab) %in% exclude), ]
  tab <- tab[tab$FDR < 0.001 & abs(tab$logFC) > 2, ]
  up_tab <- tab[tab$logFC > 0, ]
  down_tab <- tab[tab$logFC < 0, ]
  c(rownames(up_tab[order(-up_tab$logFC), ])[seq_len(min(n, nrow(up_tab)))],
    rownames(down_tab[order(down_tab$logFC), ])[seq_len(min(n, nrow(down_tab)))])
}

build_logfc_matrix <- function(genes, contrast_names) {
  m <- sapply(contrast_names, function(cn) res_list[[cn]]$table[genes, "logFC"])
  m <- matrix(m, nrow = length(genes), dimnames = list(genes, NULL))
  colnames(m) <- sapply(contrast_names, function(cn) contrast_title(contrast_specs_4h[[cn]]))
  m
}

build_sig_matrix <- function(genes, contrast_names) {
  m <- sapply(contrast_names, function(cn) {
    tab <- res_list[[cn]]$table[genes, ]
    ifelse(tab$FDR < 0.001 & abs(tab$logFC) > 2, "*", "")
  })
  matrix(m, nrow = length(genes), dimnames = list(genes, NULL))
}

build_fdr_matrix <- function(genes, contrast_names) {
  m <- sapply(contrast_names, function(cn) res_list[[cn]]$table[genes, "FDR"])
  matrix(m, nrow = length(genes), dimnames = list(genes, NULL))
}

mark_sig <- function(sig_mat) {
  function(j, i, x, y, width, height, fill) {
    grid.text(pindex(sig_mat, i, j), x, y)
  }
}

save_heatmap_pdf <- function(ht, file_name, width, height, col_labels, fontsize = 9) {
  grDevices::pdf(NULL)
  longest <- col_labels[which.max(nchar(col_labels))]
  label_reach_mm <- convertWidth(grobWidth(textGrob(longest, gp = gpar(fontsize = fontsize))),
                                 "mm", valueOnly = TRUE) * cos(pi / 4)
  dev.off()

  cairo_pdf(paste0("analysis_4h/heatmaps/", file_name, ".pdf"), width = width, height = height)
  draw(ht, padding = unit(c(2, max(25, label_reach_mm + 5), 2, 2), "mm"))
  dev.off()
}

missing_clb <- setdiff(clb_genes$locus_tag, rownames(res_list[[1]]$table))
if (length(missing_clb) > 0) {
  warning("clb genes filtered out as lowly expressed, excluded from heatmap: ",
          paste(missing_clb, collapse = ", "))
}
clb_present <- clb_genes[clb_genes$locus_tag %in% rownames(res_list[[1]]$table), ]

top20_genes <- top_updown_genes(res_list$PPNA_clbR_vs_H2O$table, exclude = clb_genes$locus_tag, n = 10)

make_clb_only_heatmap <- function(contrast_names, file_name, plot_title) {
  hm_clb <- build_logfc_matrix(clb_present$locus_tag, contrast_names)
  sig_clb <- build_sig_matrix(rownames(hm_clb), contrast_names)
  rownames(hm_clb) <- gene_name_lookup(rownames(hm_clb))

  col_fun <- colorRamp2(c(-12,-4, 0, 4), c("darkblue","steelblue", "white", "darkorange"))
  ht_clb <- Heatmap(hm_clb, name = "log2FC", col = col_fun,
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    row_title_side = "right", row_title_rot = 0,
                    row_title_gp = gpar(fontsize = 9),
                    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
                    column_names_gp = gpar(fontsize = 9), column_names_rot = 45,
                    column_title = plot_title, border = TRUE, layer_fun = mark_sig(sig_clb))

  save_heatmap_pdf(ht_clb, file_name, width = 0.55 * ncol(hm_clb) + 2,
                   height = 0.22 * nrow(hm_clb) + 2, col_labels = colnames(hm_clb))
  ht_clb
}

make_top20_heatmap <- function(contrast_names, genes, file_name, plot_title) {
  hm_top20 <- build_logfc_matrix(genes, contrast_names)
  sig_top20 <- build_sig_matrix(rownames(hm_top20), contrast_names)
  rownames(hm_top20) <- gene_name_lookup(rownames(hm_top20))

  col_fun <- colorRamp2(c(-4, 0, 4), c("steelblue", "white", "darkorange"))
  ht_top20 <- Heatmap(hm_top20, name = "log2FC", col = col_fun,
                      cluster_rows = FALSE, cluster_columns = FALSE,
                       row_title_side = "right", row_title_rot = 0,
                      row_title_gp = gpar(fontsize = 9),
                      row_names_gp = gpar(fontsize = 9, fontface = "italic"),
                      column_names_gp = gpar(fontsize = 9), column_names_rot = 45,
                      column_title = plot_title, border = TRUE, layer_fun = mark_sig(sig_top20))

  save_heatmap_pdf(ht_top20, file_name, width = 0.55 * ncol(hm_top20) + 2,
                   height = 0.22 * nrow(hm_top20) + 2, col_labels = colnames(hm_top20))
  ht_top20
}

make_clb_only_heatmap(heatmap_contrasts, "clb_cluster_heatmap_4h", "Colibactin cluster (4 hr)")
make_top20_heatmap(heatmap_contrasts, top20_genes, "top20_genes_heatmap_4h", "Top 10 up / top 10 down DE genes (4 hr)")

write_clb_top20_supplement <- function(contrast_names, genes_clb, genes_top20, path) {
  contrast_cols <- sapply(contrast_names, function(cn) contrast_title(contrast_specs_4h[[cn]]))
  logfc_names <- paste0(contrast_cols, " log2FC")
  fdr_names <- paste0(contrast_cols, " FDR")

  build_block <- function(genes, gene_set_label) {
    hm <- build_logfc_matrix(genes, contrast_names)
    fdr <- build_fdr_matrix(genes, contrast_names)
    colnames(hm) <- logfc_names
    colnames(fdr) <- fdr_names
    cbind(data.frame(locus_tag = genes, gene_name = gene_name_lookup(genes),
                     gene_set = gene_set_label, stringsAsFactors = FALSE),
         as.data.frame(hm), as.data.frame(fdr))
  }

  out <- rbind(build_block(genes_clb, "pks_cluster"), build_block(genes_top20, "top_DE_genes"))
  write.csv(out, path, row.names = FALSE)
  out
}

write_clb_top20_supplement(heatmap_contrasts, clb_present$locus_tag, top20_genes,
                           "analysis_4h/heatmaps/clb_and_top20_supplement_4h.csv")
