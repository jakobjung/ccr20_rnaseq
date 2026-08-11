## Title: RNA-seq analysis of Fusobacterium nucleatum
## Author: Jakob Jung
## Date: 2023-12-15
## Description: RNA-seq analysis of Fusobacterium nucleatum after PNA targeting, at 2 different timepoints

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
library('RUVSeq')
library(circlize)
library(ComplexHeatmap)
library(xlsx)
library(KEGGREST)
library("rtracklayer")

# Mapping statistics

# percentages of mapped read types of the 2 strains
# start with Fusobacterium nucleatum ATCC 23726
cds_counts_FNN <- unlist(read.delim("data/rna_align/counttable_fnn23.txt.summary", row.names = 1)[1,])
rrna_counts_FNN <- unlist(read.delim("data/rna_align/counttable_fnn23_rRNAs.txt.summary", row.names = 1)[1,])
tRNA_counts_FNN <- unlist(read.delim("data/rna_align/counttable_fnn23_tRNAs.txt.summary", row.names = 1)[1,])

# do same for Fusobacterium nucleatum vincentii
cds_counts_FNV <- unlist(read.delim("data/rna_align/counttable_fnv.txt.summary", row.names = 1)[1,])
rrna_counts_FNV <- unlist(read.delim("data/rna_align/counttable_fnv_rRNAs.txt.summary", row.names = 1)[1,])
tRNA_counts_FNV <- unlist(read.delim("data/rna_align/counttable_fnv_tRNAs.txt.summary", row.names = 1)[1,])

# create a dataframe with all the data
mapping_stats <- data.frame(
  counts = c(cds_counts_FNN, rrna_counts_FNN, tRNA_counts_FNN, cds_counts_FNV, rrna_counts_FNV, tRNA_counts_FNV),
  rna_type = c(rep("CDS+sRNA", length(cds_counts_FNN)), rep("rRNA", length(rrna_counts_FNN)),
               rep("tRNA", length(tRNA_counts_FNN)), rep("CDS+sRNA", length(cds_counts_FNV)),
               rep("rRNA", length(rrna_counts_FNV)), rep("tRNA", length(tRNA_counts_FNV))),
  strain = c(rep("FNN", length(cds_counts_FNN) + length(rrna_counts_FNN) + length(tRNA_counts_FNN)),
               rep("FNV", length(cds_counts_FNV) + length(rrna_counts_FNV) + length(tRNA_counts_FNV))),
  sample = c(names(cds_counts_FNN), names(rrna_counts_FNN), names(tRNA_counts_FNN),
             names(cds_counts_FNV), names(rrna_counts_FNV), names(tRNA_counts_FNV))
) %>%
  mutate(sample = gsub(".+rna_align\\.ID_\\d+_[^_]+_([^_]+_\\d+_[^_]+_[^_]+)\\.fq\\.gz\\.bam","\\1",  sample))

mapping_stats$rna_type <- factor(mapping_stats$rna_type, levels = c("tRNA", "rRNA", "CDS+sRNA"))

# now create 2 plots, 1 for each strain
g_FNN <- mapping_stats %>%
  filter(strain == "FNN") %>%
  ggplot(aes(x = sample, y = counts, fill = rna_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = met.brewer("Hokusai3")[c(1,3,5)]) +
  theme_pubr()+
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6), name = "reads (in million)",
                     breaks = seq(0,10000000, by = 1000000), limits = c(0, 10000000))+
  theme(axis.text.x = element_text(angle = 90,  size = 9, vjust = 0.5, hjust = 1)) +
  labs(x = "Sample", y = "Percentage of reads", title = "Fusobacterium nucleatum") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"),
        # remove x axis title
        axis.title.x = element_blank(),
  # remove legend title
    legend.title = element_blank())

g_FNN

g_FNV <- mapping_stats %>%
  filter(strain == "FNV") %>%
  ggplot(aes(x = sample, y = counts, fill = rna_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = met.brewer("Hokusai3")[c(1,3,5)]) +
  theme_pubr()+
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6), name = "reads (in million)",
                     breaks = seq(0,10000000, by = 1000000), limits = c(0, 10000000))+
  theme(axis.text.x = element_text(angle = 90,  size = 9, vjust = 0.5, hjust = 1)) +
  labs(x = "Sample", y = "Percentage of reads", title = "Fusobacterium vincentii") +
  # make title italic
  theme(plot.title = element_text(hjust = 0.5, face = "italic"),
        # remove x axis title
        axis.title.x = element_blank(),
  # remove legend title
    legend.title = element_blank())

g_FNV

# save the plots using cowplot
svg("analysis/RNA_mapping_stats.svg", width = 10, height = 5)
plot_grid(g_FNN, g_FNV, ncol = 2, labels = c("A", "B"), label_size = 14, label_x = 0.05, label_y = 0.95, label_fontface = "bold")
dev.off()

pdf("analysis/RNA_mapping_stats.pdf", width = 10, height = 5)
plot_grid(g_FNN, g_FNV, ncol = 2, labels = c("A", "B"), label_size = 14, label_x = 0.05, label_y = 0.95, label_fontface = "bold",
          # zoom out
          scale = 0.9)
dev.off()

# Differential gene expression analysis
count_table_FNN <- read.delim("data/rna_align/counttable_fnn23.txt", skip = 1)
count_table_FNV <- read.delim("data/rna_align/counttable_fnv.txt", skip = 1)

# get gene wise counts
gwc_FNN <- as.data.frame(count_table_FNN[,6:length(count_table_FNN)])
gwc_FNV <- as.data.frame(count_table_FNV[,6:length(count_table_FNV)])

# change colnames:
pnapat <- ".+rna_align\\.ID_\\d+_[^_]+_([^_]+_\\d+_[^_]+_[^_]+)\\.fq\\.gz\\.bam"
colnames(gwc_FNN) <- gsub(pnapat, "\\1", colnames(gwc_FNN))
colnames(gwc_FNV) <- gsub(pnapat, "\\1", colnames(gwc_FNV))

# change colnames, make from _min_ to m and _h to h
colnames(gwc_FNN) <- gsub("_min", "m", colnames(gwc_FNN))
colnames(gwc_FNV) <- gsub("_min", "m", colnames(gwc_FNV))
colnames(gwc_FNN) <- gsub("_h", "h", colnames(gwc_FNN))
colnames(gwc_FNV) <- gsub("_h", "h", colnames(gwc_FNV))


rownames(gwc_FNN) <- count_table_FNN$Geneid
rownames(gwc_FNV) <- count_table_FNV$Geneid


# get mapping stats
mapstats <- read_tsv("./scripts/mapping_stats.txt", col_names = T)
mapstats$Sample <- gsub("ID_\\d+_(.+).fq.gz", "\\1", mapstats$Sample)

# create stacked barplot for mapping stats
mstats_plot <- mapstats %>% mutate(Total_input_reads=Total_input_reads-Total_mapped_reads) %>%
  pivot_longer(2:3, names_to = "type", values_to = "counts") %>%
  ggplot(aes(x = Sample, y = counts, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c( "steelblue", "lightskyblue")) +
    theme_pubr() +
    labs(x = "Sample", y = "Reads", title = "Mapping statistics") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          axis.text.x = element_text(angle = 90,  size = 9, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          legend.title = element_blank()) +
  # change scale of y axis
  #ylim(0, 10000000) +
  # change names and ticks of y axis
    scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6), name = "reads (in million)",
                         breaks = seq(0,100000000, by = 1000000))
ggsave("analysis/mapping_stats.pdf", mstats_plot, width = 9, height = 7)

# get gene lengths
gene_lengths_FNN <- gwc_FNN$Length
gene_lengths_FNV <- gwc_FNV$Length

# Gene length vs. expression:
# normalized (TPM)
gwc_lex <- gwc_FNN[, grepl("(H2O)|(Length)", colnames(gwc_FNN))]
lex <- data.frame(length=gwc_lex$Length, counts=rowMeans(gwc_lex[,-1]))

gwcnorm_length <- data.frame(sapply(gwc_lex[,-1], function(x) x / (gwc_FNN[,1]/1000)))
gwc_tpm <- data.frame(Length = gwc_lex$Length,
                      sapply(gwcnorm_length, function(x) x * 1e6 / sum(x)), row.names = rownames(gwc_lex))
tpm <- data.frame(length=gwc_tpm$Length, counts=rowMeans(gwc_tpm[,-1]))

# create plots

lve_raw <- lex %>% ggplot(aes(x=length, y=log10(counts+1))) + geom_point() +
  scale_x_continuous(limits = c(0,1000)) + scale_y_continuous(limits = c(-0.5,5)) + theme_minimal() +
  geom_text_repel(aes( label=ifelse(length<80, rownames(lex), "")), size=2.5, max.overlaps = 15)+
  # rename y and x axis
    labs(x = "Gene length (bp)", y = "log10 counts (raw)")

lwe_tpm <- tpm %>% ggplot(aes(x=length, y=log10(counts+1))) + geom_point() +
  scale_x_continuous(limits = c(0,1000)) + scale_y_continuous(limits = c(0,5)) + theme_minimal() +
  geom_text_repel(aes( label=ifelse(length<80, rownames(lex), "")), size=2.5, max.overlaps = 15)+
  # rename y and x axis
    labs(x = "Gene length (bp)", y = "log10 counts (TPM normalized)")

pdf("analysis/gene_length_vs_expression.pdf", width = 10, height = 5)
plot_grid(lve_raw, lwe_tpm, ncol = 2, labels = c("Raw", "TPM"), label_size = 14, scale = 0.9, label_x = 0.5,
          label_fontface = "bold")
dev.off()

# do histogram of counts
distr_fnn <- gwc_FNN %>% select(-Length) %>% gather(sample, counts) %>% ggplot(aes(x=counts)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(0,1000)) + scale_y_continuous(limits = c(0,200)) + theme_minimal() +
  # rename y and x axis
    labs(x = "Counts", y = "Frequency")

distr_fnv <- gwc_FNV %>% select(-Length) %>% gather(sample, counts) %>% ggplot(aes(x=counts)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(0,1000)) + scale_y_continuous(limits = c(0,200)) + theme_minimal() +
  # rename y and x axis
    labs(x = "Counts", y = "Frequency")

pdf("analysis/histograms_raw_counts.pdf", width = 10, height = 5)
plot_grid(distr_fnn, distr_fnv, ncol = 2, labels = c("FNN", "FNV"), label_size = 14, scale = 0.9, label_x = 0.5,
          label_fontface = "bold")
dev.off()

# remove H2O_16h_B fom gwc_fnv
#gwc_FNV <- gwc_FNV[,!grepl("H2O_16h_B", colnames(gwc_FNV))]

# get test conditions usin gsub to extract the sample names e.g. from
# "...data.rna_align.ID_007569_Fnn23_PNA79_16_h_B.fq.gz.bam" to "PNA79_16_h"
test_FNN <- as.factor(gsub("_[ABC]$", "",colnames(gwc_FNN)[-1]))
test_FNN
test_FNV <- as.factor(gsub("_[ABC]$", "",colnames(gwc_FNV)[-1]))
test_FNV

# save gwc_FNV and gwc_FNN as csv. preserve rownames
# change colnames (add Fnv_ to start of each)
colnames(gwc_FNN) <- c("Length", paste0("Fnn23_", colnames(gwc_FNN)[-1]))
colnames(gwc_FNV) <- c("Length", paste0("Fnv_", colnames(gwc_FNV)[-1]))


write.csv(gwc_FNN, "./data/GEO_SUBM_2024_12/raw_counts_FNN.csv", row.names = T)
write.csv(gwc_FNV, "./data/GEO_SUBM_2024_12/raw_counts_FNV.csv", row.names = T)


# create DGEList objects
y_FNN <- DGEList(counts = gwc_FNN[-1], group = test_FNN, genes = gwc_FNN[,1,drop=FALSE])
y_FNV <- DGEList(counts = gwc_FNV[-1], group = test_FNV, genes = gwc_FNV[,1,drop=FALSE])

# filter out genes with low counts with filterByExpr
keep_FNN <- filterByExpr(y_FNN)
table(keep_FNN)
y_FNN <- y_FNN[keep_FNN,,keep.lib.sizes=FALSE]
keep_FNV <- filterByExpr(y_FNV)
table(keep_FNV)
y_FNV <- y_FNV[keep_FNV,,keep.lib.sizes=FALSE]

# create a design matrix that can be used for both strains
design_matrix <- model.matrix(~0+test_FNN)
colnames(design_matrix) <- levels(test_FNN)
rownames(design_matrix) <- colnames(y_FNN)
design_matrix

# create design matrix for FNV
design_matrix_FNV <- model.matrix(~0+test_FNV)
colnames(design_matrix_FNV) <- levels(test_FNV)
rownames(design_matrix_FNV) <- colnames(y_FNV)
design_matrix_FNV

# do tmm normalization
y_FNN <- calcNormFactors(y_FNN, method = "TMM")
y_FNV <- calcNormFactors(y_FNV, method = "TMM")

# estimate dispersion
y_FNN <- estimateDisp(y_FNN, design_matrix, robust = TRUE)
y_FNV <- estimateDisp(y_FNV, design_matrix_FNV, robust = TRUE)


# create PCA and RLE plots in a nicer way:

# get logcpms for both organisms
logCPMs <- list(FNN = cpm(y_FNN, log = TRUE, prior.count = 2),
                FNV = cpm(y_FNV, log = TRUE, prior.count = 2))

# make theme for plots:
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA, size = 1.5),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             title=element_text(colour="black", size=23),
             axis.text=element_text(colour="black", size=18),axis.ticks=element_line(colour="black"),
             axis.title=element_text(colour="black", size=21),
             plot.margin=unit(c(1,1,1,1),"line"),
              # make legend upper left
              legend.position = "top",
              # decrease legend font size
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 15),
             plot.title = element_text(size = 23, face = "bold", hjust=0.5),
              # make line of axis thicker
                axis.line = element_line(size = 1)
              )

# define colors
mycols <- c("black", "red", "grey")


# creaate PCA plots and save them as pdf:
names(logCPMs) <- c("FNN23", "FNV")
tests <- list(test_FNN, test_FNV)
ggs_pca <- lapply(1:2, function(x) {
  i <- names(logCPMs)[x]
  pca <- prcomp(t(logCPMs[[i]]))
  pca_df <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0( colnames(pca_df), " (", paste( as.character(percentage), "%", ")", sep="") )
  pca_df$group <-tests[[x]]
  pca_df$treatment <- gsub("_.+", "", pca_df$group)
  pca_df$time <- factor(gsub(".+_", "", pca_df$group), levels = c("30m", "16h"))
  p<-ggplot(pca_df,aes(x=PC1,y=PC2,label=rownames(pca_df),group=treatment, fill=treatment, shape=time,
                       colour=treatment))
  p<-p+geom_point(size=10, alpha=0.7,aes(colour = treatment, shape = time))+
    theme +
    xlab(percentage[1]) +
    ylab(percentage[2])+
    scale_fill_manual(values = mycols) +
    scale_color_manual(values = mycols) +
    scale_shape_manual(values = c(22, 24)) +
    ggtitle(i) +
    # make title bold and italic
    theme(plot.title = element_text(face = "bold", size = 40))
  p

  # gerate file name for pdf. substitute all . and spaces with _. omit __
  file_name <- paste0("PCA_", gsub("\\.+", "_", gsub(" +", "_", i)))
  file_name <- gsub("_+", "_", file_name)
  # save plot as pdf
  pdf(paste0("analysis/", file_name), width = 11, height = 10)
  print(p)
  dev.off()
  p
})

# do same for PC 1 vs 3 and 2 vs 3
ggs_pca_1_vs_3 <- lapply(1:2, function(x) {
  i <- names(logCPMs)[x]
  pca <- prcomp(t(logCPMs[[i]]))
  pca_df <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0( colnames(pca_df), " (", paste( as.character(percentage), "%", ")", sep="") )
  pca_df$group <-tests[[x]]
  pca_df$treatment <- gsub("_.+", "", pca_df$group)
  pca_df$time <- factor(gsub(".+_", "", pca_df$group), levels = c("30m", "16h"))
  p<-ggplot(pca_df,aes(x=PC1,y=PC3,label=rownames(pca_df),group=treatment, fill=treatment, shape=time,
                       colour=treatment))
  p<-p+geom_point(size=10, alpha=0.7,aes(colour = treatment, shape = time))+
    theme +
    xlab(percentage[1]) +
    ylab(percentage[3])+
    scale_fill_manual(values = mycols) +
    scale_color_manual(values = mycols) +
    scale_shape_manual(values = c(22, 24)) +
    ggtitle(i) +
    # make title bold and italic
    theme(plot.title = element_text(face = "bold", size = 40))
  p

  # gerate file name for pdf. substitute all . and spaces with _. omit __
  file_name <- paste0("PCA_", gsub("\\.+", "_", gsub(" +", "_", i)))
  file_name <- gsub("_+", "_", file_name)
  # save plot as pdf
  pdf(paste0("analysis/", file_name, "_PC1_vs_PC3"), width = 11, height = 10)
  print(p)
  dev.off()
  p
})

ggs_pca_2_vs_3 <- lapply(1:2, function(x) {
  i <- names(logCPMs)[x]
  pca <- prcomp(t(logCPMs[[i]]))
  pca_df <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0( colnames(pca_df), " (", paste( as.character(percentage), "%", ")", sep="") )
  pca_df$group <-tests[[x]]
  pca_df$treatment <- gsub("_.+", "", pca_df$group)
  pca_df$time <- factor(gsub(".+_", "", pca_df$group), levels = c("30m", "16h"))
  p<-ggplot(pca_df,aes(x=PC2,y=PC3,label=rownames(pca_df),group=treatment, fill=treatment, shape=time,
                       colour=treatment))
  p<-p+geom_point(size=10, alpha=0.7,aes(colour = treatment, shape = time))+
    theme +
    xlab(percentage[2]) +
    ylab(percentage[3])+
    scale_fill_manual(values = mycols) +
    scale_color_manual(values = mycols) +
    scale_shape_manual(values = c(22, 24)) +
    ggtitle(i) +
    # make title bold and italic
    theme(plot.title = element_text(face = "bold", size = 40))
  p

  # gerate file name for pdf. substitute all . and spaces with _. omit __
  file_name <- paste0("PCA_", gsub("\\.+", "_", gsub(" +", "_", i)))
  file_name <- gsub("_+", "_", file_name)
  # save plot as pdf
  pdf(paste0("analysis/", file_name, "_PC2_vs_PC3"), width = 11, height = 10)
  print(p)
  dev.off()
  p
})

# use cowplot to combine the plots
pdf("analysis/PCA_both_strains_raw.pdf", width = 20, height = 10)
plot_grid(ggs_pca[[1]], ggs_pca[[2]], ncol = 2,  scale = 0.95)
dev.off()

svg("analysis/PCA_both_strains_raw.svg", width = 20, height = 10)
plot_grid(ggs_pca[[1]], ggs_pca[[2]], ncol = 2,  scale = 0.95)
dev.off()

pdf("analysis/PCA_both_strains_raw_PC1_till_3.pdf", width = 20, height = 30)
plot_grid(ggs_pca[[1]], ggs_pca[[2]], ggs_pca_1_vs_3[[1]], ggs_pca_1_vs_3[[2]], ggs_pca_2_vs_3[[1]], ggs_pca_2_vs_3[[2]], ncol = 2,
          scale = 0.95, labels = c("A", "", "B", "", "C", ""), label_size = 40)
dev.off()

svg("analysis/PCA_both_strains_raw_PC1_till_3.svg", width = 20, height = 30)
plot_grid(ggs_pca[[1]], ggs_pca[[2]], ggs_pca_1_vs_3[[1]], ggs_pca_1_vs_3[[2]], ggs_pca_2_vs_3[[1]], ggs_pca_2_vs_3[[2]], ncol = 2,
          scale = 0.95, labels = c("A", "", "B", "", "C", ""), label_size = 40)
dev.off()

# save obnly 1 vs 3 and 2 vs 3 plots as pdf and tiff
pdf("analysis/PCA_both_strains_raw_PC1_vs_PC3_2_v3.pdf", width = 20, height = 20)
plot_grid(ggs_pca_1_vs_3[[1]], ggs_pca_1_vs_3[[2]], ggs_pca_2_vs_3[[1]], ggs_pca_2_vs_3[[2]], ncol = 2,
          scale = 0.95, labels = c("A", "", "B", ""), label_size = 60)
dev.off()

# save as tiff
png("analysis/PCA_both_strains_raw_PC1_vs_PC3_2_v3.png", width = 20, height = 20, units = "in", res = 500)
plot_grid(ggs_pca_1_vs_3[[1]], ggs_pca_1_vs_3[[2]], ggs_pca_2_vs_3[[1]], ggs_pca_2_vs_3[[2]], ncol = 2,
          scale = 0.95, labels = c("A", "", "B", ""), label_size = 60)
dev.off()

# save as svg
svg("analysis/PCA_both_strains_raw_PC1_vs_PC3_2_v3.svg", width = 20, height = 20)
plot_grid(ggs_pca_1_vs_3[[1]], ggs_pca_1_vs_3[[2]], ggs_pca_2_vs_3[[1]], ggs_pca_2_vs_3[[2]], ncol = 2,
          scale = 0.95, labels = c("A", "", "B", ""), label_size = 60)
dev.off()

xddfewfwe
# # do ruvseq normalization, RUVs analysis:
# matrix_y_fnv <- as.matrix(sapply(as.data.frame(cpm(y_FNV)), as.integer))
# rownames(matrix_y_fnv) <- rownames(y_FNV$counts)
# set_FNV <- newSeqExpressionSet(matrix_y_fnv,
#                                phenoData = data.frame(test_FNV, row.names = colnames(y_FNV$counts)))
#
# #set_FNN <- betweenLaneNormalization(set_FNN, which = "median")
# set_FNV <- betweenLaneNormalization(set_FNV, which = "median")
#
# plotRLE(set_FNV, outline = FALSE, col = mycols[test_FNV])
# plotPCA(set_FNV, col = mycols[test_FNV], cex = 1.2)
#
# # try ruvs for 5 ks:
# for (k in 1:5) {
#   par(mfrow = c(1, 2))
#   set_RUVs <- RUVs(set_FNV, k = k, rownames(y_FNV), makeGroups(test_FNV))
#   plotRLE(set_RUVs, outline = FALSE, ylim = c(-1, 1), col = mycols[test_FNV],
#           main = paste(c("k = ", k)))
#   plotPCA(set_RUVs, col = mycols[test_FNV], cex = 1.2, main = paste(c("k = ", k)))
# }
#
# # try ruvg with 5 ks. first determine genes that do not change much
# # between conditions. do a
#
# # k=3 looks decent. change y_FNN to it
# set_RUVs <- RUVs(set_FNV, k = 2, rownames(y_FNV), makeGroups(test_FNV))
#
# # create PCA and RLE plots in a nicer way
# logCPMs <- list(FNN = cpm(y_FNN, log = TRUE, prior.count = 2),
#                 FNV = log(normCounts(set_RUVs) + 1))
#
# names(logCPMs) <- c("F. nucleatum", "F. nucleatum vincentii")
# ggs_pca <- lapply(names(logCPMs), function(i) {
#   pca <- prcomp(t(logCPMs[[i]]))
#   pca_df <- as.data.frame(pca$x)
#   percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
#   percentage <- paste( colnames(pca_df), "(", paste( as.character(percentage), "%", ")", sep="") )
#   pca_df$group <-test_FNV
#   p<-ggplot(pca_df,aes(x=PC1,y=PC2,group=group,label=rownames(pca_df), colour=group))
#   p<-p+geom_point(size=3)+ scale_shape_identity()+
#     geom_text_repel(size=7, min.segment.length = 0, seed = 1, box.padding = 0.6, max.overlaps = 20)+
#     theme + xlab(percentage[1]) +
#     ylab(percentage[2])+ scale_color_manual(values = mycols) +
#     ggtitle(bquote(bold(bolditalic(.(as.character(i)))~"after normalization")))
#
#   # gerate file name for pdf. substitute all . and spaces with _. omit __
#   file_name <- paste0("PCA_",gsub("\\.+", "_", gsub(" +", "_", i)))
#   file_name <- gsub("_+", "_", file_name)
#   # save plot as pdf
#   pdf(paste0("analysis/", file_name), width = 11, height = 10)
#   print(p)
#   dev.off()
#   p
# })
#
# # use cowplot to combine the plots
# pdf("analysis/PCA_both_strains_norm.pdf", width = 20, height = 10)
# plot_grid(ggs_pca[[1]], ggs_pca[[2]], ncol = 2, labels = c("A", "B"), label_size = 25, scale = 0.95)
# dev.off()
#
#
# # change y_FNV to set_RUVs
# y_FNV <- DGEList(counts = counts(set_RUVs), group = test_FNV)
# options(digits = 3)
# design_matrix_FNV <- model.matrix(~0+test_FNV+W_1+W_2, data=pData(set_RUVs))
# colnames(design_matrix_FNV) <- c(levels(test_FNV), "W_1", "W_2")
# rownames(design_matrix_FNV) <- colnames(y_FNV$counts)
#
# y_FNV <- calcNormFactors(y_FNV)
# y_FNV <- estimateDisp(y_FNV, design_matrix_FNV, robust = TRUE)

# do DE analysis for all organisms. start by making contrasts betweeen test conditions that
# are of interest. First focus only on PNA, so compare PNA79_16h to H2O_16h, PNA79_16h to PNAscr_16h, and
# PNAscr_16h to H2O_16h. For both strains.
contrast_FNN <- makeContrasts(PNA79_16h_vs_H2O_16h = PNA79_16h - H2O_16h,
                              PNA79_16h_vs_PNAscr_16h = PNA79_16h - PNAscr_16h,
                              PNAscr_16h_vs_H2O_16h = PNAscr_16h - H2O_16h,
                              levels = design_matrix)

contrast_FNV <- makeContrasts(PNA79_16h_vs_H2O_16h = PNA79_16h - H2O_16h,
                                PNA79_16h_vs_PNAscr_16h = PNA79_16h - PNAscr_16h,
                                PNAscr_16h_vs_H2O_16h = PNAscr_16h - H2O_16h,
                                levels = design_matrix_FNV)

# do DE analysis
fit_FNN <- glmQLFit(y_FNN, design_matrix, robust = TRUE)
fit_FNV <- glmQLFit(y_FNV, design_matrix_FNV, robust = TRUE)

# get results
res_FNN <- list(# PNA79_16h_vs_H2O_16h
                PNA79_16h_vs_H2O_16h = glmQLFTest(fit_FNN, contrast = contrast_FNN[,1]),
                # PNA79_16h_vs_PNAscr_16h
                PNA79_16h_vs_PNAscr_16h = glmQLFTest(fit_FNN, contrast = contrast_FNN[,2]),
                # PNAscr_16h_vs_H2O_16h
                PNAscr_16h_vs_H2O_16h = glmQLFTest(fit_FNN, contrast = contrast_FNN[,3]))

res_FNV <- list(# PNA79_16h_vs_H2O_16h
                PNA79_16h_vs_H2O_16h = glmQLFTest(fit_FNV, contrast = contrast_FNV[,1]),
                # PNA79_16h_vs_PNAscr_16h
                PNA79_16h_vs_PNAscr_16h = glmQLFTest(fit_FNV, contrast = contrast_FNV[,2]),
                # PNAscr_16h_vs_H2O_16h
                PNAscr_16h_vs_H2O_16h = glmQLFTest(fit_FNV, contrast = contrast_FNV[,3]))

#

# apply FDR to all tables
res_FNN <- lapply(res_FNN, function(r) {
  r$table$FDR <- p.adjust(r$table$PValue, method = "fdr")
  r
})

res_FNV <- lapply(res_FNV, function(r) {
  r$table$FDR <- p.adjust(r$table$PValue, method = "fdr")
  r
})


# create volcano plot. use own function:
do_volcano <- function(restab, targetgene = NULL, pointsize = 2, x_limit = F, y_limit = F, show_sig = F,
                       alpha = 0.05, color_sig = T, marked_gene_names = NULL, marked_gene_title = NULL,
                       minlogfc = 1, title = "Volcano", marked_genes = NULL, add_labels = T, gene_names = NULL,
                       color_threshold_lines = "black",
                       cols = c("target" = "darkorange", "marked_genes" = "darkred", "up" = "darkorange",
                                "down" = "steelblue", "sRNA" = "darkred",
                                "other" = "darkgrey")) {
  # change rownames of restab to gene_names if gen_names is not NULL:
    if (!is.null(gene_names)) {
        rownames(restab) <- gene_names
    }
  #rownames(restab) <- gsub("^([^S][^A].+)" , "italic('\\1')" , rownames(restab))


  g <- ggplot(restab) +
    geom_point(
      data = restab,
      aes(x = logFC, y = -log10(FDR), fill = "other"),shape = 21, color="darkgrey",
      cex = pointsize, alpha = 0.4
    ) +
    theme_bw() +
    # change theme to standard black&wite.
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
          panel.grid.major.x = element_blank(), #element_line(colour="lightgrey", size=0.3),
          panel.grid.major.y = element_blank(), #element_line(colour="lightgrey", size=0.3),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
    ggtitle(title) +
    xlab(expression("log"[2] * " fold change")) +
    ylab(expression("- log"[10] * " P-value (FDR)"))

  if (x_limit) {
    g <- g + scale_x_continuous(expand = c(0, 0), breaks = seq(-6, 6, 2), limits = c(-x_limit, x_limit))
  }
  if (y_limit) {
    g <- g + scale_y_continuous(expand = c(0, 0), breaks = seq(0, 40, 5), limits = c(0, y_limit))
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
    marked_genes <- gsub("^([^S][^A].+)" , "italic('\\1')" , marked_genes)
    g <- g + geom_point(
      data = restab[marked_genes,],
      aes(x = logFC, y = -log10(FDR), fill = "marked_genes"),
      cex = pointsize , shape=21,alpha = 0.5)
  }

  # show the sign. genes:
  # show the sigficantest genes:
  if (show_sig) {
    range01 <- function(x) { (x - min(x)) / (max(x) - min(x)) }
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


    g_labels <- c(targetgene,"C4N14_09830", #marked_genes,
                  rownames(top_peaks)[!grepl("SA101588_", rownames(top_peaks))])
  }

  # add labels:
  if (add_labels == T) {
    g_labels <- unique(g_labels)
    rup <- restab[g_labels,][restab[g_labels,]$logFC > 0,]
    rownames(rup) <- gsub("(HMPREF0946_)", "", rownames(rup))
    rownames(rup) <- gsub("4\\.5S_RNA", "4.5_srRNA", rownames(rup))
    rownames(rup) <- gsub("cds-", "", rownames(rup))
#    rownames(rup) <- gsub("C4N14_", "", rownames(rup))
    print("Done")
    print(rownames(rup))
    print(rup)
    g <- g + geom_text_repel(
      data = rup,
      aes(x = logFC, y = -log10(FDR), label = rownames(rup)),
      nudge_x = x_limit - rup$logFC,
      min.segment.length = 0,
      direction = "y",
      hjust = 1,
      size = 3,
      segment.color = "black",
      segment.alpha = 0.5,
      parse = F
    )
    # add labels for downregulated genes:
    rdown <- restab[g_labels,][restab[g_labels,]$logFC < 0,]
    rownames(rdown) <- gsub("(HMPREF0946_)", "", rownames(rdown))
    rownames(rdown) <- gsub("4\\.5S_RNA", "4.5_sRNA", rownames(rdown))
    rownames(rdown) <- gsub("cds-", "", rownames(rdown))
    #rownames(rdown) <- gsub("C4N14_", "", rownames(rdown))
    print("Done")
    g <- g + geom_text_repel(
      data = rdown,
      aes(x = logFC, y = -log10(FDR), label = rownames(rdown)),
      nudge_x       = -x_limit - rdown$logFC,
      min.segment.length = 0,
      direction = "y",
      hjust = 0,
      size = 3,
      segment.color = "black",
      segment.alpha = 0.5,
      parse = F)
  }
  g + scale_color_manual(values = cols) + scale_fill_manual(values = cols)

}

volc_list <- list()
# create volcano plots for both strains. start with fnn:
for (i in names(res_FNN)){
  restab <- res_FNN[[i]]$table
  ttle <- gsub("_", " ", i)
  ttle <- gsub("vs", "vs.", ttle)
  # make volcano plot
  volc <- do_volcano(restab, pointsize = 2, x_limit = 8, y_limit = 13, show_sig = T,
                     alpha = 0.01, color_sig = T, title = ttle,
                     minlogfc = 1.5, add_labels = T,
                     color_threshold_lines = "black") +
    geom_point(data = restab[2077:2104,], aes(x = logFC, y = -log10(FDR), fill = "sRNA"),
               cex = 2, shape = 21, alpha = 0.6)
  volc_list <- c(volc_list, list(volc))
  volc
  restab <- restab[order(restab$PValue),]
  # save the results table as excel file. make first line thicker
  write.xlsx(restab, paste0("./analysis/DE_raw_data/", i, "_DE_result_FNN_16h.xlsx"), row.names = T, col.names = T)
}

volc_grid <- plot_grid(plotlist = volca_list[c(1,3,2)], ncol = 3, nrow = 1,
                       labels = c("A", "B", "C"), label_size = 20, scale=0.95)

# save the grid as pdf
ggsave("analysis/volcano_plots_FNN_16h.pdf", volc_grid, width = 35, height = 12, units = "cm")

# do same for FNV:
volc_list <- list()
for (i in names(res_FNV)){
  restab <- res_FNV[[i]]$table
  ttle <- gsub("_", " ", i)
  ttle <- gsub("vs", "vs.", ttle)
  # make volcano plot
  volc <- do_volcano(restab, pointsize = 2, x_limit = 8, y_limit = 13, show_sig = T,
                     alpha = 0.01, color_sig = T, title = ttle,
                     minlogfc = 1.5, add_labels = T,
                     color_threshold_lines = "black") +
    geom_point(data = restab[2006:2038,], aes(x = logFC, y = -log10(FDR), fill = "sRNA"),
                cex = 2, shape = 21, alpha = 0.6)

  volc_list <- c(volc_list, list(volc))

  # sort restab by p-value starting with the smallest
  restab <- restab[order(restab$PValue),]
  # save the results table as excel file. make first line thicker

  write.xlsx(restab, paste0("./analysis/DE_raw_data/", i, "_DE_result_FNV_16h.xlsx"), row.names = T, col.names = T)
}

volc_grid <- plot_grid(plotlist = volc_list[c(1,3,2)], ncol = 3, nrow = 1,
                       labels = c("A", "B", "C"), label_size = 20, scale=0.95)

# save the grid as pdf
ggsave("analysis/volcano_plots_FNV_16h.pdf", volc_grid, width = 35, height = 12, units = "cm")


# get all significantly regulated genes (log2FC > 1.5 and FDR < 0.01)
sig_genes_FNN <- lapply(res_FNN, function(x) dim(x$table[x$table$FDR < 0.01 & abs(x$table$logFC) > 1.5,]))


# save tables from res_fnn and res_fnv as excel files. first, use rtracklayer to get
# first get gff file names:
gff_FNN <- "./data/reference_sequences/ATCC23726_updated2023.gff"
gff_FNV <- "./data/reference_sequences/FNV3_1_36A.gff"

# get sgff tags that are interesting, for fnn:
gff_tags_FNN <- c("id", "Name", "product")
gff_tags_FNV <- c("id", "locus_tag", "product")
# get cols that are interesting for both gff files:
gff_cols <- c("seqid","type","start","end", "strand")

# get filter, for only CDS and ncRNA
gff_filter <- list(type=c("CDS","ncRNA"))

# load gff_FNN and gff_FNV
attributes_FNN <- readGFF(gff_FNN, filter = gff_filter, tags = gff_tags_FNN, columns = gff_cols) %>% as_tibble() %>%
  # change product and description to "sRNA" if the type is ncRNA
    mutate(product = ifelse(type == "ncRNA", "sRNA", product))

# load gff_FNV
attributes_FNV <- readGFF(gff_FNV, filter = gff_filter, tags = gff_tags_FNV, columns = gff_cols) %>% as_tibble() %>%
  # change product and description to "sRNA" if the type is ncRNA
  mutate(product = ifelse(type == "ncRNA", "sRNA", product),
  # change locus tag of ncrna to id
    locus_tag = ifelse(type == "ncRNA", id, locus_tag))



# get table of res_FNN fot pna97 vs h2o
res_FNN_PNA79_16h_vs_H2O_16h <- res_FNN$PNA79_16h_vs_H2O_16h$table %>%
  # change rownames to id column
    rownames_to_column(var = "id") %>% as_tibble()
# same for pna79 vs pnascr
res_FNN_PNA79_16h_vs_PNAscr_16h <- res_FNN$PNA79_16h_vs_PNAscr_16h$table %>%
  rownames_to_column(var = "id") %>% as_tibble()


# combine attributes_FNN with res_FNN_PNA79_16h_vs_H2O_16h
res_FNN_PNA79_16h_vs_H2O_16h_excel_ready <- inner_join( attributes_FNN, res_FNN_PNA79_16h_vs_H2O_16h, by = "id") %>%
  # sort by P-value
    arrange(PValue) %>%
    # select only id, type, product, logCPM, logFC, PValue, FDR
    select(id, type, product, logCPM, logFC, PValue, FDR) %>%
    # remove duplicates
    distinct()

# same for res_FNN_PNA79_16h_vs_PNAscr_16h
res_FNN_PNA79_16h_vs_PNAscr_16h_excel_ready <- inner_join( attributes_FNN, res_FNN_PNA79_16h_vs_PNAscr_16h, by = "id") %>%
  arrange(PValue) %>%
    # select only id, type, product, logCPM, logFC, PValue, FDR
        select(id, type, product, logCPM, logFC, PValue, FDR) %>%
    # remove duplicates
        distinct()


# ok now get the same for FNV
res_FNV_PNA79_16h_vs_H2O_16h <- res_FNV$PNA79_16h_vs_H2O_16h$table %>%
  rownames_to_column(var = "id") %>% as_tibble()

res_FNV_PNA79_16h_vs_PNAscr_16h <- res_FNV$PNA79_16h_vs_PNAscr_16h$table %>%
    rownames_to_column(var = "id") %>% as_tibble()

# combine attributes_FNV with res_FNV_PNA79_16h_vs_H2O_16h
res_FNV_PNA79_16h_vs_H2O_16h_excel_ready <- inner_join( attributes_FNV, res_FNV_PNA79_16h_vs_H2O_16h, by = "id") %>%
  arrange(PValue) %>%
    # select only id, type, product, logCPM, logFC, PValue, FDR
        select(id, type, product, logCPM, logFC, PValue, FDR) %>%
    # remove duplicates
        distinct()

# same for res_FNV_PNA79_16h_vs_PNAscr_16h
res_FNV_PNA79_16h_vs_PNAscr_16h_excel_ready <- inner_join( attributes_FNV, res_FNV_PNA79_16h_vs_PNAscr_16h, by = "id") %>%
  arrange(PValue) %>%
  # select only id, type, product, logCPM, logFC, PValue, FDR
    select(id, type, product, logCPM, logFC, PValue, FDR) %>%
  # remove duplicates
    distinct()


# ok now create an excel file, with four sheets, one for each comparison. first show FNN23 FUS79 vs PNAscr 30min,
# then FNN23 FUS79 vs PNAscr 16h, then FNV FUS79 vs PNAscr 30min, then FNV FUS79 vs PNAscr 16h

library(baizer)

write_excel(res_FNN_PNA79_16h_vs_H2O_16h_excel_ready,"analysis/Supplementary_Data_SD1.xlsx",
            sheetname = "FNN23 FUS79 vs H2O 16h")


## Now do kegg pathway analysis. First get the gene names of the DE genes:

# get kegg and check all organisms with pathways
org <- keggList("organism")

# filter for all with "Fuso" in species column
org <- org[grepl("Fuso", org[,3]),]


# get kegg list of fusobacterium nucleatum
link_kegg_FNN <- keggLink("pathway", "fnu")

# load mappings for entries in link_kegg_FNN
mappings_FNN <- read_tsv("data/reference_sequences/FNN_mapping.proteinortho.tsv", col_names = T)

# get the gene names of the DE genes:
names(link_kegg_FNN) <- gsub("fnu:", "", names(link_kegg_FNN))

# get the gene names of the portho mappings:
link_kegg_FNN <- link_kegg_FNN[names(link_kegg_FNN) %in% mappings_FNN$CDS_FNN_25586]
# rename names to the portho mappings:
names(link_kegg_FNN) <- mappings_FNN$CDS_FNN_falk[match(names(link_kegg_FNN), mappings_FNN$CDS_FNN_25586)]

# add rpoE regulon
rpoe_genes <- read.table("./data/rpoe_regulon.tsv", header = F)$V1
# make rpoe_genes the names of characters which all have the value "rpoE"
rpoe_gs <- rep("rpoE", length(rpoe_genes))
names(rpoe_gs) <- rpoe_genes
link_kegg_FNN <- c(link_kegg_FNN, rpoe_gs)


# now load for fusobacterium vincentii
link_kegg_FNV <- keggLink("pathway", "fnc")

# load mappings for entries in link_kegg_FNV
mappings_FNV <- read_tsv("data/reference_sequences/mappings_FNV.tsv", col_names = F)

# get the gene names of the DE genes:
names(link_kegg_FNV) <- gsub("fnc:", "", names(link_kegg_FNV))

# get the gene names of the portho mappings:
link_kegg_FNV <- link_kegg_FNV[names(link_kegg_FNV) %in% mappings_FNV$X2]

# rename names to the portho mappings:
names(link_kegg_FNV) <- mappings_FNV$X1[match(names(link_kegg_FNV), mappings_FNV$X2)]

# add rpoE regulon
mapings_fnn_fnv <- read.table("./data/reference_sequences/FNV_FNN_mapping.proteinortho.tsv", header = T)
rpoe_genes_fnv <- mapings_fnn_fnv$FNV_A2_CDS[mapings_fnn_fnv$CDS_FNN_Falk %in% rpoe_genes]
# make rpoe_genes the names of characters which all have the value "rpoE"
rpoe_gs_fnv <- rep("rpoE", length(rpoe_genes_fnv))
names(rpoe_gs_fnv) <- rpoe_genes_fnv
link_kegg_FNV <- c(link_kegg_FNV, rpoe_gs_fnv)


list_kegg_FNN <- keggList("pathway", "fnu")
list_kegg_FNV <- keggList("pathway", "fnc")

# add rpoE regulon to list_kegg_FNN
list_kegg_FNN <- c(list_kegg_FNN, "rpoE" = "rpoE regulon")
list_kegg_FNV <- c(list_kegg_FNV, "rpoE" = "rpoE regulon")



kegg_pw_ids_FNN <- names(list_kegg_FNN)
kegg_pw_ids_FNV <- names(list_kegg_FNV)


# filter only lts that are in DE expression results
link_kegg_FNN <- link_kegg_FNN[names(link_kegg_FNN) %in% gsub("cds-", "",
                                                              rownames(res_FNN$PNA79_16h_vs_H2O_16h$table))]
link_kegg_FNV <- link_kegg_FNV[names(link_kegg_FNV) %in% rownames(res_FNV$PNA79_16h_vs_H2O_16h$table)]

# remove path: from link_kegg:
link_kegg_FNN <- gsub("path:", "", link_kegg_FNN)
link_kegg_FNV <- gsub("path:", "", link_kegg_FNV)

# index
idx_kegg_FNN <- sapply(kegg_pw_ids_FNN, function(x){
  y <- unique(names(link_kegg_FNN[link_kegg_FNN == x])) # choose all genes, except duplucates
  y <- gsub("(.*)", "cds-\\1", y)
})


idx_kegg_FNV <- sapply(kegg_pw_ids_FNV, function(x){
  y <- unique(names(link_kegg_FNV[link_kegg_FNV == x])) # choose all genes, except duplucates
})

# OK, now we need the length of contrasts
l <- length(colnames(contrast_FNN))

# now run FRY gene set test for all contrasts:
kegg_fry_FNN <- lapply(1:l, function(x) fry(y_FNN, idx_kegg_FNN, design_matrix, contrast_FNN[,x]))
names(kegg_fry_FNN) <- colnames(contrast_FNN)
kegg_fry_FNV <- lapply(1:l, function(x) fry(y_FNV, idx_kegg_FNV, design_matrix_FNV, contrast_FNV[,x]))
names(kegg_fry_FNV) <- colnames(contrast_FNV)

# reorder kegg_fry_FNN
kegg_fry_FNN <- kegg_fry_FNN[c(1,3,2)]
kegg_fry_FNV <- kegg_fry_FNV[c(1,3,2)]

# add KEGG terms
for (fryres in names(kegg_fry_FNN)){
  kegg_fry_FNN[[fryres]][["TERM"]] <- ifelse(grepl("fnu", rownames(kegg_fry_FNN[[fryres]])),
                                             list_kegg_FNN[rownames(kegg_fry_FNN[[fryres]])],
                                             rownames(kegg_fry_FNN[[fryres]]))
  kegg_fry_FNN[[fryres]][["TERM"]] <- gsub(" - Fusobacterium nucleatum subsp. nucleatum ATCC 25586", "",
                                           kegg_fry_FNN[[fryres]][["TERM"]])
  # add a column containing all locus tags
  kegg_fry_FNN[[fryres]][["GENES"]] <- as.character(unlist(sapply(idx_kegg_FNN[rownames(kegg_fry_FNN[[fryres]])], function(x) paste(x,collapse=";"))))
  write.csv(kegg_fry_FNN[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNN_16h.csv"))
  # save as excel file
  write.xlsx(kegg_fry_FNN[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNN_16h.xlsx"), row.names = T,
             col.names = T)
}

for (fryres in names(kegg_fry_FNV)){
  kegg_fry_FNV[[fryres]][["TERM"]] <- ifelse(grepl("fnc", rownames(kegg_fry_FNV[[fryres]])),
                                             list_kegg_FNV[rownames(kegg_fry_FNV[[fryres]])],
                                             rownames(kegg_fry_FNV[[fryres]]))
  kegg_fry_FNV[[fryres]][["TERM"]] <- gsub(" - Fusobacterium vincentii 3_1_36A2", "",
                                           kegg_fry_FNV[[fryres]][["TERM"]])
  # add a column containing all locus tags
  kegg_fry_FNV[[fryres]][["GENES"]] <- as.character(unlist(sapply(idx_kegg_FNV[rownames(kegg_fry_FNV[[fryres]])], function(x) paste(x,collapse=";"))))
  write.csv(kegg_fry_FNV[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNV_16h.csv"))
    # save as excel file
    write.xlsx(kegg_fry_FNV[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNV_16h.xlsx"), row.names = T
                , col.names = T)
}

# filter for pws with FDR < 0.01 & more than 4 genes
kegg_frysig_FNN <- lapply(kegg_fry_FNN, function(x) x[x[["FDR"]]<0.01 & x[["NGenes"]]>4,])
kegg_frysig_FNV <- lapply(kegg_fry_FNV, function(x) x[x[["FDR"]]<0.1 & x[["NGenes"]]>4,])

# select siggos
kegg_siggos_FNN <- c()
kegg_siggos_FNV <- c()

for (i in names(kegg_frysig_FNN)) {
  print(i)
  print(dim(kegg_frysig_FNN[[i]]))
  print(kegg_frysig_FNN[[i]][,c(1,2,4,7)])
  kegg_siggos_FNN <- c(kegg_siggos_FNN, rownames(kegg_frysig_FNN[[i]][1:10,]))  # can be modified
}
# remove duplicates
kegg_siggos_FNN <- unique(kegg_siggos_FNN[!grepl("NA", kegg_siggos_FNN)])

for (i in names(kegg_frysig_FNV)) {
  print(i)
  print(dim(kegg_frysig_FNV[[i]]))
  print(kegg_frysig_FNV[[i]][,c(1,2,4,7)])
  kegg_siggos_FNV <- c(kegg_siggos_FNV, rownames(kegg_frysig_FNV[[i]][1:10,]))  # can be modified
}
# remove duplicates
kegg_siggos_FNV <- unique(kegg_siggos_FNV[!grepl("NA", kegg_siggos_FNV)])

# Create a heatmap-df  for KEGG:
idx_kegg_char_FNN <- lapply(idx_kegg_FNN, as.character)
idx_kegg_char_FNV <- lapply(idx_kegg_FNV, as.character)

idx_kegg_char_FNN$rpoE <- gsub("cds-FoxI", "FoxI", idx_kegg_char_FNN$rpoE)
idx_kegg_char_FNN$rpoE <- gsub("cds-FunR7", "FunR7", idx_kegg_char_FNN$rpoE)

# I create a dataframe with mean logFC values for each significant kegg-term:
hm_kegg_FNN <- t(as.data.frame(lapply(idx_kegg_char_FNN[kegg_siggos_FNN], function(x){
  sapply(names(res_FNN), function(y){
    mean(res_FNN[[y]]$table[unique(x),]$logFC)
  })
})))
hm_kegg_FNN <- as.data.frame(hm_kegg_FNN)

hm_kegg_FNV <- t(as.data.frame(lapply(idx_kegg_char_FNV[kegg_siggos_FNV], function(x){
  sapply(names(res_FNV), function(y){
    mean(res_FNV[[y]]$table[unique(x),]$logFC)
  })
})))
hm_kegg_FNV <- as.data.frame(hm_kegg_FNV)

# reorder columns
hm_kegg_FNN <- hm_kegg_FNN[,c(1,3,2)]
hm_kegg_FNV <- hm_kegg_FNV[,c(1,3,2)]

# make heatmap:
hm_kegg_FNN <- hm_kegg_FNN[order(hm_kegg_FNN[,1], decreasing = T),]
hm_kegg_FNV <- hm_kegg_FNV[order(hm_kegg_FNV[,1], decreasing = T),]

# add size of pws:
kegg_sizes_FNN <- sapply(idx_kegg_char_FNN[rownames(hm_kegg_FNN)], function(x) length(x))
kegg_sizes_FNV <- sapply(idx_kegg_char_FNV[rownames(hm_kegg_FNV)], function(x) length(x))

# add p values
pvals_FNN <- data.frame(sapply(names(kegg_fry_FNN),
                           function(x) kegg_fry_FNN[[x]][rownames(hm_kegg_FNN),"FDR"]),
                    row.names = rownames(hm_kegg_FNN))
pvals_FNV <- data.frame(sapply(names(kegg_fry_FNV),
                               function(x) kegg_fry_FNV[[x]][rownames(hm_kegg_FNV),"FDR"]),
                        row.names = rownames(hm_kegg_FNV))

#select only significant ones:
pvals_FNN <-sapply(pvals_FNN, function(x) ifelse(x<0.05, x <- "*", x<-"") )
pvals_FNV <-sapply(pvals_FNV, function(x) ifelse(x<0.05, x <- "*", x<-"") )

# add term
keggpws_FNN <- kegg_fry_FNN$PNA79_16h_vs_PNAscr_16h[rownames(hm_kegg_FNN),] [["TERM"]]
keggpws_FNV <- kegg_fry_FNV$PNA79_16h_vs_PNAscr_16h[rownames(hm_kegg_FNV),] [["TERM"]]

# add names to hm_kegg
rownames(hm_kegg_FNN) <- ifelse(!is.na(keggpws_FNN),keggpws_FNN, rownames(hm_kegg_FNN) )
rownames(hm_kegg_FNV) <- ifelse(!is.na(keggpws_FNV),keggpws_FNV, rownames(hm_kegg_FNV) )

# add coloring
col_fun <- colorRamp2(c(-1,0, 1), c("steelblue", "white", "darkred"))

colnames(hm_kegg_FNN) <- c("PNA79_vs_H20", "PNAscr_vs_H20", "PNA79_vs_PNAscr")

ht_vert_FNN <- Heatmap(hm_kegg_FNN, cluster_rows = F, cluster_columns = F,
               name = "GO-analysis", col = col_fun,
               show_heatmap_legend = F,
               #column_split = rep(c("PNA","scrambled","control"), each=4),
               row_title_side = "right", row_title_rot = 0,
               border = TRUE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.1s", pvals_FNN[i, j]), x, y)
               },
               column_names_gp = gpar(fontsize = 11),
               column_title_gp = gpar(fontsize = 15),
               row_names_gp = gpar(fontsize = 10),
               # add column title
               column_title = "Time point 16h FNN",
               row_title = NULL,
               width = unit(4, "cm"), height = unit(15, "cm"),

               right_annotation = rowAnnotation(genes = anno_barplot(kegg_sizes_FNN)))

ht_vert_FNN

lgd <- Legend(col_fun = col_fun, title = expression("mean log"[2]*" FC"), #direction = "horizontal",
             title_gp = gpar(fontsize = 20), labels = c("-1", " 0"," 1"), legend_height = unit(8, "cm"),
              grid_width = unit(1, "cm"),
              labels_gp = gpar(fontsize = 20),
             at = c(-1, 0, 1), border = "black",
             title_position = "leftcenter-rot")

svg("analysis/KEGG_heatmap_FNN_16h.svg", width = 10, height = 10)
ht_vert_FNN
draw(lgd, x = unit(2, "cm"), y = unit(15, "cm"))
dev.off()

pdf("analysis/KEGG_heatmap_FNN_16h.pdf", width = 10, height = 10)
ht_vert_FNN
draw(lgd, x = unit(2, "cm"), y = unit(15, "cm"))
dev.off()

colnames(hm_kegg_FNV) <- c("PNA79_vs_H20", "PNAscr_vs_H20", "PNA79_vs_PNAscr")
ht_vert_FNV <- Heatmap(hm_kegg_FNV, cluster_rows = F, cluster_columns = F,
                       name = "GO-analysis", col = col_fun,
                       show_heatmap_legend = F,
                       #column_split = rep(c("PNA","scrambled","control"), each=4),
                       row_title_side = "right", row_title_rot = 0,
                       border = TRUE,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.1s", pvals_FNV[i, j]), x, y)
                       },
                       column_names_gp = gpar(fontsize = 11),
                       column_title_gp = gpar(fontsize = 15),
                       row_names_gp = gpar(fontsize = 10),
                       # add column title
                       column_title = "Time point 16h FNV",
                       row_title = NULL,
                       width = unit(4, "cm"), height = unit(13, "cm"),

                       right_annotation = rowAnnotation(genes = anno_barplot(kegg_sizes_FNV)))

ht_vert_FNV

svg("analysis/KEGG_heatmap_FNV_16h.svg", width = 12, height = 10)
ht_vert_FNV
draw(lgd, x = unit(2, "cm"), y = unit(15, "cm"))
dev.off()

pdf("analysis/KEGG_heatmap_FNV_16h.pdf", width = 12, height = 10)
ht_vert_FNV
draw(lgd, x = unit(2, "cm"), y = unit(15, "cm"))
dev.off()





# Okay, now do the whole analysis for time point 30 min:
# create new contrasts
contrasts_FNN_30 <- makeContrasts(PNA79_30m_vs_H2O_30m = PNA79_30m - H2O_30m,
                                  PNA79_30m_vs_PNAscr_30m = PNA79_30m - PNAscr_30m,
                                  PNAscr_30m_vs_H2O_30m = PNAscr_30m - H2O_30m,
                                  levels = design_matrix)

contrasts_FNV_30 <- makeContrasts(PNA79_30m_vs_H2O_30m = PNA79_30m - H2O_30m,
                                    PNA79_30m_vs_PNAscr_30m = PNA79_30m - PNAscr_30m,
                                    PNAscr_30m_vs_H2O_30m = PNAscr_30m - H2O_30m,
                                    levels = design_matrix_FNV)


# do DE analysis
fit_FNN_30 <- glmQLFit(y_FNN, design_matrix, robust = TRUE)
fit_FNV_30 <- glmQLFit(y_FNV, design_matrix_FNV, robust = TRUE)

# get results
res_FNN_30 <- list(# PNA79_30m_vs_H2O_30m
                  PNA79_30m_vs_H2O_30m = glmQLFTest(fit_FNN_30, contrast = contrasts_FNN_30[,1]),
                  # PNA79_30m_vs_PNAscr_30m
                  PNA79_30m_vs_PNAscr_30m = glmQLFTest(fit_FNN_30, contrast = contrasts_FNN_30[,2]),
                  # PNAscr_30m_vs_H2O_30m
                  PNAscr_30m_vs_H2O_30m = glmQLFTest(fit_FNN_30, contrast = contrasts_FNN_30[,3]))

res_FNV_30 <- list(# PNA79_30m_vs_H2O_30m
                    PNA79_30m_vs_H2O_30m = glmQLFTest(fit_FNV_30, contrast = contrasts_FNV_30[,1]),
                    # PNA79_30m_vs_PNAscr_30m
                    PNA79_30m_vs_PNAscr_30m = glmQLFTest(fit_FNV_30, contrast = contrasts_FNV_30[,2]),
                    # PNAscr_30m_vs_H2O_30m
                    PNAscr_30m_vs_H2O_30m = glmQLFTest(fit_FNV_30, contrast = contrasts_FNV_30[,3]))


# apply FDR to all tables
res_FNN_30 <- lapply(res_FNN_30, function(r) {
  r$table$FDR <- p.adjust(r$table$PValue, method = "fdr")
  r
})

res_FNV_30 <- lapply(res_FNV_30, function(r) {
  r$table$FDR <- p.adjust(r$table$PValue, method = "fdr")
  r
})


# go through all contrasts and create volcano plots:
volc_list <- list()
for (i in names(res_FNN_30)){
  restab <- res_FNN_30[[i]]$table
  ttle <- gsub("_", " ", i)
  ttle <- gsub("vs", "vs.", ttle)
  # make volcano plot
  volc <- do_volcano(restab, pointsize = 2, x_limit = 5, y_limit = 6, show_sig = T,
                     alpha = 0.01, color_sig = T, title = ttle,
                     minlogfc = 1.5, add_labels = T,
                     color_threshold_lines = "black") +
    geom_point(data = restab[2077:2104,], aes(x = logFC, y = -log10(FDR), fill = "sRNA"),
               cex = 2, shape = 21, alpha = 0.6)
  volc_list <- c(volc_list, list(volc))
  volc
  restab <- restab[order(restab$PValue),]
  # save the results table as excel file. make first line thicker
  write.xlsx(restab, paste0("./analysis/DE_raw_data/", i, "_DE_result_FNN_30.xlsx"), row.names = T, col.names = T)
}

volc_grid <- plot_grid(plotlist = volc_list[c(1,3,2)], ncol = 3, nrow = 1,
                       labels = c("A", "B", "C"), label_size = 20, scale=0.95)

# save the grid as pdf
ggsave("analysis/volcano_plots_FNN_30.pdf", volc_grid, width = 35, height = 12, units = "cm")

# do same for FNV:
volc_list <- list()
for (i in names(res_FNV_30)){
  restab <- res_FNV_30[[i]]$table
  ttle <- gsub("_", " ", i)
  ttle <- gsub("vs", "vs.", ttle)
  # make volcano plot
  volc <- do_volcano(restab, pointsize = 2, x_limit = 5, y_limit = 6, show_sig = T,
                     alpha = 0.01, color_sig = T, title = ttle,
                     minlogfc = 1.5, add_labels = T,
                     color_threshold_lines = "black") +
    geom_point(data = restab[2006:2038,], aes(x = logFC, y = -log10(FDR), fill = "sRNA"),
                cex = 2, shape = 21, alpha = 0.6)

  volc_list <- c(volc_list, list(volc))

  # sort restab by p-value starting with the smallest
  restab <- restab[order(restab$PValue),]
  # save the results table as excel file. make first line thicker
  write.xlsx(restab, paste0("./analysis/DE_raw_data/", i, "_DE_result_FNV_30.xlsx"), row.names = T, col.names = T)
}

volc_grid <- plot_grid(plotlist = volc_list[c(1,3,2)], ncol = 3, nrow = 1,
                       labels = c("A", "B", "C"), label_size = 20, scale=0.95)

# save the grid as pdf
ggsave("analysis/volcano_plots_FNV_30.pdf", volc_grid, width = 35, height = 12, units = "cm")




# save tables from res_FNN_30 and res_FNV_30 as excel files
# start by

# get table of res_FNN_30 fot pna97 vs h2o
res_FNN_PNA79_30m_vs_H2O_30m <- res_FNN_30$PNA79_30m_vs_H2O_30m$table %>%
  # change rownames to id column
  rownames_to_column(var = "id") %>% as_tibble()

# get table of res_FNN_30 fot pna97 vs pnascr
res_FNN_PNA79_30m_vs_PNAscr_30m <- res_FNN_30$PNA79_30m_vs_PNAscr_30m$table %>%
  # change rownames to id column
  rownames_to_column(var = "id") %>% as_tibble()

# combine attributes_FNN with res_FNN_PNA79_30m_vs_H2O_30m
res_FNN_PNA79_30m_vs_H2O_30m_excel_ready <- inner_join( attributes_FNN, res_FNN_PNA79_30m_vs_H2O_30m, by = "id") %>%
  # sort by P-value
  arrange(PValue) %>%
  # select only id, type, product, logCPM, logFC, PValue, FDR
  select(id, type, product, logCPM, logFC, PValue, FDR) %>%
  # remove duplicates
  distinct()

# combine attributes_FNN with res_FNN_PNA79_30m_vs_PNAscr_30m
res_FNN_PNA79_30m_vs_PNAscr_30m_excel_ready <- inner_join( attributes_FNN, res_FNN_PNA79_30m_vs_PNAscr_30m, by = "id") %>%
  # sort by P-value
  arrange(PValue) %>%
  # select only id, type, product, logCPM, logFC, PValue, FDR
  select(id, type, product, logCPM, logFC, PValue, FDR) %>%
  # remove duplicates
  distinct()


# ok now get the same for FNV
res_FNV_PNA79_30m_vs_H2O_30m <- res_FNV_30$PNA79_30m_vs_H2O_30m$table %>%
  rownames_to_column(var = "id") %>% as_tibble()

res_FNV_PNA79_30m_vs_PNAscr_30m <- res_FNV_30$PNA79_30m_vs_PNAscr_30m$table %>%
    rownames_to_column(var = "id") %>% as_tibble()

# combine attributes_FNV with res_FNV_PNA79_30m_vs_H2O_30m
res_FNV_PNA79_30m_vs_H2O_30m_excel_ready <- inner_join( attributes_FNV, res_FNV_PNA79_30m_vs_H2O_30m, by = "id") %>%
  arrange(PValue) %>%
  # select only id, type, product, logCPM, logFC, PValue, FDR
  select(id, type, product, logCPM, logFC, PValue, FDR) %>%
  # remove duplicates
  distinct()

# combine attributes_FNV with res_FNV_PNA79_30m_vs_PNAscr_30m
res_FNV_PNA79_30m_vs_PNAscr_30m_excel_ready <- inner_join( attributes_FNV, res_FNV_PNA79_30m_vs_PNAscr_30m, by = "id") %>%
  arrange(PValue) %>%
  # select only id, type, product, logCPM, logFC, PValue, FDR
  select(id, type, product, logCPM, logFC, PValue, FDR) %>%
  # remove duplicates
  distinct()


# write excel files. start with FNN23 PNA79 vs H2O 30min, then FNN23 PNA79 vs PNAscr 30min, then FNV PNA79 vs H2O 30min,
# then FNV PNA79 vs PNAscr 30min, then FNN23 PNA79 vs H2O 16h, then FNN23 PNA79 vs PNAscr 16h, then FNV PNA79 vs H2O 16h,
# then FNV PNA79 vs PNAscr 16h. all in one file, differeent sheets.
# start with sheet FNV PNA79 vs PNAscr 30min:
write.xlsx(as.data.frame(res_FNN_PNA79_30m_vs_H2O_30m_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
                  sheetName = "FNN23_FUS79_vs_H2O_30min", row.names = F, col.names = T)
write.xlsx(as.data.frame(res_FNN_PNA79_30m_vs_PNAscr_30m_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
              sheetName = "FNN23_FUS79_vs_PNAscr_30min", row.names = F, col.names = T, append = T)
write.xlsx(as.data.frame(res_FNV_PNA79_30m_vs_H2O_30m_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
                sheetName = "FNV_FUS79_vs_H2O_30min", row.names = F, col.names = T, append = T)
write.xlsx(as.data.frame(res_FNV_PNA79_30m_vs_PNAscr_30m_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
                sheetName = "FNV_FUS79_vs_PNAscr_30min", row.names = F, col.names = T, append = T)
write.xlsx(as.data.frame(res_FNN_PNA79_16h_vs_H2O_16h_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
                sheetName = "FNN23_FUS79_vs_H2O_16h", row.names = F, col.names = T, append = T)
write.xlsx(as.data.frame(res_FNN_PNA79_16h_vs_PNAscr_16h_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
                sheetName = "FNN23_FUS79_vs_PNAscr_16h", row.names = F, col.names = T, append = T)
write.xlsx(as.data.frame(res_FNV_PNA79_16h_vs_H2O_16h_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
                sheetName = "FNV_FUS79_vs_H2O_16h", row.names = F, col.names = T, append = T)
write.xlsx(as.data.frame(res_FNV_PNA79_16h_vs_PNAscr_16h_excel_ready), "analysis/Supplementary_Data_SD1.xlsx",
                sheetName = "FNV_FUS79_vs_PNAscr_16h", row.names = F, col.names = T, append = T)






a














# Now do kegg pathway analysis. First get the gene names of the DE genes:

# OK, now we need the length of contrasts
l <- length(colnames(contrasts_FNN_30))

# now run FRY gene set test for all contrasts:
kegg_fry_FNN_30 <- lapply(1:l, function(x) fry(y_FNN, idx_kegg_FNN, design_matrix, contrasts_FNN_30[,x]))
names(kegg_fry_FNN_30) <- colnames(contrasts_FNN_30)
kegg_fry_FNV_30 <- lapply(1:l, function(x) fry(y_FNV, idx_kegg_FNV, design_matrix_FNV, contrasts_FNV_30[,x]))
names(kegg_fry_FNV_30) <- colnames(contrasts_FNV_30)

# reorder kegg_fry_FNN_30
kegg_fry_FNN_30 <- kegg_fry_FNN_30[c(1,3,2)]
kegg_fry_FNV_30 <- kegg_fry_FNV_30[c(1,3,2)]

# add KEGG terms
for (fryres in names(kegg_fry_FNN_30)){
  kegg_fry_FNN_30[[fryres]][["TERM"]] <- ifelse(grepl("fnu", rownames(kegg_fry_FNN_30[[fryres]])),
                                             list_kegg_FNN[rownames(kegg_fry_FNN_30[[fryres]])],
                                             rownames(kegg_fry_FNN_30[[fryres]]))
  kegg_fry_FNN_30[[fryres]][["TERM"]] <- gsub(" - Fusobacterium nucleatum subsp. nucleatum ATCC 25586", "",
                                           kegg_fry_FNN_30[[fryres]][["TERM"]])
  # add a column containing all locus tags
  kegg_fry_FNN_30[[fryres]][["GENES"]] <- as.character(unlist(sapply(idx_kegg_FNN[rownames(kegg_fry_FNN_30[[fryres]])], function(x) paste(x,collapse=";"))))
  write.csv(kegg_fry_FNN_30[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNN_30.csv"))
  # write as excel file
  write.xlsx(kegg_fry_FNN_30[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNN_30.xlsx"), row.names = T
                , col.names = T)
}

for (fryres in names(kegg_fry_FNV_30)){
  kegg_fry_FNV_30[[fryres]][["TERM"]] <- ifelse(grepl("fnc", rownames(kegg_fry_FNV_30[[fryres]])),
                                             list_kegg_FNV[rownames(kegg_fry_FNV_30[[fryres]])],
                                             rownames(kegg_fry_FNV_30[[fryres]]))
  kegg_fry_FNV_30[[fryres]][["TERM"]] <- gsub(" - Fusobacterium vincentii 3_1_36A2", "",
                                           kegg_fry_FNV_30[[fryres]][["TERM"]])
  # add a column containing all locus tags
  kegg_fry_FNV_30[[fryres]][["GENES"]] <- as.character(unlist(sapply(idx_kegg_FNV[rownames(kegg_fry_FNV_30[[fryres]])], function(x) paste(x,collapse=";"))))
  write.csv(kegg_fry_FNV_30[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNV_30.csv"))
  # write as excel file
    write.xlsx(kegg_fry_FNV_30[[fryres]], paste0("./analysis/KEGG_raw_data/", fryres, "_KEGG_FNV_30.xlsx"), row.names = T
                , col.names = T)
}

# filter for pws with FDR < 0.01 & more than 4 genes
kegg_frysig_FNN_30 <- lapply(kegg_fry_FNN_30, function(x) x[x[["FDR"]]<0.05 & x[["NGenes"]]>4,])
kegg_frysig_FNV_30 <- lapply(kegg_fry_FNV_30, function(x) x[x[["FDR"]]<0.05 & x[["NGenes"]]>4,])

# select siggos
kegg_siggos_FNN_30 <- c()
kegg_siggos_FNV_30 <- c()

for (i in names(kegg_frysig_FNN_30)) {
  print(i)
  print(dim(kegg_frysig_FNN_30[[i]]))
  print(kegg_frysig_FNN_30[[i]][,c(1,2,4,7)])
  kegg_siggos_FNN_30 <- c(kegg_siggos_FNN_30, rownames(kegg_frysig_FNN_30[[i]][,]))  # can be modified
}
# remove duplicates
kegg_siggos_FNN_30 <- unique(kegg_siggos_FNN_30[!grepl("NA", kegg_siggos_FNN_30)])

for (i in names(kegg_frysig_FNV_30)) {
  print(i)
  print(dim(kegg_frysig_FNV_30[[i]]))
  print(kegg_frysig_FNV_30[[i]][,c(1,2,4,7)])
  kegg_siggos_FNV_30 <- c(kegg_siggos_FNV_30, rownames(kegg_frysig_FNV_30[[i]]))  # can be modified
}
# remove duplicates
kegg_siggos_FNV_30 <- unique(kegg_siggos_FNV_30[!grepl("NA", kegg_siggos_FNV_30)])

# I create a dataframe with mean logFC values for each significant kegg-term:
hm_kegg_FNN_30 <- t(as.data.frame(lapply(idx_kegg_char_FNN[kegg_siggos_FNN_30], function(x){
  sapply(names(res_FNN_30), function(y){
    mean(res_FNN_30[[y]]$table[x,]$logFC)
  })
}))
)

hm_kegg_FNV_30 <- t(as.data.frame(lapply(idx_kegg_char_FNV[kegg_siggos_FNV_30], function(x){
  sapply(names(res_FNV_30), function(y){
    mean(res_FNV_30[[y]]$table[x,]$logFC)
  })
}))
)

# reorder columns
hm_kegg_FNN_30 <- hm_kegg_FNN_30[,c(1,3,2)]
hm_kegg_FNV_30 <- hm_kegg_FNV_30[,c(1,3,2)]

# make heatmap:
hm_kegg_FNN_30 <- hm_kegg_FNN_30[order(hm_kegg_FNN_30[,1], decreasing = T),]
hm_kegg_FNV_30 <- hm_kegg_FNV_30[order(hm_kegg_FNV_30[,1], decreasing = T),]

# add size of pws:
kegg_sizes_FNN_30 <- sapply(idx_kegg_char_FNN[rownames(hm_kegg_FNN_30)], function(x) length(x))
kegg_sizes_FNV_30 <- sapply(idx_kegg_char_FNV[rownames(hm_kegg_FNV_30)], function(x) length(x))

# add p values
pvals_FNN_30 <- data.frame(sapply(names(kegg_fry_FNN_30),
                               function(x) kegg_fry_FNN_30[[x]][rownames(hm_kegg_FNN_30),"FDR"]),
                        row.names = rownames(hm_kegg_FNN_30))
pvals_FNV_30 <- data.frame(sapply(names(kegg_fry_FNV_30),
                               function(x) kegg_fry_FNV_30[[x]][rownames(hm_kegg_FNV_30),"FDR"]),
                        row.names = rownames(hm_kegg_FNV_30))

#select only significant ones:
pvals_FNN_30 <-sapply(pvals_FNN_30, function(x) ifelse(x<0.05, x <- "*", x<-"") )
pvals_FNV_30 <-sapply(pvals_FNV_30, function(x) ifelse(x<0.05, x <- "*", x<-"") )

# add term
keggpws_FNN_30 <- kegg_fry_FNN_30$PNA79_30m_vs_H2O_30m[rownames(hm_kegg_FNN_30),] [["TERM"]]
keggpws_FNV_30 <- kegg_fry_FNV_30$PNA79_30m_vs_H2O_30m[rownames(hm_kegg_FNV_30),] [["TERM"]]

# add names to hm_kegg
rownames(hm_kegg_FNN_30) <- ifelse(!is.na(keggpws_FNN_30),keggpws_FNN_30, rownames(hm_kegg_FNN_30) )
rownames(hm_kegg_FNV_30) <- ifelse(!is.na(keggpws_FNV_30),keggpws_FNV_30, rownames(hm_kegg_FNV_30) )

colnames(hm_kegg_FNN_30) <- c("PNA79_vs_H20", "PNAscr_vs_H20", "PNA79_vs_PNAscr")

ht_vert_FNN_30 <- Heatmap(hm_kegg_FNN_30, cluster_columns = F,cluster_rows = F,
                          name = "GO-analysis", col = col_fun,
                          show_heatmap_legend = F,
                          #column_split = rep(c("PNA","scrambled","control"), each=4),
                          row_title_side = "right", row_title_rot = 0,
                          border = TRUE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.1s", pvals_FNN_30[i, j]), x, y)
                          },
                          column_names_gp = gpar(fontsize = 11),
                          column_title_gp = gpar(fontsize = 15),
                          row_names_gp = gpar(fontsize = 10),
                          # add column title
                          column_title = "Time point 30m FNN",
                          row_title = NULL,
                          width = unit(4, "cm"), height = unit(5, "cm"),

                          right_annotation = rowAnnotation(genes = anno_barplot(kegg_sizes_FNN_30)))

ht_vert_FNN_30

lgd <- Legend(col_fun = col_fun, title = expression("mean log"[2]*" FC"), #direction = "horizontal",
              title_gp = gpar(fontsize = 20), labels = c("-1", " 0"," 1"), legend_height = unit(8, "cm"),
              grid_width = unit(1, "cm"),
              labels_gp = gpar(fontsize = 20),
              at = c(-1, 0, 1), border = "black",
              title_position = "leftcenter-rot")

svg("analysis/KEGG_heatmap_FNN_30.svg", width = 9, height = 7)
ht_vert_FNN_30
draw(lgd, x = unit(2, "cm"), y = unit(10, "cm"))
dev.off()

pdf("analysis/KEGG_heatmap_FNN_30.pdf", width = 9, height = 7)
ht_vert_FNN_30
draw(lgd, x = unit(2, "cm"), y = unit(10, "cm"))
dev.off()

colnames(hm_kegg_FNV_30) <- c("PNA79_vs_H20", "PNAscr_vs_H20", "PNA79_vs_PNAscr")

ht_vert_FNV_30 <- Heatmap(hm_kegg_FNV_30, cluster_columns = F,cluster_rows = F,
                          name = "GO-analysis", col = col_fun,
                          show_heatmap_legend = F,
                          #column_split = rep(c("PNA","scrambled","control"), each=4),
                          row_title_side = "right", row_title_rot = 0,
                          border = TRUE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.1s", pvals_FNV_30[i, j]), x, y)
                          },
                          column_names_gp = gpar(fontsize = 11),
                          column_title_gp = gpar(fontsize = 15),
                          row_names_gp = gpar(fontsize = 10),
                          # add column title
                          column_title = "Time point 30m FNV",
                          row_title = NULL,
                          width = unit(4, "cm"), height = unit(8, "cm"),

                          right_annotation = rowAnnotation(genes = anno_barplot(kegg_sizes_FNV_30)))

ht_vert_FNV_30

svg("analysis/KEGG_heatmap_FNV_30.svg", width = 9, height = 7)
ht_vert_FNV_30
draw(lgd, x = unit(2, "cm"), y = unit(10, "cm"))
dev.off()

pdf("analysis/KEGG_heatmap_FNV_30.pdf", width = 9, height = 7)
ht_vert_FNV_30
draw(lgd, x = unit(2, "cm"), y = unit(10, "cm"))
dev.off()


# get fatty acid metabolism genes:
fatty_acid_genes_FNV <- unlist(idx_kegg_FNV[names(list_kegg_FNV[grep("Fatty acid metabolism", list_kegg_FNV)])])

# create volcano plot for fnv 30 pnascr vs ctrl
restab <- res_FNV_30$PNAscr_30m_vs_H2O_30m$table

volc <- do_volcano(restab, pointsize = 2, x_limit = 3, y_limit = 4, show_sig = T,
                   alpha = 0.01, color_sig = T, title = "PNAscr 30m vs H2O 30m",
                   minlogfc = 1.5, add_labels = T,
                   color_threshold_lines = "black") +
  geom_point(data = restab[which(rownames(restab) %in% fatty_acid_genes_FNV),], aes(x = logFC, y = -log10(FDR), fill = "sRNA"),
             cex = 2, shape = 21, alpha = 0.6)
svg("~/Downloads/fatty_acifd_met_volcano.svg", width = 5, height = 5)
volc
dev.off()




# get all significant genes up and downreg for FNv 30 min (log2fc > 1.5 & FDR < 0.05)
t_fnv_30 <- as_tibble(t(sapply(res_FNV_30, function(x) {
    up <- rownames(x$table[x$table$logFC > 1.5 & x$table$FDR < 0.01,])
    down <- rownames(x$table[x$table$logFC < -1.5 & x$table$FDR < 0.01,])
    c(up = length(up), down = length(down))
})), rownames = "contrast") %>% mutate(organism = "FNV", time = "30m")

# same at 16h
t_fnv_16 <- as_tibble(t(sapply(res_FNV, function(x) {
    up <- rownames(x$table[x$table$logFC > 1.5 & x$table$FDR < 0.01,])
    down <- rownames(x$table[x$table$logFC < -1.5 & x$table$FDR < 0.01,])
    c(up = length(up), down = length(down))
})), rownames = "contrast") %>% mutate(organism = "FNV", time = "16h")

# same for fnn
t_fnn_30 <- as_tibble(t(sapply(res_FNN_30, function(x) {
    up <- rownames(x$table[x$table$logFC > 1.5 & x$table$FDR < 0.01,])
    down <- rownames(x$table[x$table$logFC < -1.5 & x$table$FDR < 0.01,])
    c(up = length(up), down = length(down))
})), rownames = "contrast") %>% mutate(organism = "FNN", time = "30m")

t_fnn_16 <- as_tibble(t(sapply(res_FNN, function(x) {
    up <- rownames(x$table[x$table$logFC > 1.5 & x$table$FDR < 0.01,])
    down <- rownames(x$table[x$table$logFC < -1.5 & x$table$FDR < 0.01,])
    c(up = length(up), down = length(down))
})), rownames = "contrast") %>% mutate(organism = "FNN", time = "16h")

# make one table
t <- rbind(t_fnv_30, t_fnv_16, t_fnn_30, t_fnn_16)

# save as excel...
write.xlsx(t, "analysis/significant_genes.xlsx", row.names = T, col.names = T)



# for kegg, make a table with TERM, ngenes, FDR (-log10), Direction, strain for T16
table_kegg_fnn_pna79_vs_scr <- kegg_fry_FNN$PNA79_16h_vs_PNAscr_16h %>%
    select(TERM, NGenes, FDR, Direction) %>%
    mutate(Strain = "FNN23") %>% mutate(FDR = -log10(FDR))%>%
  # make fdr negative if direction is Down
    mutate(FDR = ifelse(Direction == "Down", -FDR, FDR))

# same for FNV
table_kegg_fnv_pna79_vs_scr <- kegg_fry_FNV$PNA79_16h_vs_PNAscr_16h %>%
    select(TERM, NGenes, FDR, Direction) %>%
    mutate(Strain = "FNV") %>% mutate(FDR = -log10(FDR)) %>%
  # make fdr negative if direction is Down
    mutate(FDR = ifelse(Direction == "Down", -FDR, FDR))

# combine
table_kegg <- rbind(table_kegg_fnn_pna79_vs_scr, table_kegg_fnv_pna79_vs_scr)

# get sgnificant pathways significant <0.01
sig_pws <- na.omit(unique(table_kegg$TERM[abs(table_kegg$FDR) > 2 & table_kegg$NGenes > 5]))

table_kegg <- table_kegg[table_kegg$TERM %in% sig_pws,]

# make term a factor ordered by abs(FDR)
table_kegg$TERM <- factor(table_kegg$TERM, levels = table_kegg[table_kegg$Strain == "FNN23",]$TERM[order(abs(table_kegg[table_kegg$Strain == "FNN23",]$FDR))])


# keep only significant pathways

# make a ggplot barplot
kegg_barplot_16h <- ggplot(table_kegg, aes(x = TERM, y = FDR, fill = Strain)) +
  geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    theme_bw() +
  ylim(-2.5, 7.5)+
  # remove background
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # add a line at 0.05
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  geom_hline(yintercept = log10(0.05), linetype = "dotted", color = "black") +
  scale_fill_manual(values = c("FNN23" = viridis(10, option="plasma")[3], "FNV" = viridis(10, option="plasma")[8])) +
  scale_y_continuous(breaks = c(log10(0.05),0, -log10(0.05)), labels = c("0.05","0", "0.05")) +
  # change x axis title
    ylab("FDR-corrected P-value") +
  # add vertical lines to all bars
    geom_vline(xintercept = 1:dim(table_kegg)[1], linetype = "dotted", color = "grey", alpha=0.5) +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          # remove y axis title
            axis.title.y=element_blank(),
          # increase font xaxis title
            axis.title.x = element_text(size = 12),
    # put legend into lower right of plot
    legend.position = c(0.85, 0.1),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size = 10))

pdf("analysis/KEGG_barplot_16h.pdf", width = 7, height = 7)
kegg_barplot_16h
dev.off()

