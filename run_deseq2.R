#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2); library(tidyverse); library(EnhancedVolcano); library(pheatmap); library(apeglm)
})
args <- commandArgs(trailingOnly = TRUE)
counts_file <- ifelse(length(args) >= 1, args[1], 'data/counts.csv')
metadata_file <- ifelse(length(args) >= 2, args[2], 'data/samples.csv')
output_dir <- ifelse(length(args) >= 3, args[3], 'results')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1)
metadata <- metadata[colnames(counts), , drop = FALSE]
keep <- rowSums(counts >= 10) >= ncol(counts) * 0.5
counts_filtered <- counts[keep, ]
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ condition)
dds <- DESeq(dds)
res_shrunk <- lfcShrink(dds, coef = 2, type = 'apeglm')
sig_genes <- as.data.frame(res_shrunk) %>% rownames_to_column('gene') %>% filter(padj < 0.05) %>% arrange(padj)
write.csv(sig_genes, file.path(output_dir, 'significant_DEGs.csv'), row.names = FALSE)
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = 'condition', returnData = TRUE)
ggsave(file.path(output_dir, 'pca_plot.png'), ggplot(pca_data, aes(PC1, PC2, colour = condition)) + geom_point(size = 3) + theme_minimal())
png(file.path(output_dir, 'volcano_plot.png'), width = 800, height = 600)
EnhancedVolcano(res_shrunk, lab = rownames(res_shrunk), x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 1)
dev.off()
cat('Pipeline complete.\n')
