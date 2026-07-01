#!/usr/bin/env Rscript

## Script to calculate gene-level counts between SCZ case and control using DESeq2

# Load libraries
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is not installed. Please install it using BiocManager before running this script.")
}
library(DESeq2)
library(ggplot2)


## LOAD COUNT AND METADATA ## 
# Load counts
counts <- read.table(
  "/projects/e6e/RNASeq/Exeter_matched/featureCounts/gene_counts_primary_bam.txt",
  header=TRUE,
  row.names=1,
  comment.char="#"
)
# Remove annotation columns
counts <- counts[,6:ncol(counts)]


# Load metadata
coldata <- read.table(
  "/projects/e6e/RNASeq/Exeter_matched/exp_design.tsv",
  header=TRUE,
  row.names=1
)

# Check sample matching
if (!all(colnames(counts) == rownames(coldata))) {
  stop("Sample names in count and metadata do not match!")
}

cat("Sample names in count and metadata successfully matched.\n")


## RUN DESEQ2 AND SAVE THE RESULTS ## 
# Create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2 and results 
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), "DE_results.csv")

cat("DE analysis complete.\n")



## GENERATE EXPLORATORY PLOTS ## 
# PCA 
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  labs(
    title = "PCA Plot",
    x = paste0("PC1: ", round(attr(pcaData, "percentVar")[1]*100), "% variance"),
    y = paste0("PC2: ", round(attr(pcaData, "percentVar")[2]*100), "% variance")
  ) +
  theme_bw()

ggsave("PCA_plot.png", pca_plot, width=6, height=5)


# MA plot 
res_df <- as.data.frame(res)

ma_plot <- ggplot(res_df, aes(baseMean, log2FoldChange)) +
  geom_point(alpha=0.3) +
  scale_x_log10() +
  geom_hline(yintercept=0, color="red") +
  labs(title="MA Plot", x="Mean expression", y="Log2 fold change") +
  theme_bw()

ggsave("MA_plot.png", ma_plot, width=6, height=5)


# Sample distance heatmap
"""
library(pheatmap)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

png(""Sample_distance_heatmap.png"")
pheatmap(sampleDistMatrix)
dev.off()
"""

# Volcano plot 
res_df$padj[is.na(res_df$padj)] <- 1

volcano_plot <- ggplot(res_df,
  aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.4) +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(adj p-value)") +
  theme_bw()

ggsave("Volcano_plot.png", volcano_plot, width=6, height=5)


# Significant genes 
res_df$significant <- res_df$padj < 0.05

sig_genes <- ggplot(res_df, 
  aes(log2FoldChange, -log10(padj), color=significant)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("grey", "red")) +
  theme_bw()

ggsave("Volcano_plot_sig_genes.png", sig_genes, width=6, height=5)
