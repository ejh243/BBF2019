#!/usr/bin/env Rscript

## Script to calculate gene-level counts between SCZ case and control using DESeq2

# Load libraries
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is not installed. Please install it using BiocManager before running this script.")
}
library(DESeq2)
library(ggplot2)


## LOAD COUNT AND METADATA ## 
# Load count data - the output .txt file from running FeatureCounts 
counts <- read.table(
  "/projects/e6e/RNASeq/Exeter_matched/featureCounts/gene_counts_primary_bam_s2.txt",
  header=TRUE,
  row.names=1,
  comment.char="#"
)
# Cleanup column names 
colnames(counts) <- sub(
  ".*\\.([0-9]+_[A-Z]+_L[0-9]{3})Aligned.*",
  "\\1",
  colnames(counts)
)
# Remove annotation columns 
counts <- counts[, 6:ncol(counts)]

# Print log
cat("Reading Gene Count data...", "\n")
cat("Genes:", nrow(counts), "\n")
cat("Samples:", ncol(counts), "\n")
cat("Reading Sample metadata...", "\n")

# Load metadata - a .tsv file containing sample IDs and experimental design
metadata <- read.table(
  "/projects/e6e/RNASeq/Exeter_matched/exp_design.tsv",
  header=TRUE,
  row.names=1
)
# Reorder samples in metadata to match count data
metadata <- metadata[match(colnames(counts), rownames(metadata)), ]


# Check sample matching
if (!all(colnames(counts) == rownames(metadata))) {
  stop("Sample names in count and metadata do not match!")
}
cat("Sample names in count and metadata successfully matched.\n")
cat("\n")


# Set output dir
outdir <- "/projects/e6e/RNASeq/Exeter_matched/gene_level_analysis"


## RUN DESEQ2 AND SAVE THE RESULTS ## 
cat("Starting DE analysis...", "\n")
cat("\n")

# Create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2 and results 
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), file.path(outdir, "DE_results.csv"), row.names = TRUE)
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

ggsave(filename = file.path(outdir, "PCA_plot.png"), plot = pca_plot, width = 6,height = 5)


# MA plot 
res_df <- as.data.frame(res)

ma_plot <- ggplot(res_df, aes(baseMean, log2FoldChange)) +
  geom_point(alpha=0.3) +
  scale_x_log10() +
  geom_hline(yintercept=0, color="red") +
  labs(title="MA Plot", x="Mean expression", y="Log2 fold change") +
  theme_bw()

ggsave(filename = file.path(outdir, "MA_plot.png"), plot = ma_plot, width = 6,height = 5)


# Volcano plot 
res_df$padj[is.na(res_df$padj)] <- 1

volcano_plot <- ggplot(res_df,
  aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.4) +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(adj p-value)") +
  theme_bw()

ggsave(filename = file.path(outdir, "Volcano_plot.png"), plot = volcano_plot, width = 6,height = 5)


# Significant genes 
res_df$significant <- res_df$padj < 0.05

sig_genes <- ggplot(res_df, 
  aes(log2FoldChange, -log10(padj), color=significant)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("grey", "red")) +
  theme_bw()

ggsave(filename = file.path(outdir, "Volcano_plot_sig_genes.png"), plot = sig_genes, width = 6,height = 5)


# End of script
cat("Exploratory plots generated.\n")
