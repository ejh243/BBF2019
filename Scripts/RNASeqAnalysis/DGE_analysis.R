#!/usr/bin/env Rscript

## Script to run Differential Gene Expression (DGE) analysis 
## And explore expressed genes in the transcriptome, matching short reads data, etc. 

## LOAD REQUIRED LIBRARIES 
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is not installed. Please install it using BiocManager before running this script.")
}
library(DESeq2)
library(ggplot2)


## Need to source Config File and load Conda environment here since there is no wrapper bash script, right? 


## LOAD AND PARSE INPUT DATA 
# Load count data - the count matrix .txt file (FeatureCounts output)
cat("Reading in gene count matrix...\n")

counts <- read.table(
  "/lfs1i3/projects/e6e/RNASeq/Exeter_matched/03_gene_quantification/gene_counts_transcriptome.txt",
  header=TRUE,
  row.names=1,
  comment.char="#"
)

# Cleanup column names (optional)
colnames(counts) <- sub(
  ".*\\.([0-9]+_[A-Z]+_L[0-9]{3})Aligned.*",
  "\\1",
  colnames(counts))

# Remove annotation columns 
counts <- counts[, 6:ncol(counts)]

# Print log
cat("Number of gene features:", nrow(counts), "\n")
cat("Number of RNA-Seq samples:", ncol(counts), "\n")


# Load metadata - a .tsv file containing sample IDs and experimental design
cat("Reading in metadata...\n")

metadata <- read.table(
  "/lfs1i3/projects/e6e/RNASeq/Exeter_matched/metadata/exp_design.tsv",
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


## RUN DESEQ2 
# Create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ batch + sex + age + condition
)

# Filter low counts 
# minimim 10 read counts in at least 3 samples # try with 4 etc and see how the hist changes
cat("Applying set pre-filters:", "\n")
cat("   minimim 10 read counts in at least 3 samples", "\n")

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2 and results 
cat("Running DESeq2...", "\n")
dds <- DESeq(dds)
res <- results(dds, name="condition_SCZ_vs_Control")


# Save output file
outdir <- ""

write.csv(as.data.frame(res), file.path(outdir, "DGE_results.csv"), row.names = TRUE)
cat("DGE analysis complete.\n")


# Count number of features that passed prefilters: 
passed <- round((nrow(dds)/nrow(counts)*100), 2)
cat("Number of gene features that passed prefilters:",nrow(dds), "(",passed,"% )", "\n")

# Get summaries 
summary(res)
summary(res$padj)

# Check no of genes in padj value below 0.05 and 0.01 
sum(res$padj < 0.05, na.rm=TRUE)
sum(res$padj < 0.01, na.rm=TRUE)

res_sig <- res[which(res$padj < 0.1), ]
res_sig[order(res_sig$padj), ]



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
