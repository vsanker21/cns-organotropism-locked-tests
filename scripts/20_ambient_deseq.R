# Repeat the mo-Mac contrast on ambient-corrected pseudobulks.
# Gene selection uses only these corrected counts.

suppressPackageStartupMessages(library(DESeq2))

out_dir <- "g:/Cancer Metastasis/results/macrophage"
clinical <- read.csv("g:/Cancer Metastasis/results/locked/gse131907_momac_samples.csv", stringsAsFactors = FALSE)
raw <- read.csv(file.path(out_dir, "momac_ambient_pseudobulk.csv"), check.names = FALSE)
genes <- raw$gene
counts <- as.matrix(raw[, clinical$sample])
storage.mode(counts) <- "integer"
rownames(counts) <- genes
coldata <- data.frame(
  origin = factor(clinical$origin, levels = c("mLN", "mBrain")),
  row.names = clinical$sample
)

keep <- rowSums(counts) >= 10
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = coldata, design = ~ origin)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("origin", "mBrain", "mLN"))
table <- as.data.frame(res)
table$gene <- rownames(table)
table <- table[order(table$pvalue), ]
write.csv(table, file.path(out_dir, "momac_ambient_brain_vs_ln.csv"), row.names = FALSE)

vsd <- assay(vst(dds, blind = TRUE))
group <- as.character(coldata[colnames(vsd), "origin"])
stat <- function(g) {
  brain <- rowMeans(vsd[, g == "mBrain", drop = FALSE])
  node <- rowMeans(vsd[, g == "mLN", drop = FALSE])
  sum((brain - node)^2)
}
observed <- stat(group)
set.seed(20261112)
null <- replicate(9999, stat(sample(group)))
perm_p <- (sum(null >= observed) + 1) / 10000
n_fdr <- sum(table$padj < 0.05, na.rm = TRUE)
n_up <- sum(table$padj < 0.05 & table$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(table$padj < 0.05 & table$log2FoldChange < 0, na.rm = TRUE)
separates <- perm_p < 0.05 || n_fdr > 0

previous <- read.csv(file.path(out_dir, "momac_brain_vs_ln_primary.csv"), stringsAsFactors = FALSE)
previous_sig <- previous$gene[which(previous$padj < 0.05)]
residual <- table$gene[which(table$padj < 0.05)]
gate <- c("TMEM119", "P2RY12", "CX3CR1", "SALL1", "GFAP", "PLP1", "APOD", "GLDN", "MS4A1", "CD79A", "CD68", "AIF1")
gate_df <- table[match(gate, table$gene), c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")]
write.csv(gate_df, file.path(out_dir, "momac_ambient_markers.csv"), row.names = FALSE)

if (separates) {
  frozen <- table[table$padj < 0.05 & !is.na(table$padj), c("gene", "log2FoldChange", "padj")]
  frozen$direction <- ifelse(frozen$log2FoldChange > 0, "brain_high", "lymph_node_high")
  write.csv(frozen, "g:/Cancer Metastasis/results/locked/momac_ambient_genes.csv", row.names = FALSE)
}

lines <- c(
  paste("N_GENES_TESTED", nrow(table)),
  paste("N_FDR05", n_fdr),
  paste("N_FDR05_BRAIN_HIGH", n_up),
  paste("N_FDR05_LN_HIGH", n_down),
  paste("PERM_P", signif(perm_p, 6)),
  paste("SEPARATES", if (separates) "yes" else "no"),
  paste("OVERLAP_WITH_UNCORRECTED_280", sum(residual %in% previous_sig)),
  paste("UNCORRECTED_280", length(previous_sig))
)
writeLines(lines, file.path(out_dir, "momac_ambient_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
print(gate_df)
