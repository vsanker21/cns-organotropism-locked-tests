# Discovery only. mo-Mac pseudobulk, brain metastasis versus lymph-node metastasis.
# Rule: results/locked/momac_discovery_rule.txt

suppressPackageStartupMessages(library(DESeq2))

out_dir <- "g:/Cancer Metastasis/results/macrophage"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read.csv("g:/Cancer Metastasis/results/locked/gse131907_momac_samples.csv", stringsAsFactors = FALSE)
ann <- read.delim(
  gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
  stringsAsFactors = FALSE, check.names = FALSE
)
mac <- ann[ann$Cell_subtype == "mo-Mac" & ann$Sample_Origin %in% c("mBrain", "mLN"), ]
cat("MOMAC_CELLS", nrow(mac), "\n")
cat("SUBTYPES_KEPT", paste(unique(mac$Cell_subtype), collapse = ","), "\n")

obj <- readRDS(gzcon(gzfile(
  "g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz",
  open = "rb"
)))
symbols <- toupper(rownames(obj))
keep <- intersect(mac$Index, colnames(obj))
cat("MATCHED", length(keep), "OF", nrow(mac), "\n")
sub <- as.matrix(obj[, keep, drop = FALSE])
rm(obj)
gc()

sample_of <- mac$Sample[match(keep, mac$Index)]
# Collapse duplicate symbols, then sum mo-Mac cells within each sample.
sub <- rowsum(sub, symbols)
mode(sub) <- "integer"
samples <- clinical$sample
missing <- setdiff(samples, unique(sample_of))
if (length(missing)) stop("Samples with no mo-Mac cells: ", paste(missing, collapse = ","))

counts <- sapply(samples, function(s) {
  use <- sample_of == s
  if (sum(use) == 1) sub[, use] else rowSums(sub[, use, drop = FALSE])
})
counts <- round(counts)
mode(counts) <- "integer"
rownames(counts) <- rownames(sub)
n_cells <- as.integer(sapply(samples, function(s) sum(sample_of == s)))
coldata <- data.frame(
  sample = samples,
  origin = factor(clinical$origin, levels = c("mLN", "mBrain")),
  sex = factor(clinical$sex),
  n_momac = n_cells,
  row.names = samples,
  stringsAsFactors = FALSE
)
write.csv(coldata, file.path(out_dir, "momac_sample_qc.csv"), row.names = FALSE)
cat("N_CELLS", paste(samples, n_cells, sep = "=", collapse = ","), "\n")

run_deseq <- function(count_mat, meta, design_formula) {
  keep_genes <- rowSums(count_mat) >= 10
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat[keep_genes, , drop = FALSE],
    colData = meta,
    design = design_formula
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("origin", "mBrain", "mLN"))
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[order(df$pvalue), ]
  list(dds = dds, table = df)
}

perm_p <- function(vsd_mat, group, n_perm = 9999) {
  stat <- function(g) {
    brain <- rowMeans(vsd_mat[, g == "mBrain", drop = FALSE])
    node <- rowMeans(vsd_mat[, g == "mLN", drop = FALSE])
    sum((brain - node)^2)
  }
  observed <- stat(group)
  set.seed(20261112)
  null <- replicate(n_perm, stat(sample(group)))
  list(stat = observed, p = (sum(null >= observed) + 1) / (n_perm + 1))
}

primary <- run_deseq(counts, coldata, ~ origin)
write.csv(primary$table, file.path(out_dir, "momac_brain_vs_ln_primary.csv"), row.names = FALSE)
vsd <- assay(vst(primary$dds, blind = TRUE))
group <- as.character(coldata[colnames(vsd), "origin"])
glob <- perm_p(vsd, group)
n_fdr <- sum(primary$table$padj < 0.05, na.rm = TRUE)
n_up <- sum(primary$table$padj < 0.05 & primary$table$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(primary$table$padj < 0.05 & primary$table$log2FoldChange < 0, na.rm = TRUE)

gate <- c("TMEM119", "P2RY12", "CX3CR1", "SALL1")
gate_df <- primary$table[primary$table$gene %in% gate, c("gene", "log2FoldChange", "pvalue", "padj")]

# Sensitivity: at least 50 mo-Mac cells.
big <- coldata$n_momac >= 50
sens_n <- run_deseq(counts[, big, drop = FALSE], coldata[big, , drop = FALSE], ~ origin)
write.csv(sens_n$table, file.path(out_dir, "momac_brain_vs_ln_min50.csv"), row.names = FALSE)

# Sensitivity: sex adjustment. Primary remains the unadjusted contrast.
sens_sex <- run_deseq(counts, coldata, ~ sex + origin)
write.csv(sens_sex$table, file.path(out_dir, "momac_brain_vs_ln_sex_adjusted.csv"), row.names = FALSE)

summary_lines <- c(
  paste("N_BRAIN", sum(coldata$origin == "mBrain")),
  paste("N_LN", sum(coldata$origin == "mLN")),
  paste("N_GENES_TESTED", nrow(primary$table)),
  paste("N_FDR05", n_fdr),
  paste("N_FDR05_BRAIN_HIGH", n_up),
  paste("N_FDR05_LN_HIGH", n_down),
  paste("PERM_STAT", signif(glob$stat, 6)),
  paste("PERM_P", signif(glob$p, 6)),
  paste("SEPARATES", if (glob$p < 0.05 || n_fdr > 0) "yes" else "no"),
  paste("MIN50_SAMPLES", sum(big)),
  paste("MIN50_FDR05", sum(sens_n$table$padj < 0.05, na.rm = TRUE)),
  paste("SEXADJ_FDR05", sum(sens_sex$table$padj < 0.05, na.rm = TRUE)),
  paste("TOP_BRAIN", paste(head(primary$table$gene[primary$table$log2FoldChange > 0], 15), collapse = ",")),
  paste("TOP_LN", paste(head(primary$table$gene[primary$table$log2FoldChange < 0], 15), collapse = ","))
)
writeLines(summary_lines, file.path(out_dir, "momac_summary.txt"))
write.csv(gate_df, file.path(out_dir, "momac_microglia_gate.csv"), row.names = FALSE)
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("GATE\n")
print(gate_df)
