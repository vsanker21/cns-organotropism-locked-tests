# Proliferation control on the same malignant pseudobulks.
# The Hallmark G2M list was locked before any compartment score was calculated.

suppressPackageStartupMessages(library(Matrix))

prolif <- readLines("g:/Cancer Metastasis/results/locked/proliferation_genes.txt")
prolif <- unique(prolif[nzchar(prolif)])
out_dir <- "g:/Cancer Metastasis/results/compartment"

ann <- read.delim(
  gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
  stringsAsFactors = FALSE, check.names = FALSE
)
mb <- ann[ann$Sample_Origin == "mBrain", ]
obj <- readRDS(gzcon(gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz", open = "rb")))
symbols <- toupper(rownames(obj))
keep_cells <- intersect(mb$Index, colnames(obj))
sub <- as.matrix(obj[, keep_cells, drop = FALSE])
rm(obj)
gc()

sum_genes <- function(mat, wanted) {
  rows <- which(symbols %in% wanted)
  if (!length(rows)) return(rep(0, ncol(mat)))
  as.numeric(Matrix::colSums(mat[rows, , drop = FALSE]))
}
epithelial <- sum_genes(sub, c("KRT8", "KRT18", "KRT19", "EPCAM", "NAPSA", "NKX2-1"))
immune <- sum_genes(sub, c("PTPRC", "CD3D", "CD3E", "CD68", "AIF1", "C1QA"))
tcell <- sum_genes(sub, c("CD3D", "CD3E"))
macrophage <- sum_genes(sub, c("CD68", "AIF1", "C1QA"))
malignant <- epithelial >= 2 & epithelial > immune & tcell == 0 & macrophage == 0
cell_sample <- mb$Sample[match(keep_cells, mb$Index)]

rows <- list()
for (sample_id in sort(unique(cell_sample))) {
  use <- malignant & cell_sample == sample_id
  total_umi <- sum(sub[, use, drop = FALSE])
  present <- symbols %in% prolif
  counts <- tapply(as.numeric(Matrix::rowSums(sub[present, use, drop = FALSE])), symbols[present], sum)
  log2cpm <- log2(counts / as.numeric(total_umi) * 1e6 + 1)
  rec <- as.list(log2cpm)
  rec$sample <- paste0("BrM_", sample_id)
  rows[[sample_id]] <- rec
}
genes <- sort(unique(unlist(lapply(rows, function(x) setdiff(names(x), "sample")))))
out <- data.frame(sample = vapply(rows, function(x) x$sample, character(1)), stringsAsFactors = FALSE)
for (gene in genes) {
  out[[gene]] <- vapply(rows, function(x) {
    value <- x[[gene]]
    if (is.null(value)) NA_real_ else as.numeric(value)
  }, numeric(1))
}
write.csv(out, file.path(out_dir, "parenchyma_proliferation_log2cpm.csv"), row.names = FALSE)
cat("PROLIF_GENES", length(genes), "SAMPLES", nrow(out), "\n")
