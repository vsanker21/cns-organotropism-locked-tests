# Parenchymal brain-metastasis pseudobulk from GSE131907.
# Uses the gene list and malignant-cell rule already frozen in results/locked.
# Author labels do not admit cells. They are recorded only as a secondary count.

suppressPackageStartupMessages(library(Matrix))

locked <- "g:/Cancer Metastasis/results/locked"
out_dir <- "g:/Cancer Metastasis/results/compartment"
lepto <- readLines(file.path(locked, "leptomeningeal_genes_human.txt"))
paren <- readLines(file.path(locked, "parenchyma_genes_human.txt"))
genes <- unique(c(lepto[nzchar(lepto)], paren[nzchar(paren)]))

ann <- read.delim(
  gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
  stringsAsFactors = FALSE, check.names = FALSE
)
mb <- ann[ann$Sample_Origin == "mBrain", ]
cat("MBRAIN_CELLS", nrow(mb), "SAMPLES", paste(sort(unique(mb$Sample)), collapse = ","), "\n")

# GEO compressed this RDS twice. One gzip layer still leaves a gzip stream.
obj <- readRDS(gzcon(gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz", open = "rb")))
cat("CLASS", paste(class(obj), collapse = ","), "\n")
cat("DIM", paste(dim(obj), collapse = "x"), "\n")

symbols <- toupper(rownames(obj))
keep_cells <- intersect(mb$Index, colnames(obj))
cat("MATCHED_CELLS", length(keep_cells), "\n")
if (length(keep_cells) < 1000) stop("Cell names in the matrix do not match the annotation index")

# Subset columns before any full-matrix copy. 29,634 genes by 208,506 cells
# does not fit a second dense copy.
sub <- as.matrix(obj[, keep_cells, drop = FALSE])
rm(obj)
gc()
cat("SUBSET", paste(dim(sub), collapse = "x"), "STORAGE", typeof(sub), "\n")

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
cat("MALIGNANT", sum(malignant), "\n")

cell_sample <- mb$Sample[match(keep_cells, mb$Index)]
cell_author <- mb$Cell_type[match(keep_cells, mb$Index)]
samples <- sort(unique(cell_sample))

qc <- list()
rows <- list()
for (sample_id in samples) {
  use <- malignant & cell_sample == sample_id
  author_epi <- sum(use & cell_author == "Epithelial cells")
  n_mal <- sum(use)
  total_umi <- if (n_mal) sum(sub[, use, drop = FALSE]) else 0
  qc[[sample_id]] <- data.frame(
    sample = paste0("BrM_", sample_id),
    compartment = "Parenchyma",
    n_cells = sum(cell_sample == sample_id),
    n_malignant = n_mal,
    n_malignant_also_author_epithelial = author_epi,
    malignant_umi = as.numeric(total_umi)
  )
  cat(sample_id, "malignant", n_mal, "author_epi", author_epi, "umi", total_umi, "\n")
  if (n_mal < 1 || total_umi == 0) next
  present <- symbols %in% genes
  counts <- as.numeric(Matrix::rowSums(sub[present, use, drop = FALSE]))
  names(counts) <- symbols[present]
  # Collapse duplicate symbols by summing.
  counts <- tapply(counts, names(counts), sum)
  log2cpm <- log2(counts / as.numeric(total_umi) * 1e6 + 1)
  rec <- as.list(log2cpm)
  rec$sample <- paste0("BrM_", sample_id)
  rec$compartment <- "Parenchyma"
  rows[[sample_id]] <- rec
}

qc_df <- do.call(rbind, qc)
write.csv(qc_df, file.path(out_dir, "parenchyma_malignant_qc.csv"), row.names = FALSE)
# Bind rows that may miss genes by going through a character matrix of names.
all_names <- sort(unique(unlist(lapply(rows, names))))
all_names <- all_names[!all_names %in% c("sample", "compartment")]
out <- data.frame(sample = vapply(rows, function(x) x$sample, character(1)), stringsAsFactors = FALSE)
out$compartment <- "Parenchyma"
for (gene in all_names) {
  out[[gene]] <- vapply(rows, function(x) {
    value <- x[[gene]]
    if (is.null(value)) NA_real_ else as.numeric(value)
  }, numeric(1))
}
write.csv(out, file.path(out_dir, "parenchyma_pseudobulk_log2cpm.csv"), row.names = FALSE)
cat("PARENCHYMA_SAMPLES_WRITTEN", nrow(out), "\n")
