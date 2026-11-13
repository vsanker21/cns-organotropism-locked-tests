# Within-study control. Frozen gene lists, GSE131907 only.
# Origins are pleural effusion, brain metastasis, and lymph-node metastasis.

suppressPackageStartupMessages(library(Matrix))

locked <- "g:/Cancer Metastasis/results/locked"
out_dir <- "g:/Cancer Metastasis/results/compartment"
lepto <- readLines(file.path(locked, "leptomeningeal_genes_human.txt"))
paren <- readLines(file.path(locked, "parenchyma_genes_human.txt"))
prolif <- readLines(file.path(locked, "proliferation_genes.txt"))
genes <- unique(c(lepto[nzchar(lepto)], paren[nzchar(paren)], prolif[nzchar(prolif)]))

ann <- read.delim(
  gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
  stringsAsFactors = FALSE, check.names = FALSE
)
keep_origin <- c("PE", "mBrain", "mLN")
ann <- ann[ann$Sample_Origin %in% keep_origin, ]
cat("CELLS", nrow(ann), "BY", paste(sapply(keep_origin, function(x) sum(ann$Sample_Origin == x)), collapse = ","), "\n")

obj <- readRDS(gzcon(gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz", open = "rb")))
symbols <- toupper(rownames(obj))
cat("DIM", paste(dim(obj), collapse = "x"), "\n")

sum_genes <- function(mat, symbol_names, wanted) {
  rows <- which(symbol_names %in% wanted)
  if (!length(rows)) return(rep(0, ncol(mat)))
  as.numeric(Matrix::colSums(mat[rows, , drop = FALSE]))
}

qc <- list()
rows <- list()
for (origin in keep_origin) {
  part <- ann[ann$Sample_Origin == origin, ]
  keep_cells <- intersect(part$Index, colnames(obj))
  cat(origin, "MATCHED", length(keep_cells), "OF", nrow(part), "\n")
  if (length(keep_cells) < 100) stop("Too few matched cells for ", origin)
  sub <- as.matrix(obj[, keep_cells, drop = FALSE])
  epithelial <- sum_genes(sub, symbols, c("KRT8", "KRT18", "KRT19", "EPCAM", "NAPSA", "NKX2-1"))
  immune <- sum_genes(sub, symbols, c("PTPRC", "CD3D", "CD3E", "CD68", "AIF1", "C1QA"))
  tcell <- sum_genes(sub, symbols, c("CD3D", "CD3E"))
  macrophage <- sum_genes(sub, symbols, c("CD68", "AIF1", "C1QA"))
  malignant <- epithelial >= 2 & epithelial > immune & tcell == 0 & macrophage == 0
  cell_sample <- part$Sample[match(keep_cells, part$Index)]
  cell_author <- part$Cell_type[match(keep_cells, part$Index)]
  for (sample_id in sort(unique(part$Sample))) {
    use <- malignant & cell_sample == sample_id
    n_mal <- sum(use)
    author_epi <- sum(use & cell_author == "Epithelial cells")
    total_umi <- if (n_mal) sum(sub[, use, drop = FALSE]) else 0
    qc[[sample_id]] <- data.frame(
      sample = sample_id,
      origin = origin,
      n_cells = sum(cell_sample == sample_id),
      n_malignant = n_mal,
      n_malignant_also_author_epithelial = author_epi,
      malignant_umi = as.numeric(total_umi),
      stringsAsFactors = FALSE
    )
    cat(origin, sample_id, "malignant", n_mal, "author_epi", author_epi, "umi", total_umi, "\n")
    if (n_mal < 1 || total_umi == 0) next
    present <- symbols %in% genes
    counts <- as.numeric(Matrix::rowSums(sub[present, use, drop = FALSE]))
    names(counts) <- symbols[present]
    counts <- tapply(counts, names(counts), sum)
    log2cpm <- log2(as.numeric(counts) / as.numeric(total_umi) * 1e6 + 1)
    names(log2cpm) <- names(counts)
    rec <- as.list(log2cpm)
    rec$sample <- sample_id
    rec$origin <- origin
    rows[[sample_id]] <- rec
  }
  rm(sub)
  gc()
}
rm(obj)
gc()

qc_df <- do.call(rbind, qc)
write.csv(qc_df, file.path(out_dir, "gse131907_within_qc.csv"), row.names = FALSE)

all_names <- sort(unique(unlist(lapply(rows, function(x) setdiff(names(x), c("sample", "origin"))))))
out <- data.frame(
  sample = vapply(rows, function(x) x$sample, character(1)),
  origin = vapply(rows, function(x) x$origin, character(1)),
  stringsAsFactors = FALSE
)
for (gene in all_names) {
  out[[gene]] <- vapply(rows, function(x) {
    value <- x[[gene]]
    if (is.null(value)) NA_real_ else as.numeric(value)
  }, numeric(1))
}
write.csv(out, file.path(out_dir, "gse131907_within_log2cpm.csv"), row.names = FALSE)
cat("SAMPLES_WRITTEN", nrow(out), "GENES", length(all_names), "\n")
