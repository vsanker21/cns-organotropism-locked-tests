# Export one raw count matrix per sample for the ambient-RNA control.
# Clusters are assigned here. Expression is not tested here.

suppressPackageStartupMessages(library(Matrix))

out_dir <- "g:/Cancer Metastasis/results/macrophage/ambient_input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
clinical <- read.csv("g:/Cancer Metastasis/results/locked/gse131907_momac_samples.csv", stringsAsFactors = FALSE)

ann <- read.delim(
  gzfile("g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
  stringsAsFactors = FALSE, check.names = FALSE
)
label <- ann$Cell_subtype
missing <- is.na(label) | !nzchar(label)
label[missing] <- ann$Cell_type[missing]
ann$ambient_label <- label

obj <- readRDS(gzcon(gzfile(
  "g:/Cancer Metastasis/data/geo/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz",
  open = "rb"
)))
symbols <- toupper(rownames(obj))
writeLines(symbols, file.path(out_dir, "genes.txt"))

for (sample_id in clinical$sample) {
  rows <- which(ann$Sample == sample_id)
  keep <- intersect(ann$Index[rows], colnames(obj))
  sub <- obj[, keep, drop = FALSE]
  sub <- rowsum(as.matrix(sub), symbols)
  mode(sub) <- "integer"
  present <- rowSums(sub) > 0
  sub <- sub[present, , drop = FALSE]
  lab <- ann$ambient_label[match(keep, ann$Index)]
  freq <- table(lab)
  lab[lab != "mo-Mac" & freq[lab] < 10] <- "Rare"
  sm <- Matrix(sub, sparse = TRUE)
  mtx_path <- file.path(out_dir, paste0(sample_id, ".mtx"))
  writeMM(sm, mtx_path)
  write.table(
    data.frame(cell = keep, label = lab, momac = lab == "mo-Mac"),
    file.path(out_dir, paste0(sample_id, ".labels.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  writeLines(rownames(sub), file.path(out_dir, paste0(sample_id, ".genes.txt")))
  cat(sample_id, "cells", ncol(sub), "genes", nrow(sub), "momac", sum(lab == "mo-Mac"), "clusters", length(unique(lab)), "\n")
  rm(sub, sm)
  gc()
}
rm(obj)
cat("EXPORT_DONE\n")
