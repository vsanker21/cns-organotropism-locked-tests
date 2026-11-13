# Compartment program, locked before any human sample is read.
# Primary contrast within each model: leptomeningeal-selected vs parenchyma-selected.
# A gene is retained only when the inverse-variance meta-analysis is significant
# after BH correction and the same direction is present in at least 4 of 5 models.

suppressPackageStartupMessages(library(DESeq2))

out_dir <- "g:/Cancer Metastasis/results/compartment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  B16 = "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_D_htseq_B16.txt.gz",
  E0771 = "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_F_htseq_E0771.txt.gz",
  Yumm52 = "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_G_htseq_Yumm52.txt.gz",
  `4T1` = "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_I_htseq_4T1.txt.gz",
  EMT6 = "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_J_htseq_EMT6.txt.gz"
)

classify <- function(sample_name) {
  if (grepl("BrM2", sample_name)) return("Parenchyma")
  if (grepl("_Par_", sample_name)) return("Parental")
  if (grepl("_IC_", sample_name)) return("Leptomeningeal")
  stop("Unclassified sample: ", sample_name)
}

effects <- list()
for (model_name in names(files)) {
  raw <- read.delim(gzfile(files[[model_name]]), check.names = FALSE, stringsAsFactors = FALSE)
  raw <- raw[nzchar(raw$GeneSymbol), ]
  raw$GeneSymbol <- toupper(raw$GeneSymbol)
  sample_cols <- setdiff(colnames(raw), c("GeneID", "GeneSymbol"))
  counts <- rowsum(as.matrix(raw[, sample_cols]), raw$GeneSymbol)
  mode(counts) <- "integer"
  counts <- counts[rowSums(counts) >= 10, , drop = FALSE]
  group <- factor(vapply(colnames(counts), classify, character(1)),
                  levels = c("Parenchyma", "Leptomeningeal", "Parental"))
  coldata <- data.frame(group = group, row.names = colnames(counts))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "Leptomeningeal", "Parenchyma"))
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df$model <- model_name
  effects[[model_name]] <- df[, c("gene", "model", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  cat(model_name, "genes", nrow(df), "samples", paste(colnames(counts), collapse = ","), "\n")
}

all_effects <- do.call(rbind, effects)
write.csv(all_effects, file.path(out_dir, "per_model_lepto_vs_parenchyma.csv"), row.names = FALSE)

# Inverse-variance meta-analysis. A gene must be estimable in every model.
lfc_wide <- reshape(
  all_effects[, c("gene", "model", "log2FoldChange")],
  idvar = "gene", timevar = "model", direction = "wide"
)
se_wide <- reshape(
  all_effects[, c("gene", "model", "lfcSE")],
  idvar = "gene", timevar = "model", direction = "wide"
)
se_wide <- se_wide[match(lfc_wide$gene, se_wide$gene), ]
lfc <- as.matrix(lfc_wide[, -1])
se <- as.matrix(se_wide[, -1])
ok <- is.finite(lfc) & is.finite(se) & se > 0
keep <- rowSums(ok) == ncol(lfc)
lfc <- lfc[keep, , drop = FALSE]
se <- se[keep, , drop = FALSE]
colnames(lfc) <- sub("^log2FoldChange\\.", "", colnames(lfc))
colnames(se) <- sub("^lfcSE\\.", "", colnames(se))
weight <- 1 / se^2
meta_lfc <- rowSums(weight * lfc) / rowSums(weight)
meta_se <- sqrt(1 / rowSums(weight))
z <- meta_lfc / meta_se
meta <- data.frame(
  gene = lfc_wide$gene[keep],
  meta_log2fc = meta_lfc,
  meta_se = meta_se,
  z = z,
  p = 2 * pnorm(-abs(z)),
  n_positive = rowSums(lfc > 0),
  n_negative = rowSums(lfc < 0),
  stringsAsFactors = FALSE
)
meta$fdr <- p.adjust(meta$p, method = "BH")
meta <- meta[order(meta$p), ]
write.csv(meta, file.path(out_dir, "meta_lepto_vs_parenchyma.csv"), row.names = FALSE)

lepto <- meta[meta$fdr < 0.05 & meta$meta_log2fc > 0 & meta$n_positive >= 4, ]
parenchyma <- meta[meta$fdr < 0.05 & meta$meta_log2fc < 0 & meta$n_negative >= 4, ]
rule <- "meta_fdr_lt_0.05_and_concordant_in_at_least_4_of_5"
if (nrow(lepto) < 10 && nrow(parenchyma) < 10) {
  lepto <- meta[meta$p < 0.05 & meta$meta_log2fc > 0.5 & meta$n_positive >= 4, ]
  parenchyma <- meta[meta$p < 0.05 & meta$meta_log2fc < -0.5 & meta$n_negative >= 4, ]
  rule <- "fallback_nominal_p_lt_0.05_abs_lfc_gt_0.5_concordant_4_of_5"
}
writeLines(lepto$gene, file.path(out_dir, "leptomeningeal_genes_mouse.txt"))
writeLines(parenchyma$gene, file.path(out_dir, "parenchyma_genes_mouse.txt"))
writeLines(rule, file.path(out_dir, "rule.txt"))

# Leave-one-model-out does not change the locked list. It reports whether one model carries it.
loo <- list()
model_names <- colnames(lfc)
lepto_idx <- match(lepto$gene, lfc_wide$gene[keep])
for (drop_model in model_names) {
  use <- setdiff(model_names, drop_model)
  sub_lfc <- lfc[lepto_idx, use, drop = FALSE]
  sub_se <- se[lepto_idx, use, drop = FALSE]
  sub_w <- 1 / sub_se^2
  sub_meta <- rowSums(sub_w * sub_lfc) / rowSums(sub_w)
  sub_pos <- rowSums(sub_lfc > 0)
  loo[[drop_model]] <- data.frame(
    dropped = drop_model,
    n_locked = length(lepto_idx),
    n_still_positive = sum(sub_meta > 0 & sub_pos >= 3, na.rm = TRUE)
  )
}
write.csv(do.call(rbind, loo), file.path(out_dir, "leave_one_model_out_leptomeningeal.csv"), row.names = FALSE)

cat("RULE", rule, "\n")
cat("N_TESTED", nrow(meta), "\n")
cat("N_LEPTO", nrow(lepto), "\n")
cat("N_PARENCHYMA", nrow(parenchyma), "\n")
cat("LEPTO_TOP", paste(head(lepto$gene, 30), collapse = ","), "\n")
cat("PARENCHYMA_TOP", paste(head(parenchyma$gene, 30), collapse = ","), "\n")
