# Lock signatures before any patient cohort is scored.
# Step A comes only from GSE138122: MDA-MB-231 brain-tropic (BrM2) versus
# lung-tropic (LM2) cells profiled in culture. Patient datasets are not read.

suppressPackageStartupMessages({
  library(DESeq2)
  library(msigdbr)
  library(babelgene)
})

out_dir <- "g:/Cancer Metastasis/results/locked"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- "g:/Cancer Metastasis/data/geo/GSE138122/suppl/GSE138122_STAR.counts.geo.txt.gz"
raw <- read.delim(gzfile(counts_path), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("ensembl", "hgnc_symbol", "brm2-1", "lm2-1", "par-1") %in% colnames(raw)))

sample_cols <- c(
  "brm2-1", "brm2-2", "brm2-3",
  "lm2-1", "lm2-2", "lm2-3",
  "par-1", "par-2", "par-3"
)
raw <- raw[nzchar(raw$hgnc_symbol), ]
raw$hgnc_symbol <- toupper(raw$hgnc_symbol)
# Sum Ensembl rows that share a symbol so inference is at gene level.
agg <- rowsum(as.matrix(raw[, sample_cols]), group = raw$hgnc_symbol)
mode(agg) <- "integer"
agg <- agg[rowSums(agg) >= 10, , drop = FALSE]

coldata <- data.frame(
  sample = colnames(agg),
  group = factor(
    c("BrM", "BrM", "BrM", "LM", "LM", "LM", "Par", "Par", "Par"),
    levels = c("LM", "BrM", "Par")
  ),
  row.names = colnames(agg)
)

dds <- DESeqDataSetFromMatrix(countData = agg, colData = coldata, design = ~ group)
dds <- DESeq(dds, quiet = TRUE)

# Primary definition: genes whose BrM-versus-LM fold change exceeds 2
# (lfcThreshold = 1 on the log2 scale), alternative = greater.
res_strict <- results(
  dds,
  contrast = c("group", "BrM", "LM"),
  lfcThreshold = 1,
  altHypothesis = "greater",
  alpha = 0.05
)
res_strict_df <- as.data.frame(res_strict)
res_strict_df$gene <- rownames(res_strict_df)
res_wald <- as.data.frame(results(dds, contrast = c("group", "BrM", "LM")))
res_wald$gene <- rownames(res_wald)
res_par <- as.data.frame(results(dds, contrast = c("group", "BrM", "Par")))
res_par$gene <- rownames(res_par)

primary <- res_strict_df[!is.na(res_strict_df$padj) & res_strict_df$padj < 0.05, ]
rule <- "lfcThreshold_gt1_padj_lt_0.05"
if (nrow(primary) < 10) {
  eligible <- res_wald[!is.na(res_wald$stat) & res_wald$log2FoldChange > 0, ]
  eligible <- eligible[order(-eligible$stat), ]
  primary <- head(eligible, 50)
  rule <- "fallback_top50_positive_wald_stat"
}

write.csv(res_strict_df, file.path(out_dir, "gse138122_brm_vs_lm_lfcThreshold1.csv"), row.names = FALSE)
write.csv(res_wald, file.path(out_dir, "gse138122_brm_vs_lm_wald.csv"), row.names = FALSE)
write.csv(res_par, file.path(out_dir, "gse138122_brm_vs_par_wald.csv"), row.names = FALSE)
writeLines(primary$gene, file.path(out_dir, "stepA_genes.txt"))
writeLines(rule, file.path(out_dir, "stepA_rule.txt"))

# Public benchmark sets. The choice rule is name-based and does not use Step A.
c4 <- msigdbr(species = "Homo sapiens", collection = "C4", subcollection = "3CA")
c4_names <- sort(unique(c4$gs_name))
writeLines(c4_names, file.path(out_dir, "c4_3ca_names.txt"))
neural_names <- c4_names[grepl("neural progenitor", c4_names, ignore.case = TRUE)]
if (length(neural_names) == 1) {
  neural_name <- neural_names
} else if (length(neural_names) > 1) {
  ranked <- neural_names[order(nchar(neural_names))]
  neural_name <- ranked[1]
} else {
  neural_name <- "GOBP_SYNAPTIC_SIGNALING"
}

if (neural_name == "GOBP_SYNAPTIC_SIGNALING") {
  go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  neural_genes <- unique(go$gene_symbol[go$gs_name == neural_name])
} else {
  neural_genes <- unique(c4$gene_symbol[c4$gs_name == neural_name])
}
writeLines(neural_genes, file.path(out_dir, "neural_genes.txt"))
writeLines(neural_name, file.path(out_dir, "neural_set_name.txt"))

c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
bos_names <- sort(unique(c2$gs_name[grepl("BOS", c2$gs_name) & grepl("BRAIN", c2$gs_name)]))
writeLines(bos_names, file.path(out_dir, "bos_candidate_names.txt"))
if (length(bos_names) >= 1) {
  # The up-regulated set is the benchmark. Prefer a name ending in UP.
  up <- bos_names[grepl("_UP$", bos_names)]
  bos_name <- if (length(up) >= 1) up[1] else bos_names[1]
  bos_genes <- unique(c2$gene_symbol[c2$gs_name == bos_name])
} else {
  bos_name <- "PARTIAL_BOS_THREE_GENES_FROM_GSE12237_RECORD"
  bos_genes <- c("PTGS2", "HBEGF", "ST6GALNAC5")
}
writeLines(bos_genes, file.path(out_dir, "bos_genes.txt"))
writeLines(bos_name, file.path(out_dir, "bos_set_name.txt"))

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
prolif <- unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_G2M_CHECKPOINT"])
writeLines(prolif, file.path(out_dir, "proliferation_genes.txt"))

ortho <- babelgene::orthologs(primary$gene, species = "mouse")
# One mouse symbol per human gene: keep the row babelgene ranks first.
ortho <- ortho[!duplicated(ortho$human_symbol), c("human_symbol", "symbol")]
colnames(ortho) <- c("human_symbol", "mouse_symbol")
write.csv(ortho, file.path(out_dir, "stepA_mouse_orthologs.csv"), row.names = FALSE)

cat("RULE", rule, "\n")
cat("N_STEP_A", nrow(primary), "\n")
cat("NEURAL_SET", neural_name, "N", length(neural_genes), "\n")
cat("BOS_SET", bos_name, "N", length(bos_genes), "\n")
cat("PROLIF_N", length(prolif), "\n")
cat("ORTHOLOGS", nrow(ortho), "\n")
cat("STEP_A_GENES", paste(primary$gene, collapse = ","), "\n")
