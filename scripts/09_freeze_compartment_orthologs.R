# Map the locked mouse compartment genes to human symbols before any human matrix is read.
# One human gene keeps the mouse gene with the larger absolute meta-analytic log2 fold change.

suppressPackageStartupMessages(library(babelgene))

meta <- read.csv("g:/Cancer Metastasis/results/compartment/meta_lepto_vs_parenchyma.csv", stringsAsFactors = FALSE)
lepto <- readLines("g:/Cancer Metastasis/results/compartment/leptomeningeal_genes_mouse.txt")
paren <- readLines("g:/Cancer Metastasis/results/compartment/parenchyma_genes_mouse.txt")
lepto <- lepto[nzchar(lepto)]
paren <- paren[nzchar(paren)]

# DESeq2 used uppercase symbols. babelgene expects the original mouse capitalization.
htseq <- c(
  "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_D_htseq_B16.txt.gz",
  "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_F_htseq_E0771.txt.gz",
  "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_G_htseq_Yumm52.txt.gz",
  "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_I_htseq_4T1.txt.gz",
  "g:/Cancer Metastasis/data/geo/GSE221358/suppl/GSE221358_Proj_09317_J_htseq_EMT6.txt.gz"
)
original <- character()
for (path in htseq) {
  symbols <- read.delim(gzfile(path), stringsAsFactors = FALSE)$GeneSymbol
  key <- toupper(symbols)
  original[key] <- symbols
}

query_mouse <- function(symbols) {
  symbols <- unique(symbols[nzchar(symbols)])
  if (!length(symbols)) return(data.frame())
  chunks <- split(symbols, ceiling(seq_along(symbols) / 200))
  rows <- lapply(chunks, function(chunk) {
    hit <- tryCatch(
      babelgene::orthologs(chunk, species = "mouse", human = FALSE, top = TRUE),
      error = function(e) NULL
    )
    if (is.null(hit)) {
      one <- lapply(chunk, function(gene) {
        tryCatch(
          babelgene::orthologs(gene, species = "mouse", human = FALSE, top = TRUE),
          error = function(e) NULL
        )
      })
      hit <- do.call(rbind, one)
    }
    hit
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab) || !nrow(tab)) {
    return(data.frame(mouse_symbol = character(), human_symbol = character(), support_n = integer()))
  }
  tab[, c("symbol", "human_symbol", "support_n")]
}

# Second pass: some current human symbols match the uppercase mouse symbol
# even when the mouse-side query misses them (Gsta4 is the example).
query_as_human <- function(upper_symbols) {
  upper_symbols <- unique(upper_symbols[nzchar(upper_symbols)])
  if (!length(upper_symbols)) return(data.frame())
  chunks <- split(upper_symbols, ceiling(seq_along(upper_symbols) / 200))
  rows <- lapply(chunks, function(chunk) {
    hit <- tryCatch(
      babelgene::orthologs(chunk, species = "mouse", top = TRUE),
      error = function(e) NULL
    )
    if (is.null(hit)) return(NULL)
    hit <- hit[toupper(hit$symbol) == hit$human_symbol, ]
    hit
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab) || !nrow(tab)) {
    return(data.frame(mouse_symbol = character(), human_symbol = character(), support_n = integer()))
  }
  data.frame(
    symbol = original[tab$human_symbol],
    human_symbol = tab$human_symbol,
    support_n = tab$support_n,
    stringsAsFactors = FALSE
  )
}

map_one <- function(genes, direction) {
  proper <- original[genes]
  missing_case <- genes[is.na(proper)]
  if (length(missing_case)) {
    cat("NO_ORIGINAL_SYMBOL", paste(head(missing_case, 10), collapse = ","), "\n")
  }
  proper <- proper[!is.na(proper)]
  ortho <- query_mouse(proper)
  found_upper <- toupper(ortho$symbol)
  missed <- genes[!genes %in% found_upper]
  extra <- query_as_human(missed)
  if (nrow(extra)) {
    names(extra)[1] <- "symbol"
    ortho <- rbind(ortho[, c("symbol", "human_symbol", "support_n")], extra[, c("symbol", "human_symbol", "support_n")])
  }
  ortho$mouse_key <- toupper(ortho$symbol)
  sub <- meta[meta$gene %in% genes, ]
  sub <- sub[order(-abs(sub$meta_log2fc)), ]
  ortho <- ortho[order(match(ortho$mouse_key, sub$gene)), ]
  ortho <- ortho[!duplicated(ortho$human_symbol), ]
  data.frame(
    mouse_symbol = ortho$mouse_key,
    mouse_symbol_original = ortho$symbol,
    human_symbol = ortho$human_symbol,
    support_n = ortho$support_n,
    direction = direction,
    stringsAsFactors = FALSE
  )
}

lepto_map <- map_one(lepto, "leptomeningeal_high")
paren_map <- map_one(paren, "parenchyma_high")
both <- rbind(lepto_map, paren_map)
# A human gene cannot sit in both programs. Keep the direction with the larger mouse effect.
effect <- setNames(abs(meta$meta_log2fc), meta$gene)
both$abs_effect <- effect[both$mouse_symbol]
both <- both[order(-both$abs_effect), ]
both <- both[!duplicated(both$human_symbol), ]

locked <- "g:/Cancer Metastasis/results/locked"
write.csv(both, file.path(locked, "compartment_human_orthologs.csv"), row.names = FALSE)
writeLines(both$human_symbol[both$direction == "leptomeningeal_high"], file.path(locked, "leptomeningeal_genes_human.txt"))
writeLines(both$human_symbol[both$direction == "parenchyma_high"], file.path(locked, "parenchyma_genes_human.txt"))
writeLines(
  "Frozen before human expression was scored. Mouse rule: inverse-variance meta-analysis of per-model DESeq2 log2 fold changes, leptomeningeal-selected versus parenchyma-selected, BH FDR < 0.05, and the same sign in at least 4 of 5 models. One human symbol per gene; the mouse ortholog with the larger absolute meta-analytic log2 fold change is kept. A human symbol present on both sides is assigned to the side with the larger absolute effect.",
  file.path(locked, "compartment_rule.txt")
)

qs <- quantile(abs(meta$meta_log2fc[meta$gene %in% c(lepto, paren)]), probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
cat("N_MOUSE_LEPTO", length(lepto), "\n")
cat("N_MOUSE_PAREN", length(paren), "\n")
cat("N_HUMAN_LEPTO", sum(both$direction == "leptomeningeal_high"), "\n")
cat("N_HUMAN_PAREN", sum(both$direction == "parenchyma_high"), "\n")
cat("ABS_LFC_QUANTILES", paste(names(qs), round(qs, 3), sep = "=", collapse = ","), "\n")
cat("LEPTO_MEDIAN_LFC", round(median(meta$meta_log2fc[meta$gene %in% lepto]), 3), "\n")
cat("PAREN_MEDIAN_LFC", round(median(meta$meta_log2fc[meta$gene %in% paren]), 3), "\n")
