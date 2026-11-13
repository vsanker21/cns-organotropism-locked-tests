suppressPackageStartupMessages(library(babelgene))
genes <- c("Gsta4", "Eno3", "Lama3", "Wisp1", "Krt80")
ortho <- babelgene::orthologs(genes, species = "mouse", human = FALSE, top = TRUE)
cat("CLASS", class(ortho), "\n")
cat("NROW", nrow(ortho), "\n")
cat("COLS", paste(colnames(ortho), collapse = ","), "\n")
print(utils::head(ortho, 10))
if (is.null(ortho) || nrow(ortho) == 0) {
  ortho2 <- babelgene::orthologs(genes, species = "mouse")
  cat("DEFAULT_NROW", nrow(ortho2), "\n")
  cat("DEFAULT_COLS", paste(colnames(ortho2), collapse = ","), "\n")
  print(utils::head(ortho2, 10))
}
