suppressPackageStartupMessages(library(msigdbr))
c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
nm <- sort(unique(c2$gs_name))
hit <- nm[grepl("BRAIN|SMID|MASSAGUE|BOS", nm)]
writeLines(hit, "g:/Cancer Metastasis/results/locked/bos_name_search.txt")
cat("HITS", length(hit), "\n")
