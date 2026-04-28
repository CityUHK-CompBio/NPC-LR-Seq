# Package loader for tumor specificity analysis modules

load_tumor_specificity_packages <- function() {
  pkgs <- c(
    "dplyr", "parallel", "rtracklayer", "svglite", "tidyr", "treeio", "seqinr",
    "stringr", "msigdbr", "tidyverse", "reshape2", "clusterProfiler", "plyr",
    "matrixStats", "ggplot2", "ggpubr", "ggstatsplot", "ggsci", "pheatmap",
    "AnnotationDbi", "org.Hs.eg.db", "DESeq2", "EnhancedVolcano", "fastmatch"
  )
 
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required packages: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  invisible(lapply(pkgs, library, character.only = TRUE))
}
