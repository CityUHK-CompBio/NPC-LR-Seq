
library(dplyr)
library(tidyr)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(matrixStats)

TPM_exp     <- "/path/to/short_read_to_long_read_transcriptome/isoform-count-matrix.txt"
sqanti_file <- "/path/to/output/sqanti3_all/SQANTI3_QC_ALL/classification.txt"
out_prefix  <- "path/to/output/pathway_enrichment"

NOVEL_CATS <- c("novel_in_catalog", "novel_not_in_catalog")

BG_MODE       <- "median"  
BG_TPM_CUTOFF <- 0.1
BG_PROP       <- 0

B      <- 200   
TOPN   <- 50   
N_BINS <- 10

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

TPM_matrix <- read.delim(TPM_exp, header = TRUE, row.names = 1, check.names = FALSE)
colnames(TPM_matrix) <- sub("\\.isoforms$", "", colnames(TPM_matrix))

sqanti3_qc <- read.delim(sqanti_file, header = TRUE, check.names = FALSE)

sqanti_use <- sqanti3_qc %>%   
  transmute(
    isoform = isoform,
    structural_category = structural_category,
    gene_ensembl = associated_gene
  ) %>%
  mutate(gene_ensembl = gsub("\\.\\d+", "", gene_ensembl)) %>%      
  distinct(isoform, .keep_all = TRUE)


ens <- unique(sqanti_use$gene_ensembl)
ens2sym <- bitr(ens, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db, drop = TRUE)
ens_map <- setNames(ens2sym$SYMBOL, ens2sym$ENSEMBL)

sqanti_use <- sqanti_use %>%
  mutate(gene_symbol = ens_map[gene_ensembl]) %>%
  filter(!is.na(gene_symbol), gene_symbol != "")


common <- intersect(rownames(TPM_matrix), sqanti_use$isoform)
length(common) 

TPM_matrix  <- TPM_matrix[common, , drop = FALSE]
sqanti_use  <- sqanti_use %>% filter(isoform %in% common)
dim(sqanti_use)

sqanti_use <- sqanti_use %>% arrange(match(sqanti_use$isoform, rownames(TPM_matrix)))
stopifnot(identical(sqanti_use$isoform, rownames(TPM_matrix)))

expr_mat <- as.matrix(apply(TPM_matrix, 2, as.numeric))
rownames(expr_mat) <- rownames(TPM_matrix)

gene_TPM <- rowsum(expr_mat, group = sqanti_use$gene_symbol, reorder = FALSE) 
dim(gene_TPM)

novel_genes <- sqanti_use %>%
  filter(structural_category %in% NOVEL_CATS) %>%
  pull(gene_symbol) %>%
  unique()

cat("Total genes (gene_TPM):", nrow(gene_TPM), "\n")
cat("Novel genes (before BG):", length(novel_genes), "\n")

if (BG_MODE == "median") {
  med <- matrixStats::rowMedians(gene_TPM, na.rm = TRUE)
  bg_genes <- names(med)[med >= BG_TPM_CUTOFF]
} else if (BG_MODE == "prop") {
  prop <- rowMeans(gene_TPM >= BG_TPM_CUTOFF, na.rm = TRUE)
  bg_genes <- names(prop)[prop >= BG_PROP]
} else {
  stop("BG_MODE must be 'median' or 'prop'")
}

novel_genes <- intersect(novel_genes, bg_genes)

cat("Expressed universe genes:", length(bg_genes), "\n")
cat("Novel genes after intersect(bg):", length(novel_genes), "\n\n")

make_bins <- function(expr_vec, n_bins = 10){
  qs <- quantile(expr_vec, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  qs <- unique(qs)
  if(length(qs) < 3) stop("Expression quantiles collapsed (too many ties). Try fewer bins.")
  cut(expr_vec, breaks = qs, include.lowest = TRUE)
}

draw_matched_set <- function(expr_df, novel_set, need_by_bin, seed = NULL, plot = TRUE){
  
  if(!is.null(seed)) set.seed(seed)
  
  pool_df <- expr_df %>% dplyr::filter(!(gene %in% novel_set))
  out <- c()
  
  for(b in names(need_by_bin)){
    k <- as.integer(need_by_bin[[b]])
    candidates <- pool_df$gene[pool_df$bin == b]
    
    if(length(candidates) < k){
      candidates <- pool_df$gene
    }
    
    out <- c(out, sample(candidates, k, replace = FALSE))
  }
  
  rand_set <- unique(out)
  
  if(plot){
    
    dev.off()
    novel_expr <- as.numeric(as.character(expr_df$medTPM[expr_df$gene %in% novel_set]))
    
    rand_expr  <- as.numeric(as.character(expr_df$medTPM[expr_df$gene %in% rand_set]))
    
    novel_expr <- na.omit(novel_expr)
    rand_expr  <- na.omit(rand_expr)

    plot(
      density(novel_expr),
      col="steelblue",
      lwd=2,
      main="Expression Distribution",
      xlab="Expression",
      xlim=c(0,1000),
      ylim=c(0, 0.1) 
    )
    
    lines(
      density(rand_expr),
      col="orange",
      lwd=2
    )
    
    legend(
      "topright",
      legend=c("Novel set","Random set"),
      col=c("steelblue","orange"),
      lwd=2
    )
    
    }
  return(rand_set)
  
  }
  

ora_with_empirical <- function(label, term2gene, gene_list, bg_genes, gene_TPM,
                               B = 1000, topN = 50, n_bins = 10, out_prefix = "ORA"){
  universe <- intersect(bg_genes, unique(term2gene$gene_symbol))
  length(universe)

  gene_use <- intersect(gene_list, universe)
  length(gene_use)
  
  if(length(gene_use) < 10){
    warning(label, ": too few genes after intersect with universe.")
    return(NULL)
  }
  
  expr_med <- matrixStats::rowMedians(gene_TPM[universe, , drop = FALSE], na.rm = TRUE)
  
  expr_df <- data.frame(gene = universe, medTPM = expr_med) %>% filter(!is.na(medTPM))
  dim(expr_df)
  expr_df$bin <- make_bins(expr_df$medTPM, n_bins = n_bins)
  
  novel_df <- expr_df %>% filter(gene %in% gene_use)
  dim(novel_df)
  need_by_bin <- table(novel_df$bin)
  
  real_df <- enricher(gene = gene_use, TERM2GENE = term2gene, universe = universe,
                      pvalueCutoff = 1, qvalueCutoff = 1) %>% as.data.frame()
  
  if(nrow(real_df) == 0){
    warning(label, ": real ORA returns empty.")
    return(NULL)
  }
  
  real_df <- real_df %>% arrange(p.adjust)
  real_top <- real_df %>% slice_head(n = min(topN, nrow(real_df))) %>%
    dplyr::select(ID, Description, pvalue, p.adjust, Count)
  
  rand_p <- matrix(NA_real_, nrow = nrow(real_top), ncol = B,
                   dimnames = list(real_top$ID, paste0("rep", seq_len(B))))
  
  for(i in seq_len(B)){
    rand_set <- draw_matched_set(expr_df, gene_use, need_by_bin, seed = 1000 + i)
    length(rand_set) 
    rand_res <- enricher(gene = rand_set, TERM2GENE = term2gene, universe = universe,
                         pvalueCutoff = 1, qvalueCutoff = 1) %>% as.data.frame()
    if(nrow(rand_res) > 0){
      m <- setNames(rand_res$pvalue, rand_res$ID)
      rand_p[, i] <- m[real_top$ID]
    }
  }
  
  real_map <- setNames(real_top$pvalue, real_top$ID)
  empirical_p <- sapply(real_top$ID, function(pid){
    rp <- rand_p[pid, ]
    rp <- rp[!is.na(rp)]
    if(length(rp) == 0) return(NA_real_)
    pr <- real_map[[pid]]
    (sum(rp <= pr) + 1) / (length(rp) + 1)
  })
  
  out_check <- real_top %>%
    mutate(empirical_p = as.numeric(empirical_p),
           rand_N_used = sapply(real_top$ID, function(pid) sum(!is.na(rand_p[pid, ]))),
           label = label)
  
  write.csv(real_df,  paste0(out_prefix, "_", label, "_ORA_full.csv"), row.names = FALSE)
  write.csv(out_check, paste0(out_prefix, "_", label, "_ORA_top", topN, "_empirical.csv"), row.names = FALSE)
  
  cat(label, ": universe=", length(universe),
      " gene=", length(gene_use),
      " real_terms=", nrow(real_df),
      " wrote CSV\n")
  
  invisible(list(full = real_df, top = out_check))
}

kegg_with_empirical <- function(label, gene_list, bg_genes, gene_TPM,
                                B = 1000, topN = 50, n_bins = 10, out_prefix = "ORA"){
  
  sym_all <- unique(c(gene_list, bg_genes))
  sym2ent <- bitr(sym_all, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  sym2ent_map <- split(sym2ent$ENTREZID, sym2ent$SYMBOL)
  
  universe_sym <- bg_genes
  gene_use_sym <- intersect(gene_list, universe_sym)
  
  expr_med <- matrixStats::rowMedians(gene_TPM[universe_sym, , drop = FALSE], na.rm = TRUE)
  expr_df <- data.frame(gene = universe_sym, medTPM = expr_med) %>% filter(!is.na(medTPM))
  expr_df$bin <- make_bins(expr_df$medTPM, n_bins = n_bins)
  
  novel_df <- expr_df %>% filter(gene %in% gene_use_sym)
  need_by_bin <- table(novel_df$bin)
  
  gene_use_ent <- unique(unlist(sym2ent_map[gene_use_sym]))
  bg_ent      <- unique(unlist(sym2ent_map[universe_sym]))
  
  real_df <- enrichKEGG(gene = gene_use_ent, organism = "hsa", universe = bg_ent,
                        pvalueCutoff = 1) %>% as.data.frame()
  
  if(nrow(real_df) == 0){
    warning(label, ": real KEGG returns empty.")
    return(NULL)
  }
  
  real_df <- real_df %>% arrange(p.adjust)
  real_top <- real_df %>% slice_head(n = min(topN, nrow(real_df))) %>%dplyr::select(ID, Description, pvalue, p.adjust, Count)
  
  rand_p <- matrix(NA_real_, nrow = nrow(real_top), ncol = B,
                   dimnames = list(real_top$ID, paste0("rep", seq_len(B))))
  
  for(i in seq_len(B)){
    rand_set_sym <- draw_matched_set(expr_df, gene_use_sym, need_by_bin, seed = 2000 + i)
    rand_set_ent <- unique(unlist(sym2ent_map[rand_set_sym]))
    rand_res <- enrichKEGG(gene = rand_set_ent, organism = "hsa", universe = bg_ent,
                           pvalueCutoff = 1) %>% as.data.frame()
    if(nrow(rand_res) > 0){
      m <- setNames(rand_res$pvalue, rand_res$ID)
      rand_p[, i] <- m[real_top$ID]
    }
  }
  
  real_map <- setNames(real_top$pvalue, real_top$ID)
  empirical_p <- sapply(real_top$ID, function(pid){
    rp <- rand_p[pid, ]
    rp <- rp[!is.na(rp)]
    if(length(rp) == 0) return(NA_real_)
    pr <- real_map[[pid]]
    (sum(rp <= pr) + 1) / (length(rp) + 1)
  })
  
  out_check <- real_top %>%
    mutate(empirical_p = as.numeric(empirical_p),
           rand_N_used = sapply(real_top$ID, function(pid) sum(!is.na(rand_p[pid, ]))),
           label = label)
  
  write.csv(real_df,  paste0(out_prefix, "_", label, "_ORA_full.csv"), row.names = FALSE)
  write.csv(out_check, paste0(out_prefix, "_", label, "_ORA_top", topN, "_empirical.csv"), row.names = FALSE)
  
  cat(label, ": SYMBOL_universe=", length(universe_sym),
      " SYMBOL_gene=", length(gene_use_sym),
      " ENTREZ_gene=", length(gene_use_ent),
      " real_terms=", nrow(real_df),
      " wrote CSV\n")
  
  invisible(list(full = real_df, top = out_check))
}

m_h <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::rename(ID = gs_name)

m_c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::rename(ID = gs_name)

m_c5bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::rename(ID = gs_name)

m_go <- bind_rows(m_c2, m_c5bp)

to_TERM2GENE <- function(df){
  df %>% transmute(term = ID, gene_symbol=gene_symbol)
}

res_hallmark <- ora_with_empirical(
  label     = "Hallmark",
  term2gene = to_TERM2GENE(m_h),
  gene_list = novel_genes,
  bg_genes  = bg_genes,
  gene_TPM  = gene_TPM,
  B         = B,
  topN      = TOPN,
  n_bins    = N_BINS,
  out_prefix = out_prefix
)

res_GO <- ora_with_empirical(
  label     = "GO_C2_C5BP",
  term2gene = to_TERM2GENE(m_go),
  gene_list = novel_genes,
  bg_genes  = bg_genes,
  gene_TPM  = gene_TPM,
  B         = B,
  topN      = TOPN,
  n_bins    = N_BINS,
  out_prefix = out_prefix
)

res_KEGG <- kegg_with_empirical(
  label      = "KEGG",
  gene_list  = novel_genes,
  bg_genes   = bg_genes,
  gene_TPM   = gene_TPM,
  B          = B,
  topN       = TOPN,
  n_bins     = N_BINS,
  out_prefix = out_prefix  
)

cat("\nAll done.\n",
    "Outputs:\n",
    out_prefix, "_Hallmark_ORA_full.csv\n", 
    out_prefix, "_Hallmark_ORA_top", TOPN, "_empirical.csv\n",
    out_prefix, "_GO_C2_C5BP_ORA_full.csv\n",
    out_prefix, "_GO_C2_C5BP_ORA_top", TOPN, "_empirical.csv\n",
    out_prefix, "_KEGG_ORA_full.csv\n",
    out_prefix, "_KEGG_ORA_top", TOPN, "_empirical.csv\n", sep="")
