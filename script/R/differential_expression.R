# Differential transcript expression (DTE) functions

library(DESeq2)
library(fastmatch)
library(dplyr)
library(pheatmap)
library(svglite)
library(EnhancedVolcano)

df_cal<-function(count_matrix,colData){
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix,0),
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, independentFiltering = TRUE)
  res$structural_category<-sqanti3_qc[match(rownames(res),sqanti3_qc$isoform),]$structural_category
  res<-res[order(abs(res$log2FoldChange),decreasing = T),]
  res$structural_category<-sqanti3_qc[match(rownames(res),sqanti3_qc$isoform),]$structural_category
  res$predicted_NMD<-sqanti3_qc[match(rownames(res),sqanti3_qc$isoform),]$predicted_NMD
  res$associated_gene<-sqanti3_qc[match(rownames(res),sqanti3_qc$isoform),]$associated_gene
  res$associated_transcript<-sqanti3_qc[match(rownames(res),sqanti3_qc$isoform),]$associated_transcript
  idx <- fmatch(res$associated_gene, reference_gtf$gene_id)
  gene_symbol <- reference_gtf$gene_name[idx]
  res$gene_symbol<-gene_symbol
  
  res$coding <-sqanti3_qc[match(rownames(res),sqanti3_qc$isoform),]$coding
  res<-res[order(res$log2FoldChange),]
  
  sig_genes <- subset(res, padj <= 0.05 & abs(log2FoldChange) >= 2)
  return(sig_genes)
  
}
