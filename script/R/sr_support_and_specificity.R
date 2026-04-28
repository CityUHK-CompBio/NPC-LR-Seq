# Short-read expression distribution
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
library(matrixStats)

sr_support<-function(){
  TPM_matrix<-read.table(TPM_exp,header=T,row.names = 1)
  TPM_matrix[TPM_matrix<1]<-0
  sub_TPM_matrix<-TPM_matrix
  sub_TPM_matrix$structural_category <-sqanti3_qc[match(rownames(TPM_matrix),sqanti3_qc$isoform),]$structural_category
  sub_TPM_matrix$structural_category[which(sub_TPM_matrix$structural_category %in% c("antisense","genic","intergenic","fusion"))]="others"
  sub_TPM_matrix$median<-sub_TPM_matrix %>% dplyr::select(contains("isoforms")) %>% as.matrix() %>% rowMedians()
  
  count_by_median_threshold <- function(df, threshold) {
    df %>%
      group_by(structural_category) %>%
      summarize(count = sum(median > threshold)) %>%
      mutate(Threshold = paste("Median >", threshold))
  }
  

  count_data <- bind_rows(
    count_by_median_threshold(sub_TPM_matrix, 1),
    count_by_median_threshold(sub_TPM_matrix, 5),
    count_by_median_threshold(sub_TPM_matrix, 10)
  )
  
  count_data<-count_data %>%filter(structural_category %in% c("full-splice_match","novel_in_catalog","novel_not_in_catalog"))
  count_data$Threshold<-factor(count_data$Threshold,levels = c("Median > 1", "Median > 5", "Median > 10"))
  svglite("/path/to/output/figures/Medians_Expression.svg", width = 10, height = 10)
  p<-ggplot(count_data, aes(x = Threshold, y = count, fill = structural_category)) +
    geom_bar(stat = "identity", position = "dodge",width = 0.8,color="black") +
    labs(y = "Number of Transcripts", x = "", fill = "") +
    theme_minimal() +
    theme_bw()+theme(legend.position = "top",legend.title= element_blank(),axis.text = element_text(size = 30,color="black"),axis.title = element_text(size = 30,color="black"),legend.text = element_text(size = 30,color="black"),plot.title = element_text(size = 30))
  print(p)
  dev.off()
  
  
  calculate_non_zero_ratio <- function(row) {
    sum(row!= 0, na.rm = TRUE) / sum(!is.na(row), na.rm = TRUE)
  }
  
  sr_detect<-data.frame(structural_category=c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog","others"))
  for (sample in grep("isoforms", names(sub_TPM_matrix), value = TRUE)){
    print(sample)
    detect_ratio<-sub_TPM_matrix %>%
    group_by(structural_category) %>% 
    dplyr::summarize(across(starts_with(sample),calculate_non_zero_ratio, .names = "{.col}"))
  
  sr_detect<-cbind(sr_detect,detect_ratio[,2])}
  
  sr_detect<-tidyr::gather(sr_detect,key="sample",value = "percent",-structural_category)
  sr_detect$percent<-100*(sr_detect$percent)
  sr_detect<-ggplot(sr_detect, aes(x=structural_category, y=percent, fill=structural_category)) +
    geom_violin(show.legend = FALSE) +geom_boxplot(width=0.2,position=position_dodge(0.9))+
      theme_bw()+ theme(legend.position = "none") +
      theme(axis.text.x=element_text(colour="black",family="Times",size=15),
            axis.text.y=element_text(family="Times",size=12,face="plain"),
            axis.title.x=element_text(family="Times",size = 16,face="plain"),
            axis.title.y=element_text(family="Times",size = 16,face="plain"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            title=element_text(family="Times",size = 14,face="plain")
      )
  
  svglite("/path/to/output/figures/sr_detect.svg", width = 15, height = 15)
  plot(sr_detect)
  dev.off()
  
}

T_specific<-function(sqanti3_qc,sqanti3_faa,TPM_exp){
  novel_transcripts<-sqanti3_qc[sqanti3_qc$structural_category %in% c("novel_in_catalog","novel_not_in_catalog"),]
  dim(novel_transcripts) 

  merged_combined_novel<-merged_combined[merged_combined$transcript_id %in% novel_transcripts$isoform]
  merged_combined_novel$gene_id<-sqanti3_qc[match(merged_combined_novel$transcript_id,sqanti3_qc$isoform),]$associated_gene

  TPM_exp<-"/path/to/short_read_to_long_read_transcriptome/isoform-count-matrix.txt"
  TPM_matrix<-read.table(TPM_exp,header=T,row.names = 1)
  TPM_matrix[TPM_matrix<1]<-0
  
  low_exp<-TPM_matrix[rowSums(TPM_matrix)<1,]
  dim(low_exp)

  TPM_matrix<-TPM_matrix[rowSums(TPM_matrix)>=1,]
  tumor<-TPM_matrix %>% dplyr::select(!matches("^NP\\d+"))
  normal<-TPM_matrix %>% dplyr::select(matches("^NP\\d+"))

  num_T_specific<-setdiff(rownames(normal[rowSums(normal)==0,]),rownames(tumor[rowSums(tumor)==0,])) %>% length()
  T_specific<-TPM_matrix[setdiff(rownames(normal[rowSums(normal)==0,]),rownames(tumor[rowSums(tumor)==0,])),]
  dim(T_specific) 
  
  N_specific<-TPM_matrix[setdiff(rownames(tumor[rowSums(tumor)==0,]),rownames(normal[rowSums(normal)==0,])),] 
  All_exp<-TPM_matrix[intersect(rownames(tumor[rowSums(tumor)!=0,]),rownames(normal[rowSums(normal)!=0,])),]
  dim(All_exp) 

  T_specific$associated_gene<-sqanti3_qc[match(rownames(T_specific),sqanti3_qc$isoform),]$associated_gene
  T_specific$structural_category<-sqanti3_qc[match(rownames(T_specific),sqanti3_qc$isoform),]$structural_category
  T_specific$gene_symbol<-merged_track[match(rownames(T_specific),merged_track$Query_transfrag_id),]$Reference_gene_id

  
  TPM_matrix$mean<- (TPM_matrix%>% as.matrix() %>% rowSums(.,) %>% as.numeric())/dim(TPM_matrix)[2]
  TPM_matrix$`log2(mean+1)`<-log2(TPM_matrix$mean+1)
  TPM_matrix$median<-TPM_matrix %>% dplyr::select(contains(".isoforms")) %>% as.matrix()  %>% rowMedians(useNames = TRUE)
  TPM_matrix$`log2(median+1)`<-log2(TPM_matrix$median+1)
  TPM_matrix$structural_category<-sqanti3_qc[match(rownames(TPM_matrix),sqanti3_qc$isoform),]$structural_category

  density<-ggplot(TPM_matrix %>% dplyr::filter(structural_category %in% c("full-splice_match","novel_in_catalog","novel_not_in_catalog")), aes(x = `log2(median+1)`,fill=structural_category)) + geom_density(alpha=.25)+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(legend.title= element_blank(),axis.text = element_text(size = 20),axis.title = element_text(size=20),legend.position = 'bottom')
  svglite("/path/to/output/figures/TPM_density.svg", width = 10, height = 10)
  plot(density)
  dev.off()
}
