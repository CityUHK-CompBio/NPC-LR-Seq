# Isoform characterization and plotting functions
library(ggstatsplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggsci)
library(svglite)
library(ggridges)
library(reshape2)
library(scales)

sq_sub<-function(sqanti3_qc){
  sqanti3_qc$structural_category[which(sqanti3_qc$structural_category %in% c("antisense","genic","intergenic","fusion"))]="others"
  isoform_category<-ggpiestats(sqanti3_qc, structural_category,
             results.subtitle = F,
             factor.levels = c('4 Cylinders', '6 Cylinders', '8 Cylinders'),
             slice.label = 'percentage', 
             perc.k = 2, 
             direction = 1, 
             palette = 'Pastel2', 
             label="both",
             label.args = list(size = 8),digits.perc = 2,
             title = 'Percentage of full-length transcripts subcategories')+theme(legend.position = 'bottom',legend.text = element_text(size = 20),legend.title=element_text(size=5),plot.title = element_text(size = 20))
  print(isoform_category)
  table(sqanti3_qc$structural_category)
  svglite("/path/to/output/figures/isoform_category.svg", width = 10, height = 10)
  plot(isoform_category)
  dev.off()
}


fl_length<-function(sqanti3_qc){
  sqanti3_qc$structural_category[which( sqanti3_qc$structural_category %in% c("antisense","genic","intergenic","fusion"))]="others"
  fl_length<-ggplot(sqanti3_qc)+geom_histogram(aes(x=length,fill=structural_category),color="black",bins=50)+theme_bw()+theme(panel.grid=element_blank(),axis.title=element_text(size=20),axis.text=element_text(size=20))
  fl_length<-ggplot(sqanti3_qc) +ggridges::geom_density_ridges(aes(x = length, y = structural_category, fill = structural_category))+theme_bw()+theme(panel.grid=element_blank(),axis.title=element_text(size=20),axis.text=element_text(size=20),legend.position = "none")+labs(y = "")
  svglite("/path/to/output/figures/fl_length.svg", width = 10, height = 10)
  plot(fl_length)
  dev.off()
}

gene_transcript<-function(sqanti3_qc,cut_off=10,fill="#365A87"){
  num_genes<-unique(sqanti3_qc$associated_gene) %>% length() 
  novel_loci<-unique(grep("_",sqanti3_qc$associated_gene,value = T)) %>% length()
  gene_freq <- sqanti3_qc %>%
    dplyr::count(associated_gene,name = 'freq', sort = TRUE) 
  freq_data<-gene_freq %>%
    group_by(`Transcripts Per Gene` = cut(freq, breaks = c(0:cut_off, Inf), labels = c(c(1:cut_off), ">10"))) %>%
    dplyr::count(`Transcripts Per Gene`,name = 'num', sort = TRUE)  %>% as.data.frame()

  names(freq_data)<-c("Transcripts Per Gene","Number of Gene")
  freq_data$`Transcripts Per Gene`<- factor(freq_data$`Transcripts Per Gene`, levels = c(c(1:cut_off), ">10"))
  freq_data$`Number of Gene`<-as.numeric(freq_data$`Number of Gene`)
  freq_data<-ggplot(freq_data,aes(x=`Transcripts Per Gene`,y=`Number of Gene`))+
    geom_bar(fill=fill,stat="identity")+theme_bw()+theme(panel.grid=element_blank(),axis.title=element_text(size=20),axis.text=element_text(size=20))
  
  svglite("/path/to/output/figures/gene_transcript.svg", width = 10, height = 10)
  plot(freq_data)
  dev.off()
}

isoform_distribution<-function(sqanti3_qc,merged_track){
dat<- sqanti3_qc %>% filter(structural_category %in% c("full-splice_match","novel_in_catalog","novel_not_in_catalog")) 
dat<-inner_join(dat,merged_track,c("isoform"="Query_transfrag_id")) 
dat <- dat %>%
  mutate(sample_group = case_when(
    samples == 1 ~ "1",
    samples >= 2 & samples <= 5 ~ "2-5",
    samples >= 6 & samples <= 10 ~ "6-10",
    samples > 10 ~ "10+"
  ))


plot_dat <- dat %>% dplyr::group_by(structural_category,sample_group) %>% tally() %>% ungroup()
plot_dat$sample_group<-factor(plot_dat$sample_group,levels=c("1","2-5","6-10","10+"))
plot_dat<-plot_dat %>%
  dplyr::group_by(structural_category) %>%
  dplyr::mutate(
    total = sum(n),
    prop = n / total,
    percent = prop * 100
  ) %>%
  dplyr::ungroup()


p<-ggplot(plot_dat, aes(x = sample_group, y = n, fill = structural_category)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~structural_category, scales = "free") +
  theme_minimal() +
  labs(
    x = "Number of Samples",
    y = "Number of Isoforms",
  )+theme(axis.title =element_text(size=25),axis.text = element_text(size = 25),legend.text = element_text(size = 25),legend.position = 'none', strip.text = element_text(size = 25),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

svglite("/path/to/output/figures/isoform_distribution.svg", width = 10, height = 10)
plot(p)
dev.off()

TPM_exp="/path/to/short_read_to_long_read_transcriptome/isoform-count-matrix.txt"
TPM_matrix<-read.table(TPM_exp,header=T,row.names = 1)
TPM_matrix[TPM_matrix<1]<-0
sub_TPM_matrix<-TPM_matrix
sub_TPM_matrix$count<-rowSums(sub_TPM_matrix!= 0, na.rm = TRUE)
sub_TPM_matrix$structural_category<-sqanti3_qc[match(rownames(sub_TPM_matrix),sqanti3_qc$isoform),]$structural_category


df_plot <- sub_TPM_matrix %>%
  filter(structural_category %in% c("full-splice_match",
                                    "novel_in_catalog",
                                    "novel_not_in_catalog")) %>%
  filter(count > 0) %>%
  mutate(
    bin = case_when(
      count == 1 ~ "1",
      count >= 2 & count <= 5 ~ "2-5",
      count >= 6 & count <= 10 ~ "6-10",
      count >= 11 ~ "10+",
      TRUE ~ NA_character_
    ),
    bin = factor(bin, levels = c("1", "2-5", "6-10", "10+")),
    structural_category = factor(structural_category,
                                levels = c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog"),
                                labels = c("FSM", "NIC", "NNC")
    )
  ) %>%
  dplyr::group_by(structural_category, bin) %>%
  dplyr::summarise(n_isoforms = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(structural_category) %>%
  dplyr::mutate(
    total_isoforms = sum(n_isoforms),
    prop = n_isoforms / total_isoforms,
    percent = prop * 100
  ) %>%
  dplyr::ungroup()

  df_plot$bin <- factor(df_plot$bin, levels = c("1", "2-5", "6-10", "10+"))  

p <- ggplot(df_plot, aes(x = bin, y = n_isoforms)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", percent)),
            vjust = -0.3, size = 3) +
  facet_wrap(~ structural_category, nrow = 1, scales = "free_y") +
  labs(x = "Number of Samples", y = "Number of Isoforms (×10^3)") +
  theme_bw()

svglite("/path/to/output/figures/isoform_distribution_sr.svg", width = 10, height = 10)
plot(p)
dev.off()

}

per_sample_gene_isoform<-function(){
  path <- "/path/to/output/sqanti3_all/per_sample"
  dirs <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  dirs <- sort(dirs)
  
  res <- data.frame(
    sample = character(),
    isoform_number = integer(),
    gene_number = integer(),
    stringsAsFactors = FALSE
  )
  for (dir in dirs){
  sqanti <- list.files(
    file.path(path, dir),
    pattern = "classification",   
    full.names = TRUE,
    recursive = TRUE          
  )
  sqanti<-read.table(sqanti,header=T)
  isoform_number<-nrow(sqanti)
  print(isoform_number)
  gene_number<-sqanti$associated_gene %>% unique() %>% length()
  print(gene_number)
  res <- rbind(
    res,
    data.frame(sample = dir, isoform_number = isoform_number, gene_number = gene_number)
  )
  
  }
  
df_long <- res %>%
  pivot_longer(cols = c(gene_number, isoform_number),
                names_to = "metric", values_to = "count") %>%
  mutate(metric = recode(metric,
                          gene_number = "Genes detected",
                          isoform_number = "Isoforms detected"),
          sample = factor(sample, levels=c(rev(normal_input_list),rev(tumor_input_list))))


p <- ggplot(df_long, aes(x = sample, y = count, fill = metric)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_x") +
  labs(x = "Sample", y = "Count") +
  theme_bw()+coord_flip() 
svglite("/path/to/output/figures/per_sample_gene_isoform.svg", width = 10, height = 10)
plot(p)
dev.off()
  
}

stack_columns_charts<-function(){
  merged_track$structural_category<-sqanti3_qc[match(merged_track$Query_transfrag_id,sqanti3_qc$isoform),]$structural_category
  merged_track$structural_category[which(merged_track$structural_category %in% c("antisense","genic","intergenic","fusion"))]="others"

  process_sample <- function(sample) {
    tmp <- merged_track %>%
      filter(!is.na(structural_category) & .[[sample]] != "-") %>%
      group_by(structural_category) %>%
      tally()
    return(as.vector(tmp$n))
  }
  count <- do.call(rbind, lapply(input_list, process_sample))
  rownames(count)<-input_list
  colnames(count)<-names(table(merged_track$structural_category))
  count<-as.data.frame(count)
  anno<-data.frame(sample=input_list,origin=c(rep(c("tumor cell line"),2),rep(c("normal cell line"),4),rep("xenograft",10),rep("tumor tissue",2)))
  count$origin<-anno[match(rownames(count),anno$sample),]$origin
  count<-bind_rows(count %>% filter(rownames(count) %in% tumor_input_list),count %>% filter(rownames(count) %in% normal_input_list))
  count$sample<-rownames(count)
  count<-reshape2::melt(count,id.vars=c("sample","origin"),variable.name = c("subcategories"),value.name = "number")
  count$sample<-factor(count$sample,levels=unique(count$sample))
  svglite(filename = "/path/to/output/figures/stacked.svg", width = 16, height = 8)
  p<-ggplot(data=count,aes(x=sample,y=number))+geom_col(aes(fill=subcategories),position=position_fill())+ scale_fill_brewer(palette = 'Set1',name = 'Subcategories') +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw(base_size = 18) +
    coord_cartesian(expand = 0,clip = 'off') +
    xlab('') + ylab('') +
    theme(plot.margin = margin(t = 1,b = 3,unit = 'cm'),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p1 <- annoRect(object = p,
                 annoPos = 'botomn',
                 aesGroup = T,
                 aesGroName = 'origin',
                 yPosition = c(-0.08,-0.05),
                 rectWidth = 0.9,
                 pCol = c('#77ACF1','#FFC947','#B85C38','#AAAAAA','#364547','#FA26A0'),pFill = c('#77ACF1',
                                                                                                          '#FFC947','#B85C38','#AAAAAA','#364547','#FA26A0'))
                                                                                                          
  p2<-annoLegend(object = p1,
                 xPosition = 10,
                 yPosition = 0.1,
                 labels = c(unique(count$origin)),
                 pch = 15,
                 col = c('#77ACF1',
                                  '#FFC947','#B85C38','#AAAAAA','#364547','#FA26A0'),
                                  vgap = 0.1,
                 cex = 1.2)
  
  p2
  dev.off()
}

fl_features<-function(sqanti3_qc){
  sqanti3_qc$structural_category[which(sqanti3_qc$structural_category %in% c("antisense","genic","intergenic","fusion"))]="others"
  table(sqanti3_qc$structural_category)
  
  process_table <- function(data, column_name) {
    dat <- table(sqanti3_qc$structural_category, data[[column_name]]) %>% as.data.frame()
    names(dat) <- c("structural_category", column_name, "Freq")
    
    if (column_name == "all_canonical") {
      dat$all_canonical <- as.logical(dat$all_canonical == "canonical")
    }
    
    if (column_name == "coding") {
      dat$coding <- as.logical(dat$coding == "coding")
      }
    
    dat <- ddply(dat, .variables = 'structural_category', .fun = transform, percent_con = Freq / sum(Freq) * 100)
    dat[, 2] <- factor(dat[, 2], levels = c(TRUE, FALSE))
    return(dat)
  }
  
  dat<-process_table(sqanti3_qc, "all_canonical")
  dat <- process_table(sqanti3_qc, "predicted_NMD")
  dat <- process_table(sqanti3_qc, "coding")
  
  feature<-ggplot(dat,aes(x=structural_category,y=percent_con,fill=dat[,2]))+
    geom_bar(stat = 'identity',width = 0.5,colour='black') + 
    scale_fill_brewer(palette = 'Set2') +
    theme_bw()+
    theme(axis.text.x=element_text(angle=0,size=20,face ='bold'),axis.text.y=element_text(angle=0,size=20,face ='bold'),axis.title.y=element_text(size=20,face ='bold'),panel.grid = element_blank(),legend.position = 'bottom') + 
    labs(y = 'Percent') + labs(x="")+
    scale_y_continuous()
  
  svglite("/path/to/output/figures/canonical.svg", width = 10, height = 10)
  plot(feature)
  dev.off()

}

compare_FSM_NIC_NNC<-function(factor){
sub_sqanti3_qc<- sqanti3_qc %>% filter(structural_category %in% c("full-splice_match","novel_in_catalog","novel_not_in_catalog"))
comparisons <- list(
  c("full-splice_match", "novel_in_catalog"),
  c("full-splice_match", "novel_not_in_catalog"),
  c("novel_in_catalog", "novel_not_in_catalog")
)
feature<-ggplot(sub_sqanti3_qc, aes(x = structural_category, y = .data[[factor]], fill = structural_category)) +
  geom_violin(trim = TRUE)+ geom_boxplot(fill="white",width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() + 
  labs(
    x="",
    y = factor
  ) +
  theme(legend.position = "none")+stat_compare_means(
    method = "wilcox.test", 
    comparisons = comparisons, 
    label = "p.signif"
  )
svglite(file.path("/path/to/output/figures", paste0(factor, ".svg")), width = 10, height = 10)
plot(feature)
dev.off()
}

predicted_NMD<-function(){
 NMD<-sqanti3_qc %>% filter(structural_category %in% c("novel_not_in_catalog","novel_in_catalog")) %>% filter(.$predicted_NMD==TRUE)
 p<-ggpiestats(NMD, subcategory,results.subtitle = F,factor.levels = c('4 Cylinders', '6 Cylinders', '8 Cylinders'),slice.label = 'percentage', 
                                     perc.k = 2, 
                                     direction = 1, 
                                     title = '')+theme(plot.title = element_text(size = 15,hjust=0.5),text = element_text(size = 15),legend.text = element_text(size = 20),legend.position = 'right') + scale_fill_npg()
  
 svglite("/path/to/output/figures/NMD_subcategory.svg", width = 20, height = 10)
 plot(p)
 dev.off()

}

num_transcripts_exons <- function(){
  sqanti3_qc
  gene_count <- table(sqanti3_qc$associated_gene)
  gene_count_df <- as.data.frame(gene_count)
  colnames(gene_count_df) <- c("associated_gene", "count")
  
  reference_gtf_df <- as.data.frame(reference_gtf)
  reference_gtf_df$exon_number <- as.numeric(reference_gtf_df$exon_number)
  reference_gtf_df <- reference_gtf_df %>%
    filter(!is.na(gene_id) & !is.na(exon_number))
  exon_summary <- as.data.frame(reference_gtf_df) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(max_exon_number = max(exon_number, na.rm = TRUE))
  
  result <- gene_count_df %>%
    left_join(exon_summary, by = c("associated_gene" = "gene_id"))
  plot(result$count, result$max_exon_number, 
       xlab = "Number of Transcripts", 
       ylab = "Number of Exons", 
       main = "",
       pch = 16,
       col =c("#ff5e00"),
       cex.axis = 1.5,cex.lab = 1.8)
  correlation<-cor(result$count, result$max_exon_number, use = "complete.obs")
  text(x = max(result$count) * 0.6, y = max(result$max_exon_number) * 0.7,
    labels = paste("Correlation: ", round(correlation, 2)),
       col = "black", cex = 1.5) 
  abline(lm(max_exon_number ~ count, data = result), col = "black", lwd = 2)
  
}
