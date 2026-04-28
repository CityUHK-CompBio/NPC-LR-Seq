# Data loading and context initialization

library(rtracklayer)
library(seqinr)
library(dplyr)
library(stringr)

load_tumor_specificity_context <- function(
  reference_gtf_path = "/path/to/reference/gencode.v38.annotation.gtf",
  reference_pc_translations_path = "/path/to/reference/gencode.v38.pc_translations.fa",
  sample_list_path = "/path/to/long_read/sample_list.txt",
  sr_sample_path = "/path/to/short_read/sr_sample.txt",
  sqanti3_qc_path = "/path/to/output/sqanti3_all/SQANTI3_QC_ALL/classification.txt",
  sqanti3_faa_path = "/path/to/output/sqanti3_all/SQANTI3_QC_ALL/corrected.faa",
  num_fl_transcript_path = "/path/to/output/cupcake/num_full_transcript.txt",
  merged_combined_gtf_path = "/path/to/output/gffcompare/merged.combined.gtf",
  merged_tracking_path = "/path/to/output/gffcompare/merged.tracking",
  sqanti3_gtf_path = "/path/to/output/sqanti3_all/SQANTI3_QC_ALL/corrected.gtf",
  TPM_exp = "/path/to/short_read_to_long_read_transcriptome/isoform-count-matrix.txt",
  normal_input_list_path = "/path/to/long_read/normal_input_list.txt"
) {  
  reference_gtf <- rtracklayer::import(reference_gtf_path, format = "gtf")
  reference_pc_translations <- seqinr::read.fasta(reference_pc_translations_path)

  input_list <- read.table(sample_list_path, stringsAsFactors = FALSE)
  input_list <- gsub(".collapsed.gff", "", input_list$V1)

  normal_input_list <- read.table(normal_input_list_path, stringsAsFactors = FALSE, header = FALSE)
  normal_input_list <- as.character(normal_input_list[[1]])

  tumor_input_list <- setdiff(input_list, normal_input_list)
  sr_sample <- read.table(sr_sample_path, stringsAsFactors = FALSE)

  sqanti3_qc <- read.csv(sqanti3_qc_path, sep = "\t")
  sqanti3_faa <- seqinr::read.fasta(sqanti3_faa_path)

  num_fl_transcript <- read.table(num_fl_transcript_path, header = FALSE)
  merged_combined <- rtracklayer::import.gff(merged_combined_gtf_path)
  merged_combined$gene_id <- sqanti3_qc[match(merged_combined$transcript_id, sqanti3_qc$isoform), ]$associated_gene

  merged_track <- read.table(merged_tracking_path, sep = "\t")
  names(merged_track) <- c("Query_transfrag_id", "Query_locus_id", "Reference_gene_id", "Class_code", input_list)
  merged_track$samples <- rowSums(dplyr::select(merged_track, dplyr::matches("NP_")) != "-")

  peptides_id <- seqinr::getName(sqanti3_faa) %>%
    stringr::str_split("\t|\\|") %>%
    do.call(rbind.data.frame, .) %>%
    dplyr::select(c(1, 4)) %>%
    setNames(c("transcript_id", "ORF_length"))

  names(sqanti3_faa) <- peptides_id$transcript_id
  sqanti3_gtf <- rtracklayer::import.gff(sqanti3_gtf_path)

  list(
    reference_gtf = reference_gtf,
    reference_pc_translations = reference_pc_translations,
    input_list = input_list,
    normal_input_list = normal_input_list,
    tumor_input_list = tumor_input_list,
    sr_sample = sr_sample,
    sqanti3_qc = sqanti3_qc,
    sqanti3_faa = sqanti3_faa,
    num_fl_transcript = num_fl_transcript,
    merged_combined = merged_combined,
    merged_track = merged_track,
    peptides_id = peptides_id,
    sqanti3_gtf = sqanti3_gtf,
    TPM_exp = TPM_exp
  )
}
