#load the libraries needed

library(tidyverse)        # dplyr, tidyr, ggplot2, stringr, readr
library(tximport)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)

##############################################################################################################
#setup the directories
CH_RiboTag_salmon_directory=file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/Data_CHale_et_al_2021_RiboTag_salmon")
list.files(CH_RiboTag_salmon_directory)

##############################################################################################################
#code for how I made these RDS files are in the "RDS_files" folder in this github repository

mm10_Tx_final <- readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/Data_mm10gtf/mm10_Tx_final.rds"))

##############################################################################################################
#looking at subcellular localization of RiboTag categories from Hale et al 2021, RiboTag on microdissected CA1 soma and neuropils
CH_RiboTag_files <- list.files(
  path = CH_RiboTag_salmon_directory,
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

CH_RiboTag_files

CH_RiboTag_sample_names <- basename(dirname(CH_RiboTag_files))
names(CH_RiboTag_files) <- CH_RiboTag_sample_names

all(file.exists(CH_RiboTag_files))

length(CH_RiboTag_files) #8

CH_RiboTag_txi <- tximport(CH_RiboTag_files, type = "salmon", tx2gene=mm10_Tx_final %>% dplyr::select(transcript_id,gene_name))

CH_RiboTag_Metadata <- data.frame(sample = CH_RiboTag_sample_names, stringsAsFactors = FALSE) %>%
  mutate(
    SampleID   = sapply(strsplit(sample, "_"), `[`, 1),
    Experiment = sapply(strsplit(sample, "_"), `[`, 2),
    Fraction   = sapply(strsplit(sample, "_"), `[`, 3),
    Mapper     = sapply(strsplit(sample, "_"), `[`, 4)
  ) %>%
  column_to_rownames("sample")

CH_RiboTag_Metadata

CH_RiboTag_Metadata$Group <- factor(CH_RiboTag_Metadata$Fraction)

CH_RT_dds <- DESeqDataSetFromTximport(CH_RiboTag_txi, CH_RiboTag_Metadata, ~ Group)
CH_RT_dds <- DESeq(CH_RT_dds)

source(file.path("~/ruthasinger_github/OptoCLIP_April2026/R_functions/run_deseq_contrasts_v17.R"))

CH_comparisons <- data.frame(
  label   = c("NPvsCB"),
  mode    = "contrast",
  treat   = c("NP"),
  control = c("CB")
)

CH_RT_data=run_deseq_contrasts(
  dds = CH_RT_dds,                         # Your DESeqDataSet object
  my_comparisons = CH_comparisons,      # Comparisons table
  Outdirectory = RT_Outdirectory,          # Main output directory
  tag = "CH_RiboTag",
  var1 = "Group",
  ref1 = "CB",
  interaction = FALSE,
  use_lrt = FALSE,
  reduced_formula = NULL,
  anno_df = NULL ,                     
  dds_id_col = "gene_name",                   
  collapse_method = "none",              
  filter_mode = "any_group_all",
  filter_threshold = 5,
  filter_counts_type = "norm",
  min_reps = "all",
  sig_metric = "padj",
  sig_threshold = 0.05,
  n_label_up = 30,
  n_label_down = 30,
  OrgDb = org.Mm.eg.db,
  do_pca = FALSE,
  do_go = FALSE,
  make_venn = FALSE,
  save_deseq_tables = FALSE,
  save_volcano = FALSE,                  
  outfile_tag = "Salmon_Collapsedgene_5norm" 
)

CH_RT_resdata=CH_RT_data$results_list$NP_vs_CB
nrow(CH_RT_resdata) #11931

CH_RT_resdata=CH_RT_resdata %>%
  dplyr::rename(CH_NPvsCB_log2FC = log2FoldChange,
                CH_NPvsCB_pvalue = pvalue,
                CH_NPvsCB_padj   = padj)

saveRDS(CH_RT_resdata, file=file.path(CH_RiboTag_salmon_directory,"CH_RT_resdata.rds"))
################################################################################################################################