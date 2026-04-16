#load the libraries needed

library(tidyverse)        # dplyr, tidyr, ggplot2, stringr, readr
library(tximport)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(corrplot)
library(VennDetail)

##############################################################################################################
#setup the directories
Base_directory=file.path("~/ruthasinger_github/OptoCLIP_April2026/Figure4_FigureS4_RiboTag")
list.files(Base_directory)

Data_directory=file.path(Base_directory,"Data")
list.files(Data_directory)

RiboTag_salmon_directory=file.path(Data_directory,"Data_RiboTag")
list.files(RiboTag_salmon_directory)

FMRP_CLIP_STAR_Directory=file.path(Data_directory,"Data_FMRP_CLIP_Tx_STAR")
list.files(FMRP_CLIP_STAR_Directory)

Outdirectory=file.path(Base_directory,paste("Output",Sys.Date(),sep="_"))
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#code for how I made these RDS files are in the "RDS_files" folder in this github repository

mm10_Tx_final <- readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/mm10_Tx_final.rds"))
mm10_txdb_lengths <- readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/mm10_txdb_lengths.rds"))
CH_RT_resdata=readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/CH_RT_resdata.rds"))

##############################################################################################################

theme_illustrator_black <- theme_classic(base_size = 18) +
  theme(
    text              = element_text(color = "black"),
    axis.text         = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    plot.title        = element_text(color = "black"),
    plot.subtitle     = element_text(color = "black"),
    plot.caption      = element_text(color = "black"),
    legend.text       = element_text(color = "black"),
    legend.title      = element_text(color = "black"),
    axis.line         = element_line(color = "black"),
    axis.ticks        = element_line(color = "black")
  )

##############################################################################################################
#Bring in Opto-RiboTag data Input vs IP enrichment
RiboTag_files <- list.files(
  path = file.path(RiboTag_salmon_directory,"Data_InputandIP"),
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

RiboTag_files

RiboTag_sample_names <- basename(dirname(RiboTag_files))
names(RiboTag_files) <- RiboTag_sample_names

all(file.exists(RiboTag_files))

length(RiboTag_files) #26

make_metadata <- function(sample_names) {
  parts <- do.call(rbind, strsplit(sample_names, "_", fixed = TRUE))
  df <- data.frame(
    sample     = sample_names,
    Experiment = parts[,1],
    Batch      = parts[,2],
    MouseID    = parts[,3],
    Condition  = parts[,4],
    Timepoint  = parts[,5],
    Fraction   = parts[,6],
    Mapper     = parts[,7],
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$sample
  return(df)
}

# Usage
RiboTag_Metadata <- make_metadata(RiboTag_sample_names)

RiboTag_Metadata <- RiboTag_Metadata %>%
  mutate(
    MouseID   = factor(MouseID),
    Batch     = factor(Batch),
    Fraction  = factor(Fraction,  levels = c("Input","IP")),
    Condition = factor(Condition, levels = c("noCre","Control","ChR2")),
    Timepoint = factor(Timepoint, levels = c("noCre","5min","30min")),
    Group     = factor(paste(Condition, Timepoint, sep = "_"))
  )

RiboTag_txi <- tximport(RiboTag_files, type = "salmon", tx2gene=mm10_Tx_final %>% dplyr::select(transcript_id,gene_name))

RiboTag_salmon_TPM_mat <- RiboTag_txi$abundance 

RiboTag_salmon_TPM_mat <- RiboTag_salmon_TPM_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene_name")

RiboTag_TPM=RiboTag_salmon_TPM_mat %>%
  dplyr::mutate(ChR2_5min_rep1_log2FC      = log2(RiboTag_Exp03_M18_ChR2_5min_IP_salmon/RiboTag_Exp03_M18_ChR2_5min_Input_salmon),
                ChR2_5min_rep2_log2FC      = log2(RiboTag_Exp03_M19_ChR2_5min_IP_salmon/RiboTag_Exp03_M19_ChR2_5min_Input_salmon),
                ChR2_5min_rep3_log2FC      = log2(RiboTag_Exp04_M36_ChR2_5min_IP_salmon/RiboTag_Exp04_M36_ChR2_5min_Input_salmon),
                ChR2_5min_rep4_log2FC      = log2(RiboTag_Exp04_M37_ChR2_5min_IP_salmon/RiboTag_Exp04_M37_ChR2_5min_Input_salmon),
                Control_30min_rep1_log2FC  = log2(RiboTag_Exp03_M26_Control_30min_IP_salmon/RiboTag_Exp03_M26_Control_30min_Input_salmon),
                Control_30min_rep2_log2FC  = log2(RiboTag_Exp03_M27_Control_30min_IP_salmon/RiboTag_Exp03_M27_Control_30min_Input_salmon),
                Control_30min_rep3_log2FC  = log2(RiboTag_Exp04_M38_Control_30min_IP_salmon/RiboTag_Exp04_M38_Control_30min_Input_salmon),
                Control_30min_rep4_log2FC  = log2(RiboTag_Exp04_M39_Control_30min_IP_salmon/RiboTag_Exp04_M39_Control_30min_Input_salmon),
                ChR2_30min_rep1_log2FC     = log2(RiboTag_Exp03_M16_ChR2_30min_IP_salmon/RiboTag_Exp03_M16_ChR2_30min_Input_salmon),
                ChR2_30min_rep2_log2FC     = log2(RiboTag_Exp03_M17_ChR2_30min_IP_salmon/RiboTag_Exp03_M17_ChR2_30min_Input_salmon),
                ChR2_30min_rep3_log2FC     = log2(RiboTag_Exp04_M34_ChR2_30min_IP_salmon/RiboTag_Exp04_M34_ChR2_30min_Input_salmon),
                ChR2_30min_rep4_log2FC     = log2(RiboTag_Exp04_M35_ChR2_30min_IP_salmon/RiboTag_Exp04_M35_ChR2_30min_Input_salmon),
                noCre_rep1_log2FC          = log2(RiboTag_Exp01_M01_noCre_noCre_IP_salmon/RiboTag_Exp01_M01_noCre_noCre_Input_salmon))

GoI=c("Neurod6",
      "Camk2a",
      "Rbfox3",
      "Snap25",
      "Nrgn",
      "Hpca",
      "Crym",
      "Chn1",
      "Gad2",
      "Sst",
      "Calb2",
      "Pdgfra",
      "Ptprz1",
      "Mag",
      "Mal",
      "Mbp",
      "Mobp",
      "Plp1",
      "Gfap",
      "Glul",
      "Aqp4",
      "Aldh1l1",
      "Pla2g7",
      "Slc1a3",
      "Aldoc")

RiboTag_TPM_GOI=RiboTag_TPM[RiboTag_TPM$gene_name %in% GoI,]

RiboTag_TPM_GOI_log2FC=RiboTag_TPM_GOI %>% dplyr::select(c(gene_name,
                                                           noCre_rep1_log2FC,
                                                           ChR2_5min_rep1_log2FC,
                                                           ChR2_5min_rep2_log2FC,
                                                           ChR2_5min_rep3_log2FC,
                                                           ChR2_5min_rep4_log2FC,
                                                           Control_30min_rep1_log2FC,
                                                           Control_30min_rep2_log2FC,
                                                           Control_30min_rep3_log2FC,
                                                           Control_30min_rep4_log2FC,
                                                           ChR2_30min_rep1_log2FC,
                                                           ChR2_30min_rep2_log2FC,
                                                           ChR2_30min_rep3_log2FC,
                                                           ChR2_30min_rep4_log2FC))

RiboTag_TPM_GOI_long <- RiboTag_TPM_GOI_log2FC %>%
  pivot_longer(names_to = "Sample", values_to = "log2FC_TPM", cols = ends_with("log2FC")) 

#order samples
sample_order <- c(
  "noCre_rep1",
  "Control_30min_rep1", "Control_30min_rep2", "Control_30min_rep3", "Control_30min_rep4",
  "ChR2_5min_rep1", "ChR2_5min_rep2", "ChR2_5min_rep3", "ChR2_5min_rep4",
  "ChR2_30min_rep1", "ChR2_30min_rep2", "ChR2_30min_rep3", "ChR2_30min_rep4"
)

RiboTag_TPM_GOI_sample_use <- RiboTag_TPM_GOI_long %>%
  mutate(
    Sample_clean = str_remove(Sample, "_log2FC$")
  ) %>%
  filter(Sample_clean %in% sample_order) %>%
  mutate(
    Sample_clean = factor(Sample_clean, levels = sample_order),
    gene_name = factor(gene_name, levels = rev(GoI))
  )

RT_IPvsInput_sample_plot <- ggplot(
  RiboTag_TPM_GOI_sample_use,
  aes(Sample_clean, gene_name, fill = log2FC_TPM)
) +
  geom_tile() +
  scale_x_discrete(limits = sample_order) +
  scale_y_discrete(limits = rev(GoI)) +
  coord_fixed(ratio = 1/6) +
  scale_fill_viridis_c(
    direction = 1,
    limits = c(-4, 2),
    breaks = c(-4, -2, 0, 2),
    oob = scales::squish
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
RT_IPvsInput_sample_plot
ggsave(file.path(Outdirectory, "RiboTag_IPvsInput_perSample_plot.pdf"),
       plot = RT_IPvsInput_sample_plot,
       device = "pdf",
       width = 8,
       height = 6)

##############################################################################################################
#Bring in Opto-RiboTag data for batch corrected TPM 
RiboTag_files <- list.files(
  path = file.path(RiboTag_salmon_directory,"Data_InputandIP"),
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

RiboTag_files <- RiboTag_files[!grepl("noCre", RiboTag_files)]

RiboTag_sample_names <- basename(dirname(RiboTag_files))
names(RiboTag_files) <- RiboTag_sample_names

all(file.exists(RiboTag_files))

length(RiboTag_files) #24

make_metadata <- function(sample_names) {
  parts <- do.call(rbind, strsplit(sample_names, "_", fixed = TRUE))
  df <- data.frame(
    sample     = sample_names,
    Experiment = parts[,1],
    Batch      = parts[,2],
    MouseID    = parts[,3],
    Condition  = parts[,4],
    Timepoint  = parts[,5],
    Fraction   = parts[,6],
    Mapper     = parts[,7],
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$sample
  return(df)
}

# Usage
RiboTag_Metadata <- make_metadata(RiboTag_sample_names)

RiboTag_Metadata <- RiboTag_Metadata %>%
  mutate(
    MouseID   = factor(MouseID),
    Batch     = factor(Batch),
    Fraction  = factor(Fraction,  levels = c("Input","IP")),
    Condition = factor(Condition, levels = c("Control","ChR2")),
    Timepoint = factor(Timepoint, levels = c("5min","30min")),
    Group     = factor(paste(Condition, Timepoint, sep = "_"))
  )

RiboTag_txi <- tximport(RiboTag_files, type = "salmon", tx2gene=mm10_Tx_final %>% dplyr::select(transcript_id,gene_name))

RT_dds <- DESeqDataSetFromTximport(
  txi     = RiboTag_txi,
  colData = as.data.frame(RiboTag_Metadata),
  design  = ~ Group + Fraction + Group:Fraction
)

RT_dds <- DESeq(RT_dds, fitType = "local")

# VST
vsd <- vst(RT_dds, blind = TRUE)

#batch corrected PCA
mat     <- assay(vsd)
coldata <- as.data.frame(colData(vsd))

# Make sure factors are as expected
coldata$MouseID   <- factor(coldata$MouseID)
coldata$Group     <- factor(coldata$Group)
coldata$Fraction  <- factor(coldata$Fraction, levels = c("Input", "IP"))

# Model matrix for biological variables you want to preserve
design_mb <- model.matrix(~ Group + Fraction, data = coldata)

# Remove MouseID as batch
mat_bc <- limma::removeBatchEffect(
  mat,
  batch  = coldata$MouseID,
  design = design_mb
)

# Store batch-corrected matrix inside vsd for later use
assay(vsd, "batch_corrected") <- mat_bc

# PCA on top variable genes from batch-corrected matrix
ntop    <- 500
rv      <- matrixStats::rowVars(mat_bc)
select  <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca     <- prcomp(t(mat_bc[select, ]), center = TRUE, scale. = FALSE)
percentVar <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- data.frame(
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  Group     = coldata$Group,
  Fraction  = coldata$Fraction,
  MouseID   = coldata$MouseID,
  row.names = rownames(coldata)
)

plot_groups   <- c("Control_30min", "ChR2_5min", "ChR2_30min")
plot_fraction <- c("Input", "IP")

pca_df_plot <- pca_df %>%
  dplyr::filter(
    Group %in% plot_groups,
    Fraction %in% plot_fraction,
  ) %>%
  dplyr::mutate(
    Group = factor(Group, levels = plot_groups),
    Fraction = factor(Fraction, levels = plot_fraction)
  )

PCA_batch_corrected <- ggplot(pca_df_plot, aes(PC1, PC2, color = Group, shape = Fraction)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", round(percentVar[1] * 100), "%)")) +
  ylab(paste0("PC2 (", round(percentVar[2] * 100), "%)")) +
  theme_classic(base_size = 9) +
  coord_fixed() +
  ylim(-25, 25) +
  scale_color_manual(
    values = setNames(
      c("gray35","dodgerblue","red3"),
      c("Control_30min", "ChR2_5min", "ChR2_30min")
    )
  )
PCA_batch_corrected
ggsave(file.path(Outdirectory,"RiboTag_PCA_batch_corrected_MouseID.pdf"), PCA_batch_corrected, width = 5, height = 5)

##############################################################################################################
#Bring in Opto-RiboTag data for DESeq- 5min
RiboTag_files <- list.files(
  path = file.path(RiboTag_salmon_directory,"Data_5min"),
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

RiboTag_sample_names <- basename(dirname(RiboTag_files))
names(RiboTag_files) <- RiboTag_sample_names

all(file.exists(RiboTag_files))

length(RiboTag_files) #8

make_metadata <- function(sample_names) {
  parts <- do.call(rbind, strsplit(sample_names, "_", fixed = TRUE))
  df <- data.frame(
    sample     = sample_names,
    Experiment = parts[,1],
    Batch      = parts[,2],
    MouseID    = parts[,3],
    Condition  = parts[,4],
    Timepoint  = parts[,5],
    Fraction   = parts[,6],
    Mapper     = parts[,7],
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$sample
  return(df)
}

# Usage
RiboTag_Metadata <- make_metadata(RiboTag_sample_names)

RiboTag_Metadata <- RiboTag_Metadata %>%
  mutate(
    MouseID   = factor(MouseID),
    Batch     = factor(Batch),
    Condition = factor(Condition, levels = c("Control","ChR2")),
    Timepoint = factor(Timepoint, levels = c("5min")),
    Group     = factor(paste(Condition, Timepoint, sep = "_"))
  )

### Run tximport to load in the salmon files
RiboTag_txi <- tximport(RiboTag_files, type = "salmon", tx2gene=mm10_Tx_final %>% dplyr::select(transcript_id,gene_name))

# ---- calculate TPMs ----
RiboTag_salmon_TPM <- as.data.frame(RiboTag_txi$abundance) %>%
  tibble::rownames_to_column("gene_name")

# add ".TPM" suffix to all sample columns
colnames(RiboTag_salmon_TPM) <- c(
  "gene_name",
  paste0(colnames(RiboTag_salmon_TPM)[-1], ".TPM")
)

# ---- order columns:
tp_levels <- c("5min")

tpm_cols <- setdiff(colnames(RiboTag_salmon_TPM), "gene_name")

get_tp <- function(x) str_match(x, "_([0-9]+min)_.*\\.TPM$")[,2]

ctrl_cols  <- tpm_cols[grepl("_Control_", tpm_cols)]
chr2_cols  <- tpm_cols[grepl("_ChR2_", tpm_cols)]

# sort within each group by timepoint, then by full name (keeps reps stable)
ctrl_cols <- ctrl_cols[order(match(get_tp(ctrl_cols), tp_levels), ctrl_cols)]
chr2_cols <- chr2_cols[order(match(get_tp(chr2_cols), tp_levels), chr2_cols)]

RiboTag_salmon_TPM <- RiboTag_salmon_TPM %>%
  dplyr::select(gene_name,  all_of(ctrl_cols), all_of(chr2_cols))

nrow(RiboTag_salmon_TPM) #28467

RT_dds <- DESeqDataSetFromTximport(
  txi     = RiboTag_txi,
  colData = as.data.frame(RiboTag_Metadata),
  design  = ~ Batch + Group
)

RT_dds <- DESeq(RT_dds)

source(file.path("~/ruthasinger_github/OptoCLIP_April2026/R_functions/run_deseq_contrasts.R"))

my_comparisons <- data.frame(
  label = c(
    "ChR2_5min_vs_Control_5min"),
  mode = "contrast",
  treat = c(
    "ChR2_5min"),
  control = c(
    "Control_5min"),
  stringsAsFactors = FALSE
)

my_comparisons <- as.data.frame(my_comparisons, stringsAsFactors = FALSE)
my_comparisons[] <- lapply(my_comparisons, as.character)

RT_5min_data <- run_deseq_contrasts(
  dds = RT_dds,
  my_comparisons = my_comparisons,
  Outdirectory = Outdirectory,
  tag = "RiboTag_5min",
  design_formula = NULL,
  var1 = NULL,
  var2 = NULL,
  ref1 = NULL,
  ref2 = NULL,
  group_col = "Group",
  interaction = FALSE,
  use_lrt = FALSE,
  reduced_formula = NULL,
  anno_df = RiboTag_salmon_TPM,
  dds_id_col = "gene_name",
  anno_id_col = "gene_name",
  anno_gene_col = "gene_name",
  collapse_method = "none",
  filter_mode = "any_group_all",
  filter_threshold = 1,
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
  outfile_tag = "Salmon_Collapsedgene_1norm"
)

################################################################################################################################################################
#Bring in Opto-RiboTag data for 30 min vs 5 min LRT comparison
RiboTag_files <- list.files(
  path = file.path(RiboTag_salmon_directory,"Data_30min"),
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

RiboTag_sample_names <- basename(dirname(RiboTag_files))
names(RiboTag_files) <- RiboTag_sample_names

all(file.exists(RiboTag_files))

length(RiboTag_files) #16

make_metadata <- function(sample_names) {
  parts <- do.call(rbind, strsplit(sample_names, "_", fixed = TRUE))
  df <- data.frame(
    sample     = sample_names,
    Experiment = parts[,1],
    Batch      = parts[,2],
    MouseID    = parts[,3],
    Condition  = parts[,4],
    Timepoint  = parts[,5],
    Fraction   = parts[,6],
    Mapper     = parts[,7],
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$sample
  return(df)
}

# Usage
RiboTag_Metadata <- make_metadata(RiboTag_sample_names)

RiboTag_Metadata <- RiboTag_Metadata %>%
  mutate(
    MouseID   = factor(MouseID),
    Batch     = factor(Batch),
    Condition = factor(Condition, levels = c("Control","ChR2")),
    Timepoint = factor(Timepoint, levels = c("5min","30min")),
    Group     = factor(paste(Condition, Timepoint, sep = "_"))
  )

RiboTag_txi <- tximport(RiboTag_files, type = "salmon", tx2gene=mm10_Tx_final %>% dplyr::select(transcript_id,gene_name))

# ---- calculate TPMs ----
RiboTag_salmon_TPM <- as.data.frame(RiboTag_txi$abundance) %>%
  tibble::rownames_to_column("gene_name")

# add ".TPM" suffix to all sample columns
colnames(RiboTag_salmon_TPM) <- c(
  "gene_name",
  paste0(colnames(RiboTag_salmon_TPM)[-1], ".TPM")
)

# ---- order columns:
tp_levels <- c("5min","30min")

tpm_cols <- setdiff(colnames(RiboTag_salmon_TPM), "gene_name")

get_tp <- function(x) str_match(x, "_([0-9]+min)_.*\\.TPM$")[,2]

ctrl_cols  <- tpm_cols[grepl("_Control_", tpm_cols)]
chr2_cols  <- tpm_cols[grepl("_ChR2_", tpm_cols)]

# sort within each group by timepoint, then by full name (keeps reps stable)
ctrl_cols <- ctrl_cols[order(match(get_tp(ctrl_cols), tp_levels), ctrl_cols)]
chr2_cols <- chr2_cols[order(match(get_tp(chr2_cols), tp_levels), chr2_cols)]

RiboTag_salmon_TPM <- RiboTag_salmon_TPM %>%
  dplyr::select(gene_name,  all_of(ctrl_cols), all_of(chr2_cols))

nrow(RiboTag_salmon_TPM) #28467

RT_dds <- DESeqDataSetFromTximport(
  txi     = RiboTag_txi,
  colData = as.data.frame(RiboTag_Metadata),
  design  = ~ Batch + Group
)

RT_dds <- DESeq(RT_dds)

source(file.path("~/ruthasinger_github/OptoCLIP_April2026/R_functions/run_deseq_contrasts.R"))

# Run DESeq2 contrasts
RT_LRT_data=run_deseq_contrasts(
  dds = RT_dds,                 
  my_comparisons = NULL,     
  Outdirectory = Outdirectory,     
  tag = "RiboTag_Interaction",
  var1 = "Condition",
  var2 = "Timepoint",
  ref1 = "Control",
  ref2 = "5min",
  interaction = TRUE,
  use_lrt = TRUE,
  reduced_formula = ~ Condition + Timepoint,
  anno_df = RiboTag_salmon_TPM ,             
  dds_id_col = "gene_name",                 
  anno_id_col = "gene_name",                
  anno_gene_col = "gene_name",
  collapse_method = "none",           
  filter_mode = "any_group_all",     
  filter_threshold = 1,
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
  outfile_tag = "Salmon_Collapsedgene_1norm_LRT"  
)

################################################################################################################################
RT_specs <- list(
  "5v0" = list(
    label = "5v0",
    results = RT_5min_data$results_list$ChR2_5min_vs_Control_5min,
    prefix = "RiboTag",
    volcano_xlim = c(-2.5, 2.5),
    
    prio_up   = "synaptic_first",
    prio_down = "nuclear_first",
    
    go_subset_up_terms = c(
      "neurotransmitter secretion",
      "calcium ion homeostasis",
      "protein localization to synapse",
      "synapse assembly",
      "regulation of synaptic plasticity",
      "synaptic vesicle cycle",
      "synaptic transmission, glutamatergic",
      "excitatory synapse assembly"
    ),
    go_subset_down_terms = c(
      "mRNA processing",
      "nuclear transport",
      "nuclear export",
      "regulation of DNA-binding transcription factor activity",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
    ),
    volcano_label_up_terms = c(
      "neurotransmitter secretion",
      "calcium ion homeostasis",
      "protein localization to synapse",
      "synapse assembly",
      "regulation of synaptic plasticity",
      "synaptic vesicle cycle",
      "synaptic transmission, glutamatergic",
      "excitatory synapse assembly"
    ),
    volcano_label_down_terms = c(
      "mRNA processing",
      "nuclear transport",
      "nuclear export",
      "regulation of DNA-binding transcription factor activity",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
    ),
    volcano_label_top_n = 30
  ),
  
  "30v5" = list(
    label = "30v5",
    results = RT_LRT_data$lrt_df,
    prefix = "RiboTag",
    volcano_xlim = c(-2.5, 2.5),
    
    prio_up   = "nuclear_first",
    prio_down = "synaptic_first",
    
    go_subset_up_terms = c(
      "mRNA processing",
      "nuclear transport",
      "nuclear export",
      "regulation of DNA-binding transcription factor activity",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
    ),
    go_subset_down_terms = c(
      "neurotransmitter secretion",
      "calcium ion homeostasis",
      "protein localization to synapse",
      "synapse assembly",
      "regulation of synaptic plasticity",
      "synaptic vesicle cycle",
      "synaptic transmission, glutamatergic",
      "excitatory synapse assembly"
    ),
    volcano_label_up_terms = c(
      "mRNA processing",
      "nuclear transport",
      "nuclear export",
      "regulation of DNA-binding transcription factor activity",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
    ),
    volcano_label_down_terms = c(
      "neurotransmitter secretion",
      "calcium ion homeostasis",
      "protein localization to synapse",
      "synapse assembly",
      "regulation of synaptic plasticity",
      "synaptic vesicle cycle",
      "synaptic transmission, glutamatergic",
      "excitatory synapse assembly"
    ),
    volcano_label_top_n = 30
  )
)

desc_col   <- "Description"
geneid_col <- "geneID"
gene_sep   <- "/"

# Keywords
nuclear_keywords <- c(
  "nuclear","chromatin","transcription","gene expression","RNA processing",
  "mRNA processing","splicing","epigenetic","histone","DNA-binding","RNA",
  "nucleotide","ribonucleotide","nucleoside","chromatid","chromosome",
  "heterochromatin","DNA","nuclear-transcribed","nucleus","RNA localization"
)

exclude_nuclear_keywords <- c(
  "external",  "maternal","internalization"
)

synaptic_keywords <- c(
  "synapse","synaptic","neuron projection","cognition","learning",
  "postsynapse","postsynaptic","memory","neurotransmitter","neuron",
  "behavior","axon","calcium","presynapse","adult behavior",
  "transmembrane","second-messenger","receptor"
)

exclude_synaptic_keywords <- c(
  "B cell","platelet-derived growth factor","fibroblast growth factor","immune response-regulating cell surface","photoreceptor"
)

check_go_cols <- function(go_df, desc_col, geneid_col) {
  if (!all(c(desc_col, geneid_col) %in% colnames(go_df))) {
    stop(paste0(
      "Could not find expected columns in GO table.\n",
      "Expected: ", desc_col, " and ", geneid_col, "\n",
      "Found: ", paste(colnames(go_df), collapse = ", ")
    ))
  }
}

flag_terms_by_keywords <- function(go_df, desc_col, geneid_col,
                                   include_keywords, exclude_keywords = character(0)) {
  check_go_cols(go_df, desc_col, geneid_col)
  
  pattern_incl <- str_c(include_keywords, collapse = "|")
  
  terms <- go_df %>%
    filter(str_detect(str_to_lower(.data[[desc_col]]), str_to_lower(pattern_incl)))
  
  if (length(exclude_keywords) > 0) {
    pattern_excl <- str_c(exclude_keywords, collapse = "|")
    terms <- terms %>%
      filter(!str_detect(str_to_lower(.data[[desc_col]]), str_to_lower(pattern_excl)))
  }
  
  terms
}

genes_from_terms <- function(terms_df, geneid_col, gene_sep = "/") {
  if (nrow(terms_df) == 0) return(tibble(gene = character(0)))
  
  terms_df %>%
    transmute(gene = .data[[geneid_col]]) %>%
    mutate(gene = as.character(gene)) %>%
    filter(!is.na(gene), gene != "") %>%
    separate_rows(gene, sep = fixed(gene_sep)) %>%
    mutate(gene = str_trim(gene)) %>%
    filter(gene != "") %>%
    distinct(gene)
}

summarize_for_bar <- function(classified_df, bar_label) {
  classified_df %>%
    dplyr::count(category, name = "n_genes") %>%
    mutate(
      bar = bar_label,
      total_genes = nrow(classified_df),
      pct = 100 * n_genes / total_genes
    )
}

classify_genes_for_stacking <- function(all_gene_set,
                                        nuclear_gene_set,
                                        synaptic_gene_set,
                                        priority = c("synaptic_first", "nuclear_first")) {
  priority <- match.arg(priority)
  
  df <- all_gene_set %>%
    mutate(
      in_nuclear  = gene %in% nuclear_gene_set$gene,
      in_synaptic = gene %in% synaptic_gene_set$gene
    )
  
  if (priority == "nuclear_first") {
    df %>%
      mutate(category = case_when(
        in_nuclear  ~ "nuclear_terms",
        in_synaptic ~ "synaptic_terms",
        TRUE        ~ "other"
      ))
  } else {
    df %>%
      mutate(category = case_when(
        in_synaptic ~ "synaptic_terms",
        in_nuclear  ~ "nuclear_terms",
        TRUE        ~ "other"
      ))
  }
}

get_genes_from_go_terms <- function(
    go_df,
    terms,
    res_df,
    lfc_col,
    pvalue_col,
    geneid_col = "geneID",
    desc_col = "Description",
    gene_sep = "/",
    top_n = 25,
    rank_by = c("pvalue", "abs_lfc")
) {
  rank_by <- match.arg(rank_by)
  
  if (is.null(go_df) || nrow(go_df) == 0 || length(terms) == 0) {
    return(character(0))
  }
  
  keep_go <- go_df %>%
    dplyr::filter(.data[[desc_col]] %in% terms)
  
  if (nrow(keep_go) == 0) {
    return(character(0))
  }
  
  gene_df <- keep_go %>%
    dplyr::select(
      term = dplyr::all_of(desc_col),
      genes = dplyr::all_of(geneid_col)
    ) %>%
    dplyr::mutate(genes = as.character(genes)) %>%
    tidyr::separate_rows(genes, sep = stringr::fixed(gene_sep)) %>%
    dplyr::mutate(
      genes = stringr::str_trim(genes)
    ) %>%
    dplyr::filter(!is.na(genes), genes != "") %>%
    dplyr::rename(gene_name = genes)
  
  if (nrow(gene_df) == 0) {
    return(character(0))
  }
  
  ranked_df <- gene_df %>%
    dplyr::left_join(
      res_df %>%
        dplyr::select(
          gene_name,
          dplyr::all_of(c(lfc_col, pvalue_col))
        ),
      by = "gene_name"
    ) %>%
    dplyr::filter(!is.na(.data[[lfc_col]]), !is.na(.data[[pvalue_col]]))
  
  if (nrow(ranked_df) == 0) {
    return(character(0))
  }
  
  if (rank_by == "pvalue") {
    ranked_df <- ranked_df %>%
      dplyr::arrange(.data[[pvalue_col]], dplyr::desc(abs(.data[[lfc_col]])))
  } else {
    ranked_df <- ranked_df %>%
      dplyr::arrange(dplyr::desc(abs(.data[[lfc_col]])), .data[[pvalue_col]])
  }
  
  ranked_df %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(gene_name)
}

run_one_direction <- function(
    go_df,
    genes_df,
    direction_label,
    priority = c("synaptic_first", "nuclear_first"),
    comp_label = NULL,            
    Outdirectory = ".",
    desc_col = "Description",
    geneid_col = "geneID",
    gene_sep = "/",
    nuclear_keywords,
    exclude_nuclear_keywords = character(0),
    synaptic_keywords,
    exclude_synaptic_keywords = character(0)
) {
  
  priority <- match.arg(priority)
  
  # genes_df should be a 1-col tibble/dataframe with gene symbols
  if (ncol(genes_df) != 1) {
    stop(paste0(direction_label, ": genes_df must be a 1-column table of gene IDs/symbols."))
  }
  
  # clean gene set
  all_gene_set <- tibble(gene = unique(as.character(genes_df[[1]]))) %>%
    mutate(gene = str_trim(gene)) %>%
    filter(!is.na(gene), gene != "")
  
  # --- Terms flagged by keywords ---
  nuclear_terms <- flag_terms_by_keywords(
    go_df, desc_col, geneid_col,
    include_keywords = nuclear_keywords,
    exclude_keywords = exclude_nuclear_keywords
  )
  
  synaptic_terms <- flag_terms_by_keywords(
    go_df, desc_col, geneid_col,
    include_keywords = synaptic_keywords,
    exclude_keywords = exclude_synaptic_keywords
  )
  
  # --- Gene sets from those terms ---
  nuclear_gene_set  <- genes_from_terms(nuclear_terms,  geneid_col, gene_sep)
  synaptic_gene_set <- genes_from_terms(synaptic_terms, geneid_col, gene_sep)
  
  # --- Gene-level classification (mutually exclusive category for stacking) ---
  gene_classification <- classify_genes_for_stacking(
    all_gene_set,
    nuclear_gene_set,
    synaptic_gene_set,
    priority = priority
  )
  
  # --- filenames: include comp_label if provided so you don't overwrite ---
  prefix <- direction_label
  if (!is.null(comp_label) && nzchar(comp_label)) {
    prefix <- paste0(comp_label, "_", direction_label)
  }
  # Save transparency tables
  dir.create(Outdirectory, showWarnings = FALSE, recursive = TRUE)
  
  # Summarize for stacked bar
  bar_df <- summarize_for_bar(gene_classification, bar_label = prefix)
  
  list(
    comp_label = comp_label,
    direction_label = direction_label,
    priority = priority,
    nuclear_terms = nuclear_terms,
    synaptic_terms = synaptic_terms,
    nuclear_gene_set = nuclear_gene_set,
    synaptic_gene_set = synaptic_gene_set,
    gene_classification = gene_classification,
    bar_df = bar_df
  )
}

RT_out_list <- list()

for (comp in names(RT_specs)) {
  
  spec   <- RT_specs[[comp]]
  prefix <- spec$prefix
  
  message("Running RT block: ", comp)
  
  res <- spec$results %>%
    dplyr::rename(
      !!paste0(prefix, "_log2FC") := log2FoldChange,
      !!paste0(prefix, "_pvalue") := pvalue,
      !!paste0(prefix, "_padj")   := padj
    ) %>%
    dplyr::select(-dplyr::any_of(c("baseMean","lfcSE","stat")))
  
  lfc_col    <- paste0(prefix, "_log2FC")
  padj_col   <- paste0(prefix, "_padj")
  status_col <- paste0("RT_", comp, "_Status") 
  
  res <- res %>%
    dplyr::mutate(
      !!status_col := dplyr::case_when(
        .data[[padj_col]] < 0.05 & .data[[lfc_col]] > 0 ~ paste0(prefix, "_Up"),
        .data[[padj_col]] < 0.05 & .data[[lfc_col]] < 0 ~ paste0(prefix, "_Down"),
        TRUE ~ "NS"
      )
    )
  
  up_df   <- res %>% dplyr::filter(.data[[lfc_col]] > 0, .data[[padj_col]] < 0.05)
  down_df <- res %>% dplyr::filter(.data[[lfc_col]] < 0, .data[[padj_col]] < 0.05)
  
  vol <- ggplot(res, aes(x = .data[[lfc_col]], y = -log10(.data[[padj_col]]), color = .data[[status_col]])) +
    geom_point(alpha = 1, size = 1) +
    scale_color_manual(
      values = setNames(
        c("red3","dodgerblue","gray70"),
        c(paste0(prefix,"_Up"), paste0(prefix,"_Down"), "NS")
      )
    ) +
    theme_classic(base_size = 9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray50") +
    labs(y = "-log10(p-adjusted)", x = paste0("Log2 Fold Change RiboTag\n", comp)) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Decreased (", nrow(down_df), ")"),
             color = "dodgerblue", hjust = 1.2, vjust = 2, size = 3) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Increased (", nrow(up_df), ")"),
             color = "red3", hjust = 1.2, vjust = 4, size = 3) +
    theme(legend.position = "none", axis.text = element_text(color = "black")) +
    xlim(spec$volcano_xlim[1], spec$volcano_xlim[2])
  print(vol)
  ggsave(file.path(Outdirectory, paste0("VolcanoPlot_", prefix,"_",comp, "_padj_unlabeled.pdf")),
         vol, width = 4, height = 4, units = "in")
  
  #GO enrichment (Up/Down)
  go_up <- clusterProfiler::enrichGO(
    gene = up_df$gene_name,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    readable = TRUE,
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "none"
  )@result %>%
    dplyr::filter(!is.na(pvalue), pvalue < 0.05)
  
  go_down <- clusterProfiler::enrichGO(
    gene = down_df$gene_name,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    readable = TRUE,
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "none"
  )@result %>%
    dplyr::filter(!is.na(pvalue), pvalue < 0.05)
  
  up_out <- run_one_direction(
    go_df = go_up,
    genes_df = up_df %>% dplyr::select(gene_name),
    direction_label = paste0(prefix, "_Up"),
    priority = spec$prio_up,
    comp_label = comp,
    Outdirectory = Outdirectory,
    nuclear_keywords = nuclear_keywords,
    exclude_nuclear_keywords = exclude_nuclear_keywords,
    synaptic_keywords = synaptic_keywords,
    exclude_synaptic_keywords = exclude_synaptic_keywords
  )
  
  down_out <- run_one_direction(
    go_df = go_down,
    genes_df = down_df %>% dplyr::select(gene_name),
    direction_label = paste0(prefix, "_Down"),
    priority = spec$prio_down,
    comp_label = comp,
    Outdirectory = Outdirectory,
    nuclear_keywords = nuclear_keywords,
    exclude_nuclear_keywords = exclude_nuclear_keywords,
    synaptic_keywords = synaptic_keywords,
    exclude_synaptic_keywords = exclude_synaptic_keywords
  )
  
  stack_df <- dplyr::bind_rows(up_out$bar_df, down_out$bar_df) %>%
    dplyr::mutate(category = factor(category, levels = c("synaptic_terms", "nuclear_terms", "other")))
  
  CH_merge <- merge(res, CH_RT_resdata, by.x = "gene_name", by.y = "gene_name")
  
  RT_out_list[[comp]] <- list(
    res = res,
    status_col = status_col,
    up_df = up_df,
    down_df = down_df,
    go_up = go_up,
    go_down = go_down,
    nucsyn_up = up_out,
    nucsyn_down = down_out,
    stack_df = stack_df,
    CH_merge = CH_merge
  )
}

nuc_syn_fill_colors <- c(
  "synaptic_terms" = "#1B9E77",
  "nuclear_terms"  = "#D95F02",
  "other"          = "#7570B3"
)

stacked_tables_rt <- list()

for (comp in names(RT_out_list)) {
  
  stack_df <- RT_out_list[[comp]]$stack_df
  if (is.null(stack_df) || nrow(stack_df) == 0) {
    message("No stack_df found for RT comp: ", comp)
    next
  }
  
  stack_df2 <- stack_df %>%
    mutate(
      comp = comp,
      direction = case_when(
        str_detect(bar, "_Down$") | str_detect(bar, "Down") ~ "Down",
        str_detect(bar, "_Up$")   | str_detect(bar, "Up")   ~ "Up",
        TRUE ~ bar
      ),
      category  = factor(category, levels = c("synaptic_terms", "nuclear_terms", "other")),
      direction = factor(direction, levels = c("Down", "Up"))
    )
  
  # Plot
  p_stack <- ggplot(stack_df2, aes(x = direction, y = pct, fill = category)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = nuc_syn_fill_colors, drop = FALSE) +
    theme_classic(base_size = 9) +
    labs(x = NULL, y = "% of genes", fill = NULL, title = paste0("RiboTag: ", comp))
  print(p_stack)
  
  stacked_tables_rt[[comp]] <- stack_df2
}

# Faceted stacked bar
if (length(stacked_tables_rt) > 0) {
  RT_all_stack_df <- bind_rows(stacked_tables_rt) %>%
    mutate(comp = factor(comp, levels = c("5v0", "30v5")))
  
  RT_all_stack_faceted <- ggplot(RT_all_stack_df, aes(x = direction, y = pct, fill = category)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = nuc_syn_fill_colors, drop = FALSE) +
    theme_classic(base_size = 9) +
    labs(x = NULL, y = "% of genes", fill = NULL) +
    facet_wrap(~comp, nrow = 1)
  print(RT_all_stack_faceted)
  ggsave(
    filename = file.path(Outdirectory, "RT_ALL_UpDown_stackedbar_nuclear_synaptic_other_faceted.pdf"),
    plot = RT_all_stack_faceted, width = 9, height = 4, units = "in"
  )
}

print(RT_all_stack_df)

#category       n_genes bar               total_genes   pct comp  direction
#nuclear_terms       13 5min_RiboTag_Up           637  2.04 5v0  Up       
#other              329 5min_RiboTag_Up           637 51.6  5v0  Up       
#synaptic_terms     295 5min_RiboTag_Up           637 46.3  5v0  Up       
#nuclear_terms      151 5min_RiboTag_Down         458 33.0  5v0  Down     
#other              235 5min_RiboTag_Down         458 51.3  5v0  Down     
#synaptic_terms      72 5min_RiboTag_Down         458 15.7  5v0  Down     
#nuclear_terms      225 Int_RiboTag_Up            640 35.2  30v5   Up       
#other              353 Int_RiboTag_Up            640 55.2  30v5   Up       
#synaptic_terms      62 Int_RiboTag_Up            640  9.69 30v5   Up       
#nuclear_terms       14 Int_RiboTag_Down          717  1.95 30v5   Down     
#other              374 Int_RiboTag_Down          717 52.2  30v5   Down     
#synaptic_terms     329 Int_RiboTag_Down          717 45.9  30v5   Down    

################################################################################################################################
#make bubble plot
nuclear_terms <- c(
  "mRNA processing",
  "nuclear transport",
  "nuclear export",
  "regulation of DNA-binding transcription factor activity",
  "protein localization to nucleus",
  "RNA splicing",
  "chromatin remodeling",
  "miRNA transcription"
)

synaptic_terms <- c(
  "neurotransmitter secretion",
  "calcium ion homeostasis",
  "protein localization to synapse",
  "synapse assembly",
  "regulation of synaptic plasticity",
  "synaptic vesicle cycle",
  "synaptic transmission, glutamatergic",
  "excitatory synapse assembly"
)

target_terms <- c(nuclear_terms, synaptic_terms)

make_go_count_bubble_df <- function(RT_out_list, target_terms) {
  
  go_tables <- list(
    "5v0_Down" = RT_out_list[["5v0"]]$go_down,
    "5v0_Up"   = RT_out_list[["5v0"]]$go_up,
    "30v5_Down"  = RT_out_list[["30v5"]]$go_down,
    "30v5_Up"    = RT_out_list[["30v5"]]$go_up
  )
  
  out <- lapply(names(go_tables), function(cat_name) {
    
    go_df <- go_tables[[cat_name]]
    
    if (is.null(go_df) || nrow(go_df) == 0) {
      return(
        tibble::tibble(
          Description = target_terms,
          Count = 0,
          p.adjust = NA_real_,
          category = cat_name
        )
      )
    }
    
    matched <- go_df %>%
      dplyr::filter(Description %in% target_terms) %>%
      dplyr::select(Description, Count, p.adjust) %>%
      dplyr::mutate(category = cat_name)
    
    missing_terms <- setdiff(target_terms, matched$Description)
    
    if (length(missing_terms) > 0) {
      matched <- dplyr::bind_rows(
        matched,
        tibble::tibble(
          Description = missing_terms,
          Count = 0,
          p.adjust = NA_real_,
          category = cat_name
        )
      )
    }
    
    matched
  })
  
  dplyr::bind_rows(out)
}

RT_go_bubble_df <- make_go_count_bubble_df(
  RT_out_list = RT_out_list,
  target_terms = target_terms
)

x_order <- c("5v0_Down", "5v0_Up", "30v5_Down", "30v5_Up")
term_order <- c(nuclear_terms, synaptic_terms)

RT_go_bubble_df <- RT_go_bubble_df %>%
  dplyr::mutate(
    biology = dplyr::case_when(
      Description %in% nuclear_terms  ~ "nuclear",
      Description %in% synaptic_terms ~ "synaptic",
      TRUE ~ "other"
    ),
    neglog10_padj_raw = -log10(p.adjust),
    neglog10_padj = pmin(neglog10_padj_raw, 12),
    category = factor(category, levels = x_order),
    Description = factor(Description, levels = rev(term_order))
  )

RT_go_bubble_plot_df <- RT_go_bubble_df %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05,
    Count > 0
  )

RT_p_go_bubble <- ggplot(
  RT_go_bubble_plot_df,
  aes(
    x = category,
    y = Description,
    size = Count,
    fill = neglog10_padj
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.25
  ) +
  scale_size_continuous(
    name = "Gene count",
    breaks = c(5,10, 20, 40)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "adj p",
    limits = c(1.3, 12),
    breaks = c( 2, 4, 6, 8, 10, 12),
    labels = c("0.01", "1e-4", "1e-6", "1e-8", "1e-10", "≤1e-12"),
    guide = guide_colorbar(barheight = grid::unit(4, "cm"))
  ) +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )
RT_p_go_bubble
ggsave(file.path(Outdirectory, "RT_selected_GO_terms_bubbleplot.pdf"),RT_p_go_bubble,width = 5.5, height = 5.5,  units = "in",  device = grDevices::cairo_pdf)

################################################################################################################################
#label volcano plots
make_labeled_volcano <- function(
    res,
    lfc_col,
    pvalue_col,
    status_col,
    prefix,
    comp,
    volcano_xlim,
    up_df,
    down_df,
    label_up_genes = character(0),
    label_down_genes = character(0)
) {
  
  plot_df <- res %>%
    dplyr::filter(!is.na(.data[[lfc_col]]), !is.na(.data[[pvalue_col]])) %>%
    dplyr::mutate(
      neglog10_pvalue = -log10(.data[[pvalue_col]])
    )
  
  p <- ggplot(
    plot_df,
    aes(x = .data[[lfc_col]], y = neglog10_pvalue, color = .data[[status_col]])
  ) +
    geom_point(alpha = 1, size = 1) +
    scale_color_manual(
      values = setNames(
        c("red3", "dodgerblue", "gray70"),
        c(paste0(prefix, "_Up"), paste0(prefix, "_Down"), "NS")
      )
    ) +
    theme_classic(base_size = 9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray50") +
    labs(
      y = "-log10(p-adj)",
      x = paste0("Log2 Fold Change RiboTag\n", comp)
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("Decreased (", nrow(down_df), ")"),
      color = "dodgerblue", hjust = 1.2, vjust = 2, size = 3
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("Increased (", nrow(up_df), ")"),
      color = "red3", hjust = 1.2, vjust = 4, size = 3
    ) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black")
    ) +
    xlim(volcano_xlim[1], volcano_xlim[2]) +
    scale_y_continuous(trans = "log1p") +
    coord_cartesian(ylim = c(0, 20))
  
  label_df_up <- plot_df %>%
    dplyr::filter(gene_name %in% label_up_genes)
  
  label_df_down <- plot_df %>%
    dplyr::filter(gene_name %in% label_down_genes)
  
  if (nrow(label_df_down) > 0) {
    p <- p +
      ggrepel::geom_label_repel(
        data = label_df_down,
        aes(
          x = .data[[lfc_col]],
          y = neglog10_pvalue,
          label = gene_name
        ),
        nudge_x = -0.6,
        nudge_y = 0.4,
        color = "dodgerblue",
        fill = "white",
        label.size = 0.15,     # thin border
        label.padding = grid::unit(0.08, "lines"),
        size = 2.7,
        box.padding = 0.25,
        point.padding = 0.15,
        segment.color = "gray40",
        segment.size = 0.2,
        min.segment.length = 0.05,
        max.overlaps = Inf,
        show.legend = FALSE
      )
  }
  
  if (nrow(label_df_up) > 0) {
    p <- p +
      ggrepel::geom_label_repel(
        data = label_df_up,
        aes(
          x = .data[[lfc_col]],
          y = neglog10_pvalue,
          label = gene_name
        ),
        nudge_x = 0.6,
        nudge_y = 0.4,
        color = "red3",
        fill = "white",
        label.size = 0.15,
        label.padding = grid::unit(0.08, "lines"),
        size = 2.7,
        box.padding = 0.25,
        point.padding = 0.15,
        segment.color = "gray40",
        segment.size = 0.2,
        min.segment.length = 0.05,
        max.overlaps = Inf,
        show.legend = FALSE
      )
  }
  p
}

for (comp in names(RT_out_list)) {
  
  spec <- RT_specs[[comp]]
  
  res        <- RT_out_list[[comp]]$res
  go_up      <- RT_out_list[[comp]]$go_up
  go_down    <- RT_out_list[[comp]]$go_down
  up_df      <- RT_out_list[[comp]]$up_df
  down_df    <- RT_out_list[[comp]]$down_df
  status_col <- RT_out_list[[comp]]$status_col
  
  prefix   <- spec$prefix
  lfc_col  <- paste0(prefix, "_log2FC")
  padj_col <- paste0(prefix, "_padj")
  
  up_terms   <- spec$volcano_label_up_terms
  down_terms <- spec$volcano_label_down_terms
  
  up_label_genes <- get_genes_from_go_terms(
    go_df = go_up,
    terms = up_terms,
    res_df = res,
    lfc_col = lfc_col,
    pvalue_col = padj_col,
    geneid_col = geneid_col,
    gene_sep = gene_sep,
    top_n = spec$volcano_label_top_n,
    rank_by = "pvalue"
  )
  
  down_label_genes <- get_genes_from_go_terms(
    go_df = go_down,
    terms = down_terms,
    res_df = res,
    lfc_col = lfc_col,
    pvalue_col = padj_col,
    geneid_col = geneid_col,
    gene_sep = gene_sep,
    top_n = spec$volcano_label_top_n,
    rank_by = "pvalue"
  )
  
  message(comp, ": labeling up genes from GO terms: ", paste(up_terms, collapse = ", "))
  message(comp, ": found ", length(up_label_genes), " up genes to label")
  
  message(comp, ": labeling down genes from GO terms: ", paste(down_terms, collapse = ", "))
  message(comp, ": found ", length(down_label_genes), " down genes to label")
  
  vol_labeled <- make_labeled_volcano(
    res = res,
    lfc_col = lfc_col,
    pvalue_col = padj_col,
    status_col = status_col,
    prefix = prefix,
    comp = comp,
    volcano_xlim = spec$volcano_xlim,
    up_df = up_df,
    down_df = down_df,
    label_up_genes = up_label_genes,
    label_down_genes = down_label_genes
  )
  
  print(vol_labeled)
  
  ggsave(
    file.path(
      Outdirectory,
      paste0("VolcanoPlot_", prefix, "_", comp, "_padj_labeled.pdf")
    ),
    vol_labeled,
    width = 4.5,
    height = 4,
    units = "in"
  )
}

################################################################################################################################
#violin plots of my RT data vs Hale et al 2021 RT data from microdissected CA1 soma and neuropils
get_rt_group_colors <- function(groups) {
  out <- setNames(rep("gray70", length(groups)), groups)
  out[grepl("_Up$", groups)]   <- "red3"
  out[grepl("_Down$", groups)] <- "dodgerblue"
  out
}

RT_violin_tables_CH_NPvsCB <- list()

for (comp in names(RT_out_list)) {
  
  message("Making CH violin for RT comp: ", comp)
  
  df <- RT_out_list[[comp]]$CH_merge
  if (is.null(df) || nrow(df) == 0) {
    message("  - No CH_merge; skipping ", comp)
    next
  }
  
  status_col <- RT_out_list[[comp]]$status_col
  if (is.null(status_col) || !status_col %in% colnames(df)) {
    message("  - status_col missing in CH_merge; skipping ", comp)
    next
  }
  
  ycol <- "CH_NPvsCB_log2FC"
  if (!ycol %in% colnames(df)) {
    message("  - ", ycol, " not found; skipping ", comp)
    next
  }
  
  df2 <- df %>%
    filter(.data[[status_col]] != "NS") %>%
    mutate(
      Group = as.character(.data[[status_col]]),
      Group = factor(
        Group,
        levels = c(
          unique(Group[grepl("_Down$", Group)]),
          unique(Group[grepl("_Up$", Group)])
        )
      )
    ) %>%
    filter(!is.na(.data[[ycol]]))
  
  if (nrow(df2) == 0 || nlevels(df2$Group) < 2) {
    message("  - Not enough groups/data to plot for ", comp)
    next
  }
  
  pw <- df2 %>%
    rstatix::pairwise_wilcox_test(
      formula = as.formula(paste0(ycol, " ~ Group")),
      p.adjust.method = "none"
    )
  
  sig_pairs <- pw %>%
    filter(p.adj < 0.05) %>%
    mutate(p.signif = rstatix::p_format(p.adj, stars = TRUE))
  
  sig_comparisons <- sig_pairs %>%
    dplyr::select(group1, group2) %>%
    split(seq_len(nrow(.))) %>%
    lapply(function(x) c(as.character(x$group1), as.character(x$group2)))
  
  grp_levels <- levels(df2$Group)
  fill_cols <- get_rt_group_colors(grp_levels)
  
  p <- ggplot(df2, aes(x = Group, y = .data[[ycol]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
    geom_boxplot(
      width = 0.15,
      fill = "white",
      color = "black",
      outlier.shape = NA,
      alpha = 0.9
    ) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    labs(
      y = "CH RiboTag log2 Fold Change (NP / CB)",
      x = NULL,
      title = paste0("RiboTag ", comp, ": Up vs Down")
    ) +
    theme_classic(base_size = 9) +
    theme(
      text = element_text(size = 9),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "none"
    )
  if (length(sig_comparisons) > 0) {
    p_sig <- p +
      ggpubr::stat_compare_means(
        comparisons = sig_comparisons,
        method = "wilcox.test",
        label = "p.signif",
        tip.length = 0.01
      )
    print(p_sig)
  } else {
    message("  - No significant pairwise tests at p<0.05 for ", comp, " (still saving plot).")
  }
  
  # ---- store for faceted plot ----
  RT_violin_tables_CH_NPvsCB[[comp]] <-
    df2 %>%
    dplyr::mutate(
      comp = comp,
      Group_simple = dplyr::case_when(
        grepl("_Down$", Group) ~ "Down",
        grepl("_Up$",   Group) ~ "Up",
        TRUE ~ NA_character_
      ),
      Group_simple = factor(Group_simple, levels = c("Down","Up"))
    ) %>%
    dplyr::filter(!is.na(Group_simple))
  
}

if (length(RT_violin_tables_CH_NPvsCB) > 0) {
  
  CH_all_violin_df <- dplyr::bind_rows(RT_violin_tables_CH_NPvsCB) %>%
    dplyr::mutate(comp = factor(comp, levels = c("5v0","30v5")))
  
  fill_cols_all <- c("Down" = "dodgerblue", "Up" = "red3")
  
  # --- stats per facet (comp) ---
  pw_facet <- CH_all_violin_df %>%
    dplyr::group_by(comp) %>%
    rstatix::wilcox_test(CH_NPvsCB_log2FC ~ Group_simple) %>%
    rstatix::add_significance("p") %>%
    rstatix::add_xy_position(x = "Group_simple")
  
  # add a little headroom per facet for brackets
  pw_facet <- pw_facet %>%
    dplyr::group_by(comp) %>%
    dplyr::mutate(
      y.position = max(CH_all_violin_df$CH_NPvsCB_log2FC[CH_all_violin_df$comp == unique(comp)], na.rm = TRUE) + 0.5
    ) %>%
    dplyr::ungroup()
  
  CH_all_violin_plot <-
    ggplot2::ggplot(
      CH_all_violin_df,
      ggplot2::aes(x = Group_simple, y = CH_NPvsCB_log2FC, fill = Group_simple)
    ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
    ggplot2::geom_boxplot(
      width = 0.15, fill = "white", color = "black",
      outlier.shape = NA, alpha = 0.9
    ) +
    ggplot2::scale_fill_manual(values = fill_cols_all, drop = FALSE) +
    ggplot2::facet_wrap(~comp, nrow = 1) +
    ggpubr::stat_pvalue_manual(
      pw_facet,
      label = "p.signif",   # stars
      tip.length = 0.01,
      hide.ns = TRUE
    ) +
    ggplot2::labs(
      y = "CH CLIPscore log2 Fold Change (NP / CB)",
      x = NULL
    ) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 9),
      legend.position = "none"
    )
  
  print(CH_all_violin_plot)
  
  ggplot2::ggsave(file.path(Outdirectory, "ViolinPlot_RTvsCH_ALLcomps_faceted_withStats.pdf"),CH_all_violin_plot,width = 7, height = 4, units = "in"  )
}

################################################################################################################################
#comparison of Up and Down genes from each of the 3 comparisons 
comp_order <- c("5v0", "30v5")

get_gene_vec <- function(RT_out_list, comp, which = c("up", "down")) {
  which <- match.arg(which)
  
  if (!comp %in% names(RT_out_list)) return(character(0))
  
  x <- if (which == "up") RT_out_list[[comp]]$up_df else RT_out_list[[comp]]$down_df
  
  if (is.null(x)) return(character(0))
  
  if (is.data.frame(x)) {
    if (!"gene_name" %in% colnames(x)) return(character(0))
    return(unique(as.character(x$gene_name)))
  } else {
    return(unique(as.character(x)))
  }
}

venn_sets_all <- list()

for (comp in comp_order) {
  venn_sets_all[[paste0(comp, "_Up")]]   <- get_gene_vec(RT_out_list, comp, "up")
  venn_sets_all[[paste0(comp, "_Down")]] <- get_gene_vec(RT_out_list, comp, "down")
}

# Drop empty sets
venn_sets_all <- venn_sets_all[lengths(venn_sets_all) > 0]

if (length(venn_sets_all) < 2) {
  stop("Not enough non-empty gene sets for venndetail(). Check p<0.05 thresholds / comps.")
}

vd_all <- venndetail(venn_sets_all)

pdf(file.path(Outdirectory, "VennDetail_RT_AllComparisons_UpDown.pdf"), width = 10, height = 6)
plot(vd_all, type = "upset")
dev.off()

plot(vd_all, type = "upset")

################################################################################################################################
#Bring in CLIP STAR Transcriptome files
CLIP_STAR_files <- grep(
  "transcriptome.unique.bed$",
  list.files(FMRP_CLIP_STAR_Directory),
  value = TRUE
)

transcripts <- vector("list", length(CLIP_STAR_files))

length(CLIP_STAR_files) #16

for (i in seq_along(CLIP_STAR_files)) {
  
  df <- readr::read_tsv(
    file.path(FMRP_CLIP_STAR_Directory, CLIP_STAR_files[[i]]),
    col_names = FALSE,
    col_types = readr::cols(
      X1 = readr::col_character(),
      X2 = readr::col_integer(),
      X3 = readr::col_integer(),
      X4 = readr::col_character(),
      X5 = readr::col_double(),
      X6 = readr::col_character()
    )
  )
  
  colnames(df) <- c("Transcript_id", "Txstart", "Txend", "read_ID", "Txscore", "Txstrand")
  
  file_name <- sub("\\..*$", "", CLIP_STAR_files[[i]])
  
  df <- df %>%
    dplyr::mutate(
      Sample = file_name,
      Condition = {
        parts <- strsplit(file_name, "_")[[1]]
        paste(tail(parts, 2), collapse = "_")
      }
    )
  
  transcripts[[i]] <- df
}

STAR_Tx <- dplyr::bind_rows(transcripts)

# join gene_name onto transcriptome tags
FMRP_STAR_genes <- STAR_Tx %>%
  dplyr::left_join(mm10_Tx_final %>% dplyr::select(transcript_id,gene_name,width), by = c("Transcript_id" = "transcript_id")) %>%
  dplyr::filter(!is.na(gene_name))

condition_levels <- c("Control_5min","Control_30min", "ChR2_5min", "ChR2_30min")
region_levels <- c("5UTR", "CDS", "3UTR")

STAR_Tx_Crepos <- FMRP_STAR_genes %>%
  dplyr::filter(Condition %in% condition_levels) %>%  
  dplyr::mutate(
    Condition = factor(Condition, levels = condition_levels)
  ) %>%
  dplyr::arrange(Condition, Sample) %>%
  dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)))

# annotate regions on transcriptome coordinates
annotate_transcriptome_tags <- function(tags_df, tx_lengths) {
  tags_df %>%
    dplyr::left_join(
      tx_lengths %>% dplyr::select(tx_name, utr5_len, cds_len, utr3_len, tx_len),
      by = c("Transcript_id" = "tx_name")
    ) %>%
    dplyr::mutate(
      region = dplyr::case_when(
        Txend <= utr5_len ~ "5UTR",
        Txend <= utr5_len + cds_len ~ "CDS",
        Txend <= utr5_len + cds_len + utr3_len ~ "3UTR",
        TRUE ~ "other"
      )
    )
}

STAR_annot_df <- annotate_transcriptome_tags(STAR_Tx_Crepos, mm10_txdb_lengths)

STAR_annot_df_filtered <- STAR_annot_df %>%
  dplyr::filter(region %in% region_levels, !is.na(gene_name)) %>%
  dplyr::mutate(
    region    = factor(region, levels = region_levels),
    Condition = factor(Condition, levels = condition_levels),
    Sample    = factor(Sample, levels = levels(STAR_Tx_Crepos$Sample))
  )

nrow(STAR_annot_df_filtered) #1798900

df_TPM_use <- STAR_annot_df_filtered %>%
  dplyr::filter(
    !is.na(gene_name), !is.na(Sample), !is.na(Condition),
    as.character(Condition) %in% condition_levels
  )

tag_counts <- df_TPM_use %>%
  dplyr::group_by(gene_name, Sample) %>%
  dplyr::summarise(tag_count = dplyr::n(), .groups = "drop")

FMRP_STARcounts <- tag_counts %>%
  tidyr::pivot_wider(names_from = Sample, values_from = tag_count, values_fill = 0)

gene_lengths_df <- mm10_Tx_final %>% dplyr::select(transcript_id,gene_name,width) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(gene_length_max = max(width), .groups = "drop")

FMRP_STARcounts_lengths <- dplyr::left_join(FMRP_STARcounts, gene_lengths_df, by = "gene_name")

cond_cols <- lapply(condition_levels, function(cond) {
  grep(paste0(cond, "$"), colnames(FMRP_STARcounts_lengths), value = TRUE)
})

names(cond_cols) <- condition_levels

FMRP_flagged <- FMRP_STARcounts_lengths %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    keep_BC = any(sapply(condition_levels, function(cond) {
      cols <- cond_cols[[cond]]
      length(cols) > 0 && all(c_across(dplyr::all_of(cols)) >= 2)
    }))
  ) %>%
  dplyr::ungroup()

FMRP_STARcounts_BC4 <- FMRP_flagged %>%
  dplyr::filter(keep_BC) %>%
  dplyr::select(-keep_BC)

nrow(FMRP_STARcounts_BC4) #4812

FMRP_STARcounts_BC4=FMRP_STARcounts_BC4 %>% relocate("gene_name","gene_length_max",
                                                     "RAS119_FMRP_CLIP_Control_5min",
                                                     "RAS127_FMRP_CLIP_Control_5min",
                                                     "RAS130_FMRP_CLIP_Control_5min",
                                                     "RAS131_FMRP_CLIP_Control_5min",
                                                     "RAS120_FMRP_CLIP_ChR2_5min",
                                                     "RAS128_FMRP_CLIP_ChR2_5min",
                                                     "RAS132_FMRP_CLIP_ChR2_5min",
                                                     "RAS133_FMRP_CLIP_ChR2_5min",
                                                     "RAS81_FMRP_CLIP_Control_30min",
                                                     "RAS91_FMRP_CLIP_Control_30min",
                                                     "RAS101_FMRP_CLIP_Control_30min",
                                                     "RAS102_FMRP_CLIP_Control_30min",
                                                     "RAS80_FMRP_CLIP_ChR2_30min",
                                                     "RAS90_FMRP_CLIP_ChR2_30min",
                                                     "RAS98_FMRP_CLIP_ChR2_30min",
                                                     "RAS100_FMRP_CLIP_ChR2_30min")

#calculate TPM for CLIP 
FMRP_STAR_TPM <- FMRP_STARcounts_BC4 %>%
  pivot_longer(names_to = "Sample", values_to = "tags", cols = starts_with("RAS")) %>%
  distinct %>%
  mutate(rpk = tags/(gene_length_max/1000))

#get scaling factors for read normalization
scaling <- FMRP_STAR_TPM %>%
  group_by(Sample) %>%
  summarise(scaling.factor = (sum(rpk))/1000000)

FMRP_STAR_TPM <- left_join(FMRP_STAR_TPM, scaling, by = "Sample")

FMRP_STAR_TPM <- FMRP_STAR_TPM %>%
  mutate(tpm = rpk/scaling.factor) 

#this should equal 1 million for each sample
FMRP_STAR_TPM %>% group_by(Sample) %>% summarise(sum(tpm))

FMRP_STAR_TPM_wider=FMRP_STAR_TPM %>% 
  dplyr::select(gene_name,gene_length_max,Sample,tags,tpm) %>%
  pivot_wider(names_from = Sample, values_from = c(tags,tpm), values_fill = 0)  

################################################################################################################################
#enrichment of RT dataset for FMRP targets
subset_vs_rest_p <- function(subset_genes, universe_genes, target_genes) {
  subset_genes <- intersect(unique(subset_genes), universe_genes)
  target_genes <- intersect(unique(target_genes), universe_genes)
  
  a <- length(intersect(subset_genes, target_genes))                      # target in subset
  b <- length(setdiff(subset_genes, target_genes))                        # non-target in subset
  c <- length(setdiff(target_genes, subset_genes))                        # target not in subset
  d <- length(setdiff(universe_genes, union(subset_genes, target_genes))) # neither
  
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  exp_min <- suppressWarnings(min(chisq.test(mat)$expected))
  if (is.finite(exp_min) && exp_min >= 5) {
    suppressWarnings(chisq.test(mat, correct = FALSE)$p.value)
  } else {
    fisher.test(mat)$p.value
  }
}

p_to_stars <- function(p) {
  if (is.na(p)) "NA"
  else if (p < 1e-4) "****"
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 0.05) "*"
  else "ns"
}

# Targets (global)
FMRP_targets_all <- unique(FMRP_STAR_TPM_wider$gene_name)

enrich_summary <- list()

for (comp in names(RT_out_list)) {
  
  message("FMRP enrichment plot: ", comp)
  
  res   <- RT_out_list[[comp]]$res
  up_df <- RT_out_list[[comp]]$up_df
  down_df <- RT_out_list[[comp]]$down_df
  
  if (is.null(res) || !("gene_name" %in% colnames(res))) {
    message("  - Missing res/gene_name for ", comp, "; skipping.")
    next
  }
  
  # Universe: all genes present in this RT results table
  U <- unique(as.character(res$gene_name))
  U <- U[!is.na(U) & U != ""]
  
  # Subsets
  RT_up   <- if (!is.null(up_df)   && "gene_name" %in% colnames(up_df))   unique(as.character(up_df$gene_name))   else character(0)
  RT_down <- if (!is.null(down_df) && "gene_name" %in% colnames(down_df)) unique(as.character(down_df$gene_name)) else character(0)
  
  # Restrict to universe
  RT_up   <- intersect(RT_up, U)
  RT_down <- intersect(RT_down, U)
  FMRP_targets <- intersect(FMRP_targets_all, U)
  
  # p-values
  p_up   <- subset_vs_rest_p(RT_up,   U, FMRP_targets)
  p_down <- subset_vs_rest_p(RT_down, U, FMRP_targets)
  
  # For transparency, also save overlap counts
  a_up   <- length(intersect(RT_up, FMRP_targets))
  a_down <- length(intersect(RT_down, FMRP_targets))
  
  # ----------------------------
  # Build bar_df (stacked percent)
  # ----------------------------
  bar_df <- bind_rows(
    tibble(group = "All RT expressed", gene = U),
    tibble(group = "RT Up",            gene = RT_up),
    tibble(group = "RT Down",          gene = RT_down)
  ) %>%
    mutate(class = ifelse(gene %in% FMRP_targets, "FMRP targets", "non-targets")) %>%
    dplyr::count(group, class, name = "n") %>%
    group_by(group) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(
      group = factor(group, levels = c("All RT expressed", "RT Down", "RT Up")),
      class = factor(class, levels = c("non-targets", "FMRP targets")),
      group_col = case_when(
        group == "All RT expressed" ~ "grey60",
        group == "RT Up"            ~ "red3",
        group == "RT Down"          ~ "dodgerblue"
      ),
      fill_col = ifelse(class == "FMRP targets", group_col, "white")
    )
  
  outline_df <- bar_df %>%
    distinct(group, group_col) %>%
    mutate(percent = 100)
  
  label_df <- tibble(
    group = factor(c("RT Up", "RT Down"), levels = levels(bar_df$group)),
    y = 104,
    label = c(p_to_stars(p_up), p_to_stars(p_down))
  )
  
  p <- ggplot() +
    geom_col(
      data = bar_df,
      aes(x = group, y = percent, fill = fill_col),
      width = 0.75,
      color = NA
    ) +
    geom_col(
      data = outline_df,
      aes(x = group, y = percent, color = group_col),
      width = 0.75,
      fill = NA,
      linewidth = 0.9
    ) +
    geom_text(
      data = label_df,
      aes(x = group, y = y, label = label),
      inherit.aes = FALSE,
      size = 5
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_y_continuous(
      limits = c(0, 110),
      breaks = c(0, 25, 50, 75, 100),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Percent",
      title = paste0("FMRP target enrichment: RiboTag ", comp)
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 12, 5.5, 5.5)
    )
  print(p)
  out_pdf <- file.path(Outdirectory, paste0("RiboTag_", comp, "_subset_FMRPtargets_enrichment_plot.pdf"))
  ggsave(out_pdf, p, width = 4, height = 6, units = "in")
  
  # Save a compact summary row
  enrich_summary[[comp]] <- tibble(
    comp = comp,
    universe_n = length(U),
    targets_in_universe_n = length(FMRP_targets),
    RT_up_n = length(RT_up),
    RT_down_n = length(RT_down),
    targets_in_RT_up_n = a_up,
    targets_in_RT_down_n = a_down,
    p_up = p_up,
    p_down = p_down
  )
  
}

################################################################################################################################
#now merge RiboTag datasets
RT_5v0_res  <- RT_out_list[["5v0"]]$res
RT_30v5_res   <- RT_out_list[["30v5"]]$res

RT_5v0_res <- RT_5v0_res %>%
  dplyr::rename(
    RiboTag_5v0_log2FC = RiboTag_log2FC,
    RiboTag_5v0_pvalue = RiboTag_pvalue,
    RiboTag_5v0_padj   = RiboTag_padj)

RT_30v5_res <- RT_30v5_res %>%
  dplyr::rename(
    RiboTag_30v5_log2FC = RiboTag_log2FC,
    RiboTag_30v5_pvalue = RiboTag_pvalue,
    RiboTag_30v5_padj   = RiboTag_padj)

colnames(RT_5v0_res)
colnames(RT_30v5_res)

RT_results_merged_all <- RT_5v0_res %>%
  dplyr::select(gene_name,
                dplyr::starts_with("RiboTag_5v0")) %>%
  dplyr::left_join(
    RT_30v5_res %>%
      dplyr::select(gene_name,
                    dplyr::starts_with("RiboTag")),
    by = "gene_name"
  )

################################################################################################################################
colnames(RT_results_merged_all)

nrow(RT_results_merged_all) #15242
nrow(FMRP_STAR_TPM_wider) #4812

#optional to write these files out
write.csv(RT_results_merged_all,file.path(Outdirectory,"Opto_RiboTag_dataset.csv"), row.names=F)
write.csv(FMRP_STAR_TPM_wider,file.path(Outdirectory,"Opto_FMRP_CLIP_dataset.csv"), row.names=F)

FMRP_RiboTag_merged <- merge(RT_results_merged_all , FMRP_STAR_TPM_wider, by.x= "gene_name", by.y= "gene_name", all.x=TRUE, all.y=FALSE) 

nrow(FMRP_RiboTag_merged) #15242

#save this df as an rds object for use in downstream analyses
saveRDS(FMRP_RiboTag_merged, file=file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/FMRP_RiboTag_merged.rds"))





