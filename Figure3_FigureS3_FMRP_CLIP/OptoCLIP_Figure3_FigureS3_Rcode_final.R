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
Base_directory=file.path("~/ruthasinger_github/OptoCLIP_April2026/Figure3_FigureS3_FMRP_CLIP")
list.files(Base_directory)

Data_directory=file.path(Base_directory,"Data")
list.files(Data_directory)

FMRP_CLIP_STAR_Directory=file.path(Data_directory,"Data_FMRP_CLIP_Tx_STAR")
list.files(FMRP_CLIP_STAR_Directory)

Outdirectory=file.path(Base_directory,paste("Output",Sys.Date(),sep="_"))
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#code for how I made these RDS files are in the "RDS_files" folder in this github repository

mm10_Tx_final <- readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/Data_mm10gtf","mm10_Tx_final.rds"))
mm10_txdb_lengths <- readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/Data_mm10gtf","mm10_txdb_lengths.rds"))

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
#Bring in CLIP STAR Transcriptome files
CLIP_STAR_files <- grep(
  "transcriptome.unique.bed$",
  list.files(FMRP_CLIP_STAR_Directory),
  value = TRUE
)

transcripts <- vector("list", length(CLIP_STAR_files))

length(CLIP_STAR_files) #15

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

# condition order for plotting/stats
condition_levels_plot <- c(
  "noCre_noCre",
  "Control_Control",
  "ChR2_5min",
  "ChR2_30min"
)

STAR_Tx_plot <- STAR_Tx %>%
  dplyr::filter(Condition %in% condition_levels_plot) %>%
  dplyr::mutate(
    Condition_plot = factor(Condition, levels = condition_levels_plot)
  )

# Stats: choose comparisons to annotate
comparisons_to_plot <- list(
  c("ChR2_5min", "ChR2_30min"),
  c("Control_Control", "ChR2_5min"),
  c("Control_Control", "ChR2_30min")
)

df_counts_tx_sample <- STAR_Tx_plot %>%
  dplyr::count(Sample, Condition_plot, name = "n_tags")

all_pairwise_tbl <- df_counts_tx_sample %>%
  rstatix::pairwise_t_test(n_tags ~ Condition_plot, p.adjust.method = "none") %>%
  dplyr::arrange(p)

pairs_df <- tibble::tibble(
  group1 = sapply(comparisons_to_plot, `[`, 1),
  group2 = sapply(comparisons_to_plot, `[`, 2)
)

stats_to_plot <- all_pairwise_tbl %>%
  dplyr::semi_join(pairs_df, by = c("group1", "group2")) %>%
  dplyr::bind_rows(
    all_pairwise_tbl %>%
      dplyr::semi_join(pairs_df, by = c("group1" = "group2", "group2" = "group1"))
  ) %>%
  dplyr::distinct(group1, group2, .keep_all = TRUE) %>%
  rstatix::add_xy_position(x = "Condition_plot", step.increase = 0.02)

# Plot colors PER CONDITION
condition_colors <- c(
  "noCre_noCre"      = "grey50",
  "Control_Control"  = "grey35",
  "ChR2_5min"        = "dodgerblue",
  "ChR2_30min"       = "red3"
)

FMRP_CLIP_STAR_tx_per_condition_bar_plot_stats <- ggplot(
  df_counts_tx_sample,
  aes(x = Condition_plot, y = n_tags, fill = Condition_plot)
) +
  stat_summary(fun = mean, geom = "col", color = "black", width = 0.8) +
  stat_summary(
    fun.data = ggplot2::mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_point(size = 1) +
  ggpubr::stat_pvalue_manual(
    stats_to_plot,
    label = "p.signif",
    tip.length = 0.01,
    size = 6,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  labs(x = NULL, y = "FMRP CLIP transcriptome tags") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
FMRP_CLIP_STAR_tx_per_condition_bar_plot_stats
ggsave(file.path(Outdirectory, "FMRP_CLIP_STAR_tx_per_condition_bar_plot_stats.pdf"), plot = FMRP_CLIP_STAR_tx_per_condition_bar_plot_stats,width = 6, height = 6)

##############################################################################################################
# join gene_name onto transcriptome tags
FMRP_STAR_genes <- STAR_Tx %>%
  dplyr::left_join(mm10_Tx_final %>% dplyr::select(transcript_id,gene_name,width), by = c("Transcript_id" = "transcript_id")) %>%
  dplyr::filter(!is.na(gene_name))

# Cre+ conditions you want downstream
condition_levels <- c("Control_Control", "ChR2_5min", "ChR2_30min")
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

nrow(STAR_annot_df_filtered) #1204858

##############################################################################################################
#make a count matrix by collapsing by gene
gene_counts <- STAR_annot_df_filtered %>%
  filter(!is.na(gene_name)) %>%
  group_by(Sample, gene_name) %>%
  summarise(Tag_Count = n(), .groups = "drop")

CLIP_count_matrix <- gene_counts %>%
  pivot_wider(names_from = Sample, values_from = Tag_Count, values_fill = 0) %>%
  column_to_rownames("gene_name")

CLIP_count_matrix <- as.data.frame(CLIP_count_matrix) 

short_names <- sapply(colnames(CLIP_count_matrix), function(x) {
  parts <- strsplit(x, "_")[[1]]
  paste(parts[c(4,5)], collapse = " ")
})
short_names

CLIP_count_matrix_short <- CLIP_count_matrix
colnames(CLIP_count_matrix_short) <- short_names

corr_matrix=cor(CLIP_count_matrix_short)

testRes = cor.mtest(CLIP_count_matrix_short, conf.level = 0.95)

pdf(file = file.path(Outdirectory,"STAR_Corrplot_pvalues.pdf"))
corrplot::corrplot(corr_matrix, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,tl.col= "black",
         insig = 'label_sig', pch.col = 'white')
dev.off()

corrplot(corr_matrix, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,tl.col= "black",
         insig = 'label_sig', pch.col = 'white')

##############################################################################################################
# Meta-transcript coverage plot
bin_size_scaled <- 0.005
smooth_k <- 3

roll_mean <- function(x, k = 3) {
  n <- length(x)
  out <- numeric(n)
  half <- floor(k / 2)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }
  out
}

# plot order + colors 
condition_levels_plot <- c("Control_Control", "ChR2_5min", "ChR2_30min")

condition_cols_plot <- c(
  "Control_Control" = "grey70",
  "ChR2_5min"     = "dodgerblue",
  "ChR2_30min"    = "red3"
)

df_rel <- STAR_annot_df_filtered %>%
  dplyr::filter(
    region %in% c("5UTR","CDS","3UTR"),
    as.character(Condition) %in% condition_levels_plot
  ) %>%
  dplyr::mutate(
    center = (Txstart + Txend) / 2,
    region = factor(region, levels = c("5UTR","CDS","3UTR")),
    Condition = factor(as.character(Condition), levels = condition_levels_plot),
    rel_pos = dplyr::case_when(
      region == "5UTR" ~ -(center / utr5_len),
      region == "CDS"  ~ (center - utr5_len) / cds_len,
      region == "3UTR" ~ 1 + (center - utr5_len - cds_len) / utr3_len
    ),
    rel_pos = pmin(rel_pos,  1.999),
    rel_pos = pmax(rel_pos, -0.999)
  )

# Top transcripts per Condition (not pooled)
top_tx <- df_rel %>%
  dplyr::group_by(Condition, Transcript_id) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Condition) %>%
  dplyr::slice_max(n, n = 1000, with_ties = FALSE) %>%
  dplyr::ungroup()

df_rel_top <- df_rel %>%
  dplyr::semi_join(top_tx, by = c("Condition","Transcript_id"))

df_bins <- df_rel_top %>%
  dplyr::mutate(bin = floor(rel_pos / bin_size_scaled) * bin_size_scaled) %>%
  dplyr::count(Condition, bin) %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(rel_cov = n / max(n)) %>%
  dplyr::arrange(bin) %>%
  dplyr::mutate(smooth_cov = roll_mean(rel_cov, k = smooth_k)) %>%
  dplyr::ungroup()

meta_plot <- ggplot(df_bins, aes(x = bin, y = smooth_cov, color = Condition)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = condition_cols_plot, drop = FALSE) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "grey80", alpha = 0.4) +
  scale_x_continuous(breaks = c(-0.5, 0, 1, 1.5)) +
  xlim(-0.5, 1.5) +
  labs(x = "Meta-transcript position", y = "Relative CLIP coverage", color = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top")
meta_plot
ggsave(file.path(Outdirectory, "FMRP_CLIP_STAR_MetaTranscript.pdf"), meta_plot, width = 5, height = 5)

##############################################################################################################
# Region ratios boxplot + Wilcoxon tests (use Condition_plot so Control is pooled)
conds_keep <- c("Control_Control", "ChR2_5min", "ChR2_30min")
region_levels <- c("5UTR", "CDS", "3UTR")

# plot-specific subset
df_use <- STAR_annot_df_filtered %>%
  dplyr::filter(
    region %in% region_levels,
    as.character(Condition) %in% conds_keep
  ) %>%
  dplyr::mutate(
    Condition = factor(as.character(Condition), levels = conds_keep),
    region    = factor(region, levels = region_levels)
  )

cond_cols_plot <- c(
  "Control_Control" = "grey70",
  "ChR2_5min"     = "dodgerblue",
  "ChR2_30min"    = "red3"
)

#top transcripts per condition 
top_tx <- df_use %>%
  dplyr::group_by(Condition, Transcript_id) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Condition) %>%
  dplyr::slice_max(order_by = n, n = 1000, with_ties = FALSE) %>%
  dplyr::ungroup()

#region ratios per transcript (within each condition)
region_ratios_all <- df_use %>%
  dplyr::semi_join(top_tx, by = c("Condition", "Transcript_id")) %>%
  dplyr::group_by(Transcript_id, Condition, region) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Transcript_id, Condition) %>%
  dplyr::mutate(total_tags = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ratio = n / total_tags)

#stats
pairwise_stats <- region_ratios_all %>%
  rstatix::group_by(region) %>%
  rstatix::wilcox_test(ratio ~ Condition, p.adjust.method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::filter(p.adj < 0.05) %>%
  rstatix::add_y_position(fun = "max", step.increase = 0.08)

print(pairwise_stats[,c(1,3,4,5,6,7,8,9,10)])

#region group1          group2        n1    n2 statistic        p    p.adj p.adj.signif
#5UTR   Control_Control ChR2_30min   365   244    30964. 1.82e-10 5.46e-10 ****        
#5UTR   ChR2_5min       ChR2_30min   407   244    36661  2.23e- 8 3.34e- 8 ****        
#CDS    Control_Control ChR2_5min    968   975   427067  2.88e- 4 4.32e- 4 ***         
#CDS    ChR2_5min       ChR2_30min   975   947   527286. 6.84e- 8 2.05e- 7 ****        
#3UTR   Control_Control ChR2_5min    942   934   486714. 6.62e- 5 9.93e- 5 ****        
#3UTR   Control_Control ChR2_30min   942   930   399840. 1   e- 3 1   e- 3 ***         
#3UTR   ChR2_5min       ChR2_30min   934   930   349840  3.59e-13 1.08e-12 **** 

# plot
plot_regions_box_all <- ggplot(
  region_ratios_all,
  aes(x = Condition, y = ratio, fill = Condition)) +
  geom_boxplot(width = 0.55, color = "black", alpha = 1, outlier.shape = NA) +
  facet_wrap(~ region, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = cond_cols_plot, drop = FALSE) +
  ggpubr::stat_pvalue_manual(
    pairwise_stats,
    label = "p.adj.signif",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = NULL, y = "Proportion of tags per transcript", title = NULL)
plot_regions_box_all
ggsave(file.path(Outdirectory, "FMRP_CLIP_STAR_RegionProportions_boxplot_stats.pdf"),plot_regions_box_all, width = 5, height = 5)

################################################################################################################################################################
#now plotting by BC
condition_levels <- c("Control_Control", "ChR2_5min", "ChR2_30min")
region_levels <- c("5UTR", "CDS", "3UTR")

# plot/calculation-specific subset (prevents NA Condition from factor())
df_BC_use <- STAR_annot_df_filtered %>%
  dplyr::filter(
    region %in% region_levels,
    !is.na(gene_name), !is.na(Sample), !is.na(Condition),
    as.character(Condition) %in% condition_levels
  ) %>%
  dplyr::mutate(
    region    = factor(region, levels = region_levels),
    Condition = factor(as.character(Condition), levels = condition_levels)
  )

#replicate-level counts per gene x region x condition x sample
tx_counts <- df_BC_use %>%
  dplyr::group_by(gene_name, region, Condition, Sample) %>%
  dplyr::summarise(tags = dplyr::n(), .groups = "drop")

#BC per gene x region x condition = number of samples with >=1 tag
tx_BC <- tx_counts %>%
  dplyr::group_by(gene_name, region, Condition) %>%
  dplyr::summarise(BC = dplyr::n_distinct(Sample[tags > 0]), .groups = "drop")

#don't hard-code 1:4; derive max BC from your data
max_BC <- max(tx_BC$BC, na.rm = TRUE)

toPlot_longer_use <- tx_BC %>%
  dplyr::filter(BC > 0) %>%
  dplyr::mutate(
    BC      = factor(BC, levels = seq_len(max_BC)),
    region  = factor(region, levels = region_levels),
    Condition = factor(Condition, levels = condition_levels))

FMRP_CLIP_STAR_byregion_BC_plot <- ggplot(
  toPlot_longer_use,
  aes(x = BC, fill = region)) +
  geom_bar(position = "fill", stat = "count", color = "black") +
  facet_wrap(~ Condition, nrow = 1, drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Number of Biological Replicates (BC)") +
  ylab("Fraction of Genes\n(Not normalized)") +
  scale_fill_manual(
    values = c("5UTR" = "#E7B800", "CDS" = "#FC4E07", "3UTR" = "#00AFBB"),
    name = "Region"
  ) +
  theme_classic(base_size = 14)
FMRP_CLIP_STAR_byregion_BC_plot
ggsave(file.path(Outdirectory, "FMRP_CLIP_STAR_byregion_BC_plot.pdf"),plot = FMRP_CLIP_STAR_byregion_BC_plot, width = 6, height = 5)


