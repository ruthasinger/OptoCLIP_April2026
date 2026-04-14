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
Base_directory=file.path("~/Dropbox (Personal)/Darnell_lab/optoCLIP/OptoCLIP_manuscript/Neuron/Raw_data_for_paper/Code_and_data/Figure5_FigureS5_Figure6_FigureS6")
list.files(Base_directory)

Outdirectory=file.path(Base_directory,paste("Output",Sys.Date(),sep="_"))
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#code for how CH_RT_resdata was made is in the "RDS_files" folder in this github repository

CH_RT_resdata=readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/Data_CHale_et_al_2021_RiboTag_salmon/CH_RT_resdata.rds"))

#code for making FMRP_RiboTag_merged is in the following script: OptoCLIP_Figure4_FigureS4_Rcode_final.R in this github repository

FMRP_RiboTag_merged=readRDS(file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/FMRP_RiboTag_merged.rds"))

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

nrow(FMRP_RiboTag_merged) #15242

FMRP_RiboTag_merged <- FMRP_RiboTag_merged %>%
  mutate(
    across(
      c(starts_with("tags_"), starts_with("tpm_")),
      ~ replace_na(., 0)
    )
  )

FMRP_RiboTag_merged=FMRP_RiboTag_merged %>% 
  ungroup() %>%
  mutate(CLIP_Control_5min_rep1_TPM=case_when(tpm_RAS119_FMRP_CLIP_Control_5min >= 1 ~ tpm_RAS119_FMRP_CLIP_Control_5min, tpm_RAS119_FMRP_CLIP_Control_5min < 1 ~ 1),
         CLIP_Control_5min_rep2_TPM=case_when(tpm_RAS127_FMRP_CLIP_Control_5min >= 1 ~ tpm_RAS127_FMRP_CLIP_Control_5min, tpm_RAS127_FMRP_CLIP_Control_5min < 1 ~ 1),
         CLIP_Control_5min_rep3_TPM=case_when(tpm_RAS130_FMRP_CLIP_Control_5min >= 1 ~ tpm_RAS130_FMRP_CLIP_Control_5min, tpm_RAS130_FMRP_CLIP_Control_5min < 1 ~ 1),
         CLIP_Control_5min_rep4_TPM=case_when(tpm_RAS131_FMRP_CLIP_Control_5min >= 1 ~ tpm_RAS131_FMRP_CLIP_Control_5min, tpm_RAS131_FMRP_CLIP_Control_5min < 1 ~ 1),
         CLIP_Opto_5min_rep1_TPM=case_when(tpm_RAS120_FMRP_CLIP_ChR2_5min >= 1 ~ tpm_RAS120_FMRP_CLIP_ChR2_5min, tpm_RAS120_FMRP_CLIP_ChR2_5min < 1 ~ 1),
         CLIP_Opto_5min_rep2_TPM=case_when(tpm_RAS128_FMRP_CLIP_ChR2_5min >= 1 ~ tpm_RAS128_FMRP_CLIP_ChR2_5min, tpm_RAS128_FMRP_CLIP_ChR2_5min < 1 ~ 1),
         CLIP_Opto_5min_rep3_TPM=case_when(tpm_RAS132_FMRP_CLIP_ChR2_5min >= 1 ~ tpm_RAS132_FMRP_CLIP_ChR2_5min, tpm_RAS132_FMRP_CLIP_ChR2_5min < 1 ~ 1),
         CLIP_Opto_5min_rep4_TPM=case_when(tpm_RAS133_FMRP_CLIP_ChR2_5min >= 1 ~ tpm_RAS133_FMRP_CLIP_ChR2_5min, tpm_RAS133_FMRP_CLIP_ChR2_5min < 1 ~ 1),
         CLIP_Control_30min_rep1_TPM=case_when(tpm_RAS81_FMRP_CLIP_Control_30min >= 1 ~ tpm_RAS81_FMRP_CLIP_Control_30min, tpm_RAS81_FMRP_CLIP_Control_30min < 1 ~ 1),
         CLIP_Control_30min_rep2_TPM=case_when(tpm_RAS91_FMRP_CLIP_Control_30min >= 1 ~ tpm_RAS91_FMRP_CLIP_Control_30min, tpm_RAS91_FMRP_CLIP_Control_30min < 1 ~ 1),
         CLIP_Control_30min_rep3_TPM=case_when(tpm_RAS101_FMRP_CLIP_Control_30min >= 1 ~ tpm_RAS101_FMRP_CLIP_Control_30min, tpm_RAS101_FMRP_CLIP_Control_30min < 1 ~ 1),
         CLIP_Control_30min_rep4_TPM=case_when(tpm_RAS102_FMRP_CLIP_Control_30min >= 1 ~ tpm_RAS102_FMRP_CLIP_Control_30min, tpm_RAS102_FMRP_CLIP_Control_30min < 1 ~ 1),
         CLIP_Opto_30min_rep1_TPM=case_when(tpm_RAS80_FMRP_CLIP_ChR2_30min >= 1 ~ tpm_RAS80_FMRP_CLIP_ChR2_30min, tpm_RAS80_FMRP_CLIP_ChR2_30min < 1 ~ 1),
         CLIP_Opto_30min_rep2_TPM=case_when(tpm_RAS90_FMRP_CLIP_ChR2_30min >= 1 ~ tpm_RAS90_FMRP_CLIP_ChR2_30min, tpm_RAS90_FMRP_CLIP_ChR2_30min < 1 ~ 1),
         CLIP_Opto_30min_rep3_TPM=case_when(tpm_RAS98_FMRP_CLIP_ChR2_30min >= 1 ~ tpm_RAS98_FMRP_CLIP_ChR2_30min, tpm_RAS98_FMRP_CLIP_ChR2_30min < 1 ~ 1),
         CLIP_Opto_30min_rep4_TPM=case_when(tpm_RAS100_FMRP_CLIP_ChR2_30min >= 1 ~ tpm_RAS100_FMRP_CLIP_ChR2_30min, tpm_RAS100_FMRP_CLIP_ChR2_30min < 1 ~ 1))

colnames(FMRP_RiboTag_merged)

FMRP_RiboTag_merged <- FMRP_RiboTag_merged %>%
  mutate(
    avg_RT_Control_5min_TPM       = rowMeans(across(c("RiboTag_Exp03_M28_Control_5min_IP_salmon.TPM",
                                                      "RiboTag_Exp03_M29_Control_5min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M40_Control_5min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M41_Control_5min_IP_salmon.TPM")),na.rm = TRUE),
    avg_RT_Opto_5min_TPM          = rowMeans(across(c("RiboTag_Exp03_M18_ChR2_5min_IP_salmon.TPM",
                                                      "RiboTag_Exp03_M19_ChR2_5min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M36_ChR2_5min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M37_ChR2_5min_IP_salmon.TPM")),na.rm = TRUE),
    avg_RT_Control_30min_TPM      = rowMeans(across(c("RiboTag_Exp03_M26_Control_30min_IP_salmon.TPM",
                                                      "RiboTag_Exp03_M27_Control_30min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M38_Control_30min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M39_Control_30min_IP_salmon.TPM")), na.rm = TRUE ),
    avg_RT_Opto_30min_TPM         = rowMeans(across(c("RiboTag_Exp03_M16_ChR2_30min_IP_salmon.TPM",
                                                      "RiboTag_Exp03_M17_ChR2_30min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M34_ChR2_30min_IP_salmon.TPM",
                                                      "RiboTag_Exp04_M35_ChR2_30min_IP_salmon.TPM")),na.rm = TRUE),
    avg_CLIP_Control_5min_TPM     = rowMeans(across(c("CLIP_Control_5min_rep1_TPM",
                                                      "CLIP_Control_5min_rep2_TPM",
                                                      "CLIP_Control_5min_rep3_TPM",
                                                      "CLIP_Control_5min_rep4_TPM")),na.rm = TRUE),
    avg_CLIP_Opto_5min_TPM        = rowMeans(across(c("CLIP_Opto_5min_rep1_TPM",
                                                      "CLIP_Opto_5min_rep2_TPM",
                                                      "CLIP_Opto_5min_rep3_TPM",
                                                      "CLIP_Opto_5min_rep4_TPM")),na.rm = TRUE),
    avg_CLIP_Control_30min_TPM    = rowMeans(across(c("CLIP_Control_30min_rep1_TPM",
                                                      "CLIP_Control_30min_rep2_TPM",
                                                      "CLIP_Control_30min_rep3_TPM",
                                                      "CLIP_Control_30min_rep4_TPM")),na.rm = TRUE),
    avg_CLIP_Opto_30min_TPM       = rowMeans(across(c("CLIP_Opto_30min_rep1_TPM",
                                                      "CLIP_Opto_30min_rep2_TPM",
                                                      "CLIP_Opto_30min_rep3_TPM",
                                                      "CLIP_Opto_30min_rep4_TPM")),na.rm = TRUE))

FMRP_RiboTag_merged_log2= FMRP_RiboTag_merged %>% 
  mutate(
    across(
      .cols = c(avg_RT_Control_5min_TPM,
                avg_RT_Opto_5min_TPM,
                avg_RT_Control_30min_TPM,
                avg_RT_Opto_30min_TPM,
                avg_CLIP_Control_5min_TPM, 
                avg_CLIP_Opto_5min_TPM, 
                avg_CLIP_Control_30min_TPM, 
                avg_CLIP_Opto_30min_TPM,
                CLIP_Control_5min_rep1_TPM, 
                CLIP_Control_5min_rep2_TPM,
                CLIP_Control_5min_rep3_TPM, 
                CLIP_Control_5min_rep4_TPM,
                CLIP_Opto_5min_rep1_TPM,
                CLIP_Opto_5min_rep2_TPM,
                CLIP_Opto_5min_rep3_TPM,
                CLIP_Opto_5min_rep4_TPM,
                CLIP_Control_30min_rep1_TPM,
                CLIP_Control_30min_rep2_TPM,
                CLIP_Control_30min_rep3_TPM,
                CLIP_Control_30min_rep4_TPM,
                CLIP_Opto_30min_rep1_TPM,
                CLIP_Opto_30min_rep2_TPM,
                CLIP_Opto_30min_rep3_TPM,
                CLIP_Opto_30min_rep4_TPM),
      .fns = log2,
      .names = "{.col}_log2"
    ))

################################################################################################################################################################
#CLIP score generation 
FMRP_RiboTag_merged_5min=FMRP_RiboTag_merged_log2 %>% dplyr::filter(avg_RT_Control_5min_TPM_log2 > 0 & avg_RT_Opto_5min_TPM_log2 > 0) 
FMRP_RiboTag_merged_30min=FMRP_RiboTag_merged_log2 %>% dplyr::filter(avg_RT_Control_30min_TPM_log2 > 0 & avg_RT_Opto_30min_TPM_log2 > 0 ) 

CLIPscore_plots_directory=file.path(Outdirectory,"CLIPscore_plots")
dir.create(CLIPscore_plots_directory,showWarnings = TRUE)

#Control- 5 min clip score code
df <- FMRP_RiboTag_merged_5min

my_x <- "avg_RT_Control_5min_TPM_log2"

#rep1
Treatment <- "Control"
Time <- "5min"
Rep <- "rep1"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep2
Rep <- "rep2"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep3
Treatment <- "Control"
Time <- "5min"
Rep <- "rep3"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep4
Rep <- "rep4"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#ChR2- 5 min clip score code
my_x <- "avg_RT_Opto_5min_TPM_log2"

#rep1
Treatment <- "Opto"
Time <- "5min"
Rep <- "rep1"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep2
Rep <- "rep2"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep3
Rep <- "rep3"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep4
Rep <- "rep4"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

CLIPscores_5min=df

#Control- 30 min clip score code

df <- FMRP_RiboTag_merged_30min

my_x <- "avg_RT_Control_30min_TPM_log2"

#rep1
Treatment <- "Control"
Time <- "30min"
Rep <- "rep1"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0("log2 ", Treatment, " RiboTag"),
    ylab = paste0("log2 ", Treatment, " ",Rep, " CLIP")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep2
Rep <- "rep2"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0("log2 ", Treatment, " RiboTag"),
    ylab = paste0("log2 ", Treatment, " ",Rep, " CLIP")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep3
Rep <- "rep3"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0("log2 ", Treatment, " RiboTag"),
    ylab = paste0("log2 ", Treatment, " ",Rep, " CLIP")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep4
Rep <- "rep4"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0("log2 ", Treatment, " RiboTag"),
    ylab = paste0("log2 ", Treatment, " ",Rep, " CLIP")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#ChR2- 30 min clip score code
my_x <- "avg_RT_Opto_30min_TPM_log2"

#rep1
Treatment <- "Opto"
Time <- "30min"
Rep <- "rep1"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep2
Rep <- "rep2"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep3
Rep <- "rep3"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep4
Rep <- "rep4"
my_y <- paste("CLIP",Treatment,Time,Rep,"TPM_log2", sep="_")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggpubr::ggscatter(
    x = my_x, y = my_y,
    color = "black", shape = 19, size = .2,
    add = "reg.line",
    add.params = list(color = "red"),
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "pearson",
      # move label location upward for small panels
      label.x.npc = "left",
      label.y.npc = "top",
      label.sep = ",",
      size = 2.5               # smaller text
    ),
    xlab = paste0(Treatment, " RiboTag log2 TPM"),
    ylab = paste0(Treatment, " ",Rep, " CLIP log2 TPM")
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic(base_size = 6) +          # globally shrink fonts
  theme(
    axis.text  = element_text(size = 5),  # axis tick labels
    axis.title = element_text(size = 6),  # axis title
    plot.margin = margin(1, 1, 1, 1),     # remove outer whitespace
    line = element_line(linewidth = .2)   # thin axis lines
  )
myplot_pearson
ggsave(file.path(CLIPscore_plots_directory,paste0("CLIP_", Treatment, "_", Time,"_",Rep, "_CLIPscores_pearson.pdf")), plot = myplot_pearson,  width = 3, height = 3,  units = "in")

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_",Time,"_",Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

CLIPscores_30min=df

################################################################################################################################################################
#make some plots and analyze CLIP score data 

CLIPscores_final <- CLIPscores_5min %>%
  left_join(
    CLIPscores_30min %>%
      dplyr::select(gene_name, ends_with("_CLIPscore")),
    by = "gene_name"
  )

# identify which columns are the CLIPscore replicate columns
clip_cols <- grep("CLIPscore", colnames(CLIPscores_final))

# replace NA with 0 only in those columns
CLIPscores_final <- CLIPscores_final %>% mutate(across(all_of(clip_cols), ~ replace_na(., 0)))

CLIPscores_merged <- CLIPscores_final %>%
  mutate(
    avg_control_5min_CLIPscores  = rowMeans(
      across(c(
        "Control_5min_rep1_CLIPscore","Control_5min_rep2_CLIPscore",
        "Control_5min_rep3_CLIPscore","Control_5min_rep4_CLIPscore"
      )),
      na.rm = TRUE
    ),
    
    avg_opto_5min_CLIPscores     = rowMeans(
      across(c(
        "Opto_5min_rep1_CLIPscore","Opto_5min_rep2_CLIPscore",
        "Opto_5min_rep3_CLIPscore","Opto_5min_rep4_CLIPscore"
      )),
      na.rm = TRUE
    ),
    
    avg_control_30min_CLIPscores = rowMeans(
      across(c(
        "Control_30min_rep1_CLIPscore","Control_30min_rep2_CLIPscore",
        "Control_30min_rep3_CLIPscore","Control_30min_rep4_CLIPscore"
      )),
      na.rm = TRUE
    ),
    
    avg_opto_30min_CLIPscores    = rowMeans(
      across(c(
        "Opto_30min_rep1_CLIPscore","Opto_30min_rep2_CLIPscore",
        "Opto_30min_rep3_CLIPscore","Opto_30min_rep4_CLIPscore"
      )),
      na.rm = TRUE
    )
  ) %>% 
  ungroup()

Control_0min_CLIPscore_plot_avg=ggplot(CLIPscores_merged %>% filter(avg_control_5min_CLIPscores <= 0 )) +
  geom_point(aes(x=avg_RT_Control_5min_TPM_log2, y=avg_CLIP_Control_5min_TPM_log2, color="NonFMRP_target"), size =2) + 
  geom_point(data= CLIPscores_merged %>% filter(avg_control_5min_CLIPscores <= 1  & avg_control_5min_CLIPscores > 0), aes(x=avg_RT_Control_5min_TPM_log2, y=avg_CLIP_Control_5min_TPM_log2, color="Low_binding"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% filter(avg_control_5min_CLIPscores <= 2  & avg_control_5min_CLIPscores > 1), aes(x=avg_RT_Control_5min_TPM_log2, y=avg_CLIP_Control_5min_TPM_log2, color="High_binding"), size =2) +
  geom_point(data= CLIPscores_merged %>% filter(avg_control_5min_CLIPscores > 2 ), aes(x=avg_RT_Control_5min_TPM_log2, y=avg_CLIP_Control_5min_TPM_log2, color="Stringent_Target"), size =2) +
  scale_color_manual(name = "FMRP CLIP score", 
                     values = c(NonFMRP_target ="black",
                                Low_binding = "grey",
                                High_binding = "grey90",
                                Stringent_Target = "grey35"), guide = "none") + 
  xlab("Control RiboTag TPM (log2)") +
  ylab("Control FMRP CLIP TPM (log2FC)") +
  xlim(c(0,15)) +
  ylim(c(0,17)) +
  theme_classic(base_size = 18) +
  theme_illustrator_black +
  annotate("text", 
           x = Inf, 
           y = 17, size=4,
           label = paste0("NonFMRP_target (",nrow(CLIPscores_merged %>% filter(avg_control_5min_CLIPscores <= 0 & avg_CLIP_Control_5min_TPM_log2 > 0)),")"),
           color = "black",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = .5, size=4,
           label = paste0("No CLIP tags (",nrow(CLIPscores_merged %>% filter(avg_control_5min_CLIPscores <= 0 & avg_CLIP_Control_5min_TPM_log2 <= 0)),")"),
           color = "black",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 16.5, size=4,
           label = paste0("Low_binding (",nrow(CLIPscores_merged %>% filter(avg_control_5min_CLIPscores <= 1  & avg_control_5min_CLIPscores > 0)),")"),
           color = "grey",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 16, size=4,
           label = paste0("High_binding (",nrow(CLIPscores_merged %>% filter(avg_control_5min_CLIPscores <= 2  & avg_control_5min_CLIPscores > 1)),")"),
           color = "grey90",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 15.5, size=4,
           label = paste0("Stringent_Target (",nrow(CLIPscores_merged %>% filter(avg_control_5min_CLIPscores > 2)),")"),
           color = "grey35",hjust=1,vjust=0) 
Control_0min_CLIPscore_plot_avg
ggsave(filename =file.path(Outdirectory,"CLIPscore_Control_0min_plot_avg.pdf"), plot =Control_0min_CLIPscore_plot_avg,  width=6, height=6, units = "in")

Opto_5min_CLIPscore_plot_avg=ggplot(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores <= 0 )) +
  geom_point(aes(x=avg_RT_Opto_5min_TPM_log2, y=avg_CLIP_Opto_5min_TPM_log2, color="NonFMRP_target"), size =2) + 
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores <= 1  & avg_opto_5min_CLIPscores > 0), aes(x=avg_RT_Opto_5min_TPM_log2, y=avg_CLIP_Opto_5min_TPM_log2, color="Low_binding"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores <= 2  & avg_opto_5min_CLIPscores > 1), aes(x=avg_RT_Opto_5min_TPM_log2, y=avg_CLIP_Opto_5min_TPM_log2, color="High_binding"), size =2) +
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores > 2 ), aes(x=avg_RT_Opto_5min_TPM_log2, y=avg_CLIP_Opto_5min_TPM_log2, color="Stringent_Target"), size =2) +
  scale_color_manual(name = "FMRP CLIP score", 
                     values = c(NonFMRP_target ="black",
                                Low_binding = "grey",
                                High_binding = "lightskyblue1",
                                Stringent_Target = "dodgerblue"), guide = "none" ) + 
  xlab("Opto 5 min RiboTag TPM (log2)") +
  ylab("Opto 5 min FMRP CLIP TPM (log2FC)") +
  xlim(c(0,15)) +
  ylim(c(0,17)) +
  theme_classic(base_size = 18) +
  theme_illustrator_black +
  annotate("text", 
           x = Inf, 
           y = 17, size=4,
           label = paste0("NonFMRP_target (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores <= 0 & avg_CLIP_Opto_5min_TPM_log2 > 0)),")"),
           color = "black",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = .5, size=4,
           label = paste0("No CLIP tags (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores <= 0 & avg_CLIP_Opto_5min_TPM_log2 <= 0)),")"),
           color = "black",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 16.5, size=4,
           label = paste0("Low_binding (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores <= 1  & avg_opto_5min_CLIPscores > 0)),")"),
           color = "grey",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 16, size=4,
           label = paste0("High_binding (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores <= 2  & avg_opto_5min_CLIPscores > 1)),")"),
           color = "lightskyblue1",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 15.5, size=4,
           label = paste0("Stringent_Target (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores > 2)),")"),
           color = "dodgerblue",hjust=1,vjust=0) 
Opto_5min_CLIPscore_plot_avg
ggsave(filename =file.path(Outdirectory,"CLIPscore_Opto_5min_plot_avg.pdf"), plot =Opto_5min_CLIPscore_plot_avg,  width=6, height=6, units = "in")

Opto_30min_CLIPscore_plot_avg=ggplot(CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores <= 0 )) +
  geom_point(aes(x=avg_RT_Opto_30min_TPM_log2, y=avg_CLIP_Opto_30min_TPM_log2, color="NonFMRP_target"), size =2) + 
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores <= 1  & avg_opto_30min_CLIPscores > 0), aes(x=avg_RT_Opto_30min_TPM_log2, y=avg_CLIP_Opto_30min_TPM_log2, color="Low_binding"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores <= 2  & avg_opto_30min_CLIPscores > 1), aes(x=avg_RT_Opto_30min_TPM_log2, y=avg_CLIP_Opto_30min_TPM_log2, color="High_binding"), size =2) +
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores > 2 ), aes(x=avg_RT_Opto_30min_TPM_log2, y=avg_CLIP_Opto_30min_TPM_log2, color="Stringent_Target"), size =2) +
  scale_color_manual(name = "FMRP CLIP score", 
                     values = c(NonFMRP_target ="black",
                                Low_binding = "grey",
                                High_binding = "lightcoral",
                                Stringent_Target = "red3"), guide = "none") + 
  xlab("Opto 30 min RiboTag TPM (log2)") +
  ylab("Opto 30 min FMRP CLIP TPM (log2FC)") +
  xlim(c(0,15)) +
  ylim(c(0,17)) +
  theme_classic(base_size = 18) +
  theme_illustrator_black +
  annotate("text", 
           x = Inf, 
           y = 17, size=4,
           label = paste0("NonFMRP_target (",nrow(CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores <= 0 & avg_CLIP_Opto_30min_TPM_log2 >0 )),")"),
           color = "black",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = .5, size=4,
           label = paste0("No CLIP tags (",nrow(CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores <= 0 & avg_CLIP_Opto_30min_TPM_log2 <= 0)),")"),
           color = "black",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 16.5, size=4,
           label = paste0("Low_binding (",nrow(CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores <= 1  & avg_opto_30min_CLIPscores > 0)),")"),
           color = "grey",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 16, size=4,
           label = paste0("High_binding (",nrow(CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores <= 2  & avg_opto_30min_CLIPscores > 1)),")"),
           color = "lightcoral",hjust=1,vjust=0) +
  annotate("text", 
           x = Inf, 
           y = 15.5, size=4,
           label = paste0("Stringent_Target (",nrow(CLIPscores_merged %>% filter(avg_opto_30min_CLIPscores > 2)),")"),
           color = "red3",hjust=1,vjust=0) 
Opto_30min_CLIPscore_plot_avg
ggsave(filename =file.path(Outdirectory,"Opto_30min_CLIPscore_plot_avg.pdf"), plot =Opto_30min_CLIPscore_plot_avg,  width=6, height=6, units = "in")

Opto5minvsControl0min_CLIPscore_plot_avg=ggplot(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores < 2 & avg_control_5min_CLIPscores < 2 )) + 
  geom_point(aes(x=avg_control_5min_CLIPscores, y=avg_opto_5min_CLIPscores, color="NonFMRP_target"), size =2) + 
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores >= 2 &  avg_control_5min_CLIPscores >= 2 ), aes(x=avg_control_5min_CLIPscores, y=avg_opto_5min_CLIPscores, color="Stringent_Targets_Both"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores >= 2 & avg_control_5min_CLIPscores < 2 ), aes(x=avg_control_5min_CLIPscores, y=avg_opto_5min_CLIPscores, color="Stringent_Targets_5min"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores < 2 & avg_control_5min_CLIPscores >= 2 ), aes(x=avg_control_5min_CLIPscores, y=avg_opto_5min_CLIPscores, color="Stringent_Targets_Control"), size =2 ) +
  scale_color_manual(name = "FMRP CLIP score", 
                     values = c(NonFMRP_target ="black",
                                Stringent_Targets_Both = "green4",
                                Stringent_Targets_5min = "dodgerblue",
                                Stringent_Targets_Control = "grey35"), guide = "none") + 
  xlim(c(-10,16)) +
  ylim(c(-10,16)) +
  theme_classic(base_size = 18) +
  theme_illustrator_black +
  ylab("Average Opto 5 min FMRP CLIP Score") +
  xlab("Average Control FMRP CLIP Score") +
  geom_hline(yintercept = 2, colour = "grey", linewidth=.3) + 
  geom_vline(xintercept = 2, colour = "grey", linewidth=.3) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("Both (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores >= 2 &  avg_control_5min_CLIPscores >= 2)),")"),
           color = "green4",hjust=2,vjust=1) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("5 min (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores >= 2 & avg_control_5min_CLIPscores < 2)),")"),
           color = "dodgerblue",hjust=2,vjust=3) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("Control (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores < 2 & avg_control_5min_CLIPscores >= 2 )),")"),
           color = "grey35",hjust=2,vjust=5) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("Neither (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores < 2 & avg_control_5min_CLIPscores < 2 )),")"),
           color = "black",hjust=2,vjust=7)
Opto5minvsControl0min_CLIPscore_plot_avg
ggsave(filename =file.path(Outdirectory,"Opto5minvsControl0min_CLIPscore_plot_avg.pdf"), plot =Opto5minvsControl0min_CLIPscore_plot_avg, width=6, height=6, units = "in")

Opto5minvsOpto30min_CLIPscore_plot_avg=ggplot(CLIPscores_merged %>% dplyr::filter(avg_opto_5min_CLIPscores < 2 & avg_opto_30min_CLIPscores < 2 )) + 
  geom_point(aes(x=avg_opto_5min_CLIPscores, y=avg_opto_30min_CLIPscores, color="NonFMRP_target"), size =2) + 
  geom_point(data= CLIPscores_merged %>% dplyr::filter(avg_opto_5min_CLIPscores >= 2 &  avg_opto_30min_CLIPscores >= 2 ), aes(x=avg_opto_5min_CLIPscores, y=avg_opto_30min_CLIPscores, color="Stringent_Targets_Both"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% dplyr::filter(avg_opto_5min_CLIPscores >= 2 & avg_opto_30min_CLIPscores < 2 ), aes(x=avg_opto_5min_CLIPscores, y=avg_opto_30min_CLIPscores, color="Stringent_Targets_5min"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% dplyr::filter(avg_opto_5min_CLIPscores < 2 & avg_opto_30min_CLIPscores >= 2 ), aes(x=avg_opto_5min_CLIPscores, y=avg_opto_30min_CLIPscores, color="Stringent_Targets_30min"), size =2 ) +
  scale_color_manual(name = "FMRP CLIP score", 
                     values = c(NonFMRP_target ="black",
                                Stringent_Targets_Both = "green4",
                                Stringent_Targets_5min = "dodgerblue",
                                Stringent_Targets_30min = "red3"), guide = "none") + 
  xlim(c(-10,16)) +
  ylim(c(-10,16)) +
  theme_classic(base_size = 18) +
  theme_illustrator_black +
  xlab("Average Opto 5 min FMRP CLIP Score") +
  ylab("Average Opto 30 min FMRP CLIP Score") +
  geom_hline(yintercept = 2, colour = "grey", linewidth=.3) + 
  geom_vline(xintercept = 2, colour = "grey", linewidth=.3) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("Both (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores >= 2 &  avg_opto_30min_CLIPscores >= 2)),")"),
           color = "green4",hjust=2,vjust=1) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("5 min (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores >= 2 & avg_opto_30min_CLIPscores < 2)),")"),
           color = "dodgerblue",hjust=2,vjust=3) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("30 min (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores < 2 & avg_opto_30min_CLIPscores >= 2 )),")"),
           color = "red3",hjust=2,vjust=5) +
  annotate("text", 
           x = Inf, 
           y = Inf, size=5,
           label = paste0("Neither (",nrow(CLIPscores_merged %>% filter(avg_opto_5min_CLIPscores < 2 & avg_opto_30min_CLIPscores < 2 )),")"),
           color = "black",hjust=2,vjust=7)
Opto5minvsOpto30min_CLIPscore_plot_avg
ggsave(filename =file.path(Outdirectory,"Opto5minvsOpto30min_CLIPscore_plot_avg.pdf"), plot =Opto5minvsOpto30min_CLIPscore_plot_avg, width=6, height=6, units = "in")

#graph to quantify
summarize_bins <- function(df, colname) {
  df %>%
    mutate(bin = case_when(
      .data[[colname]] > 2 ~ ">2",
      .data[[colname]] > 1 ~ "1–2",
      .data[[colname]] > 0 ~ "0–1",
      TRUE ~ "≤0"
    )) %>%
    group_by(bin) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(factor(bin, levels = c(">2", "1–2", "0–1", "≤0"))) %>%
    mutate(variable = colname, .before = 1)
}

# Apply to multiple columns
cols_to_check <- c("avg_control_5min_CLIPscores", 
                   "avg_opto_5min_CLIPscores", 
                   "avg_opto_30min_CLIPscores")

results <- bind_rows(lapply(cols_to_check, function(x) summarize_bins(CLIPscores_merged, x)))

# Ensure bin has the desired order
results$bin <- factor(results$bin, levels = c("≤0", "0–1", "1–2", ">2"))

results <- results %>%
  mutate(
    # 1. Order bins (low → high)
    bin = factor(bin, levels = c("≤0", "0–1", "1–2", ">2")),
    
    # 2. Order variables on x-axis
    variable = factor(
      variable,
      levels = c(
        "avg_control_5min_CLIPscores",
        "avg_control_30min_CLIPscores",
        "avg_opto_5min_CLIPscores",
        "avg_opto_30min_CLIPscores"
      )
    ),
    
    # 3. Extract timepoint from variable
    timepoint = case_when(
      grepl("5min",  variable, ignore.case = TRUE)  ~ "5min",
      grepl("30min", variable, ignore.case = TRUE)  ~ "30min",
      TRUE ~ NA_character_
    ),
    
    # 4. Collapse controls; split opto by timepoint
    condition_group = case_when(
      grepl("control", variable, ignore.case = TRUE) ~ "Control",
      grepl("opto",    variable, ignore.case = TRUE) & timepoint == "5min"  ~ "Opto_5min",
      grepl("opto",    variable, ignore.case = TRUE) & timepoint == "30min" ~ "Opto_30min",
      TRUE ~ "Other"
    ),
    
    # 5. Combined fill group
    fill_group = paste(condition_group, bin, sep = "_"),
    
    # 6. Explicit ordering of fill_group (keeps dodge order consistent)
    fill_group = factor(
      fill_group,
      levels = c(
        # Control low → high
        "Control_≤0", "Control_0–1", "Control_1–2", "Control_>2",
        # Opto 5min low → high
        "Opto_5min_≤0", "Opto_5min_0–1", "Opto_5min_1–2", "Opto_5min_>2",
        # Opto 30min low → high
        "Opto_30min_≤0", "Opto_30min_0–1", "Opto_30min_1–2", "Opto_30min_>2"
      )
    )
  )

bin_cols_control <- c(
  "Control_>2"  = "grey35",
  "Control_1–2" = "grey90",
  "Control_0–1" = "grey",
  "Control_≤0"  = "black"
)

bin_cols_opto_5 <- c(
  "Opto_5min_>2"  = "dodgerblue3",
  "Opto_5min_1–2" = "lightskyblue2",
  "Opto_5min_0–1" = "grey",
  "Opto_5min_≤0"  = "black"
)

bin_cols_opto_30 <- c(
  "Opto_30min_>2"  = "red3",
  "Opto_30min_1–2" = "lightcoral",
  "Opto_30min_0–1" = "grey",
  "Opto_30min_≤0"  = "black"
)

bin_colors <- c(bin_cols_control, bin_cols_opto_5, bin_cols_opto_30)

CLIPscores_quant_plot_avg <- ggplot(
  results,
  aes(x = variable, y = n, fill = fill_group)
) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.3, size = 3
  ) +
  scale_fill_manual(values = bin_colors) +
  labs(
    x = NULL,
    y = "Gene count",
    fill = "Bin"  ) +
  theme_classic(base_size = 14) +
  theme_illustrator_black +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
CLIPscores_quant_plot_avg
ggsave(filename =file.path(Outdirectory,"CLIPscore_per_condition_quant_plot.pdf"), plot =CLIPscores_quant_plot_avg ,width=6, height=6, units = "in")

################################################################################################################################################################
#density plots

condition_panel <- tribble(
  ~Condition,       ~col,
  "Control_0min",   "avg_control_5min_CLIPscores",
  "Opto_5min",      "avg_opto_5min_CLIPscores",
  "Opto_30min",     "avg_opto_30min_CLIPscores"
)

# Colors
condition_cols <- c(
  "Control_0min"  = "gray35",
  "Opto_5min"     = "dodgerblue",
  "Opto_30min"    = "red3"
)

mean_long <- CLIPscores_merged %>%
  dplyr::select(gene_name, all_of(condition_panel$col)) %>%
  pivot_longer(
    cols = all_of(condition_panel$col),
    names_to = "col",
    values_to = "Mean_CLIPscore"
  ) %>%
  left_join(condition_panel, by = "col") %>%
  mutate(
    Condition = factor(Condition, levels = condition_panel$Condition)
  ) %>%
  dplyr::select(gene_name, Condition, Mean_CLIPscore)

Mean_CLIPscores_density_plot <- ggplot(mean_long, aes(Mean_CLIPscore, fill = Condition)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = condition_cols, drop = FALSE) +
  theme_classic(base_size = 9) +
  labs(x = "Mean CLIPscore", y = "Density") +
  theme_illustrator_black +
  coord_cartesian(xlim = c(-10, 6))
Mean_CLIPscores_density_plot
ggsave(file.path(Outdirectory, "Mean_CLIPscores_density_plot.pdf"), Mean_CLIPscores_density_plot, width = 9, height = 6, units = "in")

subset_genes <- CLIPscores_merged %>%
  mutate(across(all_of(condition_panel$col), ~replace_na(., 0))) %>%
  filter(
    avg_control_5min_CLIPscores  >= 2,
    avg_opto_5min_CLIPscores     >= 2,
    avg_opto_30min_CLIPscores    <  2
  ) %>%
  pull(gene_name)

length(subset_genes) #261

mean_long_subset <- mean_long %>%
  filter(gene_name %in% subset_genes)

Mean_CLIPscores_density_subset_plot <- ggplot(mean_long_subset, aes(Mean_CLIPscore, fill = Condition)) + 
  geom_density(alpha = 0.3) + 
  scale_fill_manual(values = condition_cols, drop = FALSE) +
  theme_classic(base_size = 16) + labs(x = "Mean CLIPscore", y = "Density") + 
  theme_illustrator_black 
Mean_CLIPscores_density_subset_plot
ggsave(file.path(Outdirectory, "Mean_CLIPscores_density_subset_plot.pdf"), Mean_CLIPscores_density_subset_plot, width = 9, height = 6, units = "in" )

################################################################################################################################
#Limma on CLIP scores- 5 vs 0 min

CLIPscores_forLimma <- CLIPscores_merged %>%
  dplyr::select(gene_name, ends_with("_CLIPscore")) %>%
  dplyr::select(gene_name,contains("5min")) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

expr_mat <- CLIPscores_forLimma

expr_mat[expr_mat < 0] <- 0

Limma_metadata <- data.frame(sample = colnames(expr_mat)) %>%
  mutate(
    Condition = ifelse(grepl("Opto", sample, ignore.case = TRUE), "Opto", "Control"),
    Timepoint = stringr::str_extract(sample, "[0-9]+min"),
    Group     = paste(Condition, Timepoint, sep = "_")
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("Control","Opto")),
    Timepoint = factor(Timepoint, levels = c("5min")),
    Group     = factor(Group, levels = c("Control_5min","Opto_5min"))  )

rownames(Limma_metadata) <- Limma_metadata$sample

design <- stats::model.matrix(~0 + Group , data= Limma_metadata)
colnames(design) <- make.names(colnames(design))

y <- edgeR::DGEList(expr_mat)
y <- edgeR::calcNormFactors(y)
v <- limma::voom(y, design, plot = TRUE)

fit <- limma::lmFit(v, design)

cont <- limma::makeContrasts(
  Opto_5min_vs_Control_5min  = `GroupOpto_5min`  - `GroupControl_5min`,
  levels = design
)

fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

OptoCLIPscores_limma_voom_5v0_results <- limma::topTable(fit2,coef="Opto_5min_vs_Control_5min",n=Inf, sort.by="P") %>% rownames_to_column("gene_name")

nrow(OptoCLIPscores_limma_voom_5v0_results %>% dplyr::filter(P.Value < 0.05)) #449

OptoCLIPscores_limma_voom_5v0_results <- OptoCLIPscores_limma_voom_5v0_results %>%
  mutate(CLIP_5v0_Status = case_when(
    P.Value < 0.05 & logFC > 0 ~ "CLIP_Up",
    P.Value < 0.05 & logFC < 0 ~ "CLIP_Down",
    TRUE ~ "NS"
  ))  %>%
  dplyr::select(-c(AveExpr,t,B)) %>%  
  dplyr::rename(CLIP_5v0_log2FC = logFC,
                CLIP_5v0_pvalue = P.Value,
                CLIP_5v0_adjpval = adj.P.Val) 

#Limma on CLIP scores, 30 min vs 5 min
CLIPscores_forLimma <- CLIPscores_merged %>%
  dplyr::select(gene_name, ends_with("_CLIPscore")) %>%
  column_to_rownames("gene_name") 

Limma_metadata <- data.frame(sample = colnames(CLIPscores_forLimma)) %>%
  mutate(
    Condition = ifelse(grepl("Opto", sample, ignore.case = TRUE), "Opto", "Control"),
    Timepoint = stringr::str_extract(sample, "[0-9]+min"),
    Group     = paste(Condition, Timepoint, sep = "_")  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("Control","Opto")),
    Timepoint = factor(Timepoint, levels = c("5min","30min")),
    Group     = factor(Group, levels = c("Control_5min","Opto_5min","Control_30min","Opto_30min"))  )

rownames(Limma_metadata)=Limma_metadata[,1]

CLIPscores_forLimma=CLIPscores_forLimma %>% rownames_to_column("gene_name")

source(file.path("~/ruthasinger_github/OptoCLIP_April2026/R_functions/run_limma_CLIPscore.R"))

OptoCLIPscores_limma_30v5 <- run_limma_interaction_CLIPscore(
  input_df = CLIPscores_forLimma,
  id_col = "gene_name",
  meta_df = Limma_metadata,
  condition_col = "Condition",
  timepoint_col = "Timepoint",
  handling = "zeros",
  voom = TRUE,
  min_reps = 2
)

OptoCLIPscores_limma_30v5_results <- OptoCLIPscores_limma_30v5$per_term

OptoCLIPscores_limma_30v5_results <- OptoCLIPscores_limma_30v5_results %>%
  mutate(CLIP_30v5_Status = case_when(
    P.Value < 0.05 & logFC > 0 ~ "CLIP_Up",
    P.Value < 0.05 & logFC < 0 ~ "CLIP_Down",
    TRUE ~ "NS"
  )) %>%
  dplyr::select(-c(AveExpr,t,B,contrast)) %>% 
  dplyr::rename(CLIP_30v5_log2FC = logFC,
                CLIP_30v5_pvalue = P.Value,
                CLIP_30v5_adjpval = adj.P.Val) 

OptoCLIP_limma_all <- OptoCLIPscores_limma_voom_5v0_results %>%
  full_join(OptoCLIPscores_limma_30v5_results, by = "gene_name")

All_CLIP_RT_data_joined <- CLIPscores_merged %>%
  left_join(OptoCLIP_limma_all, by = "gene_name")

################################################################################################################################
#Limma visualization and GO 

make_limma_results_df <- function(df, lfc_col, p_col) {
  df %>%
    dplyr::select(gene_name, dplyr::all_of(c(lfc_col, p_col))) %>%
    dplyr::transmute(
      gene_name = gene_name,
      log2FoldChange = .data[[lfc_col]],
      pvalue = .data[[p_col]]
    ) %>%
    dplyr::filter(!is.na(log2FoldChange), !is.na(pvalue))
}

CLIP_specs <- list(
  "5v0" = list(
    label = "5v0",
    results = make_limma_results_df(All_CLIP_RT_data_joined, "CLIP_5v0_log2FC", "CLIP_5v0_pvalue"),
    prefix = "CLIP",
    volcano_xlim = c(-3, 3),
    prio_up   = "nuclear_first",
    prio_down = "synaptic_first",
    go_subset_up_terms = c(
      "nuclear export",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
    ),
    go_subset_down_terms = c(
      "calcium ion transport",
      "regulation of synaptic plasticity",
      "receptor internalization",
      "postsynaptic neurotransmitter receptor internalization",
      "vesicle-mediated transport in synapse"
    ),
    volcano_label_up_terms = c(
      "nuclear export",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
    ),
    volcano_label_down_terms = c(
      "calcium ion transport",
      "regulation of synaptic plasticity",
      "receptor internalization",
      "postsynaptic neurotransmitter receptor internalization",
      "vesicle-mediated transport in synapse"
    ),
    volcano_label_top_n = 30
  ),
  "30v5" = list(
    label = "30v5",
    results = make_limma_results_df(All_CLIP_RT_data_joined, "CLIP_30v5_log2FC", "CLIP_30v5_pvalue"),
    prefix = "CLIP",
    volcano_xlim = c(-5, 5),
    prio_up   = "synaptic_first",
    prio_down = "nuclear_first",
    go_subset_up_terms = c(
      "calcium ion transport",
      "regulation of synaptic plasticity",
      "receptor internalization",
      "postsynaptic neurotransmitter receptor internalization",
      "vesicle-mediated transport in synapse"
    ),
    go_subset_down_terms = c(
      "nuclear export",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
    ),volcano_label_up_terms = c(
      "calcium ion transport",
      "regulation of synaptic plasticity",
      "receptor internalization",
      "postsynaptic neurotransmitter receptor internalization",
      "vesicle-mediated transport in synapse"
    ),
    volcano_label_down_terms = c(
      "nuclear export",
      "protein localization to nucleus",
      "RNA splicing",
      "chromatin remodeling",
      "miRNA transcription"
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

CLIPscore_out_list <- list()

for (comp in names(CLIP_specs)) {
  
  spec   <- CLIP_specs[[comp]]
  prefix <- spec$prefix
  
  message("Running CLIPscore block: ", comp)
  
  res <- spec$results %>%
    dplyr::rename(
      !!paste0(prefix, "_log2FC") := log2FoldChange,
      !!paste0(prefix, "_pvalue") := pvalue,
    ) 
  
  lfc_col  <- paste0(prefix, "_log2FC")
  pval_col <- paste0(prefix, "_pvalue")
  
  status_col <- paste0("CLIP_", comp, "_Status")
  
  res <- res %>%
    dplyr::mutate(
      !!status_col := dplyr::case_when(
        .data[[pval_col]] < 0.05 & .data[[lfc_col]] > 0 ~ paste0(prefix, "_Up"),
        .data[[pval_col]] < 0.05 & .data[[lfc_col]] < 0 ~ paste0(prefix, "_Down"),
        TRUE ~ "NS"
      )
    )
  
  up_df   <- res %>% dplyr::filter(.data[[lfc_col]] > 0, .data[[pval_col]] < 0.05)
  down_df <- res %>% dplyr::filter(.data[[lfc_col]] < 0, .data[[pval_col]] < 0.05)
  
  # Volcano (unlabeled)
  vol <- ggplot(res, aes(x = .data[[lfc_col]], y = -log10(.data[[pval_col]]), color = .data[[status_col]])) +
    geom_point(alpha = 1, size = 1) +
    scale_color_manual(
      values = setNames(
        c("red3","dodgerblue","gray70"),
        c(paste0(prefix,"_Up"), paste0(prefix,"_Down"), "NS")
      )
    ) +
    theme_classic(base_size = 9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray50") +
    labs(y = "-log10(p-value)", x = paste0("Log2 Fold Change CLIPscores\n", comp)) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Decreased (", nrow(down_df), ")"),
             color = "dodgerblue", hjust = 1.2, vjust = 2, size = 3) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Increased (", nrow(up_df), ")"),
             color = "red3", hjust = 1.2, vjust = 4, size = 3) +
    theme(legend.position = "none", axis.text = element_text(color = "black")) +
    xlim(spec$volcano_xlim[1], spec$volcano_xlim[2])
  print(vol)
  
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
  
  #Nuclear/synaptic classification + stacked bar
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
  
  # Merge with CH + violin
  CH_merge <- merge(res, CH_RT_resdata, by.x = "gene_name", by.y = "gene_name")
  
  CLIPscore_out_list[[comp]] <- list(
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

################################################################################################################################
#make bubble plot
nuclear_terms <- c(
  "nuclear export",
  "protein localization to nucleus",
  "RNA splicing",
  "miRNA transcription",
  "chromatin remodeling"
)

synaptic_terms <- c(
  "calcium ion transport",
  "receptor internalization",
  "postsynaptic neurotransmitter receptor internalization",
  "vesicle-mediated transport in synapse",
  "regulation of synaptic plasticity"
)

target_terms <- c(nuclear_terms, synaptic_terms)

allowed_terms_by_category <- list(
  "5v0_Down" = synaptic_terms,
  "5v0_Up"   = nuclear_terms,
  "30v5_Down"  = nuclear_terms,
  "30v5_Up"    = synaptic_terms
)

make_go_count_bubble_df <- function(CLIPscore_out_list, target_terms) {
  
  go_tables <- list(
    "5v0_Down" = CLIPscore_out_list[["5v0"]]$go_down,
    "5v0_Up"   = CLIPscore_out_list[["5v0"]]$go_up,
    "30v5_Down"  = CLIPscore_out_list[["30v5"]]$go_down,
    "30v5_Up"    = CLIPscore_out_list[["30v5"]]$go_up
  )
  
  out <- lapply(names(go_tables), function(cat_name) {
    
    go_df <- go_tables[[cat_name]]
    allowed_terms <- allowed_terms_by_category[[cat_name]]
    
    if (is.null(go_df) || nrow(go_df) == 0) {
      return(
        tibble::tibble(
          Description = allowed_terms,
          Count = 0,
          pvalue = NA_real_,
          category = cat_name
        )
      )
    }
    
    matched <- go_df %>%
      dplyr::filter(Description %in% allowed_terms) %>%
      dplyr::select(Description, Count, pvalue) %>%
      dplyr::mutate(category = cat_name)
    
    missing_terms <- setdiff(allowed_terms, matched$Description)
    
    if (length(missing_terms) > 0) {
      matched <- dplyr::bind_rows(
        matched,
        tibble::tibble(
          Description = missing_terms,
          Count = 0,
          pvalue = NA_real_,
          category = cat_name
        )
      )
    }
    
    matched
  })
  
  dplyr::bind_rows(out)
}

CLIPscore_go_bubble_df <- make_go_count_bubble_df(
  CLIPscore_out_list = CLIPscore_out_list,
  target_terms = target_terms
)

x_order <- c("5v0_Down", "5v0_Up", "30v5_Down", "30v5_Up")

term_order <- c(nuclear_terms, synaptic_terms)

CLIPscore_go_bubble_df <- CLIPscore_go_bubble_df %>%
  dplyr::mutate(
    biology = dplyr::case_when(
      Description %in% nuclear_terms  ~ "nuclear",
      Description %in% synaptic_terms ~ "synaptic",
      TRUE ~ "other"
    ),
    neglog10_pval_raw = -log10(pvalue),
    neglog10_pval = pmin(neglog10_pval_raw, 12),
    category = factor(category, levels = x_order),
    Description = factor(Description, levels = rev(term_order))
  )

CLIPscore_go_bubble_df <- CLIPscore_go_bubble_df %>%
  dplyr::filter(
    !is.na(pvalue),
    pvalue < 0.05,
    Count > 0
  )

CLIPscore_p_go_bubble <- ggplot(
  CLIPscore_go_bubble_df,
  aes(
    x = category,
    y = Description,
    size = Count,
    fill = neglog10_pval
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.25
  ) +
  scale_size_continuous(
    name = "Gene count",
    breaks = c(5, 10, 20, 40)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "pvalue",
    limits = c(1.3, 12),
    breaks = c(2, 4, 6, 8, 10, 12),
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
CLIPscore_p_go_bubble
ggsave(file.path(Outdirectory, "CLIPscore_GO_terms_bubbleplot.pdf"),CLIPscore_p_go_bubble,
       width = 5.5,
       height = 5.5,
       units = "in",
       device = grDevices::cairo_pdf
)

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
    coord_cartesian(ylim = c(0, 10))
  
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

for (comp in names(CLIPscore_out_list)) {
  
  spec <- CLIP_specs[[comp]]
  
  res        <- CLIPscore_out_list[[comp]]$res
  go_up      <- CLIPscore_out_list[[comp]]$go_up
  go_down    <- CLIPscore_out_list[[comp]]$go_down
  up_df      <- CLIPscore_out_list[[comp]]$up_df
  down_df    <- CLIPscore_out_list[[comp]]$down_df
  status_col <- CLIPscore_out_list[[comp]]$status_col
  
  prefix     <- spec$prefix
  lfc_col    <- paste0(prefix, "_log2FC")
  pvalue_col <- paste0(prefix, "_pvalue")
  
  up_terms   <- spec$volcano_label_up_terms
  down_terms <- spec$volcano_label_down_terms
  top_n      <- spec$volcano_label_top_n
  
  up_label_genes <- get_genes_from_go_terms(
    go_df = go_up,
    terms = up_terms,
    res_df = res,
    lfc_col = lfc_col,
    pvalue_col = pvalue_col,
    geneid_col = geneid_col,
    desc_col = desc_col,
    gene_sep = gene_sep,
    top_n = top_n,
    rank_by = "pvalue"
  )
  
  down_label_genes <- get_genes_from_go_terms(
    go_df = go_down,
    terms = down_terms,
    res_df = res,
    lfc_col = lfc_col,
    pvalue_col = pvalue_col,
    geneid_col = geneid_col,
    desc_col = desc_col,
    gene_sep = gene_sep,
    top_n = top_n,
    rank_by = "pvalue"
  )
  
  message(comp, ": labeling up genes from GO terms: ", paste(up_terms, collapse = "; "))
  message(comp, ": found ", length(up_label_genes), " uniquely-mapped up genes to label")
  
  message(comp, ": labeling down genes from GO terms: ", paste(down_terms, collapse = "; "))
  message(comp, ": found ", length(down_label_genes), " uniquely-mapped down genes to label")
  
  vol_labeled <- make_labeled_volcano(
    res = res,
    lfc_col = lfc_col,
    pvalue_col = pvalue_col,
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
      paste0("CLIPscore_VolcanoPlot_", comp, "_labeled.pdf")
    ),
    vol_labeled,
    width = 4.5,
    height = 4,
    units = "in"
  )
}

################################################################################################################################
#comparison of Up and Down genes from each of the 3 comparisons 
comp_order <- c("5v0", "30v5")

get_gene_vec <- function(CLIPscore_out_list, comp, which = c("up", "down")) {
  which <- match.arg(which)
  
  if (!comp %in% names(CLIPscore_out_list)) return(character(0))
  
  x <- if (which == "up") CLIPscore_out_list[[comp]]$up_df else CLIPscore_out_list[[comp]]$down_df
  
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
  venn_sets_all[[paste0(comp, "_Up")]]   <- get_gene_vec(CLIPscore_out_list, comp, "up")
  venn_sets_all[[paste0(comp, "_Down")]] <- get_gene_vec(CLIPscore_out_list, comp, "down")
}

# Drop empty sets (optional but prevents failures if a comparison has 0 sig genes)
venn_sets_all <- venn_sets_all[lengths(venn_sets_all) > 0]

if (length(venn_sets_all) < 2) {
  stop("Not enough non-empty gene sets for venndetail(). Check p<0.05 thresholds / comps.")
}

vd_all <- venndetail(venn_sets_all)

pdf(file.path(Outdirectory, "VennDetail_CLIPscore.pdf"), width = 10, height = 6)
plot(vd_all, type = "upset")
dev.off()

plot(vd_all, type = "upset")

vd_res <- vd_all@result

subset_counts <- vd_res %>%
  dplyr::group_by(Subset) %>%
  dplyr::summarise(gene_count = dplyr::n(), .groups = "drop")

subsets_keep_for_go <- subset_counts %>%
  dplyr::filter(gene_count >= 10) %>%
  dplyr::pull(Subset)

vd_res_go <- vd_res %>%
  dplyr::filter(Subset %in% subsets_keep_for_go)

#looking at syn/nuclear breakdown for all venn categories
analyze_overlap_nucsyn <- function(
    gene_vec,
    label,
    priority = c("synaptic_first", "nuclear_first"),
    Outdirectory = ".",
    OrgDb = org.Mm.eg.db,
    desc_col = "Description",
    geneid_col = "geneID",
    gene_sep = "/",
    nuclear_keywords,
    exclude_nuclear_keywords = character(0),
    synaptic_keywords,
    exclude_synaptic_keywords = character(0)
) {
  
  priority <- match.arg(priority)
  
  gene_vec <- unique(as.character(gene_vec))
  gene_vec <- stringr::str_trim(gene_vec)
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]
  
  if (length(gene_vec) == 0) {
    message("No genes found for overlap set: ", label)
    return(NULL)
  }
  
  message(label, ": running GO enrichment on ", length(gene_vec), " genes")
  
  go_obj <- clusterProfiler::enrichGO(
    gene = gene_vec,
    OrgDb = OrgDb,
    keyType = "SYMBOL",
    readable = TRUE,
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "none"
  )
  
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) {
    message("No significant GO terms for: ", label)
    
    empty_go <- tibble::tibble()
    
    empty_bar <- tibble::tibble(
      bar = label,
      category = c("synaptic_terms", "nuclear_terms", "other"),
      n_genes = c(0, 0, length(gene_vec)),
      total_genes = length(gene_vec),
      pct = c(0, 0, 100)
    )
    
    return(list(
      go_obj = NULL,
      go_df = empty_go,
      bar_df = empty_bar
    ))
  }
  
  go_df <- go_obj@result 
  
  out <- run_one_direction(
    go_df = go_df,
    genes_df = tibble::tibble(gene_name = gene_vec),
    direction_label = label,
    priority = priority,
    comp_label = NULL,
    Outdirectory = Outdirectory,
    desc_col = desc_col,
    geneid_col = geneid_col,
    gene_sep = gene_sep,
    nuclear_keywords = nuclear_keywords,
    exclude_nuclear_keywords = exclude_nuclear_keywords,
    synaptic_keywords = synaptic_keywords,
    exclude_synaptic_keywords = exclude_synaptic_keywords
  )
  
  out
}

# helper to extract genes from a venn subset
extract_subset_genes <- function(vd_result_df, subset_name) {
  x <- vd_result_df %>%
    dplyr::filter(Subset == subset_name) %>%
    dplyr::pull(Detail)
  
  x <- unique(as.character(x))
  x <- stringr::str_trim(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

# get all venn subset names
all_subsets <- vd_res_go %>%
  dplyr::pull(Subset) %>%
  unique()

# optional: keep only subsets with at least 1 gene
all_subsets <- all_subsets[!is.na(all_subsets) & all_subsets != ""]

# run analysis on every subset
all_overlap_out <- list()

for (ss in all_subsets) {
  
  gene_vec <- extract_subset_genes(vd_all@result, ss)
  
  if (length(gene_vec) == 0) {
    message("Skipping empty subset: ", ss)
    next
  }
  
  # choose priority automatically based on label
  # edit this logic if you want something different
  priority_map <- c(
    "5v0_Down_30v5_Down" = "synaptic_first",
    "5v0_Up_30v5_Down"   = "nuclear_first",
    "30v5_Down"           = "nuclear_first",
    "5v0_Down_30v5_Up"   = "synaptic_first",
    "30v5_Up"             = "synaptic_first",
    "5v0_Down"          = "synaptic_first",
    "5v0_Up"            = "nuclear_first"
  )
  
  this_priority <- priority_map[ss]
  
  # fallback if something unexpected appears
  if (is.na(this_priority)) {
    message("No priority defined for ", ss, " — defaulting to nuclear_first")
    this_priority <- "nuclear_first"
  }
  
  all_overlap_out[[ss]] <- analyze_overlap_nucsyn(
    gene_vec = gene_vec,
    label = ss,
    priority = this_priority,
    Outdirectory = file.path(Outdirectory),
    nuclear_keywords = nuclear_keywords,
    exclude_nuclear_keywords = exclude_nuclear_keywords,
    synaptic_keywords = synaptic_keywords,
    exclude_synaptic_keywords = exclude_synaptic_keywords
  )
}

# keep only successful outputs
all_overlap_out <- all_overlap_out[!vapply(all_overlap_out, is.null, logical(1))]

if (length(all_overlap_out) == 0) {
  stop("No overlap subsets produced usable output.")
}

# combine stack dfs from all subsets
overlap_stack_df <- dplyr::bind_rows(
  lapply(names(all_overlap_out), function(ss) {
    all_overlap_out[[ss]]$bar_df %>%
      dplyr::mutate(
        overlap = ss,
        x = "overlap"
      )
  })
) %>%
  dplyr::mutate(
    category = factor(category, levels = c("synaptic_terms", "nuclear_terms", "other")),
    overlap = factor(overlap, levels = all_subsets)
  )

overlap_stack_df$overlap <- factor(
  overlap_stack_df$overlap,
  levels = c(
    "30v5_Down",
    "5v0_Down",
    "5v0_Up",
    "30v5_Up",
    "5v0_Up_30v5_Down",
    "5v0_Down_30v5_Up",
    "5v0_Down_30v5_Down"
  )
)

print(overlap_stack_df)

#category       n_genes bar              total_genes   pct overlap              x      
#nuclear_terms      197 30v5_Down                 543  36.3 30v5_Down           overlap
#other              250 30v5_Down                 543  46.0 30v5_Down           overlap
#synaptic_terms      96 30v5_Down                 543  17.7 30v5_Down           overlap
#nuclear_terms       35 5v0_Down                  161  21.7 5v0_Down            overlap
#other               71 5v0_Down                  161  44.1 5v0_Down            overlap
#synaptic_terms      55 5v0_Down                  161  34.2 5v0_Down            overlap
#nuclear_terms       32 5v0_Up                    108  29.6 5v0_Up              overlap
#other               56 5v0_Up                    108  51.9 5v0_Up              overlap
#synaptic_terms      20 5v0_Up                    108  18.5 5v0_Up              overlap
#nuclear_terms        9 30v5_Up                    52  17.3 30v5_Up             overlap
#other               21 30v5_Up                    52  40.4 30v5_Up             overlap
#synaptic_terms      22 30v5_Up                    52  42.3 30v5_Up             overlap
#nuclear_terms       53 5v0_Up_30v5_Down          144  36.8 5v0_Up_30v5_Down    overlap
#other               66 5v0_Up_30v5_Down          144  45.8 5v0_Up_30v5_Down    overlap
#synaptic_terms      25 5v0_Up_30v5_Down          144  17.4 5v0_Up_30v5_Down    overlap
#nuclear_terms        5 5v0_Down_30v5_Up          31  16.1 5v0_Down_30v5_Up     overlap
#other               10 5v0_Down_30v5_Up          31  32.3 5v0_Down_30v5_Up     overlap
#synaptic_terms      16 5v0_Down_30v5_Up          31  51.6 5v0_Down_30v5_Up     overlap

nuc_syn_fill_colors <- c(
  "synaptic_terms" = "#1B9E77",
  "nuclear_terms"  = "#D95F02",
  "other"          = "#7570B3"
)

# plot all subsets
p_overlap <- ggplot(
  overlap_stack_df,
  aes(x = x, y = pct, fill = category)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = nuc_syn_fill_colors, drop = FALSE) +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = "% of genes", fill = NULL) +
  facet_wrap(~overlap, nrow=1) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 7)
  )
print(p_overlap)
ggsave(file.path(Outdirectory, "CLIPscore_VennOverlap_NucSyn_stackedbar.pdf"),p_overlap,width = 10, height = 8, units = "in")

################################################################################################################################

tp_tbl <- tibble::tribble(
  ~comp,  ~clip_lfc,            ~clip_p,            ~rt_lfc,              ~rt_p,
  "5v0", "CLIP_5v0_log2FC",   "CLIP_5v0_pvalue", "RiboTag_5v0_log2FC",  "RiboTag_5v0_pvalue",
  "30v5",  "CLIP_30v5_log2FC",    "CLIP_30v5_pvalue",  "RiboTag_30v5_log2FC",   "RiboTag_30v5_pvalue"
)

get_tp <- function(comp) {
  x <- tp_tbl[tp_tbl$comp == comp, ]
  stopifnot(nrow(x) == 1)
  x
}

clip_dir_colors <- c(
  "CLIP_Up"   = "red3",
  "CLIP_Down" = "dodgerblue3",
  "NS"        = "grey70"
)

p_to_stars <- function(p) {
  if (is.na(p)) "NA"
  else if (p < 1e-4) "****"
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 0.05) "*"
  else "ns"
}

#Per-comp ridges + Wilcoxon
run_ridges_and_wilcox_onecomp <- function(df, comp, Outdirectory,
                                          xlim = c(-1, 1),
                                          clip_alpha = 0.05,
                                          rt_alpha = NULL) {
  # rt_alpha:
  #  - NULL means "do not filter RT by significance" (plot all; stats on all)
  #  - numeric means "filter to RT_sig_tmp == TRUE" before plotting + stats
  
  m <- get_tp(comp)
  clip_lfc <- m$clip_lfc[[1]]
  clip_p   <- m$clip_p[[1]]
  rt_lfc   <- m$rt_lfc[[1]]
  rt_p     <- m$rt_p[[1]]
  
  # sanity checks
  needed <- c(clip_lfc, clip_p, rt_lfc)
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) stop("Missing columns in df for comp ", comp, ": ", paste(missing, collapse=", "))
  
  # RT p col optional
  has_rt_p <- rt_p %in% colnames(df)
  
  df2 <- df %>%
    mutate(
      CLIP_sig_tmp = !is.na(.data[[clip_p]]) & .data[[clip_p]] < clip_alpha,
      CLIP_direction_tmp = case_when(
        !CLIP_sig_tmp ~ "NS",
        .data[[clip_lfc]] > 0 ~ "CLIP_Up",
        .data[[clip_lfc]] < 0 ~ "CLIP_Down",
        TRUE ~ "NS"
      ),
      CLIP_direction_tmp = factor(CLIP_direction_tmp, levels = c("CLIP_Down","CLIP_Up","NS")),
      RT_sig_tmp = if (has_rt_p) (!is.na(.data[[rt_p]]) & .data[[rt_p]] < (rt_alpha %||% Inf)) else NA
    )
  
  # optional filter to RT significant only
  if (!is.null(rt_alpha) && has_rt_p) {
    df2 <- df2 %>% filter(!is.na(.data[[rt_p]]), .data[[rt_p]] < rt_alpha)
  }
  
  # plot
  p <- ggplot(df2, aes(x = .data[[rt_lfc]], y = CLIP_direction_tmp, fill = CLIP_direction_tmp)) +
    ggridges::geom_density_ridges(alpha = 0.6) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 9) +
    xlim(xlim) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.02))) +
    scale_fill_manual(values = clip_dir_colors, drop = FALSE) +
    labs(
      x = paste0("RiboTag log2FC (", comp, ")"),
      y = NULL,
      title = paste0(comp, if (!is.null(rt_alpha) && has_rt_p) paste0(" (RT p<", rt_alpha, ")") else "")
    ) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = .3)
  print(p)
  ggsave(file.path(Outdirectory, paste0("CLIPscore_RiboTag_densityridge_", comp, if (!is.null(rt_alpha) && has_rt_p) paste0("_RTp", rt_alpha) else "",".pdf")),plot = p, width = 6, height = 4)
  
  # stats helper that won't error if group is empty
  safe_wilcox <- function(x, g1, g2, alternative) {
    x1 <- x[g1]; x2 <- x[g2]
    x1 <- x1[is.finite(x1)]; x2 <- x2[is.finite(x2)]
    if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
    wilcox.test(x1, x2, alternative = alternative)$p.value
  }
  
  x <- df2[[rt_lfc]]
  g_down <- df2$CLIP_direction_tmp == "CLIP_Down"
  g_up   <- df2$CLIP_direction_tmp == "CLIP_Up"
  g_ns   <- df2$CLIP_direction_tmp == "NS"
  
  p_down_vs_ns <- safe_wilcox(x, g_down, g_ns, alternative = "greater")
  p_up_vs_ns   <- safe_wilcox(x, g_up,   g_ns, alternative = "less")
  p_down_vs_up <- safe_wilcox(x, g_down, g_up, alternative = "greater")
  
  stats_tbl <- tibble::tibble(
    comp = comp,
    rt_filter = if (!is.null(rt_alpha) && has_rt_p) paste0("<", rt_alpha) else "none",
    test = c("Down > NS", "Up < NS", "Down > Up"),
    p_value = c(p_down_vs_ns, p_up_vs_ns, p_down_vs_up),
    stars  = sapply(c(p_down_vs_ns, p_up_vs_ns, p_down_vs_up), p_to_stars),
    n_down = sum(g_down, na.rm = TRUE),
    n_up   = sum(g_up,   na.rm = TRUE),
    n_ns   = sum(g_ns,   na.rm = TRUE)
  )
  list(plot = p, stats = stats_tbl, df_used = df2)
}

#Run all comps
ridges_out <- list()

for (comp in tp_tbl$comp) {  
  ridges_out[[comp]] <- run_ridges_and_wilcox_onecomp(
    df = All_CLIP_RT_data_joined,
    comp = comp,
    Outdirectory = Outdirectory,
    xlim = c(-1, 1),
    clip_alpha = 0.05,
    rt_alpha = NULL      # set to 0.05 if you want to restrict to RT sig only
  )
}

all_stats <- dplyr::bind_rows(lapply(ridges_out, `[[`, "stats"))

print(all_stats)

#comp   rt_filter test       p_value stars   n_down  n_up  n_ns
#5v0    none      Down > NS 7.03e- 2 ns       197   252 12128
#5v0    none      Up < NS   3.39e- 7 ****     197   252 12128
#5v0    none      Down > Up 4.22e- 7 ****     197   252 12128
#30v5   none      Down > NS 2.53e-20 ****     692    83 11802
#30v5   none      Up < NS   5.78e- 3 **       692    83 11802
#30v5   none      Down > Up 8.96e-10 ****     692    83 11802

################################################################################################################################
run_clip_rt_overlap_onecomp <- function(df, comp, Outdirectory) {
  
  m <- get_tp(comp)
  clip_lfc <- m$clip_lfc[[1]]
  clip_p   <- m$clip_p[[1]]
  rt_lfc   <- m$rt_lfc[[1]]
  rt_p     <- m$rt_p[[1]]
  
  # universe
  universe <- df %>%
    dplyr::filter(
      !is.na(.data[[clip_p]]), !is.na(.data[[clip_lfc]]),
      !is.na(.data[[rt_p]]),   !is.na(.data[[rt_lfc]])
    ) %>%
    dplyr::pull(gene_name) %>% unique()
  
  # define sets
  sets <- list(
    CLIP_Up = df %>%
      dplyr::filter(.data[[clip_p]] < 0.05, .data[[clip_lfc]] > 0) %>%
      dplyr::pull(gene_name) %>% unique(),
    
    CLIP_Down = df %>%
      dplyr::filter(.data[[clip_p]] < 0.05, .data[[clip_lfc]] < 0) %>%
      dplyr::pull(gene_name) %>% unique(),
    
    RT_pval_Up = df %>%
      dplyr::filter(.data[[rt_p]] < 0.05, .data[[rt_lfc]] > 0) %>%
      dplyr::pull(gene_name) %>% unique(),
    
    RT_pval_Down = df %>%
      dplyr::filter(.data[[rt_p]] < 0.05, .data[[rt_lfc]] < 0) %>%
      dplyr::pull(gene_name) %>% unique()
  )
  
  RT_sig <- union(sets$RT_pval_Up, sets$RT_pval_Down)
  
  CLIP_order <- c("CLIP_Down","CLIP_Up")
  CLIP_RT_order <- c("CLIP_Down_RT_Up","CLIP_Down_RT_Down","CLIP_Up_RT_Down","CLIP_Up_RT_Up")
  
  overlap_df <- dplyr::bind_rows(
    tibble::tibble(CLIP="CLIP_Up",   gene=sets$CLIP_Up),
    tibble::tibble(CLIP="CLIP_Down", gene=sets$CLIP_Down)
  ) %>%
    dplyr::filter(gene %in% RT_sig) %>%
    dplyr::mutate(
      RT_status_pval = dplyr::case_when(
        gene %in% sets$RT_pval_Up   ~ "RT_Up",
        gene %in% sets$RT_pval_Down ~ "RT_Down"
      ),
      RT_status_pval = factor(RT_status_pval, levels = c("RT_Down","RT_Up")),
      CLIP = factor(CLIP, levels = CLIP_order),
      CLIP_RT = factor(paste(CLIP, RT_status_pval, sep="_"), levels = CLIP_RT_order)
    ) %>%
    dplyr::count(CLIP, RT_status_pval, CLIP_RT) %>%
    dplyr::group_by(CLIP) %>%
    dplyr::mutate(fraction = n / sum(n)) %>%
    dplyr::ungroup()
  
  CLIP_RT_cat_colors <- c(
    "CLIP_Up_RT_Down"   = "#223F99",
    "CLIP_Up_RT_Up"     = "green3",
    "CLIP_Down_RT_Down" = "gold",
    "CLIP_Down_RT_Up"   = "violet"
  )
  
  p <- ggplot2::ggplot(overlap_df, ggplot2::aes(x = CLIP, y = fraction, fill = CLIP_RT)) +
    ggplot2::geom_col(width = 0.7, color = "black", size = 0.2, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(values = CLIP_RT_cat_colors) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      x = NULL,
      y = "% of differentially bound\nFMRP targets by RiboTag category",
      fill = NULL,
      title = comp
    )
  print(p)
  
  # hypergeom pvals (same as you do)
  overlap_pval <- function(setA, setB, universe) {
    A <- intersect(unique(setA), universe)
    B <- intersect(unique(setB), universe)
    U <- unique(universe)
    
    k <- length(intersect(A, B))
    m <- length(A)
    n <- length(U) - m
    K <- length(B)
    
    stats::phyper(q = k - 1, m = m, n = n, k = K, lower.tail = FALSE)
  }
  
  pairwise_pvals <- expand.grid(
    A = names(sets),
    B = names(sets),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(A < B) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pval = overlap_pval(sets[[A]], sets[[B]], universe),
      overlap_n = length(intersect(intersect(sets[[A]], universe), intersect(sets[[B]], universe))),
      size_A = length(intersect(sets[[A]], universe)),
      size_B = length(intersect(sets[[B]], universe)),
      universe_n = length(universe)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(padj = p.adjust(pval, method = "BH"))
  print(pairwise_pvals)
  
  list(sets = sets, overlap_df = overlap_df, plot = p, pairwise_pvals = pairwise_pvals)
}

CLIP_RT_cat_colors <- c(
  "CLIP_Up_RT_Down"   = "#223F99",
  "CLIP_Up_RT_Up"     = "green3",
  "CLIP_Down_RT_Down" = "gold",
  "CLIP_Down_RT_Up"   = "violet"
)

clip_rt_out <- list()
clip_rt_overlap_tables <- list()

for (comp in c("5v0","30v5")) {
  clip_rt_out[[comp]] <- run_clip_rt_overlap_onecomp(All_CLIP_RT_data_joined, comp, Outdirectory)
  clip_rt_overlap_tables[[comp]] <- clip_rt_out[[comp]]$overlap_df %>%
    dplyr::mutate(comp = comp)
}

if (length(clip_rt_overlap_tables) > 0) {
  
  CLIP_RT_overlap_all_df <- dplyr::bind_rows(clip_rt_overlap_tables) %>%
    dplyr::mutate(
      comp = factor(comp, levels = c("5v0", "30v5")),
      CLIP = factor(CLIP, levels = c("CLIP_Down", "CLIP_Up")),
      CLIP_RT = factor(CLIP_RT, levels = c(
        "CLIP_Down_RT_Up",
        "CLIP_Down_RT_Down",
        "CLIP_Up_RT_Down",
        "CLIP_Up_RT_Up"
      ))
    )
  
  CLIP_RT_overlap_stackedbar_faceted <- ggplot(CLIP_RT_overlap_all_df,
                                               aes(x = CLIP, y = fraction, fill = CLIP_RT)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2,
             position = position_stack(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = CLIP_RT_cat_colors, drop = FALSE) +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = "% of differentially bound\nFMRP targets by RiboTag category",
      fill = NULL
    ) +
    facet_wrap(~comp, nrow = 1)
  
  print(CLIP_RT_overlap_stackedbar_faceted)
  
  ggsave(
    file.path(Outdirectory, "CLIPscore_RiboTag_overlap_stackedbar_faceted.pdf"),
    plot = CLIP_RT_overlap_stackedbar_faceted,
    width = 9, height = 4
  )
}

print(CLIP_RT_overlap_all_df)

#CLIP      RT_status_pval CLIP_RT               n fraction comp 
#CLIP_Down RT_Down        CLIP_Down_RT_Down    20    0.385 5v0 
#CLIP_Down RT_Up          CLIP_Down_RT_Up      32    0.615 5v0 
#CLIP_Up   RT_Down        CLIP_Up_RT_Down      48    0.738 5v0 
#CLIP_Up   RT_Up          CLIP_Up_RT_Up        17    0.262 5v0 
#CLIP_Down RT_Down        CLIP_Down_RT_Down    34    0.173 30v5  
#CLIP_Down RT_Up          CLIP_Down_RT_Up     163    0.827 30v5  
#CLIP_Up   RT_Down        CLIP_Up_RT_Down      17    0.773 30v5  
#CLIP_Up   RT_Up          CLIP_Up_RT_Up         5    0.227 30v5  

clip_rt_out$'5v0'$pairwise_pvals

#A            B                 pval overlap_n size_A size_B universe_n      padj
#CLIP_Down    CLIP_Up      1                 0    197    252      12577 1        
#CLIP_Up      RT_pval_Up   0.980            17    252   1288      12577 1        
#CLIP_Down    RT_pval_Up   0.00567          32    197   1288      12577 0.0170   
#RT_pval_Down RT_pval_Up   1                 0   1289   1288      12577 1        
#CLIP_Up      RT_pval_Down 0.0000160        48    252   1289      12577 0.0000962
#CLIP_Down    RT_pval_Down 0.553            20    197   1289      12577 1     

clip_rt_out$'30v5'$pairwise_pvals

#A            B                 pval overlap_n size_A size_B universe_n     padj
#CLIP_Down    CLIP_Up      1    e+ 0         0    692     83      12577 1   e+ 0
#CLIP_Up      RT_pval_Up   9.75 e- 1         5     83   1489      12577 1   e+ 0
#CLIP_Down    RT_pval_Up   4.90 e-19       163    692   1489      12577 2.94e-18
#RT_pval_Down RT_pval_Up   1    e+ 0         0   1376   1489      12577 1   e+ 0
#CLIP_Up      RT_pval_Down 7.68 e- 3        17     83   1376      12577 2.30e- 2
#CLIP_Down    RT_pval_Down 1.000e+ 0        34    692   1376      12577 1   e+ 0

################################################################################################################################

tp_tbl <- tibble::tribble(
  ~comp,  ~clip_lfc,            ~clip_p,            ~rt_lfc,              ~rt_p,
  "5v0", "CLIP_5v0_log2FC",   "CLIP_5v0_pvalue", "RiboTag_5v0_log2FC",  "RiboTag_5v0_pvalue",
  "30v5",  "CLIP_30v5_log2FC",    "CLIP_30v5_pvalue",  "RiboTag_30v5_log2FC",   "RiboTag_30v5_pvalue"
)

get_set_from_cliprt <- function(clip_rt_out, comp, set_name) {
  if (!comp %in% names(clip_rt_out)) return(character(0))
  s <- clip_rt_out[[comp]]$sets[[set_name]]
  unique(as.character(s))
}

get_tp <- function(comp) {
  x <- tp_tbl[tp_tbl$comp == comp, ]
  stopifnot(nrow(x) == 1)
  x
}

CLIP_RT_cat_colors <- c(
  "CLIP_Up_RT_Down"   = "#223F99",
  "CLIP_Up_RT_Up"     = "green3",
  "CLIP_Down_RT_Down" = "gold",
  "CLIP_Down_RT_Up"   = "violet"
)

CLIP_RT_fill_colors <- c(
  "synaptic_terms" = "#1B9E77",
  "nuclear_terms"  = "#D95F02",
  "other"          = "#7570B3"
)

priority_lookup_cliprt <- tibble::tribble(
  ~comp,  ~quadrant,           ~priority,
  "5v0", "CLIP_Down_RT_Up",   "synaptic_first",
  "5v0", "CLIP_Up_RT_Down",   "nuclear_first",
  "30v5",  "CLIP_Down_RT_Up",   "nuclear_first",
  "30v5",  "CLIP_Up_RT_Down",   "synaptic_first"
)

nuclear_subset_terms <- c(
  "nucleus organization",
  "nuclear transport",
  "nuclear export",
  "histone modification",
  "RNA localization",
  "protein localization to nucleus",
  "RNA splicing",
  "DNA-templated transcription initiation"
)

synaptic_subset_terms <- c(
  "chemical synaptic transmission, postsynaptic",
  "synapse assembly",
  "receptor internalization",
  "vesicle-mediated transport in synapse",  
  "receptor-mediated endocytosis",  
  "calcium ion homeostasis",  
  "regulation of postsynaptic membrane potential",
  "ionotropic glutamate receptor signaling pathway",
  "regulation of NMDA receptor activity"
)

get_priority_cliprt <- function(comp, quadrant) {
  out <- priority_lookup_cliprt %>%
    dplyr::filter(comp == !!comp, quadrant == !!quadrant) %>%
    dplyr::pull(priority)
  if (length(out) != 1) stop("No unique priority for ", comp, " / ", quadrant)
  out
}

get_subset_terms_for_priority <- function(priority,
                                          nuclear_terms,
                                          synaptic_terms) {
  if (priority == "nuclear_first") return(nuclear_terms)
  if (priority == "synaptic_first") return(synaptic_terms)
  stop("Unknown priority: ", priority)
}

priority_label <- function(priority) {
  if (priority == "nuclear_first") return("Nuclear-prioritized")
  if (priority == "synaptic_first") return("Synaptic-prioritized")
  as.character(priority)
}

make_clip_rt_quadrant_sets <- function(df, comp, alpha = 0.05) {
  m <- get_tp(comp)
  clip_lfc <- m$clip_lfc[[1]]
  clip_p   <- m$clip_p[[1]]
  rt_lfc   <- m$rt_lfc[[1]]
  rt_p     <- m$rt_p[[1]]
  
  df_sig <- df %>%
    dplyr::filter(!is.na(.data[[clip_p]]), !is.na(.data[[rt_p]]),
                  !is.na(.data[[clip_lfc]]), !is.na(.data[[rt_lfc]])) %>%
    dplyr::filter(.data[[clip_p]] < alpha, .data[[rt_p]] < alpha)
  
  sets <- list(
    CLIP_Down_RT_Up = df_sig %>%
      dplyr::filter(.data[[clip_lfc]] < 0, .data[[rt_lfc]] > 0) %>%
      dplyr::pull(gene_name) %>% unique(),
    CLIP_Down_RT_Down = df_sig %>%
      dplyr::filter(.data[[clip_lfc]] < 0, .data[[rt_lfc]] < 0) %>%
      dplyr::pull(gene_name) %>% unique(),
    CLIP_Up_RT_Up = df_sig %>%
      dplyr::filter(.data[[clip_lfc]] > 0, .data[[rt_lfc]] > 0) %>%
      dplyr::pull(gene_name) %>% unique(),
    CLIP_Up_RT_Down = df_sig %>%
      dplyr::filter(.data[[clip_lfc]] > 0, .data[[rt_lfc]] < 0) %>%
      dplyr::pull(gene_name) %>% unique()
  )
  
  list(
    df_sig = df_sig,
    sets = sets,
    clip_lfc = clip_lfc,
    clip_p = clip_p,
    rt_lfc = rt_lfc,
    rt_p = rt_p
  )
}

safe_enrichGO <- function(genes, ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "none") {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) < 5) return(NULL)
  
  clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    readable = TRUE,
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod
  )
}

plot_clip_rt_scatter_unlabeled <- function(df, comp, xlim = c(-4,4), ylim = c(-2,2), alpha = 0.05) {
  q <- make_clip_rt_quadrant_sets(df, comp, alpha = alpha)
  clip_lfc <- q$clip_lfc
  rt_lfc   <- q$rt_lfc
  
  ggplot(df) +
    ylab(paste0("RiboTag log2FC\n", comp)) +
    xlab(paste0("FMRP CLIP log2FC\n", comp)) +
    xlim(xlim) +
    ylim(ylim) +
    theme_classic(base_size = 9) +
    geom_vline(xintercept = 0, linewidth = .3) +
    geom_hline(yintercept = 0, linewidth = .3) +
    theme(legend.position = "none") +
    geom_point(
      data = q$df_sig %>% filter(.data[[rt_lfc]] > 0, .data[[clip_lfc]] < 0),
      aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Down_RT_Up"), size = 2
    ) +
    geom_point(
      data = q$df_sig %>% filter(.data[[rt_lfc]] < 0, .data[[clip_lfc]] < 0),
      aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Down_RT_Down"), size = 2
    ) +
    geom_point(
      data = q$df_sig %>% filter(.data[[rt_lfc]] > 0, .data[[clip_lfc]] > 0),
      aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Up_RT_Up"), size = 2
    ) +
    geom_point(
      data = q$df_sig %>% filter(.data[[rt_lfc]] < 0, .data[[clip_lfc]] > 0),
      aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Up_RT_Down"), size = 2
    ) +
    scale_color_manual(values = CLIP_RT_cat_colors, drop = FALSE)
}

clip_rt_comp_out <- list()
cliprt_nucsyn_tables <- list()

for (comp in c("5v0","30v5")) {
  
  message("=== CLIP/RT quadrant GO + classification: ", comp, " ===")
  
  p_unlab <- plot_clip_rt_scatter_unlabeled(
    df = All_CLIP_RT_data_joined,
    comp = comp,
    xlim = c(-4,4),
    ylim = c(-2,2),
    alpha = 0.05
  )
  print(p_unlab)
  
  q <- make_clip_rt_quadrant_sets(All_CLIP_RT_data_joined, comp, alpha = 0.05)
  
  ego_down_up <- safe_enrichGO(q$sets$CLIP_Down_RT_Up, ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "none")
  ego_up_down <- safe_enrichGO(q$sets$CLIP_Up_RT_Down, ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "none")
  
  go_down_up <- if (is.null(ego_down_up)) tibble() else ego_down_up@result
  go_up_down <- if (is.null(ego_up_down)) tibble() else ego_up_down@result
  
  out_kw <- list()
  
  prio_down_up <- get_priority_cliprt(comp, "CLIP_Down_RT_Up")
  prio_up_down <- get_priority_cliprt(comp, "CLIP_Up_RT_Down")
  
  if (nrow(go_down_up) > 0) {
    out_kw$Down_RT_Up <- run_one_direction(
      go_df = go_down_up,
      genes_df = tibble::tibble(gene_name = q$sets$CLIP_Down_RT_Up),
      direction_label = paste0("CLIP_", comp, "_Down_RT_Up"),
      priority = prio_down_up,
      comp_label = comp,
      Outdirectory =   Outdirectory,
      nuclear_keywords = nuclear_keywords,
      exclude_nuclear_keywords = exclude_nuclear_keywords,
      synaptic_keywords = synaptic_keywords,
      exclude_synaptic_keywords = exclude_synaptic_keywords
    )
  } else {
    message("Skipping keyword classification: no GO results for ", comp, " CLIP_Down_RT_Up")
  }
  
  if (nrow(go_up_down) > 0) {
    out_kw$Up_RT_Down <- run_one_direction(
      go_df = go_up_down,
      genes_df = tibble::tibble(gene_name = q$sets$CLIP_Up_RT_Down),
      direction_label = paste0("CLIP_", comp, "_Up_RT_Down"),
      priority = prio_up_down,
      comp_label = comp,
      Outdirectory = Outdirectory,
      nuclear_keywords = nuclear_keywords,
      exclude_nuclear_keywords = exclude_nuclear_keywords,
      synaptic_keywords = synaptic_keywords,
      exclude_synaptic_keywords = exclude_synaptic_keywords
    )
  } else {
    message("Skipping keyword classification: no GO results for ", comp, " CLIP_Up_RT_Down")
  }
  
  p_stack <- NULL
  stack_df <- tibble()
  
  if (!is.null(out_kw$Down_RT_Up) && !is.null(out_kw$Up_RT_Down)) {
    
    stack_df <- dplyr::bind_rows(
      out_kw$Down_RT_Up$bar_df,
      out_kw$Up_RT_Down$bar_df
    ) %>%
      dplyr::mutate(
        category = factor(category, levels = c("synaptic_terms","nuclear_terms","other")),
        bar = factor(bar, levels = c(
          paste0(comp, "_CLIP_", comp, "_Down_RT_Up"),
          paste0(comp, "_CLIP_", comp, "_Up_RT_Down")
        ))
      )
    
    stack_df <- dplyr::bind_rows(
      out_kw$Down_RT_Up$bar_df %>% dplyr::mutate(direction = "CLIP_Down_RT_Up"),
      out_kw$Up_RT_Down$bar_df %>% dplyr::mutate(direction = "CLIP_Up_RT_Down")
    ) %>%
      dplyr::mutate(
        category = factor(category, levels = c("synaptic_terms","nuclear_terms","other")),
        direction = factor(direction, levels = c("CLIP_Down_RT_Up","CLIP_Up_RT_Down"))
      )
    
    cliprt_nucsyn_tables[[comp]] <- stack_df %>%
      dplyr::mutate(comp = comp)
    
    p_stack <- ggplot(stack_df, aes(x = direction, y = pct, fill = category)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = CLIP_RT_fill_colors, drop = FALSE) +
      theme_classic(base_size = 9) +
      labs(x = NULL, y = "% of genes", fill = NULL, title = comp)
    
    print(p_stack)
  }
  
  prio_down_up <- get_priority_cliprt(comp, "CLIP_Down_RT_Up")
  prio_up_down <- get_priority_cliprt(comp, "CLIP_Up_RT_Down")
  
  # choose which term list to use based on priority
  terms_for_down_up <- if (prio_down_up == "nuclear_first") nuclear_subset_terms else synaptic_subset_terms
  fill_for_down_up  <- if (prio_down_up == "nuclear_first") "#D95F02" else "#1B9E77"
  label_for_down_up <- if (prio_down_up == "nuclear_first") "nuclear" else "synaptic"
  
  terms_for_up_down <- if (prio_up_down == "nuclear_first") nuclear_subset_terms else synaptic_subset_terms
  fill_for_up_down  <- if (prio_up_down == "nuclear_first") "#D95F02" else "#1B9E77"
  label_for_up_down <- if (prio_up_down == "nuclear_first") "nuclear" else "synaptic"
  
  p_labeled <- NULL
  
  if (!is.null(out_kw$Down_RT_Up) && !is.null(out_kw$Up_RT_Down)) {
    
    m <- get_tp(comp)
    clip_lfc <- m$clip_lfc[[1]]
    clip_p   <- m$clip_p[[1]]
    rt_lfc   <- m$rt_lfc[[1]]
    rt_p     <- m$rt_p[[1]]
    
    prio_down_up <- get_priority_cliprt(comp, "CLIP_Down_RT_Up")
    prio_up_down <- get_priority_cliprt(comp, "CLIP_Up_RT_Down")
    
    get_gene_set <- function(one_dir_out, which = c("synaptic", "nuclear")) {
      which <- match.arg(which)
      if (which == "synaptic") {
        if (!is.null(one_dir_out$synaptic_gene_set) && nrow(one_dir_out$synaptic_gene_set) > 0) {
          return(unique(as.character(one_dir_out$synaptic_gene_set$gene)))
        } else return(character(0))
      } else {
        if (!is.null(one_dir_out$nuclear_gene_set) && nrow(one_dir_out$nuclear_gene_set) > 0) {
          return(unique(as.character(one_dir_out$nuclear_gene_set$gene)))
        } else return(character(0))
      }
    }
    
    label_genes_down_up <- if (prio_down_up == "synaptic_first") {
      get_gene_set(out_kw$Down_RT_Up, "synaptic")
    } else {
      get_gene_set(out_kw$Down_RT_Up, "nuclear")
    }
    
    label_genes_up_down <- if (prio_up_down == "synaptic_first") {
      get_gene_set(out_kw$Up_RT_Down, "synaptic")
    } else {
      get_gene_set(out_kw$Up_RT_Down, "nuclear")
    }
    
    label_color_down_up <- if (prio_down_up == "synaptic_first") "#1B9E77" else "#D95F02"
    label_color_up_down <- if (prio_up_down == "synaptic_first") "#1B9E77" else "#D95F02"
    
    label_down_up <- All_CLIP_RT_data_joined %>%
      dplyr::filter(
        gene_name %in% label_genes_down_up,
        !is.na(.data[[clip_p]]), !is.na(.data[[rt_p]]),
        .data[[clip_p]] < 0.05, .data[[rt_p]] < 0.05,
        .data[[clip_lfc]] < 0,  .data[[rt_lfc]] > 0
      ) %>%
      dplyr::arrange(
        dplyr::desc(-log10(.data[[rt_p]]) + -log10(.data[[clip_p]]))
      ) %>%
      dplyr::slice_head(n = 30)
    
    label_up_down <- All_CLIP_RT_data_joined %>%
      dplyr::filter(
        gene_name %in% label_genes_up_down,
        !is.na(.data[[clip_p]]), !is.na(.data[[rt_p]]),
        .data[[clip_p]] < 0.05, .data[[rt_p]] < 0.05,
        .data[[clip_lfc]] > 0, .data[[rt_lfc]] < 0
      ) %>%
      dplyr::arrange(
        dplyr::desc(-log10(.data[[rt_p]]) + -log10(.data[[clip_p]]))
      ) %>%
      dplyr::slice_head(n = 30)
    
    nudge_down_up <- c(x = -0.5, y =  0.5)
    nudge_up_down <- c(x =  0.5, y = -0.5)
    
    p_labeled <- ggplot(All_CLIP_RT_data_joined) +
      ylab(paste0("RiboTag log2FC\n", comp)) +
      xlab(paste0("FMRP CLIP log2FC\n", comp)) +
      xlim(-4, 4) +
      ylim(-2, 2) +
      theme_classic(base_size = 9) +
      geom_vline(xintercept = 0, linewidth = .3) +
      geom_hline(yintercept = 0, linewidth = .3) +
      theme(legend.position = "none") +
      # points (sig in both)
      geom_point(
        data = q$df_sig %>% dplyr::filter(.data[[rt_lfc]] > 0, .data[[clip_lfc]] < 0),
        aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Down_RT_Up"),
        size = 2
      ) +
      geom_point(
        data = q$df_sig %>% dplyr::filter(.data[[rt_lfc]] < 0, .data[[clip_lfc]] < 0),
        aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Down_RT_Down"),
        size = 2
      ) +
      geom_point(
        data = q$df_sig %>% dplyr::filter(.data[[rt_lfc]] > 0, .data[[clip_lfc]] > 0),
        aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Up_RT_Up"),
        size = 2
      ) +
      geom_point(
        data = q$df_sig %>% dplyr::filter(.data[[rt_lfc]] < 0, .data[[clip_lfc]] > 0),
        aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], color = "CLIP_Up_RT_Down"),
        size = 2
      ) +
      scale_color_manual(values = CLIP_RT_cat_colors, drop = FALSE) +
      ggrepel::geom_label_repel(
        data = label_down_up,
        aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], label = gene_name),
        nudge_x = nudge_down_up["x"], nudge_y = nudge_down_up["y"],
        color = label_color_down_up, fill = "white", size = 2.8,
        box.padding = 0.25, point.padding = 0.15,
        min.segment.length = 0.1,
        segment.color = "black", segment.size = 0.25,
        max.overlaps = Inf, show.legend = FALSE
      ) +
      ggrepel::geom_label_repel(
        data = label_up_down,
        aes(x = .data[[clip_lfc]], y = .data[[rt_lfc]], label = gene_name),
        nudge_x = nudge_up_down["x"], nudge_y = nudge_up_down["y"],
        color = label_color_up_down, fill = "white", size = 2.8,
        box.padding = 0.25, point.padding = 0.15,
        min.segment.length = 0.1,
        segment.color = "black", segment.size = 0.25,
        max.overlaps = Inf, show.legend = FALSE
      )
    
    print(p_labeled)
    ggsave(
      file.path(Outdirectory, paste0("CLIPscore_RiboTag_scatterplot_labeled_", comp, ".pdf")),
      p_labeled, width = 6, height = 6, units = "in"
    )
  }
  
  clip_rt_comp_out[[comp]] <- list(
    scatter_unlabeled = p_unlab,
    sets = q$sets,
    go_down_up = go_down_up,
    go_up_down = go_up_down,
    kw = out_kw,
    stack_df = stack_df,
    stack_plot = p_stack,
    scatter_labeled = p_labeled
  )
}

if (length(cliprt_nucsyn_tables) > 0) {
  
  CLIP_RT_all_stack_df <- dplyr::bind_rows(cliprt_nucsyn_tables) %>%
    dplyr::mutate(
      comp = factor(comp, levels = c("5v0", "30v5")),
      direction = factor(direction, levels = c("CLIP_Down_RT_Up", "CLIP_Up_RT_Down")),
      category = factor(category, levels = c("synaptic_terms", "nuclear_terms", "other"))
    )
  
  CLIP_RT_all_stack_faceted <- ggplot(CLIP_RT_all_stack_df,
                                      aes(x = direction, y = pct, fill = category)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = CLIP_RT_fill_colors, drop = FALSE) +
    theme_classic(base_size = 9) +
    labs(x = NULL, y = "% of genes", fill = NULL) +
    facet_wrap(~comp, nrow = 1) +
    scale_x_discrete(labels = c(
      "CLIP_Down_RT_Up" = "CLIP Down /\nRT Up",
      "CLIP_Up_RT_Down" = "CLIP Up /\nRT Down"
    ))
  print(CLIP_RT_all_stack_faceted)
  ggsave(filename = file.path(Outdirectory, "CLIPscore_RiboTag_nuclear_synaptic_other_faceted.pdf"), plot = CLIP_RT_all_stack_faceted, width = 9, height = 4, units = "in")
}

print(CLIP_RT_all_stack_df)

#category       n_genes bar                       total_genes   pct direction       comp 
#other               12 5v0_CLIP_5v0_Down_RT_Up          32  37.5 CLIP_Down_RT_Up 5v0 
#synaptic_terms      20 5v0_CLIP_5v0_Down_RT_Up          32  62.5 CLIP_Down_RT_Up 5v0 
#nuclear_terms       18 5v0_CLIP_5v0_Up_RT_Down          48  37.5 CLIP_Up_RT_Down 5v0 
#other               22 5v0_CLIP_5v0_Up_RT_Down          48  45.8 CLIP_Up_RT_Down 5v0 
#synaptic_terms       8 5v0_CLIP_5v0_Up_RT_Down          48  16.7 CLIP_Up_RT_Down 5v0 
#nuclear_terms       67 30v5_CLIP_30v5_Down_RT_Up           163  41.1 CLIP_Down_RT_Up 30v5  
#other               76 30v5_CLIP_30v5_Down_RT_Up           163  46.6 CLIP_Down_RT_Up 30v5  
#synaptic_terms      20 30v5_CLIP_30v5_Down_RT_Up           163  12.3 CLIP_Down_RT_Up 30v5  
#other                6 30v5_CLIP_30v5_Up_RT_Down            17  35.3 CLIP_Up_RT_Down 30v5  
#synaptic_terms      11 30v5_CLIP_30v5_Up_RT_Down            17  64.7 CLIP_Up_RT_Down 30v5 

################################################################################################################################
# make bubble plot for CLIP/RT quadrant GO terms

nuclear_subset_terms <- c(
  "nuclear export",
  "histone modification",
  "protein localization to nucleus",
  "chromatin remodeling",
  "miRNA transcription",
  "nucleus organization",
  "nuclear transport",
  "RNA localization",
  "RNA splicing",
  "DNA-templated transcription initiation"
)

synaptic_subset_terms <- c(
  "chemical synaptic transmission, postsynaptic",
  "synapse assembly",
  "receptor internalization",
  "vesicle-mediated transport in synapse",  
  "receptor-mediated endocytosis",  
  "regulation of postsynaptic membrane potential",
  "ionotropic glutamate receptor signaling pathway",
  "regulation of NMDA receptor activity",
  "protein localization to postsynapse",
  "monoatomic ion transmembrane transport",
  "calcium ion homeostasis"
)

target_terms <- c(nuclear_subset_terms, synaptic_subset_terms)

make_cliprt_go_bubble_df <- function(clip_rt_comp_out, target_terms) {
  
  go_tables <- list(
    "5v0_Down_RT_Up" = clip_rt_comp_out[["5v0"]]$go_down_up,
    "5v0_Up_RT_Down" = clip_rt_comp_out[["5v0"]]$go_up_down,
    "30v5_Down_RT_Up"  = clip_rt_comp_out[["30v5"]]$go_down_up,
    "30v5_Up_RT_Down"  = clip_rt_comp_out[["30v5"]]$go_up_down
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

CLIPRT_go_bubble_df <- make_cliprt_go_bubble_df(
  clip_rt_comp_out = clip_rt_comp_out,
  target_terms = target_terms
)

x_order <- c(
  "5v0_Down_RT_Up",
  "5v0_Up_RT_Down",
  "30v5_Down_RT_Up",
  "30v5_Up_RT_Down"
)

term_order <- c(nuclear_subset_terms, synaptic_subset_terms)

CLIPRT_go_bubble_df <- CLIPRT_go_bubble_df %>%
  dplyr::mutate(
    biology = dplyr::case_when(
      Description %in% nuclear_subset_terms  ~ "nuclear",
      Description %in% synaptic_subset_terms ~ "synaptic",
      TRUE ~ "other"
    ),
    neglog10_padj_raw = -log10(p.adjust),
    neglog10_padj = pmin(neglog10_padj_raw, 12),
    category = factor(category, levels = x_order),
    Description = factor(Description, levels = rev(term_order))
  )

CLIPRT_go_bubble_plot_df <- CLIPRT_go_bubble_df %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05,
    Count > 0
  )

CLIPRT_p_go_bubble <- ggplot(
  CLIPRT_go_bubble_plot_df,
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
    breaks = c( 2, 4, 8, 10, 12)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "P.value",
    limits = c(1.3, 8),
    breaks = c(1.3, 2, 4, 6, 8),
    labels = c("0.05","0.1", "1e-4", "1e-6", "≤1e-8"),
    guide = guide_colorbar(barheight = grid::unit(4, "cm"))
  ) +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )
CLIPRT_p_go_bubble
ggsave(file.path(Outdirectory, "CLIPscore_RiboTag_GO_terms_bubbleplot.pdf"), CLIPRT_p_go_bubble, width = 6, height = 5.5, units = "in", device = grDevices::cairo_pdf)

################################################################################################################################
comp_order <- c("5v0","30v5")
quadrants <- c("CLIP_Down_RT_Up","CLIP_Down_RT_Down","CLIP_Up_RT_Up","CLIP_Up_RT_Down")

get_quad_set <- function(clip_rt_comp_out, comp, quad) {
  if (!comp %in% names(clip_rt_comp_out)) return(character(0))
  s <- clip_rt_comp_out[[comp]]$sets[[quad]]
  unique(as.character(s))
}

venn_sets_quads <- list()
for (quad in quadrants) {
  for (comp in comp_order) {
    venn_sets_quads[[paste0(comp, "_", quad)]] <- get_quad_set(clip_rt_comp_out, comp, quad)
  }
}

vd_quads <- venndetail(venn_sets_quads)

pdf(file.path(Outdirectory, "VennDetail_CLIPscore_RiboTag.pdf"), width = 14, height = 7)
plot(vd_quads, type = "upset")
dev.off()

plot(vd_quads, type = "upset")

##############################################################################################################
#merge with Hale et al 2021 RT data from microdissected CA1 soma and neuropils

CH_RAS_merged_resdata=merge(All_CLIP_RT_data_joined,CH_RT_resdata, by.x= "gene_name" , by.y="gene_name")
nrow(CH_RAS_merged_resdata) #12117

#CLIP plus Ribo- 5v0
CLIP_RT_5v0 <- list(
  CLIP_Bound_5v0        = CH_RAS_merged_resdata %>% filter(CLIP_5v0_pvalue < 0.05, RiboTag_5v0_pvalue >= 0.05) %>% pull(gene_name),
  CLIP_Down_RT_Up_5v0   = CH_RAS_merged_resdata %>% filter(CLIP_5v0_pvalue < 0.05, CLIP_5v0_log2FC < 0, RiboTag_5v0_pvalue < 0.05, RiboTag_5v0_log2FC > 0) %>% pull(gene_name),
  CLIP_Down_RT_Down_5v0 = CH_RAS_merged_resdata %>% filter(CLIP_5v0_pvalue < 0.05, CLIP_5v0_log2FC < 0, RiboTag_5v0_pvalue < 0.05, RiboTag_5v0_log2FC < 0) %>% pull(gene_name),
  CLIP_Up_RT_Up_5v0     = CH_RAS_merged_resdata %>% filter(CLIP_5v0_pvalue < 0.05, CLIP_5v0_log2FC > 0, RiboTag_5v0_pvalue < 0.05, RiboTag_5v0_log2FC > 0) %>% pull(gene_name),
  CLIP_Up_RT_Down_5v0   = CH_RAS_merged_resdata %>% filter(CLIP_5v0_pvalue < 0.05, CLIP_5v0_log2FC > 0, RiboTag_5v0_pvalue < 0.05, RiboTag_5v0_log2FC < 0) %>% pull(gene_name))

# Define color palette
CLIP_RT_CH_cat_colors <- c(
  "CLIP_Bound_5v0"   = "grey",
  "CLIP_Down_RT_Up_5v0"     = "violet",
  "CLIP_Down_RT_Down_5v0"   = "gold",
  "CLIP_Up_RT_Up_5v0"       = "green3",
  "CLIP_Up_RT_Down_5v0"     = "#223F99"
)

CLIP_RT_interaction_genes <- CLIP_RT_5v0

all_gene_data <- list()

for (group_name in names(CLIP_RT_interaction_genes)) {
  gene_list <- CLIP_RT_interaction_genes[[group_name]]
  
  gene_data <- CH_RAS_merged_resdata %>%
    filter(gene_name %in% gene_list) %>%
    mutate(Group = group_name)
  
  all_gene_data[[group_name]] <- gene_data
}

CLIP_RT_CH_combined_df_5v0 <- bind_rows(all_gene_data)

CLIP_RT_CH_combined_df_5v0$Group <- factor(CLIP_RT_CH_combined_df_5v0$Group, levels = c(
  "CLIP_Bound_5v0",
  "CLIP_Down_RT_Up_5v0",
  "CLIP_Down_RT_Down_5v0",
  "CLIP_Up_RT_Up_5v0",
  "CLIP_Up_RT_Down_5v0"))

pairwise_stats <- CLIP_RT_CH_combined_df_5v0 %>%
  rstatix::pairwise_wilcox_test(CH_NPvsCB_log2FC ~ Group, p.adjust.method = "none")

sig_pairs <- pairwise_stats %>%
  filter(p.adj < 0.05) %>%
  mutate(p.signif = rstatix::p_format(p.adj, stars = TRUE))

my_sig_comparisons <- sig_pairs %>%
  dplyr::select(group1, group2) %>%
  split(., seq(nrow(.))) %>%
  lapply(as.character)

CLIP_RT_CH_violin_plot_5v0 <- ggplot(
  CLIP_RT_CH_combined_df_5v0,
  aes(x = Group, y = CH_NPvsCB_log2FC, fill = Group)
) +
  # violin
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  
  # add boxplot on top (like ggviolin)
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    alpha = 0.9
  ) +
  scale_fill_manual(values = CLIP_RT_CH_cat_colors) +
  ggpubr::stat_compare_means(
    comparisons = my_sig_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01
  ) +
  ylab("CH RiboTag log2 Fold Change NP/CB") +
  xlab(NULL) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(size = 9),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  )
CLIP_RT_CH_violin_plot_5v0
ggsave(file.path(Outdirectory,"CLIPscore_RiboTag_subcellular_ViolinPlot_PairwiseWilcoxon_5v0.pdf"), plot = CLIP_RT_CH_violin_plot_5v0, width = 10, height = 6)

#30v5
CLIP_RT_30v5 <- list(
  CLIP_Bound_30v5        = CH_RAS_merged_resdata %>% filter(CLIP_30v5_pvalue < 0.05, RiboTag_30v5_pvalue >= 0.05) %>% pull(gene_name),
  CLIP_Down_RT_Up_30v5   = CH_RAS_merged_resdata %>% filter(CLIP_30v5_pvalue < 0.05, CLIP_30v5_log2FC < 0, RiboTag_30v5_pvalue < 0.05, RiboTag_30v5_log2FC > 0) %>% pull(gene_name),
  CLIP_Down_RT_Down_30v5 = CH_RAS_merged_resdata %>% filter(CLIP_30v5_pvalue < 0.05, CLIP_30v5_log2FC < 0, RiboTag_30v5_pvalue < 0.05, RiboTag_30v5_log2FC < 0) %>% pull(gene_name),
  CLIP_Up_RT_Up_30v5     = CH_RAS_merged_resdata %>% filter(CLIP_30v5_pvalue < 0.05, CLIP_30v5_log2FC > 0, RiboTag_30v5_pvalue < 0.05, RiboTag_30v5_log2FC > 0) %>% pull(gene_name),
  CLIP_Up_RT_Down_30v5   = CH_RAS_merged_resdata %>% filter(CLIP_30v5_pvalue < 0.05, CLIP_30v5_log2FC > 0, RiboTag_30v5_pvalue < 0.05, RiboTag_30v5_log2FC < 0) %>% pull(gene_name))

# Define color palette
CLIP_RT_CH_cat_colors <- c(
  "CLIP_Bound_30v5"   = "grey",
  "CLIP_Down_RT_Up_30v5"     = "violet",
  "CLIP_Down_RT_Down_30v5"   = "gold",
  "CLIP_Up_RT_Up_30v5"       = "green3",
  "CLIP_Up_RT_Down_30v5"     = "#223F99"
)

CLIP_RT_interaction_genes <- CLIP_RT_30v5

all_gene_data <- list()

for (group_name in names(CLIP_RT_interaction_genes)) {
  gene_list <- CLIP_RT_interaction_genes[[group_name]]
  
  gene_data <- CH_RAS_merged_resdata %>%
    filter(gene_name %in% gene_list) %>%
    mutate(Group = group_name)
  
  all_gene_data[[group_name]] <- gene_data
}

CLIP_RT_CH_combined_df_30v5 <- bind_rows(all_gene_data)

CLIP_RT_CH_combined_df_30v5$Group <- factor(CLIP_RT_CH_combined_df_30v5$Group, levels = c(
  "CLIP_Bound_30v5",
  "CLIP_Down_RT_Up_30v5",
  "CLIP_Down_RT_Down_30v5",
  "CLIP_Up_RT_Up_30v5",
  "CLIP_Up_RT_Down_30v5"))

pairwise_stats <- CLIP_RT_CH_combined_df_30v5 %>%
  rstatix::pairwise_wilcox_test(CH_NPvsCB_log2FC ~ Group, p.adjust.method = "none")

sig_pairs <- pairwise_stats %>%
  filter(p.adj < 0.05) %>%
  mutate(p.signif = rstatix::p_format(p.adj, stars = TRUE))

my_sig_comparisons <- sig_pairs %>%
  dplyr::select(group1, group2) %>%
  split(., seq(nrow(.))) %>%
  lapply(as.character)

CLIP_RT_CH_violin_plot_30v5 <- ggplot(
  CLIP_RT_CH_combined_df_30v5,
  aes(x = Group, y = CH_NPvsCB_log2FC, fill = Group)
) +
  # violin
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  
  # add boxplot on top (like ggviolin)
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    alpha = 0.9
  ) +
  scale_fill_manual(values = CLIP_RT_CH_cat_colors) +
  ggpubr::stat_compare_means(
    comparisons = my_sig_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01
  ) +
  ylab("CH RiboTag log2 Fold Change NP/CB") +
  xlab(NULL) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(size = 9),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  )
CLIP_RT_CH_violin_plot_30v5
ggsave(file.path(Outdirectory,"CLIPscore_RiboTag_subcellular_ViolinPlot_PairwiseWilcoxon_30v5.pdf"), plot = CLIP_RT_CH_violin_plot_30v5, width = 10, height = 6)

################################################################################################################################
#optional write out dfs

#write_csv(All_CLIP_RT_data_joined,file.path(Outdirectory, "All_CLIP_RT_data_joined.csv"))

#write_csv(CH_RAS_merged_resdata,file.path(Outdirectory, "All_CLIP_RT_CH_data_joined.csv"))




