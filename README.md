# Opto-CLIP Analysis Pipeline



## Overview

This repository includes all scripts used to generate the bioinformatic components main and supplementary figures for the following manuscript:

**Neuronal activation triggers rapid biphasic FMRP-mediated translational control linking synaptic and nuclear regulation**

This study implements a multi-modal analysis pipeline integrating FMRP CLIP and RiboTag RNA-seq to investigate activity-dependent RNA regulation in CA1 neurons following optogenetic stimulation.



## Repository Structure
```
├── Figure2_FigureS2_Ephys/
  └── OptoCLIP_Figure2_FigureS2_Rcode_final.R
  └── Data
      └── 5trains5Hz
      └── CurrentClamp
      └── DifferentHz
├── Figure3_FigureS3_FMRP_CLIP/
  └── OptoCLIP_Figure3_FigureS3_Rcode_final.R
  └── Data
      └── Data_FMRP_CLIP_Tx_STAR
├── Figure4_FigureS4_RiboTag/
  └── OptoCLIP_Figure4_FigureS4_Rcode_final.R
  └── Data
      └── Data_FMRP_CLIP_Tx_STAR
      └── Data_RiboTag
├── Figure5_FigureS5_Figure6_FigureS6/
  └── OptoCLIP_Figure5_FigureS5_Figure6_FigureS6_Rcode_final.R
├── R_functions/
├── RDS_files/
```



## **Figures**

### Figure 2 and Figure S2 — Electrophysiology

#### **Script:**

- OptoCLIP_Figure2_FigureS2_Rcode_final.R

#### **Input data:**

- **CurrentClamp**

  This folder contains 3 raw abf files corresponding to 3 different patched cells that underwent stepwise injection of -50 to 200 pA current in 25 pA   increments.

- **DifferentHz**

  This folder contains 6 raw abf files. The files correspond to either ChR2:mCherry+ cells or 3 Control (ChR2 negative) cells that were patched during exposure to 3 different blue LED frequencies (1, 5, and 10 Hz)

- **5trains5Hz**

  This folder contains 5 raw abf files corresponding to 3 ChR2:mCherry+ cells or 2 Control (ChR2 negative) cells that were patched during exposure to 5 trains of 30" 5 Hz blue LED light followed by 30" light off

#### **Output (11 plots):**

- **Figures 2C-2E:** Whole cell patch clamp traces of ChR2+ and ChR2- cells exposed to 10 seconds of 1Hz, 5Hz, or 10Hz blue LED light pulses

- **Figures 2C'-2E':** Zoomed in 1 second traces of 2C-2E

- **Figures 2F and 2G:** Whole cell patch clamp traces of ChR2+ and ChR2- cells exposed to trains of 30" 5Hz blue light followed by 30" dark repeated for a total of 5 trains. 

- **Figure S2A and S2B:** Current clamp traces of ChR2+ cells that underwent stepwise injection of -50 to 200 pA current in 25 pA increments. Rheobase current (+50 pA) is highlighted in green

- **Figure S2C:** Quantification of number of induced action potentials by pA of injected current. Shading indicates +/- SEM from 3 different patched cells. 



### Figure 3 and Figure S3 — Opto FMRP-CLIP

#### **Script:**

- OptoCLIP_Figure3_FigureS3_Rcode_final.R

#### **Input data:**

- **Data_FMRP_CLIP_Tx_STAR**

  This folder contains 15 bed files (*.transcriptome.unique.bed). 3 Cre negative samples, 4 Control samples, 4 Opto_5minute samples, and 4 Opto_30minute samples. Code for generating these files from raw, demultiplexed fastq files (found on GEO: GSE286379) can be found in this data folder entitled "RAS_bashscript_FMRPCLIP_fastq2transcriptome"

#### **Output (5 plots):** 

- **Figure 3D:** Quantification of uniquely mapped FMRP CLIP tags by condition. 

- **Figure 3E:** Bar plots showing genomic annotations of uniquely mapped FMRP CLIP tags separated by transcripts and regions found in increasing biological replicates.

- **Figure 3F:** Meta transcript coverage of 1000 transcripts with most CLIP tags separated by condition.

- **Figure 3G:** Boxplots depicting the distribution of transcript-level CLIP tag proportions in 5′UTR, CDS, and 3′UTR regions separated by condition.

- **Figure S3A:** Heatmap of pairwise pearson correlation analysis of each sample compared to every other sample.



### Figure 4 and Figure S4 — Opto RiboTag

#### **Script:**

- OptoCLIP_Figure4_FigureS4_Rcode_final.R

#### **Input data:**

- **Data_FMRP_CLIP_Tx_STAR**

  This folder contains 16 bed files (*.transcriptome.unique.bed). 4 Control_5minute samples, 4 Control_30minute samples, 4 Opto_5minute samples, and 4 Opto_30minute samples. 

- **Data_RiboTag**

  - **Data_InputandIP**

    This folder contains 29 folders with nested quant.sf files output from salmon transcriptome mapping. Each sample has an Input and an IP folder.

  - **Data_5min**

    This folder contains 8 folders with nested quant.sf files output from salmon transcriptome mapping of just the IP samples for 5 minute timepoint.

  - **Data_30min**

    This folder contains 16 folders with nested quant.sf files output from salmon transcriptome mapping of just the IP samples for all timepoints. Code for generating these quant.sf files from raw, paired-end fastq files (found on GEO: GSE286381) can be found in Data_RiboTag/RAS_fastq2salmon.R

- **RDS files**

  There are 3 RDS files will need to be read in as detailed in the beginning of the OptoCLIP_Figure4_FigureS4_Rcode_final.R script. These files are located in "~/ruthasinger_github/OptoCLIP_April2026/RDS_files"

#### **Output (9 plots):**

- **Figure 4E and 4F:** Volcano plots from DESeq comparing control 0 min to opto_5min and opto_30min to chr2_5min

- **Figure 4G:** Stacked bar plots of DEGs mapping to synaptic, nuclear, or other Gene Ontology (GO) terms 

- **Figure 4H:** Bubble plot showing GO analysis on Opto RiboTag DEGs  

- **Figure 4I:** Violin and box plots showing subcellular localization of Opto RiboTag DEGs 

- **Figure 4J:** Stacked bar plot showing the percentage of FMRP target and non-target transcripts among subsets of Opto RiboTag transcripts 

- **Figure S4A:** PCA on OptoRiboTag samples

- **Figure S4B:** Heatmap of enrichment of CA1 excitatory marker genes and de-enrichment of other cell type markers in RiboTag IP compared to Input samples

- **Figure S4C:** Upset plot showing overlap of RiboTag DEGs across different comparisons

- **RDS object FMRP_RiboTag_merged:** This RDS object will be used in Figures 5, S5, 6, and S6



### Figures 5, S5, 6, and S6- FMRP CLIP score and comparison to RiboTag data

#### **Script:**

- OptoCLIP_Figure5_FigureS5_Figure6_FigureS6_Rcode_final.R

#### **Input data:**

- **RDS files**

  There are 2 RDS files that will need to be read in as detailed in the beginning of the OptoCLIP_Figure5_FigureS5_Figure6_FigureS6_Rcode_final.R script. These files are located in "~/ruthasinger_github/OptoCLIP_April2026/RDS_files"

#### **Output (37 plots):**

- **Figures 5B-5D:** Scatter plots of FMRP CLIP binding vs RiboTag for each condition

- **Figures 5E:** Density plots showing the distribution of average CLIP scores for transcripts each condition. 

- **Figures 5F:** Subset of the 5E density plot showing only transcripts that are not stringent in Opto 30 minutes but are stringent

- **Figures 5G and 5H:** Scatter plots of per transcript FMRP CLIP scores comparing different conditions

- **Figures 5I and 5J:** Volcano plots comparing FMRP CLIP scores across different conditions

- **Figures 5K:** UpSet plot showing overlap between transcripts differentially bound by FMRP across different conditions. 

- **Figures 5L:** Functional classification of transcripts within each intersection category shown in (5K). 

- **Figures 6A and 6B:** Scatter plots differentially bound FMRP CLIP transcripts versus differentially ribosome associated transcripts 

- **Figures 6C-6F:** Stacked bar plots showing classification of transcripts within each regulatory category
shown in 6A and 6B. Output graphs were split in Adobe.Illustrator and aligned according to category.

- **Figures S5A-S5C:** Scatter plots of FMRP CLIP binding from individual samples compared to average RiboTag across the matching. These plots were used to calculate the linear relationship between those variables and calculate a CLIP score per transcript per sample.

- **Figures S5D:** Bar plot showing quantification of transcripts in each category of CLIP score stringency across conditions

- **Figures S5E:** Bubble plot showing Gene Ontology (GO) analysis on transcripts differentially bound by FMRP in different conditions.

- **Figures S6A and S6B:** Density ridge plots showing the distribution of RiboTag log2 fold-change subset by transcripts differentially bound by FMRP 

- **Figures S6C:** Bubble plot showing Gene Ontology (GO) analysis on transcripts differentially bound by FMRP in different conditions.

- **Figures S6D and S6E:** Violin and box plots showing subcellular localization of transcripts differentially bound by FMRP in different conditions

- **Figures S6F:** Upset plot showing overlap of transcripts differentially bound by FMRP in different



## **R_functions**

Reusable analysis functions:

- **run_deseq_contrasts.R**
    - Flexible DESeq2 pipeline with:
        - filtering strategies
        - GO enrichment
        - volcano plots
        - Venn/UpSet analysis

- **run_limma_CLIPscore.R**
    - limma-based comparison of CLIPscore data across conditions
    


## **RDS_files**

Precomputed reference objects and intermediate datasets:

- mm10_Tx_final.rds and mm10_txdb_lengths.rds
  - mm10 transcript annotations
  - Script showing how this object was made: mm10_gtf_RDS_code.R
  - Input data to make this object: RDS_files/Data_mm10gtf 
  - Used in:
    - OptoCLIP_Figure3_FigureS3_Rcode_final.R
    - OptoCLIP_Figure4_FigureS4_Rcode_final.R
- CH_RT_resdata.rds
  - External dataset from subcellular RiboTag on CA1 microdissected neuropils and cell bodies (Hale et al 2021; PMID: 34939924)
  - Script showing how this object was made: CH_RiboTag_RDS_code.R
  - Input data to make this object: RDS_files/Data_CHale_et_al_2021_RiboTag_salmon
  - Used in:
    - OptoCLIP_Figure4_FigureS4_Rcode_final.R
    - OptoCLIP_Figure5_FigureS5_Figure6_FigureS6_Rcode_final.R
- FMRP_RiboTag_merged.rds
  - Intermediate dataframe created in OptoCLIP_Figure4_FigureS4_Rcode_final.R for simplicity of downstream analyses
  - Used in OptoCLIP_Figure5_FigureS5_Figure6_FigureS6_Rcode_final.RUsed



## **Contact**

Ruth A. Singer, PhD

Postdoctoral Research Associate

Laboratory of Robert B. Darnell

The Rockefeller University

email: <rsinger01@rockefeller.edu>

