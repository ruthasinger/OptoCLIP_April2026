#load the libraries needed

library(tidyverse)   
library(rtracklayer)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GenomicFeatures)

##############################################################################################################
#setup the directories

mm10_gtf_directory=file.path("~/ruthasinger_github/OptoCLIP_April2026/RDS_files/Data_mm10gtf")
list.files(mm10_gtf_directory)

##############################################################################################################
#Bring in GTF for transcript information
mm10_gtf=rtracklayer::import(file.path(mm10_gtf_directory,"mm10.ensGene.gtf.gz"))
mm10_gtf_df=mm10_gtf %>% as_tibble
mm10_Tx=mm10_gtf_df %>% dplyr::filter(type=="transcript")
mm10_Tx_use=mm10_Tx %>% dplyr::select(seqnames,start,end,width,strand,gene_id,transcript_id)
mm10_Tx_use$gene_name<- AnnotationDbi::mapIds(org.Mm.eg.db,
                                 keys = mm10_Tx_use$gene_id,
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")

mm10_Tx_final_noNAs=subset(mm10_Tx_use, is.na(gene_name) == FALSE)

mm10_Tx2genename=mm10_Tx_final_noNAs %>% dplyr::select(transcript_id,gene_name)

mm10_Tx_genename_width=mm10_Tx_final_noNAs %>% dplyr::select(transcript_id,gene_name,width)

mm10_txdb <- GenomicFeatures::makeTxDbFromGFF(file.path(mm10_gtf_directory,"mm10.ensGene.gtf.gz"))

mm10_txdb_lengths=GenomicFeatures::transcriptLengths(mm10_txdb, with.cds_len=TRUE,
                                    with.utr5_len=TRUE, with.utr3_len=TRUE)

mm10_Tx_final=merge(mm10_Tx_final_noNAs,mm10_txdb_lengths %>% dplyr::select(tx_name,tx_len), by.x="transcript_id", by.y="tx_name")

##############################################################################################################
#save RDS files for future use
saveRDS(mm10_Tx_final, file=file.path(mm10_gtf_directory,"mm10_Tx_final.rds"))
saveRDS(mm10_txdb_lengths, file=file.path(mm10_gtf_directory,"mm10_txdb_lengths.rds"))
