if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Herper")

library(dplyr)
library(Herper)
library(GEOquery)
library(Rsamtools)
library(GenomicFeatures)
library(stringr)

output_dir="~/ruthasinger_github/OptoCLIP_April2026/Figure4_FigureS4_RiboTag/Data/Data_RiboTag"
setwd(output_dir)

myMiniconda <- file.path(output_dir,"RAS_BRC_conda")
myMiniconda
dir.create(myMiniconda,showWarnings = TRUE)

#myMiniconda equal to general_conda in BRC workflow description

pathToConda <- install_CondaTools(c("samtools","salmon", "kallisto"), "herper", updateEnv = TRUE, pathToMiniConda = myMiniconda)

pathToConda$pathToConda

pathToSalmon <- file.path(pathToConda$pathToEnvBin,"salmon")
Res <- system2(command=pathToSalmon, args = "-h",stdout = TRUE)
Res

yml_name <- paste0("herper_", format(Sys.Date(), "%Y%m%d"), ".yml")
export_CondaEnv("herper", yml_name, pathToMiniConda = myMiniconda)

general_env <- {
  Herper::import_CondaEnv(yml_import = "RAS_BRC_conda.yml", pathToMiniConda = myMiniconda)

}

####################################################################################
#Create whitelisted fasta file
#this section describes how we create our reference fasta file

mm10_rnaseq_references_whitelist <- {
  options(timeout = 600 * 60)
  download.file("https://rubioinformatics.s3.amazonaws.com/Reference/mm10whiteListChr.txt", file.path(output_dir, "mm10whiteListChr.txt"))
  mm10_rnaseq_references_whitelist_MD5 <- tools::md5sum(file.path(output_dir, "mm10whiteListChr.txt"))
  if (mm10_rnaseq_references_whitelist_MD5 == "40b5c1dfd9ceaf9983f21c75e6c0204e") {
    write.table(mm10_rnaseq_references_whitelist_MD5, file.path(output_dir, "mm10_rnaseq_references_whitelist_MD5.txt"))
  } else {
    unlink(file.path(output_dir, "mm10whiteListChr.txt"))
    stop("Calculated MD5 didn't match expected")
  }
  file.path(output_dir, "mm10whiteListChr.txt")
}

mm10_rnaseq_references_fasta <- {
  options(timeout = 600 * 60)
  download.file("https://rubioinformatics.s3.amazonaws.com/Reference/mm10.fa.gz", file.path(output_dir, "mm10.fa.gz"))
  mm10_rnaseq_references_fasta_MD5 <- tools::md5sum(file.path(output_dir, "mm10.fa.gz"))
  if (mm10_rnaseq_references_fasta_MD5 == "db005b65828db31735f384e4c5787be5") {
    write.table(mm10_rnaseq_references_fasta_MD5, file.path(output_dir, "mm10_rnaseq_references_fasta_MD5.txt"))
  } else {
    unlink(file.path(output_dir, "mm10.fa.gz"))
    stop("Calculated MD5 didn't match expected")
  }
  file.path(output_dir, "mm10.fa.gz")
}

getwd()

gunzip(mm10_rnaseq_references_fasta, destname = gsub("[.]gz$", "",mm10_rnaseq_references_fasta), remove = FALSE)
bgzip(gsub("[.]gz$", "",mm10_rnaseq_references_fasta), dest = file.path(output_dir, "mm10_Bzipped.fa.bgz"), overwrite = TRUE)
indexFa(file.path(output_dir, "mm10_Bzipped.fa.bgz"))

mm10_BzippedWhite <- {
  FaFile <- FaFile(file = file.path(output_dir, "mm10_Bzipped.fa.bgz"), index = file.path(output_dir, "mm10_Bzipped.fa.bgz.fai"))
  WhiteListChrs <- read.table(mm10_rnaseq_references_whitelist, header = TRUE, row.names = 1) %>% rownames()
  dnastringset <- scanFa(FaFile, as = c("DNAStringSet"))
  fastawhite_File <- tempfile()
  writeXStringSet(dnastringset[WhiteListChrs], file = fastawhite_File)
  bgzip(fastawhite_File, dest = file.path(output_dir, "mm10_BzippedWhite.fa.bgz"))
  file.path(output_dir, "mm10_BzippedWhite.fa.bgz")
}

mm10_BzippedWhiteIndex <- {
  indexFa(mm10_BzippedWhite)
  file.path(output_dir, "mm10_BzippedWhite.fa.bgz.fai")
}

####################################################################################
#Create whitelisted gene models as GTF
#this section descibes how we create the GTF files containing gene models to count over.

mm10_rnaseq_references_gtf <- {
  options(timeout = 600 * 60)
  download.file("https://rubioinformatics.s3.amazonaws.com/Reference/mm10.ensGene.gtf.gz", file.path(output_dir, "mm10.ensGene.gtf.gz"))
  mm10_rnaseq_references_gtf_MD5 <- tools::md5sum(file.path(output_dir, "mm10.ensGene.gtf.gz"))
  if (mm10_rnaseq_references_gtf_MD5 == "3c961801204c9f37c4c4d86096c1f825") {
    write.table(mm10_rnaseq_references_gtf_MD5, file.path(output_dir, "mm10_rnaseq_references_gtf_MD5.txt"))
  } else {
    unlink(file.path(output_dir, "mm10.ensGene.gtf.gz"))
    stop("Calculated MD5 didn't match expected")
  }
  file.path(output_dir, "mm10.ensGene.gtf.gz")
}

mm10_BzippedGTF <- {
  WhiteListChrs <- read.table(mm10_rnaseq_references_whitelist, header = TRUE, row.names = 1) %>% rownames()
  gtfGR <- rtracklayer::import.gff(mm10_rnaseq_references_gtf)
  gtfGR <- gtfGR[seqnames(gtfGR) %in% WhiteListChrs]
  tempgtf <- paste0(tempfile(), ".gtf")
  rtracklayer::export.gff(gtfGR, con = tempgtf)
  bgzip(tempgtf, dest = file.path(output_dir, "mm10_gBzipped.gtf.bgz"), overwrite = TRUE)
}

####################################################################################
#Build transcript fasta
#this section describes how we build the fasta file of transcript sequences.
mm10_WhiteTransFastaBzipped <- {
  WhiteListChrs <- read.table(mm10_rnaseq_references_whitelist, header = TRUE, row.names = 1) %>% rownames()
  outFile <- paste0(tempfile(), ".fa")
  txdb <- makeTxDbFromGFF(mm10_BzippedGTF)
  allTranscripts <- exonsBy(txdb, use.names = T, by = "tx")
  fasta_file <- FaFile(mm10_BzippedWhite, index = mm10_BzippedWhiteIndex)
  transSeqs <- extractTranscriptSeqs(fasta_file, allTranscripts[seqnames(allTranscripts) %in%
                                                                  WhiteListChrs])
  writeXStringSet(transSeqs, file = outFile)
  bgzip(outFile, dest = file.path(output_dir, "mm10_WhiteTransFasta.fa.bgz"))
  file.path(output_dir, "mm10_WhiteTransFasta.fa.bgz")
}

mm10_gentrome <- {
  transscriptStringset <- readDNAStringSet(mm10_WhiteTransFastaBzipped)
  genomeStringset <- readDNAStringSet(mm10_BzippedWhite)
  gentromeStringset <- c(transscriptStringset, genomeStringset)
  outFile <- paste0(tempfile(), ".fa")
  writeXStringSet(gentromeStringset, file = outFile)
  bgzip(outFile, file.path(output_dir, "mm10_gentrome.fa.bgz"))
  file.path(output_dir, "mm10_gentrome.fa.bgz")
}

####################################################################################
#Build Salmon index
#this section describes how we create the salmon index files.
mm10_decoy4salmon <- {
  FastaFile <- FaFile(file = mm10_BzippedWhite, index = mm10_BzippedWhiteIndex)
  chroms <- seqinfo(FastaFile)
  write.table(as.data.frame(chroms) %>% rownames(), file = file.path(output_dir, "mm10_decoy4salmon.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  file.path(output_dir, "mm10_decoy4salmon.txt")
}

pathToSalmon <- file.path(pathToConda$pathToEnvBin,"salmon")
mm10_salmonIndex <- system2(command=pathToSalmon, 
                            args <- c("index", "-t", mm10_gentrome, "-d", mm10_decoy4salmon, "-p", 4, "-i", file.path(output_dir, "mm10_salmonIndex")),
                            stdout = TRUE)

####################################################################################
#Salmon quantification of transcript abundance
#salmon can quantify the transcript abundance for genes directly from the fastqs.
#here we calculate the number of reads per transcript using the salmon quant function.

fq_dir="~/fq"
salmon_dir=file.path(output_dir,"salmon_output")

pathToSalmon <- file.path(pathToConda$pathToEnvBin,"salmon")

RiboTag_fq_files=grep("R1_001.fastq.gz$",list.files(fq_dir),value=TRUE)

for (file in RiboTag_fq_files) {
  R2=paste0(str_extract(file, "^([^_]+_[^_]+_[^_]+_[^_]+_[^_]+)"),"_R2_001.fastq.gz")
  salmon_fname=paste0(str_extract(file, "^([^_]+_[^_]+_[^_]+_[^_]+)"),"_salmon")
  print(paste0("R1 filename is: ",file))
  print(paste0("R2 filename is: ",R2))
  print(paste0("salmon output name is:  ",salmon_fname))
  salmon_quant=system2(command=pathToSalmon, 
                       args <- c(paste0("quant"), 
                                 "-i", file.path(output_dir, "mm10_salmonIndex"), 
                                 "-l", "A", 
                                 "-1", file.path(fq_dir, file), 
                                 "-2", file.path(fq_dir, R2),
                                 "--validateMappings", 
                                 "-p", 9, 
                                 "-o", file.path(salmon_dir, salmon_fname)),
                       stdout = TRUE)
}




