# RNAseq data: PE
# AML
# Dr. Jerome T
# Mapping using Salmon (salmon 1.1.0)
# Data import using tximport
# Reference: RNA-seq workflow: gene-level exploratory analysis and differential expression
# Prasad Chaskar

####################################################
# Analysis setup

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path 

dirpath<-dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# Load libraries
source("./scripts/libraries.R")

# Load raw data

load("./RDS/raw_data.RData")

# Removed samples 206 and 223 sample samples run twice
# Removed samples and 208
raw_counts<-raw_counts[,-c(19, 21)]

colnames(raw_counts)

####################################################
# Meta Data Import
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

# Sample sheet/ Meta data file
sample_info<-read.table(file = "./metadata/sample_sheet_revised_Oct2020.txt",
                        row.names = 1, sep = "\t", header = TRUE, check.names = FALSE)

dim(sample_info)

#Rename samples 1>J1
rownames(sample_info) <- gsub("^", "J", rownames(sample_info))

#Add sample names to the meta data
sample_info$ID<-rownames(sample_info)

#Reorder, select the samples based on the rawdata
sample_info<-sample_info[match(colnames(raw_counts), sample_info$ID),]

dim(sample_info)

#Set the desired variables/ columns as factor
cols <- c("EED","ELN_2017","ETV6","EZH2","FLT3_ITD","FLT3_TKD","G12","GATA2",
          "ID","IDH1","IDH2","JAK2","KIT","KRAS","LDH","MLL","NF1","NF1_QPCR",
          "NF1_del","NF1_mut","NPM1","NRAS","OS_status","PFS_status","PHF6",
          "PTPN11","RASOPATHY","RUNX1","Response","SETBP1","SF3B1","SOS1",
          "SRSF2","STAG2","SUZ12","Sex","TET2","TP53","U2AF1","WT1","ZRSR2",
          "allo","aneuploidy_TP53","chromatin_spliceosome","followup_PFS",
          "followup_month_OS","inv3","splenomegaly","treatment")

sample_info[cols] <- lapply(sample_info[cols], factor)  

#Updated Deseq object 
coldata<-sample_info

dds <- DESeqDataSetFromMatrix(raw_counts, coldata, ~RASOPATHY)

meta_info<-sample_info

# Save multiple objects
save(dds, meta_info, file = "./RDS/deseq2_object.RData")
