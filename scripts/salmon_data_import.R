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

####################################################
# Meta Data Import
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

# Sample sheet/ Meta data file
sample_info<-read.table(file = "./metadata/sample_sheet.csv",
                        row.names = 1, sep = ",", header = TRUE, check.names = FALSE)

dim(sample_info)

sample_names<-rownames(sample_info)

####################################################
# Data Import
# Tximport salmon data 

# Salmon files

files <- file.path("/media/NAS/nas-hug/ongoing_projects/aml_unige_Jerome_20-20/new_GC/mapping", "salmon", 
                   sample_names, "quant.sf")
files

# Transcript to gene annotation
# Creating tx2gene data.frame 

txdb <- makeTxDbFromGFF("/media/NAS/nas-hug/ngs_data/cdna_homo_sapiens/release-101/Homo_sapiens.GRCh38.101.gtf")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
head(tx2gene)

# Transcript level count matrix
# Gene level count matrix based on genecode annotation (Homo_sapiens.GRCh38.99.gtf)

txi <- tximport(files, type = "salmon", tx2gene=tx2gene[,c("TXNAME", "GENEID")],  ignoreTxVersion = TRUE,
                countsFromAbundance="lengthScaledTPM")

####################################################
# Analyzing RNA-seq data with DESeq2
# Source: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#how-to-get-help-for-deseq2

# Sample information table, which we will name coldata

coldata <- sample_info

#####################
# Construct a DESeqDataSet

dds <- DESeqDataSetFromTximport(txi, coldata, ~RASOPATHY)
nrow(dds) #<- 35425

#View(counts(dds))

# Collapsed replicates (technical)
# merge counts from samples sequenced twice 

ddsColl <- collapseReplicates(dds, dds$ID)
#View(counts(ddsColl))
dim(counts(ddsColl))

#####################
# Raw counts for petros
raw_counts<- counts(ddsColl)

##########################
# Data frame of gene ID conversion
##########################
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annots <- getBM(mart=mart,
                attributes=c("ensembl_gene_id", "entrezgene_id"),
                filter="ensembl_gene_id",
                values=rownames(raw_counts),
                uniqueRows=TRUE)
annots_ent <- annots[!duplicated(annots[,1]),]
rownames(annots_ent)<-annots_ent$ensembl_gene_id

annots <- getBM(mart=mart,
                attributes=c("ensembl_gene_id", "hgnc_symbol"),
                filter="ensembl_gene_id",
                values=rownames(raw_counts),
                uniqueRows=TRUE)
annots_hgnc <- annots[!duplicated(annots[,1]),]
rownames(annots_hgnc)<-annots_hgnc$ensembl_gene_id

geneid<-annots_ent
geneid$hgnc_symbol<- annots_hgnc$hgnc_symbol
  
geneid<-apply(geneid, 2, function(x) gsub("^$|^ $", NA, x))

geneid<-as.data.frame(geneid)
geneid<-geneid[complete.cases(geneid$entrezgene_id), ]
geneid<-geneid[complete.cases(geneid$hgnc_symbol), ]

# Save multiple objects
save(raw_counts, geneid, file = "./RDS/raw_data.RData")
