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
# load functions
source("./scripts/functions.R")

# Load DESEQ2 Data
load("./RDS/DGE_edgeR.RData")

####################################
# Volcano plot
library(dplyr)
library(ggrepel)
library(EnhancedVolcano)

de_valcano<-DGE_edgeR
colnames(de_valcano)<-c("gene.id","log2FoldChange","logCPM",
                        "LR","pvalue","padj","genes")

#retain genes with valid padj values
de_valcano<-de_valcano[complete.cases(de_valcano$padj), ]

res_tableOE_tb <- de_valcano %>%
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Create an empty column to indicate which genes to label
res_tableOE_tb <- res_tableOE_tb %>% mutate(genelabels = "")

## Sort by padj values
res_tableOE_tb <- res_tableOE_tb %>% arrange(padj)

## Populate the genelabels column with contents of the gene symbols column for the first 10 rows, i.e. the top 10 most significantly expressed genes
res_tableOE_tb$genelabels <- as.character(res_tableOE_tb$genes)

View(res_tableOE_tb)

# Custom colors
# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  (res_tableOE_tb$log2FoldChange < -1.5 & res_tableOE_tb$pvalue < 0.05), 'red3',
  ifelse((res_tableOE_tb$log2FoldChange > 1.5 & res_tableOE_tb$pvalue < 0.05), 'blue','gray33'))
keyvals[is.na(keyvals)] <- 'gray33'
names(keyvals)[keyvals == 'blue'] <- 'Up in Rasopathy'
names(keyvals)[keyvals == 'gray33'] <- 'Not Significant'
names(keyvals)[keyvals == 'red3'] <- 'Down in Rasopathy'

selected_deg<-c("FREM2", "CLIP3", "RET", "GRM7", "PLK2", 
                "HOXB7", "ARHGEF10")

v<-EnhancedVolcano(res_tableOE_tb,
                   lab = res_tableOE_tb$genelabels,
                   selectLab = selected_deg,
                   x = 'log2FoldChange',
                   y = 'pvalue',
                   xlim = c(-8, 8),
                   title='Rasopathy yes versus Rasopathy no',
                   titleLabSize = 12,
                   subtitle = "Differential expression analysis",
                   caption = bquote(~Log[2]~ "fold change cutoff, 1.5; p-value cutoff, 0.05"),
                   captionLabSize = 10,
                   subtitleLabSize = 12,
                   pCutoff = 0.05,
                   FCcutoff = 1.5,
                   #pointSize = 1.0,
                   pointSize = c(ifelse(res_tableOE_tb$genelabels %in% selected_deg , 3, 1)),
                   labSize = 2.5,
                   labCol = 'black',
                   labFace = 'bold',
                   boxedLabels = TRUE,
                   colCustom = keyvals,
                   colAlpha = 4/5,
                   legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                  'p-value & Log (base 2) FC'),
                   legendPosition = 'right',
                   legendLabSize = 10,
                   legendIconSize = 4.0,
                   axisLabSize = 14,
                   drawConnectors = TRUE,
                   widthConnectors = 0.65,
                   colConnectors = 'black',
                   gridlines.major = FALSE, gridlines.minor = FALSE
                   #col=c('gray', 'gray', 'gray', 'blue'),
)

pdf(paste("./images/DEG_EdgeR",Sys.Date(), "Volcano.pdf", sep = "_"), width = 10, height = 8)
#    units = 'in', res = 400)
print(v)
dev.off()



