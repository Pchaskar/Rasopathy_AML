# Figure Manuscript AML

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

library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(1)

####################################################
# Meta Data Import
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

# Sample sheet/ Meta data file
sample_info<-read.table(file = "./metadata/sample_sheet_revised_Oct2020.txt",
                        row.names = 1, sep = "\t", header = TRUE, check.names = FALSE)

dim(sample_info)

####################################################
# selected parameters

cols<-c("RASOPATHY", "Adverse_K", "Complex_K", "ELN_2017", "NF1_del", "NF1_mut", 
        "KRAS", "NRAS", "PTPN11", "CBL", "BRAF", "SOS1", "EZH2", "SUZ12", 
        "EED", "ASXL1", "BCOR", "CEBPA_mono", "DNMT3A", 
        "FLT3_ITD", "FLT3_TKD", "IDH1", "IDH2", "PHF6", "RUNX1", "SF3B1", 
        "SRSF2", "STAG2", "TET2", "TP53")

sample_info_subset<-sample_info[,cols]

#order according to rasopathy
sample_info_subset<-sample_info_subset[order(sample_info_subset$RASOPATHY, decreasing = TRUE),]

# reassign certain values
# as character
sample_info_subset[cols] <- lapply(sample_info_subset[cols], as.character)

sample_info_subset$RASOPATHY[sample_info_subset$RASOPATHY == "1"] <- "5"
sample_info_subset$RASOPATHY[sample_info_subset$RASOPATHY == "0"] <- "6"

sample_info_subset$ELN_2017[sample_info_subset$ELN_2017 == "0"] <- "7"
sample_info_subset$ELN_2017[sample_info_subset$ELN_2017 == "1"] <- "8"
sample_info_subset$ELN_2017[sample_info_subset$ELN_2017 == "2"] <- "9"


sample_info_subset$NF1_del[sample_info_subset$NF1_del == "0"] <- "3"
sample_info_subset$NF1_del[sample_info_subset$NF1_del == "1"] <- "4"

#as factor
sample_info_subset[cols] <- lapply(sample_info_subset[cols], factor)
#######################################

# Transpose

sample_info_subset_tp<-t(sample_info_subset)
colnames(sample_info_subset_tp) <- paste0("J", colnames(sample_info_subset_tp), sep = "")

sample_info_subset_tp<-as.data.frame(sample_info_subset_tp)
sample_info_subset_tp$gene<-rownames(sample_info_subset_tp)
sample_info_subset_tp$gene<-factor(sample_info_subset_tp$gene)

# reorder gene
sample_info_subset_tp<-sample_info_subset_tp[match(cols, sample_info_subset_tp$gene),]

#sample_info_subset_tp<-as.matrix(sample_info_subset_tp)

# LOng format

library("reshape")

sample_info_melted<-melt(sample_info_subset_tp, id.vars="gene")

colnames(sample_info_melted)<-c("gene", "sample", "mutated")

sample_info_melted$mutated<-factor(sample_info_melted$mutated)
sample_info_melted$gene<-factor(sample_info_melted$gene)

sample_info_melted$gene<- factor(sample_info_melted$gene, levels=sample_info_subset_tp$gene)
sample_info_melted<-sample_info_melted[order(sample_info_melted$gene),]

f1a<-ggplot(sample_info_melted, aes(x=sample, y=gene, fill=mutated)) + geom_tile(color="white", size=0.5) +
  coord_equal() +
  labs(x=NULL, y=NULL, title="") +
  theme_tufte(base_family="Helvetica") +
  scale_fill_manual(values = c("gray", "red", "gold", "darkgreen",
                               "blue", "lightgreen", 
                               "pink", "magenta", "darkred")) +
  theme(axis.ticks=element_blank()) + 
  theme(axis.text.x=element_text(angle = 90, hjust = 1, colour = "black"),
        axis.text.y=element_text(colour = "black"))  

f1a<-f1a+theme(legend.position = "none")

png(paste("./images/Figure1",Sys.Date(), "matrix.png", sep = "_"), width = 10, height = 8,
    units = 'in', res = 300)
print(f1a)
dev.off()
