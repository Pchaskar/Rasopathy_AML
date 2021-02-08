# DEG analysis with edgeR

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
load("./RDS/deseq2_object.RData")
load("./RDS/raw_data.RData")

# Remove un-necessary data 
# (eg. Retain NRAS)

dds_nocontrol <- dds[,dds$Control == "0"]

meta_info <- meta_info[meta_info$Control == "0",]

dds_dge <- dds_nocontrol

# Construct edgeR object
se<-assay(dds_dge)
coldata<-meta_info[match(colnames(se),rownames(meta_info)), ]

library("edgeR")

genetable <- data.frame(gene.id=rownames(se))

y <- DGEList(counts=se, 
             samples=coldata, 
             genes=genetable)
names(y)

# DEsign formula same as deseq

design <- model.matrix(~ RASOPATHY, y$samples)

# Filter genes with low counts
keep <- filterByExpr(y, design)
table(keep)

# recompute library size
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalization

y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

# differential expression testing with edgeR

fit <- glmFit(y, design,robust=TRUE)
lrt <- glmLRT(fit)
tt <- topTags(lrt, n=nrow(y), p.value=0.05)[[1]]
tt.all <- topTags(lrt, n = nrow(y), sort.by = "none") # all genes

Sorted_res_diff<-as.data.frame(tt.all)

Sorted_res_diff<-getHGNC(Sorted_res_diff)

Sorted_res_diff$Gene<-rownames(Sorted_res_diff)

#write xlsx
library("openxlsx")

write.xlsx(Sorted_res_diff, file = paste("./DEG/DGE_edgeR",Sys.Date(), ".xlsx", sep="_"))

DGE_edgeR<-Sorted_res_diff

# Save multiple objects
save(DGE_edgeR, file = "./RDS/DGE_edgeR.RData")
