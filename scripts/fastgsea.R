#Fast GSEA analysis setup

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

# Load DEGs
load("./RDS/DGE_edgeR.RData")

#FGSEA, preraked gsea analysis
library(fgsea)

#GSEA Preranked
er<- DGE_edgeR

er$fcsign <- sign(er$logFC)
er$logP=-log10(er$PValue)
er$metric= er$logP/er$fcsign

er <- er[complete.cases(er), ]

er<-er[ order(-er[,"metric"]), ]

ranks<-er$metric
names(ranks) <- rownames(er)

#write.table(ranks,file="./gsea/GSEA_pre-ranked.rnk",quote=F,sep="\t",row.names=T)

#plot(ranks)
barplot(ranks)

#Load pathways:
#C6 Msigdb
pathways <- gmtPathways("./gsea/c6.all.v7.2.symbols.gmt")

# Additional signatures

sigras <- readRDS("./RDS/sigras.rds")
sigmpas <- readRDS("./RDS/sigmpas.rds")
sigmerk <- readRDS("./RDS/sigerk.rds")
oxphos<-read.table(file = "./RDS/Signature_JE_Sarry.txt",
                   header = FALSE, check.names = FALSE)

oxphos <- oxphos[!duplicated(oxphos$V1),]


pathways$JE_Sarry_OXPHOS<-oxphos
pathways$Sigras<-names(sigras$sig.hgnc)
pathways$Sigmpas<-names(sigmpas$sig.hgnc)
pathways$Sigerk_IEG<-sigmerk$siglist.hgnc$IEG
pathways$Sigerk_ILG<-sigmerk$siglist.hgnc$ILG
pathways$Sigerk_DEG<-sigmerk$siglist.hgnc$DEG

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nperm=1000)
class(fgseaRes)
dim(fgseaRes)
head(fgseaRes)  

#Plot results
library(dplyr)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  tail()


fgseaResTidy %>% top_n(20,NES)  # highest values
fgseaResTidy %>% top_n(-20,NES) # lowest values

fgseaResTidy_top<-bind_rows(fgseaResTidy %>% top_n(20,NES), 
                            fgseaResTidy %>% top_n(-20,NES))

pathway_names<-fgseaResTidy_top$pathway
pathway_names_split <- strsplit(pathway_names, "%")

fgseaResTidy_top$pathway <- sapply(pathway_names_split, "[", 1)


# only plot the top 20 pathways
plot_gsea<-ggplot(fgseaResTidy_top, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj<0.1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from fgsea") + 
  theme_classic()

svg("./images/fgsea.svg", width = 10, height = 10)
plot_gsea
dev.off()

svg("./images/fgsea_KRAS.600.LUNG.BREAST_UP.V1_DN.svg", width = 6, height = 6)
plotEnrichment(pathways[["KRAS.600.LUNG.BREAST_UP.V1_DN"]],
               ranks) + labs(title="KRAS.600.LUNG.BREAST_UP.V1_DN")
dev.off()

svg("./images/fgsea_Sigerk_IEG.svg", width = 6, height = 6)
plotEnrichment(pathways[["Sigerk_IEG"]],
               ranks) + labs(title="Sigerk_IEG")
dev.off()

svg("./images/fgsea_Sigerk_ILG.svg", width = 6, height = 6)
plotEnrichment(pathways[["Sigerk_ILG"]],
               ranks) + labs(title="Sigerk_ILG")
dev.off()

svg("./images/fgsea_Sigerk_DEG.svg", width = 6, height = 6)
plotEnrichment(pathways[["Sigerk_DEG"]],
               ranks) + labs(title="Sigerk_DEG")
dev.off()


svg("./images/fgsea_JE_Sarry_OXPHOS.svg", width = 6, height = 6)
plotEnrichment(pathways[["JE_Sarry_OXPHOS"]],
               ranks) + labs(title="JE_Sarry_OXPHOS")
dev.off()

svg("./images/fgsea_Sigras.svg", width = 6, height = 6)
plotEnrichment(pathways[["Sigras"]],
               ranks) + labs(title="Sigras")
dev.off()

svg("./images/fgsea_Sigmpas.svg", width = 6, height = 6)
plotEnrichment(pathways[["Sigmpas"]],
               ranks) + labs(title="Sigmpas")
dev.off()
