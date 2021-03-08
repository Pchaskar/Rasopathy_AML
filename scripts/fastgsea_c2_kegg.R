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
#C2 Msigdb
c2<-readRDS("./gsea/c2list.rds")

c2BroadSets<-c2$c2list.symbol

kegg <- c2BroadSets[c(grep("^K_", names(c2BroadSets)))]

pathways <- kegg


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
  labs(x="Gene Sets", y="Normalized Enrichment Score",
       title="C2: kegg gene sets NES from fgsea") + 
  theme_classic()

svg("./images/fgsea_c2_kegg.svg", width = 10, height = 10)
plot_gsea
dev.off()

# Write xlsx file
library("openxlsx")

write.xlsx(fgseaResTidy,row.names = TRUE, file = paste("./gsea/c2_kegg","fgsea.xlsx", sep="_"))
