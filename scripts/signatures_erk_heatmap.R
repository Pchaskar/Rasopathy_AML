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
library(ggpubr)

# load functions
source("./scripts/functions.R")

# Load DESEQ2 Data
load("./RDS/raw_data.RData")
load("./RDS/deseq2_object.RData")

# remove duplicate naems
meta_info<-subset(meta_info, select = -c(NPM1) )
meta_info<-meta_info[meta_info$Control == "0",]
meta_info[is.na(meta_info)] = 0


# format meta data
library(car)

recode2 <- function ( data, fields, recodes, as.factor.result = FALSE ) {
  for ( i in which(names(data) %in% fields) ) { # iterate over column indexes that are present in the passed dataframe that are also included in the fields list
    data[,i] <- car::recode( data[,i], recodes)
  }
  data
}

meta_info<-recode2(meta_info, fields = c('NRAS', 'KRAS', 'NF1','RASOPATHY','PTPN11'), 
             recodes = "'1' = 'Raso'; '0' = 'No_raso'", as.factor.result = TRUE)

# Data normalization
log_cpm<-get.cpm(dds, meta_info)

# Subset 
log_cpm<-log_cpm[,meta_info$ID]

log_cpm<-getHGNC(log_cpm)
#########################################
# Ras signatures

sigmerk <- readRDS("./RDS/sigerk.rds")

erk<-c(sigmerk$siglist.hgnc$hIEG, sigmerk$siglist.hgnc$hILG)

log_cpm<-log_cpm[rownames(log_cpm) %in% erk, ]

log_cpm<-log_cpm[match(erk, rownames(log_cpm)), ]

###########################################
# Scale and center
log_cpm_scaled<-t(scale(t(log_cpm), center = TRUE, scale = TRUE))

##################
#Heatmap color

library(circlize)

mycols <-  colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

cm = ColorMapping(name = "",
                  col_fun = colorRamp2(breaks=c(-1.5, 0, 1.5), color=c("blue", "white", "red")))

grid.newpage()
color_mapping_legend(cm, title_gp = gpar(fontsize = 16))

##################
# Annotation data frame

#reorder metadata
meta_info<- arrange(meta_info, RASOPATHY)

log_cpm_scaled<-log_cpm_scaled[ ,match(rownames(meta_info), colnames(log_cpm_scaled))]

annotation = data.frame(RASOPATHY=meta_info$RASOPATHY)

ha = HeatmapAnnotation(df = annotation,
                       col = list(RASOPATHY=c("No_raso" = "brown", "Raso" = "violet")),
                       na_col = "white",
                       height = unit(0.75, "cm"),
                       simple_anno_size_adjust = TRUE,
                       show_legend = TRUE
)

hIEG = which(rownames(log_cpm_scaled) %in% sigmerk$siglist.hgnc$hIEG)
hILG = which(rownames(log_cpm_scaled) %in% sigmerk$siglist.hgnc$hILG)

hIEG = rowAnnotation(link = anno_mark(hIEG, 
                                          labels = rownames(log_cpm_scaled[hIEG,]), 
                                          labels_gp = gpar(fontsize = 10, col="red")))

hILG = rowAnnotation(link = anno_mark(hILG, 
                                            labels = rownames(log_cpm_scaled[hILG,]), 
                                            labels_gp = gpar(fontsize = 10, col="blue")))


tiff(paste("./images/ERK_hIEG_hILG",Sys.Date(), "heatmap.tiff", sep = "_"), width = 20, height = 22, 
     units = 'in', res = 600, compression =  "lzw+p")
p_deg<- Heatmap(log_cpm_scaled,
                width = unit(18, "cm"), height = unit(20, "cm"),
                show_row_names = FALSE, show_row_dend = FALSE,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                heatmap_legend_param = list(title = "z-score", 
                                            color_bar = "continuous", labels_gp = gpar(fontsize = 6)),
                col=mycols, column_names_gp = gpar(fontsize = 10), 
                row_names_gp = gpar(fontsize = 8),
                row_title="genes",
                top_annotation = ha)+hIEG+hILG
draw(p_deg, annotation_legend_side = "bottom")
dev.off()
