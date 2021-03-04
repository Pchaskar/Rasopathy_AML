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

#########################################
# Ras signatures

#sigras <- readRDS("./RDS/sigras.rds")
#meta_info$sigras <-as.numeric(sigras$calcsig(log_cpm, gene.names="ID"))

#sigmpas <- readRDS("./RDS/sigmpas.rds")
#meta_info$sigmpas <-as.numeric(sigmpas$calcsig(log_cpm, gene.names="ID"))

sigmerk <- readRDS("./RDS/sigerk.rds")

erk <-sigmerk$calcsig(log_cpm, gene.names="ID")

meta_info$sigmerk_hIEG <-erk$hIEG

meta_info$sigmerk_hILG <-erk$hILG

###########################################################
#sigmerk_hIEG

p1 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$NRAS=="No_raso"),],
                x = "NRAS", y = "sigmerk_hIEG",
                color = "NRAS", palette = "jco", add = "jitter") + ylim(2, 7)
p1<-p1 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p2 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$KRAS=="No_raso"),],
                x = "KRAS", y = "sigmerk_hIEG",
                color = "KRAS", palette = "jco", add = "jitter") + ylim(2, 7)
p2<-p2 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p3 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$PTPN11=="No_raso"),], 
                x = "PTPN11", y = "sigmerk_hIEG",
                color = "PTPN11", palette = "jco", add = "jitter") + ylim(2, 7) 
p3<-p3 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p4 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$NF1=="No_raso"),], 
                x = "NF1", y = "sigmerk_hIEG",
                color = "NF1", palette = "jco", add = "jitter") + ylim(2, 7)
p4<-p4 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p5 <- ggboxplot(data=meta_info,
                x = "RASOPATHY", y = "sigmerk_hIEG",
                color = "RASOPATHY", palette = "jco", add = "jitter") + ylim(2, 7)
p5<-p5 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

library(multipanelfigure)

figure2 <- multi_panel_figure(
  width = 200, height = 200,
  columns = 3, rows = 2, row_spacing=5,
  column_spacing=5)

figure2 %<>% fill_panel(p1, column = 1, row = 1)
figure2 %<>% fill_panel(p2, column = 2, row = 1)
figure2 %<>% fill_panel(p3, column = 3, row = 1)
figure2 %<>% fill_panel(p4, column = 1, row = 2)
figure2 %<>% fill_panel(p5, column = 2, row = 2)


# Save figure
figure2 %>% save_multi_panel_figure(filename = "./images/sigmerk_hIEG.png", dpi = 300)

###########################################################
#sigmerk_hILG

p1 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$NRAS=="No_raso"),],
                x = "NRAS", y = "sigmerk_hILG",
                color = "NRAS", palette = "jco", add = "jitter") + ylim(2, 6)
p1<-p1 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p2 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$KRAS=="No_raso"),],
                x = "KRAS", y = "sigmerk_hILG",
                color = "KRAS", palette = "jco", add = "jitter") + ylim(2, 6)
p2<-p2 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p3 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$PTPN11=="No_raso"),], 
                x = "PTPN11", y = "sigmerk_hILG",
                color = "PTPN11", palette = "jco", add = "jitter") + ylim(2, 6) 
p3<-p3 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p4 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY=="Raso" & meta_info$NF1=="No_raso"),], 
                x = "NF1", y = "sigmerk_hILG",
                color = "NF1", palette = "jco", add = "jitter") + ylim(2, 6)
p4<-p4 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

p5 <- ggboxplot(data=meta_info,
                x = "RASOPATHY", y = "sigmerk_hILG",
                color = "RASOPATHY", palette = "jco", add = "jitter") + ylim(2, 6)
p5<-p5 + stat_compare_means(label.x = 1.2, method = "t.test")+ theme(legend.position = "none")

library(multipanelfigure)

figure2 <- multi_panel_figure(
  width = 200, height = 200,
  columns = 3, rows = 2, row_spacing=5,
  column_spacing=5)

figure2 %<>% fill_panel(p1, column = 1, row = 1)
figure2 %<>% fill_panel(p2, column = 2, row = 1)
figure2 %<>% fill_panel(p3, column = 3, row = 1)
figure2 %<>% fill_panel(p4, column = 1, row = 2)
figure2 %<>% fill_panel(p5, column = 2, row = 2)

# Save figure
figure2 %>% save_multi_panel_figure(filename = "./images/sigmerk_hILG.png", dpi = 300)

#ERK heatmap
meta_info_erk<-subset(meta_info, select = c(sigmerk_hIEG, sigmerk_hILG))

svg("./images/Sig_erk.svg", width = 6, height = 8)
heatmap.2(as.matrix(meta_info_erk))
dev.off()

