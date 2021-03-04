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

# Data normalization
log_cpm<-get.cpm(dds, meta_info)

# Subset meta data
meta_info<-subset(meta_info, select = -c(NPM1) )

meta_info[is.na(meta_info)] = 0

meta_info<-meta_info[meta_info$Control == "0",]
log_cpm<-log_cpm[,meta_info$ID]

#########################################
# Ras signatures

sigras <- readRDS("./RDS/sigras.rds")
meta_info$sigras <-as.numeric(sigras$calcsig(log_cpm, gene.names="ID"))

sigmpas <- readRDS("./RDS/sigmpas.rds")
meta_info$sigmpas <-as.numeric(sigmpas$calcsig(log_cpm, gene.names="ID"))

sigmerk <- readRDS("./RDS/sigerk.rds")

erk <-sigmerk$calcsig(log_cpm, gene.names="ID")

meta_info$sigmerk_IEG <-erk$hIEG

meta_info$sigmerk_ILG <-erk$hILG


#########################################
#Sigras

p1 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NRAS==0),],
                x = "NRAS", y = "sigras",
                color = "NRAS", palette = "jco", add = "jitter") +  ylim(4, 8)
p1<-p1 + stat_compare_means(label.x = 1.5, method = "t.test")

p2 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$KRAS==0),],
                x = "KRAS", y = "sigras",
                color = "KRAS", palette = "jco", add = "jitter") +  ylim(4, 8)
p2<-p2 + stat_compare_means(label.x = 1.5, method = "t.test")

p3 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$PTPN11==0),], 
                x = "PTPN11", y = "sigras",
                color = "PTPN11", palette = "jco", add = "jitter") +  ylim(4, 8)
p3<-p3 + stat_compare_means(label.x = 1.5, method = "t.test")

p4 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NF1==0),], 
                x = "NF1", y = "sigras",
                color = "NF1", palette = "jco", add = "jitter") +  ylim(4, 8)
p4<-p4 + stat_compare_means(label.x = 1.5, method = "t.test")

p5 <- ggboxplot(data=meta_info,
                x = "RASOPATHY", y = "sigras",
                color = "RASOPATHY", palette = "jco", add = "jitter") +  ylim(4, 8)
p5<-p5 + stat_compare_means(label.x = 1.5, method = "t.test")

library(multipanelfigure)

figure1 <- multi_panel_figure(
  width = 200, height = 200,
  columns = 3, rows = 2, row_spacing=5,
  column_spacing=5)

figure1 %<>% fill_panel(p1, column = 1, row = 1)
figure1 %<>% fill_panel(p2, column = 2, row = 1)
figure1 %<>% fill_panel(p3, column = 3, row = 1)
figure1 %<>% fill_panel(p4, column = 1, row = 2)
figure1 %<>% fill_panel(p5, column = 2, row = 2)

# Save figure
figure1 %>% save_multi_panel_figure(filename = "./images/sigras.png", dpi = 300)

#########################################################################
#Sigmpas

p1 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NRAS==0),],
                x = "NRAS", y = "sigmpas",
                color = "NRAS", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p1<-p1 + stat_compare_means(label.x = 1.5, method = "t.test")

p2 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$KRAS==0),],
                x = "KRAS", y = "sigmpas",
                color = "KRAS", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p2<-p2 + stat_compare_means(label.x = 1.5, method = "t.test")

p3 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$PTPN11==0),], 
                x = "PTPN11", y = "sigmpas",
                color = "PTPN11", palette = "jco", add = "jitter") + ylim(-0.04, 6) 
p3<-p3 + stat_compare_means(label.x = 1.5, method = "t.test")

p4 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NF1==0),], 
                x = "NF1", y = "sigmpas",
                color = "NF1", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p4<-p4 + stat_compare_means(label.x = 1.5, method = "t.test")

p5 <- ggboxplot(data=meta_info,
                x = "RASOPATHY", y = "sigmpas",
                color = "RASOPATHY", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p5<-p5 + stat_compare_means(label.x = 1.5, method = "t.test")

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
figure2 %>% save_multi_panel_figure(filename = "./images/sigmpas.png", dpi = 300)

###########################################################
#sigmerk_hIEG

p1 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NRAS==0),],
                x = "NRAS", y = "sigmerk_IEG",
                color = "NRAS", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p1<-p1 + stat_compare_means(label.x = 1.5, method = "t.test")

p2 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$KRAS==0),],
                x = "KRAS", y = "sigmerk_IEG",
                color = "KRAS", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p2<-p2 + stat_compare_means(label.x = 1.5, method = "t.test")

p3 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$PTPN11==0),], 
                x = "PTPN11", y = "sigmerk_IEG",
                color = "PTPN11", palette = "jco", add = "jitter") + ylim(-0.04, 6) 
p3<-p3 + stat_compare_means(label.x = 1.5, method = "t.test")

p4 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NF1==0),], 
                x = "NF1", y = "sigmerk_IEG",
                color = "NF1", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p4<-p4 + stat_compare_means(label.x = 1.5, method = "t.test")

p5 <- ggboxplot(data=meta_info,
                x = "RASOPATHY", y = "sigmerk_IEG",
                color = "RASOPATHY", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p5<-p5 + stat_compare_means(label.x = 1.5, method = "t.test")

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
figure2 %>% save_multi_panel_figure(filename = "./images/sigmerk_IEG.png", dpi = 300)

###################################

###########################################################
#sigmerk_hILG

p1 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NRAS==0),],
                x = "NRAS", y = "sigmerk_ILG",
                color = "NRAS", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p1<-p1 + stat_compare_means(label.x = 1.5, method = "t.test")

p2 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$KRAS==0),],
                x = "KRAS", y = "sigmerk_ILG",
                color = "KRAS", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p2<-p2 + stat_compare_means(label.x = 1.5, method = "t.test")

p3 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$PTPN11==0),], 
                x = "PTPN11", y = "sigmerk_ILG",
                color = "PTPN11", palette = "jco", add = "jitter") + ylim(-0.04, 6) 
p3<-p3 + stat_compare_means(label.x = 1.5, method = "t.test")

p4 <- ggboxplot(data=meta_info[!(meta_info$RASOPATHY==1 & meta_info$NF1==0),], 
                x = "NF1", y = "sigmerk_ILG",
                color = "NF1", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p4<-p4 + stat_compare_means(label.x = 1.5, method = "t.test")

p5 <- ggboxplot(data=meta_info,
                x = "RASOPATHY", y = "sigmerk_ILG",
                color = "RASOPATHY", palette = "jco", add = "jitter") + ylim(-0.04, 6)
p5<-p5 + stat_compare_means(label.x = 1.5, method = "t.test")

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
figure2 %>% save_multi_panel_figure(filename = "./images/sigmerk_ILG.png", dpi = 300)

