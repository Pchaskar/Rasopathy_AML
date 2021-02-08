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
load("./RDS/logCPM_tmm_normalized_data.RData")

#########################################
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

#########################################

# Boxplot Rasopathy T test

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
figure1 %>% save_multi_panel_figure(filename = "./images/Box_Signatures_sigras.png", dpi = 300)
#########################################################################

# Boxplot Rasopathy T test
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
figure2 %>% save_multi_panel_figure(filename = "./images/Box_Signatures_sigmpas.png", dpi = 300)


# For Petros
annotation<-subset(meta_info, select = c(ID, NRAS, KRAS, NF1, PTPN11, RASOPATHY, sigras, sigmpas) )

# Save multiple objects
save(annotation, file = "./RDS/signature_annotation.RData")

