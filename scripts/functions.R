# Functions
getHGNC<-function(data)
{
  geneid<-as.data.frame(apply(geneid,2,function(x)gsub('\\s+', '',x)))
  
  # add gene symbol
  data <- data[which(rownames(data) %in% rownames(geneid)),]
  geneid_data <- geneid[which(rownames(geneid) %in% rownames(data)),]
  geneid_data <- geneid_data[!duplicated(geneid_data$hgnc_symbol),]

  data <- data[match(rownames(geneid_data), rownames(data)),]
  rownames(data) <- geneid_data$hgnc_symbol
  return(data)
}

getEntrz<-function(data)
{
  geneid<-as.data.frame(apply(geneid,2,function(x)gsub('\\s+', '',x)))
  
  # add gene symbol
  data <- data[which(rownames(data) %in% rownames(geneid)),]
  geneid_data <- geneid[which(rownames(geneid) %in% rownames(data)),]
  geneid_data <- geneid_data[!duplicated(geneid_data$entrezgene_id),]

  data <- data[match(rownames(geneid_data), rownames(data)),]
  rownames(data) <- geneid_data$entrezgene_id
  return(data)
}

##############
## Simple function to write .GMT (gene matrix transposed) format
## files.
## Get a list of genesets (ie list of character vectors)
## and write them to "filename"
write.gmt <- function(geneset.list, description.list=NA, filename)
{
  geneset.names <- names(geneset.list)
  
  if (is.na(description.list)) {
    write.set <- function(x) {
      write.table(t(c(geneset.names[[x]],"na",geneset.list[[x]])),
                  file=filename, col.names=FALSE, row.names=FALSE,
                  quote=FALSE, append=TRUE, sep="\t")
      cat("\n",file=filename,append=TRUE)
      paste("Geneset [",geneset.names[[x]],"] written to file", sep="")
    }
  } else {
    write.set <- function(x) {
      write.table(t(c(geneset.names[[x]],
                      description.list[[x]],
                      geneset.list[[x]])),
                  file=filename, col.names=FALSE, row.names=FALSE,
                  quote=FALSE, append=TRUE, sep="\t")
      cat("\n",file=filename,append=TRUE)
      paste("Geneset [",geneset.names[[x]],"] written to file", sep="")
    }
  }
  lapply(seq_along(geneset.list), write.set)
}

write.gct <- function(data, annot, filename)
{
  cat("#1.2\n",file=filename)
  cat(dim(data)[1:2],sep="\t",file=filename, append=TRUE)
  cat("\n",file=filename,append=TRUE)
  
  #outdf <- data.frame(t(data))
  outdf <- data.frame(data)
  
  outdf <- cbind("Description"=annot[dimnames(outdf)[[1]],]$RASOPATHY,outdf)
  outdf <- cbind("NAME"=rownames(outdf), outdf)
  rownames(outdf) <- NULL
  write.table(outdf,file=filename,append=TRUE,quote=FALSE,
              sep="\t",row.names=FALSE, na="")
}


call.gsea <- function(output.name="GSEA.analysis", iterations=1000)
{
  GSEA(
    ## IO
    input.ds = "aml.gct",
    input.cls = "aml.cls",
    gene.ann = "",
    gs.db = "aml.gmt",
    gs.ann = "",
    output.directory = "GSEA-output",
    doc.string = output.name,
    ## Base parameteres
    #non.interactive.run = F,
    reshuffling.type = "sample.labels",
    nperm = iterations,
    weighted.score.type = 1,
    nom.p.val.threshold = -1,
    fwer.p.val.threshold = -1,
    fdr.q.val.threshold = 0.25,
    topgs = 20,
    adjust.FDR.q.val = F,
    gs.size.threshold.min = 5,
    gs.size.threshold.max = 500,
    reverse.sign = F,
    preproc.type = 0,
    random.seed = 123456,
    ## Advanced parameters
    perm.type = 0,
    fraction = 1.0,
    replace = F,
    save.intermediate.results = F,
    #OLD.GSEA = F,
    use.fast.enrichment.routine = T)
}

###################
# heatmap function
###################

# input scale d matrix

get_heat<-function (scaled, chip)
{
  # Distance calculation based on correlation
  c <- cor((scaled), method="pearson")
  d <- as.dist(1-c)
  
  #Clustering of sample based on correlation distances
  hc <- hclust(d, method = "complete", members=NULL)
  
  #Distance calculation based on correlation
  c_g <- cor(t(scaled), method="pearson")
  d_g <- as.dist(1-c_g)
  
  #Clustering of sample based on correlation distances
  hr <- hclust(d_g, method = "complete", members=NULL)
  
  #Heatmap color
  
  library(circlize)
  library(ComplexHeatmap)
  
  mycols <-  colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
  
  cm = ColorMapping(name = "",
                    col_fun = colorRamp2(breaks=c(-1.5, 0, 1.5), color=c("blue", "white", "red")))
  
  grid.newpage()
  color_mapping_legend(cm, title_gp = gpar(fontsize = 16))
  
  # Annotation data frame
  
  #reorder metadata
  meta_info<-sample_info[match(colnames(scaled),rownames(sample_info)), ]
  
  annotation = data.frame(TKI=meta_info$TKI, Mutation=meta_info$Mutation, Rep=meta_info$Rep)
  
  ha = HeatmapAnnotation(df = annotation,
                         col = list(TKI=c("AC220" = "brown", "DMSO" = "violet",
                                          "PKC" = "blue", "PONA" = "green"),
                                    Mutation=c("D835V" = "orange", "D835Y" = "red",
                                               "F69L" = "gold", "ITD" = "gray"),
                                    Rep=c("R1"="black", "R2"="gray")
                                    
                         ),
                         height = unit(2, "cm"),
                         simple_anno_size_adjust = TRUE,
                         show_legend = TRUE
  )
  
  png(paste("./images/heatmap_", chip, "_", Sys.Date(), ".png", sep = ""), width = 20, height = 25,
      units = 'in', res = 300)
  p_deg<- Heatmap(scaled,
                  width = unit(18, "cm"), height = unit(25, "cm"),
                  show_row_names = TRUE, show_row_dend = TRUE,
                  show_column_names = TRUE,
                  cluster_columns = hc,
                  cluster_rows = hr,
                  heatmap_legend_param = list(title = "z-score",
                                              color_bar = "continuous", labels_gp = gpar(fontsize = 6)),
                  col=mycols, column_names_gp = gpar(fontsize = 8),
                  row_names_gp = gpar(fontsize = 4),
                  row_title="proteins",
                  top_annotation = ha)
  draw(p_deg, annotation_legend_side = "bottom")
  dev.off()
  
}

#####################
#Redraw heatmap with nclust order
get_heat_nclust<-function (scaled, chip)
{
  
  #Heatmap color
  
  library(circlize)
  library(ComplexHeatmap)
  
  mycols <-  colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
  
  cm = ColorMapping(name = "",
                    col_fun = colorRamp2(breaks=c(-1.5, 0, 1.5), color=c("blue", "white", "red")))
  
  grid.newpage()
  color_mapping_legend(cm, title_gp = gpar(fontsize = 16))
  
  # Annotation data frame
  #reorder metadata
  meta_info<-sample_info[match(colnames(scaled),rownames(sample_info)), ]
  
  annotation = data.frame(TKI=meta_info$TKI, Mutation=meta_info$Mutation, Rep=meta_info$Rep)
  
  ha = HeatmapAnnotation(df = annotation,
                         col = list(TKI=c("AC220" = "brown", "DMSO" = "violet",
                                          "PKC" = "blue", "PONA" = "green"),
                                    Mutation=c("D835V" = "orange", "D835Y" = "red",
                                               "F69L" = "gold", "ITD" = "gray"),
                                    Rep=c("R1"="black", "R2"="gray")
                                    
                         ),
                         height = unit(2, "cm"),
                         simple_anno_size_adjust = TRUE,
                         show_legend = TRUE
  )
  
  png(paste("./images/heatmap_nclust_", chip, "_", Sys.Date(), ".png", sep = ""), width = 20, height = 25,
      units = 'in', res = 300)
  p_deg<- Heatmap(scaled,
                  width = unit(18, "cm"), height = unit(25, "cm"),
                  show_row_names = TRUE, show_row_dend = TRUE,
                  show_column_names = TRUE,
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  heatmap_legend_param = list(title = "z-score",
                                              color_bar = "continuous", labels_gp = gpar(fontsize = 6)),
                  col=mycols, column_names_gp = gpar(fontsize = 8),
                  row_names_gp = gpar(fontsize = 4),
                  row_title="proteins"
                  #,
                  #top_annotation = ha
  )
  draw(p_deg, annotation_legend_side = "bottom")
  dev.off()
  
}

#Function to clustering genes 
ICAclust<-function(data=data, cluster_method="ward.D", mojena=3.75, times = c(20,40,70,90))
{            
  X=data
  Y=X[,3:ncol(X)]
  X.ICA<-fastICA(Y,n.comp=ncol(Y), alg.typ = "parallel", fun = "exp", alpha = 1.0,
                 method = "C", row.norm = FALSE, maxit = 5000, tol = 1e-03, verbose = TRUE)
  hc.ICA<-hclust(dist(X.ICA$S), method = cluster_method, members=NULL)
  dendro<-as.dendrogram(hc.ICA)
  plot(hc.ICA,labels = NULL, hang = 0.0,
       axes = TRUE, frame.plot = FALSE, ann = TRUE,
       main = "Cluster Dendrogram",
       sub = NULL, xlab = "", ylab = "Height")
  hc.ICA$order
  mojema=mean(hc.ICA$height)+mojena*sd(hc.ICA$height)
  g=length(hc.ICA$height[hc.ICA$height>mojema]) + 1
  rect.hclust(hc.ICA, k =g)
  group_number<-g
  group<-cutree(hc.ICA, k=g)
  number_gene<-table(group)
  grupos1<-as.data.frame(X$geneID)
  comp<-cbind(group,grupos1)
  dadosnovos <-cbind(group,X[,-c(1,2)])
  dadosnovoss <-cbind(group,X[,2],X[,-c(1,2)])
  names(dadosnovoss)<-c("group",names(X)[-1])
  dadosnovos1 <-as.data.frame(dadosnovos)
  a1=by(dadosnovos1, dadosnovos1[,"group"],function(x) colMeans(x))
  x1<-do.call("rbind",as.list(a1))
  write.csv(dadosnovoss, file = "groups.csv")
  write.table(dadosnovoss, file = "groups.txt")
  write.csv(x1, file = "averagevalues.csv")
  write.table(x1, file = "averagevalues.txt")
  result<-list( "nGroups" = group_number, "gbGroup"= number_gene)
  return(result)
}

get.cpm <- function(dds, meta_info)
{
  
  # Construct edgeR object
  
  se<-assay(dds)
  
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
  
  # recompute library siz
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # normalization
  lcpm <- cpm(y, log = TRUE, prior.count = 0.25)
  
  topMatrix <- lcpm
  
  return(topMatrix)
}

  

