params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Bioconductor methods

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Bioconductor methods

---
"    
  )
  
}



## ----fig.width=7,fig.height=4-------------------------------------------------

library(DropletUtils)



## ----fig.width=7,fig.height=4-------------------------------------------------

h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/NeuroD1_filtered_feature_bc_matrix.h5"
local_h5file <- basename(h5file)
download.file(h5file,local_h5file)

sce.NeuroD1_filtered <- read10xCounts(local_h5file, col.names=TRUE,row.names = "symbol")

class(sce.NeuroD1_filtered)



## ----include=FALSE------------------------------------------------------------
unlink(h5file)


## ----fig.width=7,fig.height=4-------------------------------------------------
sce.NeuroD1_filtered


## ----loadSCE_pres222,include=TRUE,echo=TRUE,eval=TRUE-------------------------
# cell information
colnames(sce.NeuroD1_filtered)[1:2]
# gene information
rownames(sce.NeuroD1_filtered)[1:2]


## ----loadSCE_pres12,include=TRUE,echo=TRUE,eval=TRUE--------------------------
# cell information
colData(sce.NeuroD1_filtered)[1:2,]
# gene information
rowData(sce.NeuroD1_filtered)[1:2,]


## ----loadSCE_pres2,include=TRUE,echo=TRUE,eval=TRUE---------------------------
# cell information
reducedDimNames(sce.NeuroD1_filtered)
# gene information
metadata(sce.NeuroD1_filtered)


## ----loadSCE_pres3,include=TRUE,echo=TRUE,eval=TRUE---------------------------
h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/raw_feature_bc_matrix.h5"
local_h5file <-"NeuroD1_raw_feature_bc_matrix.h5"
download.file(h5file,local_h5file)
sce.NeuroD1_unfiltered <- read10xCounts(local_h5file, col.names=TRUE)
sce.NeuroD1_unfiltered


## ----include=FALSE,eval=FALSE-------------------------------------------------
## unlink(local_h5file)


## ----fig.width=7,fig.height=4-------------------------------------------------
bcrank <- barcodeRanks(counts(sce.NeuroD1_unfiltered))
bcrank


## ----fig.width=7,fig.height=4-------------------------------------------------
bcrank$filtered <- rownames(bcrank) %in% colnames(sce.NeuroD1_filtered)
bc_plot <- as.data.frame(bcrank)
bc_plot <- bc_plot[order(bc_plot$filtered,decreasing=TRUE),]
bc_plot <- bc_plot[!duplicated(bc_plot$rank),]


## ----fig.width=7,fig.height=4-------------------------------------------------
require(ggplot2)
  ggplot(bc_plot,aes(x=rank,y=total,colour=filtered,alpha=0.001))+
  geom_point()+
  scale_y_log10()+
  scale_x_log10()+
  theme_minimal()+
  geom_hline(yintercept = metadata(bcrank)$inflection,colour="darkgreen",linetype=2)+
  geom_hline(yintercept = metadata(bcrank)$knee,colour="dodgerblue",linetype=2)



## ----eval=FALSE,echo=TRUE-----------------------------------------------------
## e.out <- emptyDrops(counts(sce.NeuroD1_unfiltered))
## table(e.out)
## head(e.out)


## ----eval=TRUE,echo=FALSE-----------------------------------------------------
e.out <- readRDS(system.file("extdata/data/e.out.RData", mustWork = TRUE, package = "scRNASeq"))
head(e.out)


## ----include=FALSE,eval=FALSE-------------------------------------------------
## unlink(local_h5file)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Bioconductor  -  QC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Bioconductor  -  QC

---
"    
  )
  
}



## ----fig.width=7,fig.height=4-------------------------------------------------
sce.NeuroD1_filtered


## ----fig.width=7,fig.height=4-------------------------------------------------
is.mito <- grepl("^mt",rowData(sce.NeuroD1_filtered)$Symbol)
is.ribo <- grepl("^Rps",rowData(sce.NeuroD1_filtered)$Symbol)

table(is.mito)
table(is.ribo)



## ----fig.width=7,fig.height=4-------------------------------------------------
library(scuttle)
sce.NeuroD1_filtered <- addPerCellQCMetrics(sce.NeuroD1_filtered, 
                                            subsets=list(Mito=is.mito,Ribosomal=is.ribo))
sce.NeuroD1_filtered


## ----fig.width=7,fig.height=4-------------------------------------------------
qc_df <- colData(sce.NeuroD1_filtered)
head(qc_df,n=2)


## ----fig.width=7,fig.height=4-------------------------------------------------
library(scater)
plotColData(sce.NeuroD1_filtered,x="Sample", y="sum")


## ----fig.width=7,fig.height=4-------------------------------------------------
library(scater)
plotColData(sce.NeuroD1_filtered,x="Sample",y="detected")



## ----fig.width=7,fig.height=4-------------------------------------------------
library(scater)
p1 <- plotColData(sce.NeuroD1_filtered,x="Sample",y="subsets_Mito_percent")
p2 <- plotColData(sce.NeuroD1_filtered,x="Sample",y="subsets_Ribosomal_percent")
gridExtra::grid.arrange(p1,p2,ncol=2)


## ----fig.width=7,fig.height=4-------------------------------------------------
library(scater)
p1 <- plotColData(sce.NeuroD1_filtered,x="sum",y="detected")
p2 <- plotColData(sce.NeuroD1_filtered,x="sum",y="subsets_Mito_percent")
p3 <- plotColData(sce.NeuroD1_filtered,x="detected",y="subsets_Mito_percent")
gridExtra::grid.arrange(p1,p2,p3,ncol=3)


## ----fig.width=7,fig.height=4-------------------------------------------------
qc.high_lib_size <- colData(sce.NeuroD1_filtered)$sum > 125000
qc.min_detected <- colData(sce.NeuroD1_filtered)$detected < 200
qc.mito <- colData(sce.NeuroD1_filtered)$subsets_Mito_percent > 25
discard <- qc.high_lib_size | qc.mito | qc.min_detected
DataFrame(LibSize=sum(qc.high_lib_size), Detected=sum(qc.min_detected),MitoProp=sum(qc.mito), Total=sum(discard))


## ----fig.width=7,fig.height=4-------------------------------------------------
colData(sce.NeuroD1_filtered) <- cbind(colData(sce.NeuroD1_filtered),DataFrame(toDiscard=discard))
p1 <- plotColData(sce.NeuroD1_filtered,x="sum",y="detected",colour_by = "toDiscard")
p2 <- plotColData(sce.NeuroD1_filtered,x="sum",y="subsets_Mito_percent",colour_by = "toDiscard")
p3 <- plotColData(sce.NeuroD1_filtered,x="detected",y="subsets_Mito_percent",colour_by = "toDiscard")
gridExtra::grid.arrange(p1,p2,p3,ncol=3)


## ----fig.width=7,fig.height=4-------------------------------------------------
sce.NeuroD1_filtered_QCed <- sce.NeuroD1_filtered[,sce.NeuroD1_filtered$toDiscard %in% "FALSE"]
sce.NeuroD1_filtered_QCed


## ----fig.width=7,fig.height=4-------------------------------------------------
sce.NeuroD1_filtered_QCed <- logNormCounts(sce.NeuroD1_filtered_QCed)
assayNames(sce.NeuroD1_filtered_QCed)


## ----fig.width=7,fig.height=4-------------------------------------------------
library(scran)
clust.sce.NeuroD1_filtered_QCed <- quickCluster(sce.NeuroD1_filtered_QCed) 
sce.NeuroD1_filtered_trimmed <- computeSumFactors(sce.NeuroD1_filtered_QCed, cluster=clust.sce.NeuroD1_filtered_QCed)
sce.NeuroD1_filtered_QCed <- logNormCounts(sce.NeuroD1_filtered_QCed)
assayNames(sce.NeuroD1_filtered_QCed)



## ----fig.width=7,fig.height=4-------------------------------------------------
##
dec.NeuroD1_filtered_QCed <- modelGeneVar(sce.NeuroD1_filtered_QCed)
ggplot(dec.NeuroD1_filtered_QCed,aes(x=mean,y=total))+geom_point()


## ----fig.width=7,fig.height=4-------------------------------------------------
top.NeuroD1_filtered_QCed <- getTopHVGs(dec.NeuroD1_filtered_QCed, n=3000)



## ----fig.width=7,fig.height=4-------------------------------------------------
set.seed(100) # See below.
sce.NeuroD1_filtered_QCed <- fixedPCA(sce.NeuroD1_filtered_QCed,
                                      subset.row=top.NeuroD1_filtered_QCed) 
reducedDimNames(sce.NeuroD1_filtered_QCed)



## ----fig.width=7,fig.height=4-------------------------------------------------
plotReducedDim(sce.NeuroD1_filtered_QCed, dimred="PCA",colour_by = "subsets_Mito_percent")



## ----fig.width=7,fig.height=4-------------------------------------------------
percent.var <- attr(reducedDim(sce.NeuroD1_filtered_QCed), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
# library(PCAtools)
# chosen.elbow <- findElbowPoint(percent.var)


## ----fig.width=7,fig.height=4-------------------------------------------------

sce.NeuroD1_filtered_QCed <- runTSNE(sce.NeuroD1_filtered_QCed,n_dimred=30)
sce.NeuroD1_filtered_QCed <- runUMAP(sce.NeuroD1_filtered_QCed,n_dimred=30)
reducedDimNames(sce.NeuroD1_filtered_QCed)


## ----fig.width=7,fig.height=4-------------------------------------------------
plotUMAP(sce.NeuroD1_filtered_QCed,colour_by="subsets_Mito_percent")
plotTSNE(sce.NeuroD1_filtered_QCed,colour_by="subsets_Mito_percent")


## ----fig.width=7,fig.height=4-------------------------------------------------
require(bluster)
clust.louvain <- clusterCells(sce.NeuroD1_filtered_QCed, use.dimred="PCA", 
                                BLUSPARAM=NNGraphParam(cluster.fun="louvain",cluster.args = list(resolution=0.8)))
clust.default <- clusterCells(sce.NeuroD1_filtered_QCed, use.dimred="PCA")



## ----fig.width=7,fig.height=4-------------------------------------------------
colLabels(sce.NeuroD1_filtered_QCed) <- clust.louvain
colData(sce.NeuroD1_filtered_QCed)$DefaultLabel <- clust.default
plotUMAP(sce.NeuroD1_filtered_QCed,colour_by="label")



## ----fig.width=7,fig.height=4-------------------------------------------------
plotColData(sce.NeuroD1_filtered_QCed,x="label",y="subsets_Mito_percent",colour_by = "label")


## ----fig.width=7,fig.height=4-------------------------------------------------
tab <- table(Walktrap=clust.default, Louvain=clust.louvain)
rownames(tab) <- paste("Walktrap", rownames(tab))
colnames(tab) <- paste("louvain", colnames(tab))
library(pheatmap)
pheatmap(log10(tab+10), color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)


## ----eval=TRUE----------------------------------------------------------------
eec_paper_meta <- read.delim("https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/GSE224223_EEC_metadata.csv",sep=";")
nrd1_cells <- as.data.frame(eec_paper_meta[eec_paper_meta$Library == "Nd1",])
   
nrd1_cells$X <- gsub("_2$","",nrd1_cells$X)

newMeta <- merge(colData(sce.NeuroD1_filtered_QCed),nrd1_cells,by.x="Barcode",by.y="X",all.x=TRUE,all.y=FALSE,order=FALSE)
rownames(newMeta) <- newMeta$Barcode
newMeta$All_Cell_Types[is.na(newMeta$All_Cell_Types)] <- "Filtered"
colData(sce.NeuroD1_filtered_QCed) <- newMeta
plotUMAP(sce.NeuroD1_filtered_QCed,colour_by="All_Cell_Types")



## ----fig.width=7,fig.height=4-------------------------------------------------
tab <- table(Walktrap=sce.NeuroD1_filtered_QCed$All_Cell_Types, Louvain=sce.NeuroD1_filtered_QCed$label)
rownames(tab) <- paste("CellType", rownames(tab))
colnames(tab) <- paste("Leiden", colnames(tab))
library(pheatmap)
pheatmap(log10(tab+10), color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)


## ----fig.width=7,fig.height=4-------------------------------------------------
markers.sce.NeuroD1 <- scoreMarkers(sce.NeuroD1_filtered_QCed,sce.NeuroD1_filtered_QCed$label)
chosen <- markers.sce.NeuroD1[["4"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4]) # showing basic stats only, for brevity.


## ----fig.width=7,fig.height=4-------------------------------------------------
plotExpression(sce.NeuroD1_filtered_QCed, features=head(rownames(ordered),n=1), 
    x="label", colour_by="label")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Seurat

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Seurat

---
"    
  )
  
}



## ----fig.width=7,fig.height=4-------------------------------------------------
library(Seurat)
h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/NeuroD1_filtered_feature_bc_matrix.h5"
local_h5file <- basename(h5file)
download.file(h5file,local_h5file)
Nd1T.mat <- Read10X_h5(filename = local_h5file)
Nd1T.obj <- CreateSeuratObject(Nd1T.mat, project = "Nd1")



## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj


## ----fig.width=7,fig.height=4-------------------------------------------------
rownames(Nd1T.obj)[1:2]
colnames(Nd1T.obj)[1:2]


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj@meta.data[1:2,]


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj@active.assay
Nd1T.obj@active.ident[1:2]


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj@assays$RNA$counts[1:2,1:2]


## ----fig.width=7,fig.height=4-------------------------------------------------
mito.genes <- grep("^mt", rownames(Nd1T.obj), value = T)
mito.genes[1:10]


## ----fig.width=7,fig.height=4-------------------------------------------------
percent.mt <- PercentageFeatureSet(Nd1T.obj, features = mito.genes)
Nd1T.obj <- AddMetaData(Nd1T.obj,metadata = percent.mt,col.name = "percent.mt")
head(Nd1T.obj@meta.data)


## ----fig.width=7,fig.height=4-------------------------------------------------
VlnPlot(Nd1T.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


## ----fig.width=7,fig.height=4-------------------------------------------------
plot1 <- FeatureScatter(Nd1T.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Nd1T.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Nd1T.obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1+plot2+plot3


## ----fig.width=7,fig.height=4-------------------------------------------------

Nd1T.obj.filt <- subset(Nd1T.obj, subset = `nCount_RNA` < 125000 & percent.mt < 25 & nFeature_RNA > 200) 
Nd1T.obj.filt


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj.filt <- NormalizeData(Nd1T.obj.filt)


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj.filt <- FindVariableFeatures(Nd1T.obj.filt, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(Nd1T.obj.filt), 10)
top10


## ----fig.width=7,fig.height=4-------------------------------------------------
plot1 <- VariableFeaturePlot(Nd1T.obj.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


## ----fig.width=7,fig.height=4-------------------------------------------------
all.genes <- rownames(Nd1T.obj.filt)
Nd1T.obj.filt <- ScaleData(Nd1T.obj.filt, features = all.genes)
Nd1T.obj.filt@assays


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj.filt <- RunPCA(Nd1T.obj.filt, features = VariableFeatures(object = Nd1T.obj))
ElbowPlot(Nd1T.obj.filt,ndims = 50)


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj.filt <- RunUMAP(Nd1T.obj.filt, dims = 1:30)


## ----fig.width=7,fig.height=4-------------------------------------------------
FeaturePlot(Nd1T.obj.filt, reduction = "umap",features = "percent.mt")


## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.obj.filt <- FindNeighbors(Nd1T.obj.filt, dims = 1:30)
Nd1T.obj.filt <- FindClusters(Nd1T.obj.filt, resolution = 0.8)
Nd1T.obj.filt@active.ident


## ----fig.width=7,fig.height=4-------------------------------------------------
DimPlot(Nd1T.obj.filt, reduction = "umap")


## ----fig.width=7,fig.height=4-------------------------------------------------
bioc.clusters <- data.frame(Bioc=sce.NeuroD1_filtered_QCed$label,row.names = sce.NeuroD1_filtered_QCed$Barcode)
seurat.clusters <- as.data.frame(Nd1T.obj.filt@active.ident)
colnames(seurat.clusters) <- "Seurat"
clustersAll <- merge(bioc.clusters,seurat.clusters,by=0,all=FALSE)



## ----fig.width=7,fig.height=4-------------------------------------------------
tab <- table(Bioc=clustersAll$Bioc, Seurat=clustersAll$Seurat)
rownames(tab) <- paste("Bioc", rownames(tab))
colnames(tab) <- paste("Seurat", colnames(tab))
library(pheatmap)
pheatmap(log10(tab+10), color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)



## ----fig.width=7,fig.height=4-------------------------------------------------
Nd1T.markers <- FindAllMarkers(Nd1T.obj.filt, only.pos = TRUE)
Nd1T.markers[Nd1T.markers$cluster == 0,]


## ----fig.width=7,fig.height=4-------------------------------------------------
VlnPlot(Nd1T.obj.filt, features = head(Nd1T.markers[Nd1T.markers$cluster == 0,"gene"],n=1))


