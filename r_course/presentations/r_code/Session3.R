params <-
list(isSlides = "no")

## ----include=FALSE-------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
library(bluster)
library(decontX)
library(harmony)
library(scCustomize)
library(scater)
library(scran)
library(scDblFinder)
library(scuttle)
library(DropletUtils)
library(celda)
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------

h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/filtered_feature_bc_matrix.h5"
local_h5file <- basename(h5file)
download.file(h5file,local_h5file)

sce.NeuroD1_filtered <- read10xCounts(local_h5file, col.names=TRUE)
sce.NeuroD1_filtered



## ----loadSCE_pres3,include=TRUE,echo=TRUE,eval=TRUE----------------------------
h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/raw_feature_bc_matrix.h5"
local_h5file <-"NeuroD1_raw_feature_bc_matrix.h5"
download.file(h5file,local_h5file)
sce.NeuroD1_unfiltered <- read10xCounts(local_h5file, col.names=TRUE)
sce.NeuroD1_unfiltered


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(decontX)


## ----include=FALSE,eval=FALSE--------------------------------------------------
## sce.NeuroD1_filtered <- decontX(sce.NeuroD1_filtered, background = sce.NeuroD1_unfiltered)
## saveRDS(sce.NeuroD1_filtered,file = "decontx.RDS")
## sce.NeuroD1_filtered


## ----echo=FALSE,eval=TRUE------------------------------------------------------
rdsfile <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/decontx.RDS"
local_rdsfile <-"decontx.RDS"
download.file(rdsfile,local_rdsfile)
sce.NeuroD1_filtered <- readRDS(local_rdsfile)
sce.NeuroD1_filtered


## ----echo=FALSE,eval=TRUE------------------------------------------------------
assays(sce.NeuroD1_filtered)[["decontXcounts"]][1:2,1:4]


## ----echo=FALSE,eval=TRUE------------------------------------------------------
colData(sce.NeuroD1_filtered)


## ----echo=FALSE,eval=TRUE------------------------------------------------------
reducedDim(sce.NeuroD1_filtered)[1:2,]


## ----echo=FALSE,eval=TRUE------------------------------------------------------
umap <- reducedDim(sce.NeuroD1_filtered, "decontX_UMAP")
plotDimReduceCluster(x = sce.NeuroD1_filtered$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])


## ----echo=FALSE,eval=TRUE------------------------------------------------------
plotDecontXContamination(sce.NeuroD1_filtered)


## ----echo=FALSE,eval=TRUE------------------------------------------------------
rdsfile <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/sce.NeuroD1_filtered_QCed.RDS"
local_rdsfile <-"sce.NeuroD1_filtered_QCed.RDS"
download.file(rdsfile,local_rdsfile)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.NeuroD1_filtered_QCed <- readRDS("sce.NeuroD1_filtered_QCed.RDS")
sce.NeuroD1_filtered_QCed$decontX_contamination <- sce.NeuroD1_filtered$decontX_contamination[match(colnames(sce.NeuroD1_filtered_QCed),colnames(sce.NeuroD1_filtered),nomatch = 0)]


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotColData(sce.NeuroD1_filtered_QCed,x = "label",y="decontX_contamination")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.NeuroD1_filtered_QCed,colour_by ="label")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.NeuroD1_filtered_QCed,colour_by ="decontX_contamination")


## ----fig.width=7,fig.height=4--------------------------------------------------

h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/cellbender_v0.3.2_filtered.h5"
local_h5file <- basename(h5file)
download.file(h5file,local_h5file)

sce.NeuroD1_CellBender <- read10xCounts(local_h5file, col.names=TRUE,row.names = "symbol")
sce.NeuroD1_CellBender$Sample <- "Nd1"
class(sce.NeuroD1_CellBender)



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.NeuroD1_CellBender


## ----fig.width=7,fig.height=4--------------------------------------------------
is.mito <- grepl("^MT",rowData(sce.NeuroD1_CellBender)$Symbol)
sce.NeuroD1_CellBender <- addPerCellQCMetrics(sce.NeuroD1_CellBender, 
                                            subsets=list(Mito=is.mito))
p1 <- plotColData(sce.NeuroD1_CellBender,x="sum",y="detected")
p2 <- plotColData(sce.NeuroD1_CellBender,x="sum",y="subsets_Mito_percent")
p3 <- plotColData(sce.NeuroD1_CellBender,x="detected",y="subsets_Mito_percent")
gridExtra::grid.arrange(p1,p2,p3,ncol=3)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
qc.high_lib_size <- colData(sce.NeuroD1_CellBender)$sum > 125000
qc.min_detected <- colData(sce.NeuroD1_CellBender)$detected < 200
qc.mito <- colData(sce.NeuroD1_CellBender)$subsets_Mito_percent > 25
discard <- qc.high_lib_size | qc.mito | qc.min_detected
colData(sce.NeuroD1_CellBender) <- cbind(colData(sce.NeuroD1_CellBender),DataFrame(toDiscard=discard))
sce.NeuroD1_CellBender_QCed <- sce.NeuroD1_CellBender[,sce.NeuroD1_CellBender$toDiscard %in% "FALSE"]


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
set.seed(100)
clust.sce.NeuroD1_CellBender_QCed <- quickCluster(sce.NeuroD1_CellBender_QCed) 
sce.NeuroD1_CellBender_trimmed <- computeSumFactors(sce.NeuroD1_CellBender_QCed, cluster=clust.sce.NeuroD1_CellBender_QCed)
sce.NeuroD1_CellBender_QCed <- logNormCounts(sce.NeuroD1_CellBender_QCed)
dec.NeuroD1_CellBender_QCed <- modelGeneVar(sce.NeuroD1_CellBender_QCed)
top.NeuroD1_CellBender_QCed <- getTopHVGs(dec.NeuroD1_CellBender_QCed, n=3000)
sce.NeuroD1_CellBender_QCed <- fixedPCA(sce.NeuroD1_CellBender_QCed,
                                      subset.row=top.NeuroD1_CellBender_QCed) 
sce.NeuroD1_CellBender_QCed <- runTSNE(sce.NeuroD1_CellBender_QCed,n_dimred=30)
sce.NeuroD1_CellBender_QCed <- runUMAP(sce.NeuroD1_CellBender_QCed,n_dimred=30)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(bluster)
clust.louvain <- clusterCells(sce.NeuroD1_CellBender_QCed, use.dimred="PCA", 
                                BLUSPARAM=NNGraphParam(cluster.fun="louvain",
                                                       cluster.args = list(resolution=0.8)))
clust.default <- clusterCells(sce.NeuroD1_CellBender_QCed, use.dimred="PCA")
colLabels(sce.NeuroD1_CellBender_QCed) <- clust.louvain
colData(sce.NeuroD1_CellBender_QCed)$DefaultLabel <- clust.default


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.NeuroD1_CellBender_QCed,colour_by="label")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
eec_paper_meta <- read.delim("https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/GSE224223_EEC_metadata.csv",sep=";")
nrd1_cells <- as.data.frame(eec_paper_meta[eec_paper_meta$Library == "Nd1",])

nrd1_cells$X <- gsub("_2$","",nrd1_cells$X)

sce.NeuroD1_CellBender_QCed$All_Cell_Types <- nrd1_cells$All_Cell_Types[match(colnames(sce.NeuroD1_CellBender_QCed),nrd1_cells$X,nomatch = NA)]


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.NeuroD1_CellBender_QCed,colour_by="All_Cell_Types")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(scDblFinder)
sce.NeuroD1_CellBender_QCed <- scDblFinder(sce.NeuroD1_CellBender_QCed, clusters=colLabels(sce.NeuroD1_CellBender_QCed))


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotColData(sce.NeuroD1_CellBender_QCed,x = "label",y="scDblFinder.score")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
p1 <- plotColData(sce.NeuroD1_CellBender_QCed,x = "detected",y="scDblFinder.score")
p2 <- plotColData(sce.NeuroD1_CellBender_QCed,x = "sum",y="scDblFinder.score")
gridExtra::grid.arrange(p1,p2,ncol=2)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
p1 <- plotUMAP(sce.NeuroD1_CellBender_QCed,colour_by="label")
p2 <- plotUMAP(sce.NeuroD1_CellBender_QCed,colour_by="scDblFinder.score")
gridExtra::grid.arrange(p1,p2,ncol=2)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.NeuroD1_CellBender_QCDbled <- sce.NeuroD1_CellBender_QCed[sce.NeuroD1_CellBender_QCed$scDblFinder.class != "doublet"]
sce.NeuroD1_CellBender_QCDbled


## ----include=FALSE-------------------------------------------------------------


rdsfile <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_Ngn3/outs/sce.Ngn3_CellBender_QCDbled.RDS"
local_rdsfile <-"sce.Ngn3_CellBender_QCDbled.RDS"
download.file(rdsfile,local_rdsfile)
sce.Ngn3_CellBender_QCDbled <- readRDS(local_rdsfile)
sce.Ngn3_CellBender_QCDbled$Sample <- "Ngn3"
unlink(local_rdsfile)

rdsfile <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_Ngn3/outs/dec.Ngn3_CellBender_QCed.RDS"
local_rdsfile <-"dec.Ngn3_CellBender_QCed.RDS"
download.file(rdsfile,local_rdsfile)
dec.Ngn3_CellBender_QCed <- readRDS(local_rdsfile)
unlink(local_rdsfile)



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.Ngn3_CellBender_QCDbled,colour_by="All_Cell_Types")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
uni <- intersect(rownames(sce.NeuroD1_CellBender_QCDbled),
rownames(sce.Ngn3_CellBender_QCDbled))

sce.NeuroD1 <- sce.NeuroD1_CellBender_QCDbled[uni,]
sce.Ngn3 <- sce.Ngn3_CellBender_QCDbled[uni,]



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------

mcols(sce.NeuroD1) <- mcols(sce.NeuroD1)[,-4]
mcols(sce.Ngn3) <- mcols(sce.Ngn3)[,-4]



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.all <- cbind(sce.NeuroD1,sce.Ngn3)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
dec.NeuroD1 <- dec.NeuroD1_CellBender_QCed[uni,]
dec.Ngn3 <- dec.Ngn3_CellBender_QCed[uni,]
combined.dec <- combineVar(dec.NeuroD1, dec.Ngn3)
chosen.hvgs <- combined.dec$bio > 0


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(scater)
set.seed(0010101010)
sce.all <- runPCA(sce.all, subset_row=chosen.hvgs)
sce.all <- runTSNE(sce.all, dimred="PCA")


## ----fig.height=4, warning=FALSE, r,fig.width=7--------------------------------
plotTSNE(sce.all, colour_by="Sample")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotTSNE(sce.all, colour_by="Sample")+facet_wrap(~colour_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotTSNE(sce.all, colour_by="All_Cell_Types",shape_by="Sample")+facet_wrap(~shape_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotTSNE(sce.all, colour_by="Sample",shape_by="All_Cell_Types")+facet_wrap(~shape_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(batchelor)
set.seed(1000101001)
sce.all.mnn <- fastMNN(sce.all,batch = sce.all$Sample,d=50, k=20, subset.row=chosen.hvgs)



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
eec_paper_meta <- read.delim("https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/GSE224223_EEC_metadata.csv",sep=";")
nrd1_cells <- as.data.frame(eec_paper_meta)

nrd1_cells$X <- gsub("_2$|_1$","",nrd1_cells$X)

sce.all.mnn$All_Cell_Types <- nrd1_cells$All_Cell_Types[match(colnames(sce.all.mnn),nrd1_cells$X,nomatch = NA)]


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(bluster)
clust.louvain <- clusterCells(sce.all.mnn, use.dimred="corrected", 
                                BLUSPARAM=NNGraphParam(cluster.fun="louvain",
                                                       cluster.args = list(resolution=0.5)))
colData(sce.all.mnn)$MNNLabel <- clust.louvain


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.all.mnn <- runTSNE(sce.all.mnn, dimred="corrected",name="MNN.TSNE")
sce.all.mnn <- runUMAP(sce.all.mnn, dimred="corrected",name="MNN.UMAP")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.all.mnn,
         colour_by="All_Cell_Types",
         shape_by="batch",dimred="MNN.UMAP")+facet_wrap(~shape_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.all.mnn,
         colour_by="MNNLabel",
         shape_by="batch",dimred="MNN.UMAP")+facet_wrap(~shape_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.all.mnn,
         colour_by="batch",
         shape_by="All_Cell_Types",dimred="MNN.UMAP")+facet_wrap(~shape_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(harmony)
sce.all.harmony <- RunHarmony(sce.all,group.by.vars="Sample")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.all.harmony$All_Cell_Types <- nrd1_cells$All_Cell_Types[match(colnames(sce.all.harmony),nrd1_cells$X,nomatch = NA)]


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(bluster)
clust.louvain <- clusterCells(sce.all.harmony, use.dimred="HARMONY", 
                                BLUSPARAM=NNGraphParam(cluster.fun="louvain",
                                                       cluster.args = list(resolution=0.5)))
colData(sce.all.harmony)$HARMONYLabel <- clust.louvain


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotReducedDim(sce.all.harmony,dimred = "HARMONY",colour_by = "subsets_Mito_percent")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------

sce.all.harmony <- runUMAP(sce.all.harmony, dimred="HARMONY",name="HARMONY.UMAP")

plotUMAP(sce.all.harmony,colour_by="HARMONYLabel",dimred = "HARMONY.UMAP",
         shape_by="Sample")+facet_wrap(~shape_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------

plotUMAP(sce.all.harmony,colour_by="All_Cell_Types",
         ,dimred = "HARMONY.UMAP",
         shape_by="Sample")+facet_wrap(~shape_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
library(Seurat)
library(scCustomize)
h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/cellbender_v0.3.2_filtered.h5"
local_h5file <- basename(h5file)
download.file(h5file,local_h5file)
Neurod1.mat <- Read_CellBender_h5_Mat(file_name = local_h5file)
Neurod1.obj <- CreateSeuratObject(Neurod1.mat, project = "Nd1")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
mito.genes <- grep("^MT", rownames(Neurod1.obj), value = T)

percent.mt <- PercentageFeatureSet(Neurod1.obj, features = mito.genes)
Neurod1.obj <- AddMetaData(Neurod1.obj,metadata = percent.mt,col.name = "percent.mt")
Neurod1.obj.filt <- subset(Neurod1.obj, subset = `nCount_RNA` < 125000 & percent.mt < 25 & nFeature_RNA > 200)





## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
Neurod1.obj.filt <- NormalizeData(Neurod1.obj.filt)
Neurod1.obj.filt <- FindVariableFeatures(Neurod1.obj.filt, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Neurod1.obj.filt)
Neurod1.obj.filt <- ScaleData(Neurod1.obj.filt, features = all.genes)
Neurod1.obj.filt <- RunPCA(Neurod1.obj.filt, features = VariableFeatures(object = Neurod1.obj))
Neurod1.obj.filt <- RunUMAP(Neurod1.obj.filt, dims = 1:30)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
Neurod1.obj.filt <- FindNeighbors(Neurod1.obj.filt, dims = 1:30)
Neurod1.obj.filt <- FindClusters(Neurod1.obj.filt, resolution = 0.7)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
DimPlot(Neurod1.obj.filt, reduction = "umap")


## ----include=FALSE-------------------------------------------------------------

options(timeout = 60*10000)
rdsfile <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_Ngn3/outs/Ngn3.obj.filt.RDS"
local_rdsfile <-"Ngn3.obj.filt.RDS"
download.file(rdsfile,local_rdsfile)
Ngn3.obj.filt <- readRDS(local_rdsfile)
unlink(local_rdsfile)




## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
DimPlot(Ngn3.obj.filt, reduction = "umap")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj <- merge(Neurod1.obj.filt, y = Ngn3.obj.filt, add.cell.ids = c("Nd1", "Ngn3"), project = "All",merge.data = TRUE)
All.obj


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj <- NormalizeData(All.obj)
All.obj <- FindVariableFeatures(All.obj)
All.obj <- ScaleData(All.obj)
All.obj <- RunPCA(All.obj)
All.obj <- FindNeighbors(All.obj, dims = 1:30, reduction = "pca")
All.obj <- FindClusters(All.obj, resolution = 2, cluster.name = "unintegrated_clusters")
All.obj <- RunUMAP(All.obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
DimPlot(All.obj, reduction = "umap.unintegrated",group.by = "orig.ident")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj <- IntegrateLayers(
  object = All.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = TRUE
)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj[["RNA"]] <- JoinLayers(All.obj[["RNA"]])
All.obj <- FindNeighbors(All.obj, reduction = "integrated.harmony", dims = 1:30)
All.obj <- FindClusters(All.obj, resolution = 0.5, cluster.name = "harmony_clusters")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj <- RunUMAP(All.obj, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
DimPlot(
  All.obj,
  reduction = "umap.harmony",
  group.by = "harmony_clusters",
  combine = FALSE, label.size = 2
)



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
DimPlot(
  All.obj,
  reduction = "umap.harmony",
  group.by = "orig.ident",
  combine = FALSE, label.size = 2
)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotUMAP(sce.all.mnn,colour_by="MNNLabel",dimred = "MNN.UMAP")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.all$MNNLabel <- sce.all.mnn$MNNLabel
sce.all$batch <- sce.all.mnn$batch


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
m.out <- findMarkers(sce.all, sce.all$MNNLabel, block=sce.all$batch,
    direction="up", lfc=1)



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
m.out[["1"]][,c("summary.logFC", "Top", "p.value", "FDR")]


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotExpression(sce.all,features = c("Chgb","Reg4"),x = "MNNLabel",colour_by = "batch")+facet_wrap(Feature~colour_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
plotExpression(sce.all.mnn,features = c("Chgb","Reg4"),x = "MNNLabel",exprs_values = "reconstructed",colour_by = "batch")+facet_wrap(Feature~colour_by)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
sce.all$Condition_Cluster <- paste(sce.all$Sample,sce.all$MNNLabel,sep="_")
Nd1_vs_Ngn3_Cluster9 <- scoreMarkers(sce.all,sce.all$Condition_Cluster,pairings=c("Nd1_9","Ngn3_9"))


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
Nd1_vs_Ngn3_Cluster9_Res <- Nd1_vs_Ngn3_Cluster9$Nd1_9
Nd1_vs_Ngn3_Cluster9_Res <- Nd1_vs_Ngn3_Cluster9_Res[order(Nd1_vs_Ngn3_Cluster9_Res$mean.AUC, decreasing=TRUE),]
plotExpression(sce.all,
               features = c("Xist"),
               x = "MNNLabel",colour_by = "batch")+facet_wrap(Feature~colour_by)



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
Idents(All.obj) <- "harmony_clusters"
harmony.markers <- FindConservedMarkers(All.obj, ident.1 = 1, grouping.var = "orig.ident", verbose = FALSE)
VlnPlot(All.obj,features = "Slc38a11",split.by = "orig.ident")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj[["Condition_Cluster"]] <- paste(All.obj$orig.ident,All.obj$harmony_clusters,sep="_")
Idents(All.obj) <- "Condition_Cluster"
res <- FindMarkers(All.obj, ident.1 = "Nd1_1", ident.2="Ngn3_1", verbose = FALSE)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
Neurod1.obj.scTransform <- subset(Neurod1.obj, subset = `nCount_RNA` < 125000 & percent.mt < 25 & nFeature_RNA > 200)
Neurod1.obj.scTransform <- SCTransform(Neurod1.obj.scTransform,vars.to.regress = "percent.mt")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
Neurod1.obj.scTransform


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
Neurod1.obj.scTransform <- RunPCA(Neurod1.obj.scTransform)
Neurod1.obj.scTransform <- RunUMAP(Neurod1.obj.scTransform, dims = 1:30)
Neurod1.obj.scTransform <- FindNeighbors(Neurod1.obj.scTransform, dims = 1:30)
Neurod1.obj.scTransform <- FindClusters(Neurod1.obj.scTransform, resolution = 0.7)



## ----echo=FALSE,eval=TRUE------------------------------------------------------
rdsfile <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_Ngn3/outs/Ngn3.obj.scTranform.RDS"
local_rdsfile <-"Ngn3.obj.scTranform.RDS"
download.file(rdsfile,local_rdsfile)
Ngn3.obj.scTransform <- readRDS(local_rdsfile)
Ngn3.obj.scTransform


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj.sct <- merge(Neurod1.obj.scTransform, 
                 y = Ngn3.obj.scTransform, add.cell.ids = c("Nd1", "Ngn3"), project = "All",merge.data = TRUE)
All.obj.sct <- SCTransform(All.obj.sct,
                           vars.to.regress = "percent.mt")


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj.sct <- RunPCA(All.obj.sct)
All.obj.sct <- RunUMAP(All.obj.sct, dims = 1:30)
All.obj.sct <- IntegrateLayers(object = All.obj.sct, 
                        method = HarmonyIntegration, 
                        normalization.method = "SCT", 
                        verbose = F)



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj.sct <- FindNeighbors(All.obj.sct, reduction = "harmony", dims = 1:30)
All.obj.sct <- FindClusters(All.obj.sct, resolution = 0.3)
All.obj.sct <- RunUMAP(All.obj.sct, dims = 1:30, reduction = "harmony")

DimPlot(All.obj.sct, reduction = "umap", group.by = c("orig.ident","SCT_snn_res.0.3"))



## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
All.obj.sct <- PrepSCTFindMarkers(All.obj.sct)


## ----fig.width=7,fig.height=4,warning=FALSE------------------------------------
harmonySCT.markers <- FindMarkers(All.obj.sct, ident.1 = 3, verbose = FALSE)
VlnPlot(All.obj.sct,features = "Tph1",split.by = "orig.ident")


