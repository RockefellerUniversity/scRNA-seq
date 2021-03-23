params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(bioMart))
suppressPackageStartupMessages(require(SeuratWrappers))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("Create Seurat object and Generare QC plots

---
"    
  )
  
}



## ----loadCR2Seurat,echo=TRUE,eval=FALSE,include=TRUE--------------------------
## tar_dir <- "path to raw data"


## ----loadCR2Seurat2,echo=TRUE,eval=FALSE,include=TRUE-------------------------
## library(Seurat)
## samID <- "CTRL"
## X10_file <- Seurat::Read10X(tar_dir)
## obj <- CreateSeuratObject(counts = X10_file, project = samID, min.cells = 5,min.features = 200)
## obj$dset <- samID
## obj <- Seurat::RenameCells(obj,add.cell.id=samID)
## #
## obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
## #
## cellID <- rownames(obj@meta.data)
## cellID_sel <- sample(cellID,1000,replace = FALSE)
## obj_sel <- obj[,cellID_sel]
## saveRDS(obj_sel,file.path("data",paste0("scSeq_",samID,"_1kCell_ori.rds")))

## ----loadCR2Seurat3,echo=TRUE,eval=TRUE,include=FALSE-------------------------
#obj <- readRDS("../data/scSeq_CTRL_1kCell_ori.rds")
obj <- readRDS("data/scSeq_CTRL_1kCell_ori.rds")


## ----objInfo,echo=TRUE,eval=TRUE,include=TRUE---------------------------------
obj


## ----objInfo_cell,echo=TRUE,eval=TRUE,include=TRUE----------------------------
head(obj@meta.data)


## ----objInfo_count,echo=TRUE,eval=TRUE,include=TRUE---------------------------
message("raw counts")
head(obj@assays$RNA@counts)
message("scaled data")
head(obj@assays$RNA@data)


## ----eval_readCount_geneCount_mitoCont,echo=TRUE,eval=TRUE,include=TRUE-------
# obj <- readRDS("../data/scSeq_CTRL_1kCell_ori.rds")
VlnPlot(obj,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0.2)


## ----eval_readCount_geneCount,echo=TRUE,eval=TRUE,include=TRUE----------------
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


## ----eval_readCount_mitoCont,echo=TRUE,eval=TRUE,include=TRUE-----------------
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")


## ----rmMito,echo=TRUE,eval=TRUE,include=TRUE----------------------------------
summary(obj@meta.data$percent.mt)
mt_cutH <- 10
obj <- subset(obj,subset = percent.mt < mt_cutH)


## ----rmMito_vln,echo=TRUE,eval=TRUE,include=TRUE------------------------------
VlnPlot(obj,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0.2)


## ----rmMito_scatter,echo=TRUE,eval=TRUE,include=TRUE--------------------------
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") 


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("Cell cycle phase determination

---
"    
  )
  
}



## ----cellCycle_est1,echo=TRUE,eval=TRUE,include=TRUE--------------------------
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = x , 
                   mart = human,
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)}
#
s_gene <- convertHumanGeneList(cc.genes$s.genes)
g2m_gene <- convertHumanGeneList(cc.genes$g2m.genes)


## ----cellCycle_est2,echo=TRUE,eval=TRUE,include=TRUE--------------------------
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = x , 
                   mart = human,
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)}
#
s_gene <- convertHumanGeneList(cc.genes$s.genes)
g2m_gene <- convertHumanGeneList(cc.genes$g2m.genes)
obj <- CellCycleScoring(obj, s.features = s_gene,
                              g2m.features = g2m_gene, set.ident = TRUE)
obj@meta.data[1:2,]


## ----cellCycle_plot1,echo=TRUE,eval=TRUE,include=TRUE-------------------------
yd_dat <- as.data.frame(table(obj@meta.data$dset,obj@meta.data$Phase))
head(yd_dat)


## ----cellCycle_plot2,echo=TRUE,eval=TRUE,include=TRUE-------------------------
library(ggplot2)
ggplot(yd_dat,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity",position="stack")+labs(x="",y="Counts",fill="Phase")+theme_classic()


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("Regression and clustering

---
"    
  )
  
}



## ----regress,echo=TRUE,eval=TRUE,include=TRUE---------------------------------
obj <- ScaleData(obj,
                 vars.to.regress = c("percent.mt","S.score","G2M.score","Phase"))


## ----prcomp,echo=TRUE,eval=TRUE,include=TRUE----------------------------------
set.seed(1000)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)


## ----pcSel,echo=TRUE,eval=TRUE,include=TRUE-----------------------------------
ElbowPlot(obj,ndims=30)
#
numPC <- 15


## ----clust,echo=TRUE,eval=TRUE,include=TRUE-----------------------------------
set.seed(1000)
maxPC <- numPC
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:maxPC)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:maxPC)
obj <- RunTSNE(obj, reduction = "pca", dims = 1:maxPC)


## ----seurat_UMAP,echo=TRUE,eval=TRUE,include=TRUE-----------------------------
DimPlot(obj,reduction="umap") # defaul it "umap"


## ----seurat_tSNE,echo=TRUE,eval=TRUE,include=TRUE-----------------------------
DimPlot(obj,reduction = "tsne")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("Identify marker genes for each cluster

---
"    
  )
  
}



## ----markGene,echo=TRUE,eval=TRUE,include=TRUE--------------------------------
obj <- SetIdent(obj,value = "seurat_clusters")
clust.markers <- FindAllMarkers(obj, 
                                only.pos = TRUE,
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
head(clust.markers)


## ----markGene_sub,echo=TRUE,eval=TRUE,include=TRUE----------------------------
topG <- clust.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
head(topG)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("Advanced plots

---
"    
  )
  
}



## ----dimPlot_ptSize,echo=TRUE,eval=TRUE,include=TRUE--------------------------
DimPlot(obj,pt.size = 0.2)


## ----dimPlot_label,echo=TRUE,eval=TRUE,include=TRUE---------------------------
DimPlot(obj,pt.size = 0.2,label=TRUE)+NoLegend()


## ----markGene_heatmap,echo=TRUE,eval=TRUE,include=TRUE------------------------
DoHeatmap(obj, features = topG$gene) + NoLegend()


## ----makGene_featurePlot,echo=TRUE,eval=TRUE,include=TRUE---------------------
gene_marker <- c("Krt1","Pthlh","Krt14","Cenpa","Shh")
FeaturePlot(obj,features = gene_marker,pt.size = 0.2)


## ----makGene_RidgePlott,echo=TRUE,eval=TRUE,include=TRUE----------------------
RidgePlot(obj,features = gene_marker)


## ----cellcycle_cluster,echo=TRUE,eval=TRUE,include=TRUE-----------------------
tbl <- table(obj@meta.data$seurat_clusters,obj@meta.data$Phase)
tbl_dat <- as.data.frame(tbl)
to <- rowSums(tbl)
names(to) <- rownames(tbl)
tbl_dat$to <- to[match(names(to),tbl_dat$Var1)]
tbl_dat$prop <- tbl_dat$Freq / tbl_dat$to
tbl_dat[1:2,]


## ----cellcycle_cluster_plot1,echo=TRUE,eval=FALSE,include=TRUE----------------
## ggplot(tbl_dat,aes(x=Var1,y=prop,fill=Var2))+
##   geom_bar(stat="identity",position="stack")+
##   labs(x="Seurat_clusters",y="Proportion",fill="Phase")+
##   theme_classic()


## ----cellcycle_cluster_plot2,echo=FALSE,eval=TRUE,include=TRUE----------------
ggplot(tbl_dat,aes(x=Var1,y=prop,fill=Var2))+
  geom_bar(stat="identity",position="stack")+
  labs(x="Seurat_clusters",y="Proportion",fill="Phase")+
  theme_classic()


## ----saveRDS_surat,echo=TRUE,eval=TRUE,include=TRUE---------------------------
cellID <- rownames(obj@meta.data)[!obj@meta.data$seurat_clusters==2]
obj_sub <- obj[,cellID]
saveRDS(obj_sub,"scSeq_Seurat_clean.rds")

