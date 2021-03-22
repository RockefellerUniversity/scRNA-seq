params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(phaetmap))
suppressPackageStartupMessages(require(dplyr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# ?? name of course?? (part 2)

---
"    
  )
  
}



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Analysis with Seurat

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Set Up

---
"    
  )
  
}



## ----loadCR2Seurat------------------------------------------------------------
tar_dir <- "~/Documents/packDev_ALL/scRNA-seq/scRNASeq/inst/extdata/data/CTRL"
# tar_dir <- "data/CTRL"
dir(tar_dir)


## ----loadCR2Seurat2,eval=FALSE------------------------------------------------
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


## ----eval_readCount_geneCount_mitoCont,dpi=300--------------------------------
obj <- readRDS("data/scSeq_CTRL_1kCell_ori.rds")
# obj <- readRDS("../data/scSeq_CTRL_1kCell_ori.rds")
VlnPlot(obj,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0.2)
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
#
summary(obj@meta.data$percent.mt)
mt_cutH <- 10
obj <- subset(obj,subset = percent.mt < mt_cutH)

