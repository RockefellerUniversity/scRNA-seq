params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(DropletUtils))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(scuttle))
suppressPackageStartupMessages(require(scDblFinder))
knitr::opts_chunk$set(echo = TRUE, tidy = T)



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Reticulate and python

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("#  Reticulate and python

---
"    
  )
  
}



## ----installPythonAndPkgsA,include=FALSE--------------------------------------

if(!require(reticulate)){
  install.packages("reticulate")
  require(reticulate)
}





## ----installPythonAndPkgsB,eval=FALSE-----------------------------------------
## install.packages("reticulate")
## library(reticulate)
## 
## 
## 


## ----installPythonAndPkgsC,include=FALSE--------------------------------------

if(!dir.exists(reticulate::miniconda_path())){
  install_miniconda()
}





## ----installPythonAndPkgsD,eval=FALSE-----------------------------------------
## 
## install_miniconda()
## 
## 


## ----installPythonAndPkgsE,include=FALSE--------------------------------------

if(!dir.exists(reticulate::miniconda_path())){
  install_miniconda()
}
conda_install(
  packages=c("python-igraph","scanpy","louvain","leidenalg","loompy")
)





## ----installPythonAndPkgsF,eval=FALSE-----------------------------------------
## 
## conda_install(
##   packages=c("python-igraph","scanpy","louvain","leidenalg","loompy")
## )
## 
## 
## 


## ----usePython----------------------------------------------------------------
use_condaenv()
py_config()



## ----configPython-------------------------------------------------------------
use_condaenv()
py_config()



## ----importPython-------------------------------------------------------------
sc <- reticulate::import("scanpy")
sc


## ----installAnnData_1,include=FALSE-------------------------------------------
if(!require(anndata)){
  install.packages("anndata")
  require(anndata)
}



## ----installAnnData_2,eval=FALSE----------------------------------------------
## install.packages("anndata")
## library(anndata)
## 


## ----downloadPBMC-------------------------------------------------------------
pathToDownload <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(pathToDownload,"pbmc3k_filtered_gene_bc_matrices.tar.gz")
utils::untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")


## ----readPBMC-----------------------------------------------------------------
adata <- sc$read_10x_mtx("filtered_gene_bc_matrices/hg19/",
                         var_names="gene_ids") 

adata


## ----dimensions---------------------------------------------------------------
dim(adata)


## ----dimnames-----------------------------------------------------------------
colnames(adata)[1:3]
rownames(adata)[1:3]


## ----countsMatrix-------------------------------------------------------------
countsMatrix <- adata$X
countsMatrix[1:2,1:2]


## ----vardf--------------------------------------------------------------------
var_df <- adata$var
var_df[1:2,,drop=FALSE]


## ----obsdf--------------------------------------------------------------------
obs_df <- adata$obs
obs_df
rownames(obs_df)[1:2]


## ----addsomeQC----------------------------------------------------------------
sc$pp$filter_cells(adata, min_genes=0)
sc$pp$filter_genes(adata, min_cells=0)


## ----vardf2-------------------------------------------------------------------
var_df <- adata$var
var_df[1:2,]


## ----obsdf2-------------------------------------------------------------------
obs_df <- adata$obs
obs_df[1:2,,drop=FALSE]


## ----vardfRe------------------------------------------------------------------
var_df$mito <- grepl("MT-",var_df$gene_symbols)
var_df[var_df$mito == TRUE,]


## ----vardfRe2-----------------------------------------------------------------
adata$var <- var_df
adata$var[1:2,]


## ----mitoQC,tidy=FALSE--------------------------------------------------------
sc$pp$calculate_qc_metrics(adata,
                           qc_vars = list("mito"),
                           inplace = TRUE)
adata


## ----moreQC,tidy=FALSE--------------------------------------------------------
adata$var[5:10,]


## ----moreQC2,tidy=FALSE-------------------------------------------------------
adata$obs[1:2,]


## ----plotHiExpr_1,include=FALSE-----------------------------------------------
sc$settings$figdir = getwd()
sc$pl$highest_expr_genes(adata,gene_symbols = "gene_symbols",save=".png")


## ----plotHiExpr_2,eval=FALSE--------------------------------------------------
## sc$pl$highest_expr_genes(adata,gene_symbols = "gene_symbols")


## ----plotQCpy2_1,include=FALSE------------------------------------------------
sc$pl$violin(adata, list('n_genes_by_counts', 'total_counts', 'pct_counts_mito'),
             jitter=0.4, multi_panel=TRUE,save=".png")


## ----plotQCpy2_2,eval=FALSE---------------------------------------------------
## sc$pl$violin(adata, list('n_genes_by_counts', 'total_counts', 'pct_counts_mito'),
##              jitter=0.4, multi_panel=TRUE)


## ----plotQCpy3_1,include=FALSE------------------------------------------------
sc$pl$scatter(adata, x='total_counts', y='n_genes_by_counts',save=".png")



## ----plotQCpy3_2,eval=FALSE---------------------------------------------------
## sc$pl$scatter(adata, x='total_counts', y='n_genes_by_counts')
## 


## ----filtercells--------------------------------------------------------------
adata = adata[adata$obs$n_genes < 2500]
adata = adata[adata$obs$n_genes > 200]
adata = adata[adata$obs$pct_counts_mito < 5]
adata


## ----filtergenes--------------------------------------------------------------
sc$pp$filter_genes(adata,
                   min_cells=3)
adata


## ----normCells----------------------------------------------------------------
sc$pp$normalize_per_cell(adata, counts_per_cell_after=10000)


## ----logCells-----------------------------------------------------------------
# log transform the data.
sc$pp$log1p(adata)



## ----variableGenes------------------------------------------------------------
# identify highly variable genes.
sc$pp$highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# sc$pl$highly_variable_genes(adata)


## ----filterToVariable_2,eval=FALSE--------------------------------------------
## adata = adata[,adata$var$highly_variable]
## 


## ----filterToVariable_1,include=FALSE-----------------------------------------
# keep only highly variable genes:
adata = adata$copy()
adata = adata[,adata$var$highly_variable]
adata = adata$copy()
# regress out total counts per cell and the percentage of mitochondrial genes expressed
# adata.copy = adata
# sc$pp$regress_out(adata, list('n_counts', 'pct_counts_mito'),copy=FALSE) #, n_jobs=args.threads)


## ----regress------------------------------------------------------------------
sc$pp$regress_out(adata, list('n_counts', 'pct_counts_mito'))


## ----scale--------------------------------------------------------------------
sc$pp$scale(adata, max_value=10)


## ----pca----------------------------------------------------------------------
sc$tl$pca(adata, svd_solver='arpack')


## ----neighbours---------------------------------------------------------------
sc$pp$neighbors(adata, n_neighbors=10L, n_pcs=40L)
sc$tl$umap(adata)


## ----louvain------------------------------------------------------------------
sc$tl$louvain(adata)
sc$tl$leiden(adata)


## ----plotumapCells_1,include=FALSE--------------------------------------------
sc$pl$umap(adata, color=list('leiden'),save=".png")


## ----plotumapCells_2,eval=FALSE-----------------------------------------------
## sc$pl$umap(adata, color=list('leiden'))


## ----saveToh5ad---------------------------------------------------------------
adata$write_h5ad(filename = "PBMC_Scanpy.h5ad")


## ----saveToLoom---------------------------------------------------------------
adata$write_loom(filename = "PBMC_Scanpy.loom")


## ----readLoomb,eval=FALSE-----------------------------------------------------
## 
## BiocManager::install("LoomExperiment")
## require(LoomExperiment)
## 
## loom <- import(con = "PBMC_Scanpy.loom")
## loom


## ----readLooma,include=FALSE--------------------------------------------------
if(!require(LoomExperiment)){
  BiocManager::install("LoomExperiment")
  require(LoomExperiment)
}
loom <- LoomExperiment::import(con = "PBMC_Scanpy.loom")
loom


## ----loomToSCE----------------------------------------------------------------
sce_loom <- as(loom,"SingleCellExperiment")
sce_loom


## ----replaceUMAP_a,eval=FALSE-------------------------------------------------
## library(scater)
## reducedDim(sce_loom, "UMAP") <- adata$obsm[["X_umap"]]
## reducedDim(sce_loom, "PCA") <- adata$obsm[["X_pca"]]


## ----replaceUMAP_b,include=FALSE----------------------------------------------
if(!require(scater)){
  BiocManager::install("scater")
  require(scater)
}
reducedDim(sce_loom, "UMAP") <- adata$obsm[["X_umap"]]
reducedDim(sce_loom, "PCA") <- adata$obsm[["X_pca"]]


## ----plotUMAP,fig.width=5,fig.height=5----------------------------------------
plotUMAP(sce_loom,colour_by="leiden")

