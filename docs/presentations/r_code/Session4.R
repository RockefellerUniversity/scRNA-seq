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

# Loading data and empty droplets

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("#  Loading data and empty drops

---
"    
  )
  
}



## ----installPythonAndPkgs-----------------------------------------------------

if(!require(reticulate)){
  install.packages("reticulate")
  require(reticulate)
}
if(!dir.exists(reticulate::miniconda_path())){
  install_miniconda()
}
conda_install(
  packages=c("python-igraph","scanpy","louvain")
)





## ----usePython----------------------------------------------------------------
use_condaenv()
py_config()



## ----importPython-------------------------------------------------------------
sc <- import("scanpy")
sc


## ----installAnnData-----------------------------------------------------------
if(!require(anndata)){
  install.packages("anndata")
  require(anndata)
}



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


## ----obsdf2-------------------------------------------------------------------
var_df <- adata$var
var_df[1:2,]


## ----vardf2-------------------------------------------------------------------
obs_df <- adata$obs
obs_df[1:2,,drop=FALSE]

