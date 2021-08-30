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
  packages=c("python-igraph","scanpy","louvain","leidenalg","loompy")
)





## ----usePython----------------------------------------------------------------
use_condaenv()
py_config()



## ----importPython-------------------------------------------------------------
sc <- reticulate::import("scanpy")
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


## ----santize------------------------------------------------------------------
# sc$utils$sanitize_anndata(adata)


## ----mitoQC,tidy=FALSE--------------------------------------------------------
sc$pp$calculate_qc_metrics(adata,
                           qc_vars = list("mito"),
                           inplace = TRUE)
adata


## ----moreQC,tidy=FALSE--------------------------------------------------------
adata$var[5:10,]


## ----moreQC2,tidy=FALSE-------------------------------------------------------
adata$obs[1:2,]


## ----plotHiExpr---------------------------------------------------------------
sc$pl$highest_expr_genes(adata,gene_symbols = "gene_symbols")


## ----plotQCpy2----------------------------------------------------------------
sc$pl$violin(adata, list('n_genes_by_counts', 'total_counts', 'pct_counts_mito'),
             jitter=0.4, multi_panel=TRUE)


## ----plotQCpy3----------------------------------------------------------------
sc$pl$scatter(adata, x='total_counts', y='n_genes_by_counts')



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


## ----filterToVariable---------------------------------------------------------
# keep only highly variable genes:
adata = adata$copy()
adata = adata[,adata$var$highly_variable]
adata = adata$copy()
# regress out total counts per cell and the percentage of mitochondrial genes expressed
# adata.copy = adata
# sc$pp$regress_out(adata, list('n_counts', 'pct_counts_mito'),copy=FALSE) #, n_jobs=args.threads)


## ----regress------------------------------------------------------------------
sc$pp$regress_out(adata, list('n_counts', 'pct_counts_mito'))

# adata = adata2
# scale each gene to unit variance, clip values exceeding SD 10.


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
# 


## ----plotumapCells------------------------------------------------------------
sc$pl$umap(adata, color=list('leiden'))


## ----saveToh5ad---------------------------------------------------------------
adata$write_h5ad(filename = "PBMC_Scanpy.h5ad")


## ----saveToLoom---------------------------------------------------------------
adata$write_loom(filename = "PBMC_Scanpy.loom")


## ----readLoom-----------------------------------------------------------------
if(!require(LoomExperiment)){
  BiocManager::install("LoomExperiment")
  require(LoomExperiment)
}
loom <- LoomExperiment::import(con = "PBMC_Scanpy.loom")
loom


## ----loomToSCE----------------------------------------------------------------
sce_loom <- as(loom,"SingleCellExperiment")
sce_loom


## ----replaceUMAP--------------------------------------------------------------
if(!require(scater)){
  BiocManager::install("scater")
  require(scater)
}
reducedDim(sce_loom, "UMAP") <- adata$obsm[["X_umap"]]
reducedDim(sce_loom, "PCA") <- adata$obsm[["X_pca"]]


## ----plotUMAP-----------------------------------------------------------------
plotUMAP(sce_loom,colour_by="leiden")

