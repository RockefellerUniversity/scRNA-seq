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
sc$settings

