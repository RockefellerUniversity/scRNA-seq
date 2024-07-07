params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Cell Ranger -  Output files

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell Ranger - Output files

---
"    
  )
  
}



## -----------------------------------------------------------------------------

library(DropletUtils)



## -----------------------------------------------------------------------------

h5file <- "https://rubioinformatics.s3.amazonaws.com/scRNA_graduate/SRR_NeuroD1/outs/NeuroD1_filtered_feature_bc_matrix.h5"
local_h5file <- basename(h5file)
download.file(h5file,local_h5file)

sce.NeuroD1_filtered <- read10xCounts(local_h5file, col.names=TRUE)

class(sce.NeuroD1_filtered)



## ----include=FALSE------------------------------------------------------------
unlink(h5file)


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
bcrank <- barcodeRanks(counts(sce.NeuroD1_unfiltered))
bcrank


## -----------------------------------------------------------------------------
bcrank$filtered <- rownames(bcrank) %in% colnames(sce.NeuroD1_filtered)
bc_plot <- as.data.frame(bcrank)
bc_plot <- bc_plot[order(bc_plot$filtered,decreasing=TRUE),]
bc_plot <- bc_plot[!duplicated(bc_plot$rank),]


## -----------------------------------------------------------------------------
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


## ----eval=FALSE,echo=FALSE----------------------------------------------------
## e.out <- readRDS(system.file("extdata/data/e.out.RData", mustWork = TRUE, package = "scRNASeq"))
## )
## table(e.out)
## head(e.out)


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



## -----------------------------------------------------------------------------
sce.NeuroD1_filtered


## -----------------------------------------------------------------------------
is.mito <- grepl("^mt",rowData(sce.NeuroD1_filtered)$Symbol)
is.ribo <- grepl("^Rps",rowData(sce.NeuroD1_filtered)$Symbol)

table(is.mito)
table(is.ribo)



## -----------------------------------------------------------------------------
library(scuttle)
sce.NeuroD1_filtered <- addPerCellQCMetrics(sce.NeuroD1_filtered, 
                                            subsets=list(Mito=is.mito,Ribosomal=is.ribo))
sce.NeuroD1_filtered


## -----------------------------------------------------------------------------
qc_df <- colData(sce.NeuroD1_filtered)
head(qc_df,n=2)


## -----------------------------------------------------------------------------
library(scater)
plotColData(sce.NeuroD1_filtered,x="Sample", y="sum")


## -----------------------------------------------------------------------------
library(scater)
plotColData(sce.NeuroD1_filtered,x="Sample",y="detected")



## -----------------------------------------------------------------------------
library(scater)
p1 <- plotColData(sce.NeuroD1_filtered,x="Sample",y="subsets_Mito_percent")
p2 <- plotColData(sce.NeuroD1_filtered,x="Sample",y="subsets_Ribosomal_percent")
gridExtra::grid.arrange(p1,p2,ncol=2)


## -----------------------------------------------------------------------------
library(scater)
p1 <- plotColData(sce.NeuroD1_filtered,x="sum",y="detected")
p2 <- plotColData(sce.NeuroD1_filtered,x="sum",y="subsets_Mito_percent")
p3 <- plotColData(sce.NeuroD1_filtered,x="detected",y="subsets_Mito_percent")
gridExtra::grid.arrange(p1,p2,p3,ncol=3)


## -----------------------------------------------------------------------------
qc.high_lib_size <- colData(sce.NeuroD1_filtered)$sum > 125000
qc.min_detected <- colData(sce.NeuroD1_filtered)$detected < 200
qc.mito <- colData(sce.NeuroD1_filtered)$subsets_Mito_percent > 25
discard <- qc.high_lib_size | qc.mito | qc.min_detected
DataFrame(LibSize=sum(qc.high_lib_size), Detected=sum(qc.min_detected),MitoProp=sum(qc.mito), Total=sum(discard))


## -----------------------------------------------------------------------------
colData(sce.NeuroD1_filtered) <- cbind(colData(sce.NeuroD1_filtered),DataFrame(toDiscard=discard))
p1 <- plotColData(sce.NeuroD1_filtered,x="sum",y="detected",colour_by = "toDiscard")
p2 <- plotColData(sce.NeuroD1_filtered,x="sum",y="subsets_Mito_percent",colour_by = "toDiscard")
p3 <- plotColData(sce.NeuroD1_filtered,x="detected",y="subsets_Mito_percent",colour_by = "toDiscard")
gridExtra::grid.arrange(p1,p2,p3,ncol=3)


## -----------------------------------------------------------------------------
sce.NeuroD1_filtered_QCed <- sce.NeuroD1_filtered[,sce.NeuroD1_filtered$toDiscard %in% "FALSE"]
sce.NeuroD1_filtered_QCed


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Bioconductor  -  Normalisation and Clustering

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Bioconductor  -  Normalisation and Clustering

---
"    
  )
  
}



## -----------------------------------------------------------------------------
sce.NeuroD1_filtered_QCed <- logNormCounts(sce.NeuroD1_filtered_QCed)
assayNames(sce.NeuroD1_filtered_QCed)


## -----------------------------------------------------------------------------
sce.NeuroD1_filtered_QCed <- logNormCounts(sce.NeuroD1_filtered_QCed)
assayNames(sce.NeuroD1_filtered_QCed)


## -----------------------------------------------------------------------------
library(scran)
clust.sce.NeuroD1_filtered_QCed <- quickCluster(sce.NeuroD1_filtered_QCed) 
sce.NeuroD1_filtered_trimmed <- computeSumFactors(sce.NeuroD1_filtered_QCed, cluster=clust.sce.NeuroD1_filtered_QCed)
sce.NeuroD1_filtered_QCed <- logNormCounts(sce.NeuroD1_filtered_QCed)
assayNames(sce.NeuroD1_filtered_QCed)



## -----------------------------------------------------------------------------
##
dec.NeuroD1_filtered_QCed <- modelGeneVar(sce.NeuroD1_filtered_QCed)
ggplot(dec.NeuroD1_filtered_QCed,aes(x=mean,y=total))+geom_point()


## -----------------------------------------------------------------------------
top.NeuroD1_filtered_QCed <- getTopHVGs(dec.NeuroD1_filtered_QCed, n=3000)



## -----------------------------------------------------------------------------
set.seed(100) # See below.
sce.NeuroD1_filtered_QCed <- fixedPCA(sce.NeuroD1_filtered_QCed,
                                      subset.row=top.NeuroD1_filtered_QCed) 
reducedDimNames(sce.NeuroD1_filtered_QCed)



## -----------------------------------------------------------------------------
plotReducedDim(sce.NeuroD1_filtered_QCed, dimred="PCA",colour_by = "subsets_Mito_percent")



## -----------------------------------------------------------------------------
percent.var <- attr(reducedDim(sce.NeuroD1_filtered_QCed), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
# library(PCAtools)
# chosen.elbow <- findElbowPoint(percent.var)


## -----------------------------------------------------------------------------

sce.NeuroD1_filtered_QCed <- runTSNE(sce.NeuroD1_filtered_QCed,n_dimred=30)
sce.NeuroD1_filtered_QCed <- runUMAP(sce.NeuroD1_filtered_QCed,n_dimred=30)
reducedDimNames(sce.NeuroD1_filtered_QCed)


## -----------------------------------------------------------------------------
plotUMAP(sce.NeuroD1_filtered_QCed,colour_by="subsets_Mito_percent")
plotTSNE(sce.NeuroD1_filtered_QCed,colour_by="subsets_Mito_percent")

