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
if(params$isSlides != "yes"){
  cat("# load data and access empty drops

---
"    
  )
  
}



## ----loadSCE,eval=FALSE,include=TRUE,echo=TRUE--------------------------------
## library(DropletUtils)
## library(DropletTestFiles)
## fname <- "file path to Cell Ranger results"
## sce <- read10xCounts(fname, col.names=TRUE)
## cellID <- colData(sce)$Barcode
## cellID_sel <- sample(cellID,100000,replace = FALSE) # make subsets
## sce_sub <- sce[,cellID_sel]


## ----loadSCE_io,include=TRUE,echo=FALSE,eval=TRUE-----------------------------
sce_sub <- readRDS("data/scSeq_CTRL_sceSub.rds")


## ----loadSCE_pres,include=TRUE,echo=TRUE,eval=TRUE----------------------------
sce_sub
# cell information
colData(sce_sub)[1:2,]
# gene information
rowData(sce_sub)[1:2,]
sce <- sce_sub


## ----rankUMI,eval=TRUE--------------------------------------------------------
bcrank <- barcodeRanks(counts(sce))
bcrank[1:2,]


## ----rankUMI_knee,eval=FALSE,echo=TRUE,include=TRUE---------------------------
## uniq <- !duplicated(bcrank$rank)
## #
## plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
##     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
## 
## abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
## abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
## 
## legend("bottomleft", legend=c("Inflection", "Knee"),
##         col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)


## ----rankUMI_knee2,eval=TRUE,echo=FALSE,include=TRUE--------------------------
uniq <- !duplicated(bcrank$rank)
#
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
    xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
        col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)


## ----estDroplet,eval=TRUE,echo=TRUE,include=TRUE------------------------------
set.seed(100)
limit <- 100   
e.out <- emptyDrops(counts(sce),lower=limit, test.ambient=TRUE)
#
e.out


## ----estDroplet2,eval=TRUE,echo=TRUE,include=TRUE-----------------------------
# Testeed by FDR
summary(e.out$FDR <= 0.001)
# Concordance by testing with FDR and limited
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)


## ----estDroplet3,eval=TRUE,echo=TRUE,include=TRUE-----------------------------
hist(e.out$PValue[e.out$Total <= limit & e.out$Total > 0],
    xlab="P-value", main="", col="grey80") 


## ----estDroplet4,eval=TRUE,echo=TRUE,include=TRUE-----------------------------
sce2 <- sce[,which(e.out$FDR <= 0.001)]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Normalization and clustering

---
"    
  )
  
}


## ----countNorm,eval=TRUE,echo=TRUE,include=TRUE-------------------------------
library(scran)
library(scuttle)
library(scater)
clusters <- quickCluster(sce2)
sce2 <- computeSumFactors(sce2, cluster=clusters)
sce2 <- logNormCounts(sce2)
sce2


## ----featureIdent,eval=TRUE,echo=TRUE,include=TRUE----------------------------
set.seed(1000)
# modeling variables
dec.pbmc <- modelGeneVarByPoisson(sce2)
# calcualte top features
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)


## ----plotUMAP,eval=TRUE,echo=TRUE,include=TRUE--------------------------------
set.seed(1000)
# Evaluate PCs
sce2 <- denoisePCA(sce2, subset.row=top.pbmc, technical=dec.pbmc)
# make TSNE plot
sce2 <- runTSNE(sce2, dimred="PCA")
# make UMAP plot
sce2 <- runUMAP(sce2, dimred="PCA")


## ----clust,eval=TRUE,echo=TRUE,include=TRUE-----------------------------------
g <- buildSNNGraph(sce2, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce2) <- factor(clust)
#
colData(sce2)


## ----plotUMAP_2,eval=TRUE,echo=TRUE,include=TRUE------------------------------
plotUMAP(sce2,colour_by="label")


## ----plotTSNE,eval=TRUE,echo=TRUE,include=TRUE--------------------------------
plotTSNE(sce2,colour_by ="label")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Remove ambient RNA

---
"    
  )
  
}


## ----estAmb1,eval=TRUE,echo=TRUE,include=TRUE---------------------------------
# extrat potential ambient RNA and thee estimated score
amb <- metadata(e.out)$ambient[,1]
head(amb)


## ----estAmb2,eval=TRUE,echo=TRUE,include=TRUE---------------------------------
library(scater)
stripped <- sce2[names(amb),]
out <- removeAmbience(counts(stripped), ambient=amb,groups = colLabels(stripped))


## ----recalAmb,eval=TRUE,echo=TRUE,include=TRUE--------------------------------
counts(stripped, withDimnames=FALSE) <- out
stripped <- logNormCounts(stripped)


## ----compRMAmb,eval=FALSE,echo=TRUE,include=TRUE------------------------------
## ensmbl_id <- rowData(sce2)$ID[rowData(sce2)$Symbol=="Hba-a1"]
## plotExpression(sce2, x="label", colour_by="label", features=ensmbl_id) +
##         ggtitle("Before")
## plotExpression(stripped, x="label", colour_by="label", features=ensmbl_id) +
##         ggtitle("After")


## ----compRMAmb_pres,eval=TRUE,echo=FALSE,include=TRUE-------------------------
ensmbl_id <- rowData(sce2)$ID[rowData(sce2)$Symbol=="Hba-a1"]
plotExpression(sce2, x="label", colour_by="label", features=ensmbl_id) + 
        ggtitle("Before")
plotExpression(stripped, x="label", colour_by="label", features=ensmbl_id) + 
        ggtitle("After")


## ----compRMAmb2,eval=FALSE,echo=TRUE,include=TRUE-----------------------------
## ensmbl_id <- rowData(sce2)$ID[rowData(sce2)$Symbol=="Krt17"]
## plotExpression(sce2, x="label", colour_by="label", features=ensmbl_id) +
##         ggtitle("Before")
## plotExpression(stripped, x="label", colour_by="label", features=ensmbl_id) +
##         ggtitle("After")


## ----compRMAmb2_pres,eval=TRUE,echo=FALSE,include=TRUE------------------------
ensmbl_id <- rowData(sce2)$ID[rowData(sce2)$Symbol=="Krt17"]
plotExpression(sce2, x="label", colour_by="label", features=ensmbl_id) + 
        ggtitle("Before")
plotExpression(stripped, x="label", colour_by="label", features=ensmbl_id) + 
        ggtitle("After")


## ----rmAmb_store,eval=TRUE,echo=TRUE,include=TRUE-----------------------------
saveRDS(stripped,"scSeq_CTRL_sceSub_rmAmbRNA.rds")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Remove doublets

---
"    
  )
  
}


## ----normClust,eval=FALSE,echo=TRUE,include=TRUE------------------------------
## dec <- modelGeneVar(stripped)
## hvgs <- getTopHVGs(dec,n=1000)
## stripped <- runPCA(stripped, ncomponents=10, subset_row=hvgs)
## stripped <- runUMAP(stripped, dimred="PCA")
## g <- buildSNNGraph(stripped, k=10, use.dimred = 'PCA')
## clust <- igraph::cluster_walktrap(g)$membership
## colLabels(stripped) <- factor(clust)
## plotUMAP(stripped,colour_by="label")


## ----normClust2,eval=TRUE,echo=FALSE,include=TRUE-----------------------------
dec <- modelGeneVar(stripped)
hvgs <- getTopHVGs(dec,n=1000)
stripped <- runPCA(stripped, ncomponents=10, subset_row=hvgs)
stripped <- runUMAP(stripped, dimred="PCA")
g <- buildSNNGraph(stripped, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(stripped) <- factor(clust)
plotUMAP(stripped,colour_by="label")


## ----estDoublet,eval=TRUE,echo=TRUE,include=TRUE------------------------------
dbl.dens <- computeDoubletDensity(stripped, #subset.row=top.mam, 
    d=ncol(reducedDim(stripped)),subset.row=hvgs)
summary(dbl.dens)
stripped$DoubletScore <- dbl.dens


## ----doublerScore,eval=TRUE,echo=TRUE,include=TRUE----------------------------
plotUMAP(stripped,colour_by="DoubletScore")


## ----doublerScorebyClust,eval=TRUE,echo=TRUE,include=TRUE---------------------
plotColData(stripped, x="label", y="DoubletScore", colour_by="label")+
  geom_hline(yintercept = quantile(colData(stripped)$DoubletScore,0.95),lty="dashed",color="red")


## ----rmDoublet,eval=TRUE,echo=TRUE,include=TRUE-------------------------------
cut_off <- quantile(stripped$DoubletScore,0.95)
stripped$isDoublet <- c("no","yes")[factor(as.integer(stripped$DoubletScore>=cut_off),levels=c(0,1))]
table(stripped$isDoublet)
sce_clean <- stripped[stripped$isDoublet=="no",]
saveRDS(sce_clean,"scSeq_CTRL_sceSub_rmAmbRNA_rmDoublet.rds")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# QC plots after cleanrance

---
"    
  )
  
}


## ----evaQC_cal,eval=TRUE,echo=TRUE,include=TRUE-------------------------------
library(scater)
mtGene <- rowData(sce_clean)$ID[grepl(rowData(sce_clean)$Symbol,pattern = "mt-")]
is.mito<- names(sce_clean) %in% mtGene
sce_clean <- addPerCellQC(sce_clean, subsets=list(Mito=is.mito))


## ----qc_mrg,eval=TRUE,echo=TRUE,include=TRUE----------------------------------
plotColData(sce_clean,x="label", y="sum", colour_by="label")+ggtitle("read counts")
plotColData(sce_clean,x="label", y="detected", colour_by="label")+ggtitle("gene counts")
plotColData(sce_clean,x="label", y="subsets_Mito_percent", colour_by="label")+ggtitle("mitocondrial content")


## ----qc_complex,eval=TRUE,echo=TRUE,include=TRUE------------------------------
plotColData(sce_clean,x="sum",y="subsets_Mito_percent",colour_by="label")+ggtitle("is.mito vs read counts")
plotColData(sce_clean,x="sum",y="detected",colour_by="label")+ggtitle("gene counts vs read counts")


## ----varExp,eval=TRUE,echo=TRUE,include=TRUE----------------------------------
vars <- getVarianceExplained(sce_clean, 
    variables=c("DoubletScore","label","sum","detected","subsets_Mito_percent"))
plotExplanatoryVariables(vars)

