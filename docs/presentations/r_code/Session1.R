params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Cell Ranger results

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell Ranger results


---
"    
  )
  
}



## ----scSeq_GI_cellCCount,echo=FALSE,out.width = "50%",fig.align="center"------
knitr::include_graphics("imgs/cell_count.png")


## ----scSeq_GI_sessionInfot,echo=FALSE,out.width = "75%",fig.align="center"----
knitr::include_graphics("imgs/session_info.png")


## ----scSeq_GI_seqQC,echo=FALSE,out.width = "75%",fig.align="center"-----------
knitr::include_graphics("imgs/seq_qc.png")


## ----scSeq_GI_mapQC,echo=FALSE,out.width = "75%",fig.align="center"-----------
knitr::include_graphics("imgs/map_info.png")


## ----scSeq_GI_cellCount2,echo=FALSE,out.width = "75%",fig.align="center"------
knitr::include_graphics("imgs/cell_count.png")


## ----scSeq_GI_detCellNum,echo=FALSE,out.width = "50%",fig.align="center"------
knitr::include_graphics("imgs/knee_plot.png")


## ----scSeq_GI_geneEvl,echo=FALSE,out.width = "75%",fig.align="center"---------
knitr::include_graphics("imgs/gene_eval.png")


## ----scSeq_GI_satPlot,echo=FALSE,out.width = "50%",fig.align="center"---------
knitr::include_graphics("imgs/umi_eval.png")


## ----scSeq_GI_tsnePlot,echo=FALSE,out.width = "100%",fig.align="center"-------
knitr::include_graphics("imgs/tsne.png")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Customized Analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Customized Analysis


---
"    
  )
  
}


