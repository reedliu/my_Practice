### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-24
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-24  First version 
### From: http://www.bio-info-trainee.com/3750.html
### ---------------

## Exercise 5: Retrive genes from GEO to plot heatmap

rm(list=ls())
options(stringsAsFactors = F)

## first download GSE and load expr
GSE_expr <- function(GSE){
  library(GEOquery)
  if(!file.exists(GSE)){
    geo <- getGEO(GSE, destdir = '.', getGPL = F, AnnotGPL = F)
    gdata <- paste0(GSE,'.eSet.Rdata')
    save(geo, file = gdata)
  }
  load(gdata)
  expr <<- exprs(geo[[1]])
  return(expr)
}
GSE_expr("GSE17215")

dim(expr)
expr[1:4,1:4]

## 



