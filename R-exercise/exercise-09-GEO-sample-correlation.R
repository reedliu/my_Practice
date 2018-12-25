### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-25
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-25  First version 
### From: http://www.bio-info-trainee.com/3750.html
### ---------------

## Exercise 9: Retrive GEO data and calculate correlation among samples

rm(list=ls())
options(stringsAsFactors = F)

#################################
## FIRST download GSE and load expr
#################################
GSE_expr <- function(GSE){
  library(GEOquery)
  if(!file.exists(GSE)){
    geo <<- getGEO(GSE, destdir = '.', getGPL = F, AnnotGPL = F)
    gdata <- paste0(GSE,'.eSet.Rdata')
    save(geo, file = gdata)
  }
  message ("This is created by Reed Liu (jieandze1314@gmail.com)")
  load(gdata)
  expr <<- exprs(geo[[1]])
}
GSE_expr("GSE24673")

dim(expr)
expr[1:4,1:4]

#################################
## SECOND get sample info
#################################
pdata <- pData(geo[[1]]) 
# 自己根据pdata第八列做一个grp
grp <- c('rbc','rbc','rbc',
         'rbn','rbn','rbn',
         'rbc','rbc','rbc',
         'normal','normal')

#################################
## THIRD sample correlation heatmap
#################################
grp_df <- data.frame(group=grp)
rownames(grp_df) <- colnames(expr)
new_cor <- cor(expr)
pheatmap::pheatmap(new_cor, annotation_col = grp_df)





