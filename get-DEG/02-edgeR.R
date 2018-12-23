### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-22
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-22  First version about edgeR
###
### ---------------
rm(list = ls())
options(stringsAsFactors = F)

load('airway.expreSet.Rdata')
library(edgeR)
e <- DGEList(counts=expr,group=factor(grp))
e$samples$lib.size <- colSums(e$counts)
e <- calcNormFactors(e)
e$samples

DEG=e
design <- model.matrix(~0+factor(grp))
rownames(design)<-colnames(DEG)
colnames(design)<-levels(factor(grp))

DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

fit <- glmFit(DEG, design)

lrt <- glmLRT(fit,  contrast=c(1,0)) # accoding to design to modify
nrDEG=topTags(lrt, n=nrow(expr))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
edgeR_DEG=nrDEG
#write.csv(edgeR_DEG,"DEG_edgeR.csv",quote = F)






