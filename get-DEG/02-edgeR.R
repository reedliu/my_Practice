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
keep <- rowSums(cpm(e)>1) >= 2
e <- e[keep, , keep.lib.sizes=FALSE]
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
# https://www.biostars.org/p/110861/
lrt <- glmLRT(fit,  contrast=c(-1,1)) # accoding to design to modify
nrDEG=topTags(lrt, n=nrow(DEG))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
edgeR_DEG=nrDEG
#write.csv(edgeR_DEG,"DEG_edgeR.csv",quote = F)






