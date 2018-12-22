### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-22
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-22  First version about limma-voom
###
### ---------------
rm(list = ls())
options(stringsAsFactors = F)

load('airway.expreSet.Rdata')
suppressMessages(library(limma))
design <- model.matrix(~0+factor(grp))
colnames(design)=levels(factor(grp))
rownames(design)=colnames(expr)


DEG <- edgeR::DGEList(counts=expr)
DEG <- calcNormFactors(DEG)
# logCPM <- cpm(DEG, log=TRUE, prior.count=3)

v <- voom(DEG,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

# for one comparison (trt-untrt)
cont.matrix=makeContrasts(contrasts=c('trt-untrt'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

result = topTable(fit2, coef='trt-untrt', n=Inf)
limma_voom_DEG = na.omit(result)
write.csv(limma_voom_DEG,"limma_voom_DEG.csv",quote = F)

# for two comparisons
if(F){
  cont.matrix=makeContrasts(contrasts=c('treat_1-control','treat_2-control'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  result1 = topTable(fit2, coef='treat_1-control', n=Inf)
  trt1_limma_voom_DEG = na.omit(result1)
  write.csv(limma_voom_DEG,"limma_voom_DEG.csv",quote = F)
  
  result2 = topTable(fit2, coef='treat_2-control', n=Inf)
  trt2_limma_voom_DEG = na.omit(result2)
  write.csv(trt2_limma_voom_DEG,"trt2_limma_voom_DEG.csv",quote = F)
}
