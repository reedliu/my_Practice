### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-22
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-22  First version about limma-voom
### Update Log: 2018-12-25  Add new limma DEG function
###
### ---------------
rm(list = ls())
options(stringsAsFactors = F)

load('airway.expreSet.Rdata')
suppressMessages(library(limma))
library(edgeR)
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
## save three main DEG pkgs results
save("Deseq_DEG", "limma_voom_DEG", "edgeR_DEG", file="DEG_result.Rdata")


# Update on 2018-12-25
#################################
## DEG with limma
#################################
suppressMessages(library(limma))
#limma needs：表达矩阵（expr）、分组矩阵（design）、比较矩阵（contrast）
#先做一个分组矩阵～design，说明progres是哪几个样本，stable又是哪几个，其中1代表“是”
design <- model.matrix(~0+factor(grp))
colnames(design) <- levels(factor(grp))
rownames(design) <- colnames(expr)
design
#再做一个比较矩阵【一般是case比control】
contrast<-makeContrasts(paste0(unique(grp),collapse = "-"),levels = design)
contrast

DEG <- function(expr,design,contrast){
  ##step1
  fit <- lmFit(expr,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast) 
  fit2 <- eBayes(fit2)  
  ##step3
  mtx = topTable(fit2, coef=1, n=Inf)
  deg_mtx = na.omit(mtx) 
  return(deg_mtx)
}
DEG_mtx <- DEG(expr,design,contrast)
View(DEG_mtx)
