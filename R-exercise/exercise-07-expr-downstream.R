### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-24
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-24  First version correlation
### ---------------

## Exercise 7: expression downstream analysis
rm(list=ls())
options(stringsAsFactors = F)
suppressMessages(library(CLL))

#################################
## load data and get expr & grp
#################################
data(sCLLex)
# sCLLex
expr <- exprs(sCLLex)
pdata <- pData(sCLLex)
grp <- as.character(pdata$Disease)

# table(grp)
# dim(expr)
# expr[1:4,1:4]
# boxplot(expr)
# 单对一个基因：boxplot(expr[1,]~grp)
# 单对一个基因：t.test(expr[1,]~grp)
# 但是对成千上万基因，这样循环比较复杂；并且统计过程不是特别严谨，limma加了矫正的过程

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
