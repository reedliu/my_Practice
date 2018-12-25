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

## Exercise 3: Gene Expression
# 从内置数据集的表达矩阵中找TP53基因的表达量
rm(list=ls())
options(stringsAsFactors = F)

suppressMessages(library(CLL))
data(sCLLex)
# sCLLex
###################################
# TIPS:
# 探针的三大块内容：表达矩阵assay/exprs、探针信息featureData、样本信息phenoData
# https://upload-images.jianshu.io/upload_images/9376801-899b4c4bd2e81771.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240
###################################

exprSet <- exprs(sCLLex) 
# exprSet[1:4,1:4]
pdata <- pData(sCLLex) 

library(hgu95av2.db) # annotation info
p2s <- toTable(hgu95av2SYMBOL);head(p2s)

# boxplot [find TP53 has 3 probe IDs]
boxplot(exprSet["1939_at",] ~ pdata$Disease)
boxplot(exprSet["1974_s_at",] ~ pdata$Disease)
boxplot(exprSet["31618_at",] ~ pdata$Disease)









