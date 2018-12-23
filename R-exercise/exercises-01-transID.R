### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-23
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-23  First version 
### From: http://www.bio-info-trainee.com/3750.html
### ---------------

## Exercise 1: Trans ID
# 根据ensembl ID找对应symbol ID
# 核心：利用toTable和merge

rm(list=ls())
options(stringsAsFactors = F)

a=read.table('e1.txt')
suppressMessages(library(org.Hs.eg.db))
# ls("package:org.Hs.eg.db")
g2s <- toTable(org.Hs.egSYMBOL);head(g2s)
g2e <- toTable(org.Hs.egENSEMBL);head(g2e)
library(stringr)
a$ensembl_id = lapply(a$V1, function(x){
  str_split(x,'[.]')[[1]][1]
}) %>% unlist
tmp <- merge(a,g2e, by="ensembl_id")
tmp <- merge(tmp,g2s, by="gene_id")
# 根据org.Hs.eg.db信息，可以继续添加其他信息，比如CHR信息
g2c <- toTable(org.Hs.egCHR);head(g2c)
tmp <- merge(tmp, g2c, by="gene_id")

## 找到symbol ID对应的gene type
# HGNC官网下载：https://www.genenames.org/cgi-bin/genegroup/download-all
hgnc <- read.csv('~/Downloads/gene-group.txt', sep = "\t")
hgnc <- hgnc[,c("NCBI.Gene.ID","Locus.type")]
names(tmp)[1] <- "NCBI.Gene.ID"
tmp <- merge(tmp,hgnc, by="NCBI.Gene.ID")




