### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-24
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-24  First version 
### Update Log: 2018-12-25  Add another solution (using dplyr) from Xiaojie 
### From: http://www.bio-info-trainee.com/3750.html
### ---------------

## Exercise 5: Retrive genes from GEO to plot heatmap

rm(list=ls())
options(stringsAsFactors = F)

#################################
## FIRST download GSE and load expr
## COPYRIGHT: Reed Liu
#################################
GSE_expr <- function(GSE){
  library(GEOquery)
  if(!file.exists(GSE)){
    geo <<- getGEO(GSE, destdir = '.', getGPL = F, AnnotGPL = F)
    gdata <- paste0(GSE,'.eSet.Rdata')
    save(geo, file = gdata)
  }
  echo "This is created by Reed Liu (jieandze1314@gmail.com)"
  load(gdata)
  expr <<- exprs(geo[[1]])
}
GSE_expr("GSE17215")

dim(expr)
expr[1:4,1:4]

#################################
## SECOND get annotation 
#################################
# geo => check Annotation: GPL...
# GPL to Bioc_anno please refers to https://www.jianshu.com/p/40b560755cdf

for (pkg in c("hgu133a.db")){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

#################################
## THIRD probe2geneID in expr
#################################
p2s=toTable(hgu133aSYMBOL);head(p2s)
dim(expr); expr <- expr[p2s$probe_id,]; dim(expr)
# 方法一：【弊端：不能判断去掉什么样的重复行，只是保留第一次出现的】
# https://www.biostars.org/p/109248/
probeID <- rownames(expr)
tmp <- data.frame(Gene=unlist(mget(x = probeID,envir = hgu133aSYMBOL)))
tmp <- data.frame(tmp,probe=rownames(tmp))
tmp <- tmp[!duplicated(tmp$Gene),] 
rownames(expr) <- tmp[match(rownames(expr),tmp$probe),][,1]

# 方法二：针对重复行的问题，采用了中位数的方法来判断，将中位数大的放第一个，就可以保留

p2s$median <- apply(expr, 1, median)
p2s <- p2s[order(p2s$symbol, p2s$median, decreasing = T),]
p2s <- p2s[!duplicated(p2s$symbol),]
expr <- expr[p2s$probe_id,]
rownames(expr) <- p2s$symbol

dim(expr)
View(expr)

#################################
## FOURTH match wanted genes to expr 
#################################
# get a goup of genes called "gp"
gp <- "ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T"
gp <- strsplit(gp, ' ')[[1]]
# 芯片的年代决定了有的新基因在表达矩阵中不存在，因此要先看看有多少基因是
# 没有注释在芯片表达矩阵中的，然后再过滤
gp %in% rownames(expr) %>% table 
f_gp <- gp[gp %in% rownames(expr)] # filered gene group
gp_expr <- expr[f_gp,] 

#################################
## FIFTH plot heatmap 
#################################
# 根据gp_expr做热图（做之前先log2处理下，对数据降维）
gp_expr <- log2(gp_expr)
# 做热图主要是看样本间的差异，所以要忽略基因间的差异，用scale按行(gene)做归一化处理
pheatmap::pheatmap(gp_expr, scale = "row") 


#################################
## 另一种解法 from 小洁
#################################
# 先得到expr和p2s
tmp <- dplyr::filter(p2s, p2s$symbol %in% gp)
tmp2 <- tibble::rownames_to_column(data.frame(expr),"probe_id")
tmp3 <- merge(tmp,tmp2,by="probe_id")

tmp3$median <- apply(tmp3[,3:ncol(tmp3)], 1, median)
tmp3 <- tmp3[order(tmp3$symbol, tmp3$median, decreasing = T),]
tmp3 <- tmp3[!duplicated(tmp3$symbol),]
rownames(tmp3) <- tmp3$symbol
tmp3 <- tmp3[,-c(1,2,ncol(tmp3))]

gp_expr <- log2(tmp3)
pheatmap::pheatmap(gp_expr, scale = "row") 




