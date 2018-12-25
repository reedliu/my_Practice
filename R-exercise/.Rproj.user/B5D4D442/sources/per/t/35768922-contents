### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-24
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-24  First version correlation
### ---------------

## Exercise 6: Correlation Analysis
# 核心：cor()、reduce dimensions of expr
rm(list=ls())
options(stringsAsFactors = F)

#########################################
# TIP1: 相关性高可能有两个原因：一个是本身的技术重复；
# 另一个可能没有去除表达矩阵中的0，太多的0会忽略掉真实值
#########################################

## 用bioconductor的数据包airway包来做
#########################################
# TIP2: Bioconductor有三种包：软件包、数据包、注释包
#########################################
suppressMessages(library(airway))
data(airway)
airway
expr <- assay(airway)
dim(expr);colnames(expr)

## first look at original correlation and plot heatmap
library(pheatmap)
pheatmap(cor(expr)) 
# 可以看到，原始矩阵的各个样本建相关性都很高，这样其实
# 不对，因为毕竟有处理、对照组。为了再次确定，可以去看
# 样本信息
grp <- colData(airway)[,3]
# 在heatmap中添加样本信息，构建的grp_df的行名需要是热图数据的列名
grp_df <- data.frame(grp)
rownames(grp_df) <- colnames(expr)
pheatmap(cor(expr), annotation_col = grp_df,filename = "before-trans-cor.png")

## second remove/filter those non-representative data (e.g. 0)
#########################################
# TIP3: HOW DOES APPLY WORK?
# x means each row of expr; we want to know if every value in x > 1 => "x>1" (return T/F)
# then we calculate how many values in x >1 => "sum(x>1)" (T means 1, so it can be summed)
# then if there are five values > 1 => "sum(x>1)>5" (return T/F)
# finally we select these x rows 
#########################################
dim(expr) # previous
expr <- expr[apply(expr, 1, function(x) sum(x>1) > 5),]
dim(expr) # after filter

# If we want to see filtered genes belongs to which group (coding, non-coding...)
# we need to use "exercise-01" Rscript
if(F){
  suppressMessages(library(org.Hs.eg.db))
  # ls("package:org.Hs.eg.db")
  g2s <- toTable(org.Hs.egSYMBOL);head(g2s)
  g2e <- toTable(org.Hs.egENSEMBL);head(g2e)
  a <- as.data.frame(rownames(expr))
  colnames(a) <- "ensembl_id"
  tmp <- merge(a,g2e, by="ensembl_id")
  tmp <- merge(tmp,g2s, by="gene_id")
  hgnc <- read.csv('HGNC-gene-group.txt', sep = "\t")
  hgnc <- hgnc[,c("NCBI.Gene.ID","Locus.type")]
  names(tmp)[1] <- "NCBI.Gene.ID"
  tmp <- merge(tmp,hgnc, by="NCBI.Gene.ID")
  View(tmp)
  gtype <- dplyr::count(tmp,Locus.type,sort = T)
  
}

## third we use edgeR::cpm to remove library size effect
expr <- log(edgeR::cpm(expr)+1)
dim(expr)

## forth we choose 500 signicant-expressed genes 
expr <- expr[sort(apply(expr,1,mad),decreasing = T)[1:500] %>% names,]
dim(expr)

## fifth we do log2 transform and re-calculate the correlation
newcor <- cor(log2(expr+1))
pheatmap(newcor, annotation_col = grp_df,filename = "e6-after-trans-cor.png")





