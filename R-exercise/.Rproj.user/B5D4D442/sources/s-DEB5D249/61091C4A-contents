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

## Exercise 2: Probe to Symbol ID
# 根据探针名找对应symbol ID
rm(list=ls())
options(stringsAsFactors = F)

# install bioconductor pkgs
######################################
# BUG REPORT: 安装什么软件都提示这个信息 
# “Error in install.packages : error reading from connection”
# TROUBLE SHOTTING: 
# 修改Global options中的pkgs的安装源即可
######################################
if(T){
  if(!require("BiocManager")) 
    install.packages("BiocManager",update = F,ask = F)
  if(length(getOption("BioC_mirror"))==0) 
    options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
  for (pkg in c("hgu133a.db")){
    if (! require(pkg,character.only=T) ) {
      BiocManager::install(pkg,ask = F,update = F)
      require(pkg,character.only=T) 
    }
  }
}
library(hgu133a.db)
# ls("package:hgu133a.db")
p2s=toTable(hgu133aSYMBOL);head(p2s)

a=read.table('e2.txt')
names(a) <- names(p2s)[1]
# 方法一：利用merge
tmp1 <- merge(a,p2s, by="probe_id")
# 方法二：利用match得到第一组向量在第二组中的坐标
tmp2 <- p2s[match(a$probe_id,p2s$probe_id),]

######################################
## 附：判断得到的两组结果是否一致
######################################
# 法一：
identical(tmp1,tmp2) #返回逻辑值
# 法二：
dplyr::setdiff(tmp1,tmp2) #返回两组的差别【没差就返回空】




