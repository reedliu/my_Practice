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

## Exercise 4: BRAC1 in TCGA expression
# BRCA1基因在TCGA数据库的乳腺癌数据集(TCGA PanCancer Atlas)表达
# http://www.cbioportal.org/index.do 
rm(list=ls())
options(stringsAsFactors = F)

f <- read.csv("e4-plot.txt", sep = "\t")
## boxplot
colnames(f) <- c("id", "subtype", "expression", "mut")
da <- f
# install.packages("ggstatsplot")
library(ggstatsplot)
ggbetweenstats(data = da, 
               x = subtype,
               y = expression)
library(ggplot2)
ggsave("e4-BRCA1-TCGA.png")


