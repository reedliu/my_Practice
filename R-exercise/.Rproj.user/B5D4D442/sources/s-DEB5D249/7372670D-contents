### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-25
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-25  First version 
### From: http://www.bio-info-trainee.com/3750.html
### ---------------

## Exercise 10: Retrive GEO data and （sort + DEG analysis）
## 核心：sort/ apply/ unlist+lapply/ strsplit

rm(list=ls())
options(stringsAsFactors = F)

GSE_expr <- function(GSE){
  library(GEOquery)
  if(!file.exists(GSE)){
    geo <<- getGEO(GSE, destdir = '.', getGPL = F, AnnotGPL = F)
    gdata <- paste0(GSE,'.eSet.Rdata')
    save(geo, file = gdata)
  }
  message ("This is created by Reed Liu (jieandze1314@gmail.com)")
  load(gdata)
  expr <<- exprs(geo[[1]])
}
GSE_expr("GSE42872")

dim(expr)
expr[1:4,1:4]

#########################################
# 选出所有样本的(mean/sd/mad/)最大的探针
#########################################
sort(apply(expr,1,mean),decreasing = T)[1]
sort(apply(expr,1,sd),decreasing = T)[1]
sort(apply(expr,1,mad),decreasing = T)[1]

#########################################
# 得到表型信息，然后用limma做差异分析
#########################################
pdata <- pData(geo[[1]]) 
grp <- unlist(lapply(pdata$title, function(x){
  strsplit(x, ' ')[[1]][4]
}))

########################
# 开始limma分析
########################
suppressMessages(library(limma))
#limma needs：表达矩阵（expr：需要用log后的矩阵）、分组矩阵（design）、比较矩阵（contrast）
#先做一个分组矩阵～design，说明progres是哪几个样本，stable又是哪几个
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

############
# 火山图
############
DEG=DEG_mtx

if(T){
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  
  DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  title <- paste0('log2FoldChange cutoff: ',round(logFC_cutoff,3),
                  '\nUp-regulated genes: ',nrow(DEG[DEG$change =='UP',]) ,
                  '\nDown-regulated genes: ',nrow(DEG[DEG$change =='DOWN',])
  )
}
library(ggplot2)
vol1 = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2FoldChange") + ylab("-log10 p-value") +
  ggtitle( title ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) # according to the levels(DEG$change)
print(vol1)
