### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-23
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-23  First version about three DEG pkgs comparison
###
### ---------------
rm(list = ls())
options(stringsAsFactors = F)

load('DEG_result.Rdata')
####################################
## get all genes and FC to check correlation 
# 帮助判断不同方法是否正确（理论上两种方法得到的结果应该在一条直线上）
####################################
deseq <- data.frame(gene=rownames(Deseq_DEG),
                    Deseq=Deseq_DEG$log2FoldChange)

limma <- data.frame(gene=rownames(limma_voom_DEG),
                    limma=limma_voom_DEG$logFC)

edgeR <- data.frame(gene=rownames(edgeR_DEG),
                    edgeR=edgeR_DEG$logFC)
# e.g. deseq & limma
test1 <- merge(deseq,limma,by="gene")
plot(test1[2:3])

####################################
## get UP and DOWN genes for each method to draw Venn plot
# 分别得到不同方法的上调、下调基因，然后分别做Venn图
####################################
if(T){
  # get Deseq2 UP and DOWN genes
  DEG = Deseq_DEG
  logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
  
  DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  Deseq_up <- rownames(DEG[DEG$change =='UP',])
  Deseq_down <- rownames(DEG[DEG$change =='DOWN',])
  
  # get limma UP and DOWN genes
  DEG = limma_voom_DEG
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  
  DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  limma_up <- rownames(DEG[DEG$change =='UP',])
  limma_down <- rownames(DEG[DEG$change =='DOWN',])
  
  # get edgeR UP and DOWN genes
  DEG = edgeR_DEG
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  
  DEG$change = as.factor(ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  edgeR_up <- rownames(DEG[DEG$change =='UP',])
  edgeR_down <- rownames(DEG[DEG$change =='DOWN',])
}

## Venn plot
if(!require(VennDiagram))install.packages('VennDiagram')
library (VennDiagram)
Deseq=Deseq_up #for DOWN need to re-define
limma=limma_up
edgeR=edgeR_up
if(T){
  # first look at UP comparison
  venn.diagram(x= list(Deseq = Deseq,limma = limma,edgeR = edgeR),
               filename = "compare_UP.png",
               height = 450, width = 450,
               resolution =300,
               imagetype="png",
               col="transparent",
               fill=c("green","yellow","darkorchid1"),
               alpha = 0.50,
               cex=0.45,
               cat.cex=0.45)
  # DOWN comparison
  venn.diagram(x= list(Deseq = Deseq,limma = limma,edgeR = edgeR),
               filename = "compare_DOWN.png",
               height = 450, width = 450,
               resolution =300,
               imagetype="png",
               col="transparent",
               fill=c("green","yellow","darkorchid1"),
               alpha = 0.50,
               cex=0.45,
               cat.cex=0.45)
  
  
}
