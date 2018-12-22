### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2018-12-22
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
### Update Log: 2018-12-22  First version about Deseq2, Heatmap, Vocalno
###
### ---------------
rm(list = ls())
options(stringsAsFactors = F)

# If download original exprSet, firstly, read in data as matrix
# try fread: 
# install.packages("data.table")
# library(data.table) NOT A GOOD CHOICE HERE
expr <- read.table("expression.txt", header = T, row.names = 1)
expr <- as.matrix(expr)

# Here, we use a built-in exprSet data
if(F){
  suppressMessages(library(airway))
  data(airway) # airway is an object => type "airway" 
  expr <- assay(airway)
  grp <- colData(airway)[,3]
  save("expr","grp", file= 'airway.expreSet.Rdata')
}
rm(list = ls())
options(stringsAsFactors = F)
load('airway.expreSet.Rdata')

suppressMessages(library(DESeq2)) 
(colData <- data.frame(row.names=colnames(expr), 
                       grp=grp))
dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = colData,
                              design = ~ grp)
dds <- DESeq(dds)
res <- results(dds, 
               contrast=c("grp","trt","untrt"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG_deseq=as.data.frame(resOrdered)
DEG_deseq = na.omit(DEG_deseq)
# write.csv(DEG_deseq,"DEG_deseq.csv")

DEG=DEG_deseq
#################################
## Heatmap
#################################
library(pheatmap)
part_gene=head(rownames(DEG),50) 
part_matrix=expr[part_gene,]
part_matrix=t(scale(t(part_matrix)))

annotation_col = data.frame(
  Case = factor(rep(c("treat", "untreat"), each=4))
)
rownames(annotation_col) = c("SRR1039508", "SRR1039520", "SRR1039512","SRR1039516", 
                             "SRR1039517", "SRR1039513", "SRR1039509", "SRR1039521")
pheatmap(part_matrix, annotation_col = annotation_col)

#################################
## Volcano plot
#################################
#colnames(DEG_deseq)
#plot(DEG_deseq$log2FoldChange,-log10(DEG_deseq$pvalue))

DEG=DEG_deseq
if(T){
  logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
  
  DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  title <- paste0('log2FoldChange cutoff: ',round(logFC_cutoff,3),
                  '\nUp-regulated genes: ',nrow(DEG[DEG$change =='UP',]) ,
                  '\nDown-regulated genes: ',nrow(DEG[DEG$change =='DOWN',])
  )
}

library(ggplot2)
vol1 = ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(pvalue), color=significant)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2FoldChange") + ylab("-log10 p-value") +
  ggtitle( title ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) # according to the levels(DEG$change)
print(vol1)
ggsave(vol1,filename = 'Deseq2_volcano.png')

## DIY vocalno plot
library(ggrepel)
library(ggplot2)
data=DEG_deseq
data$gene <- rownames(data)
data$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                              ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)

vol2 = ggplot(data=data, aes(x=log2FoldChange, y =-log10(pvalue),color =change)) +
  #设置点的透明度、大小等
  geom_point(alpha=0.4, size=1.75) +
  #设置点的颜色
  scale_color_manual(values =c("blue","black","red"))+
  #设置横竖阈值线（横：p为0.05，竖：FC的阈值）
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(logFC_cutoff,-logFC_cutoff),lty=4,lwd=0.6,alpha=0.8)+
  #更换背景主题
  theme_bw()+
  #更换背景格子
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  #加标题(见底部title)
  labs(title=title, x="log2 (fold change)",y="-log10 (q-value)")+
  #标题居中
  theme(plot.title = element_text(hjust = 0.5))+
  #设置图例
  #theme(legend.position='none')+
  #标出异常点并且避免重叠，还能换文字颜色
  geom_text_repel(data=subset(data, abs(log2FoldChange) > 5), 
                  aes(label=gene),col="black",alpha = 0.8)+
  geom_text(aes(x=4,y=75,label="Reed Liu !"),col ="orange",size = 6)

print(vol2)
ggsave(vol2,filename = 'DIY_Deseq2_volcano.png')

# title：
title <- paste0('log2FoldChange cutoff: ',round(logFC_cutoff,3),
                '\nUp-regulated genes: ',nrow(DEG[DEG$change =='UP',]) ,
                '\nDown-regulated genes: ',nrow(DEG[DEG$change =='DOWN',])
)
