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
deseq <- data.frame(gene=rownames(Deseq_DEG),
                    Deseq=Deseq_DEG$log2FoldChange)

limma <- data.frame(gene=rownames(limma_voom_DEG),
                    limma=limma_voom_DEG$logFC)

edgeR <- data.frame(gene=rownames(edgeR_DEG),
                    edgeR=edgeR_DEG$logFC)

test1 <- merge(deseq,limma,by="gene")
plot(test1[2:3])

test2 <- merge(deseq,edgeR,by="gene")
plot(test2[2:3])
