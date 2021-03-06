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

## Exercise 5: Survival Analysis
# Use http://www.oncolnc.org/ to get raw csv da
rm(list=ls())
options(stringsAsFactors = F)

f <- read.csv('e5-BRCA_7157_50_50.csv')
colnames(f)
da <- f
#################################
## first boxplot
#################################
library(ggstatsplot)
ggbetweenstats(data = da,
               x = Group,
               y = Expression)
#################################
## second survival analysis
#################################
library(ggplot2)
library(survival)
library(survminer)
table(da$Status)
da$Status <- ifelse(da$Status == "Dead", 1, 0)
survf <- survfit(Surv(Days,Status)~Group, data=da)
# survf
# summary(survf)
# simple survplot
ggsurvplot(survf, conf.int = F, pval = T)

# complex survplot
ggsurvplot(survf,palette = c("orange", "navyblue"),
           risk.table = T, pval = T,
           conf.int = T, xlab = "Time of days",
           ggtheme = theme_light(),
           ncensor.plot = T)
ggsave("survival_TP53_in_BRCA_taga.png")
