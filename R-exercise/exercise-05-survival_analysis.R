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
# Use http://www.oncolnc.org/ to get raw csv data
rm(list=ls())
options(stringsAsFactors = F)

a <- read.csv('e5-BRCA_7157_50_50.csv')
colnames(a)
dat <- a
## first boxplot
library(ggstatsplot)
ggbetweenstats(data = dat,
               x = Group,
               y = Expression)
## second survival analysis
library(ggplot2)
library(survival)
library(survminer)
table(dat$Status)
dat$Status <- ifelse(dat$Status == "Dead", 1, 0)
sfit <- survfit(Surv(Days,Status)~Group, data=dat)
# sfit
# summary(sfit)
# simple survplot
ggsurvplot(sfit, conf.int = F, pval = T)

# complex survplot
ggsurvplot(sfit,palette = c("orange", "navyblue"),
           risk.table = T, pval = T,
           conf.int = T, xlab = "Time of days",
           ggtheme = theme_light(),
           ncensor.plot = T)

