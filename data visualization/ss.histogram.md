---
title: "ss.historam"
author: "KimSS"
date: "2017년 4월 12일"
output: html_document
note: This document was made for Mice CR Liver Project.
---

# 180816 THU - ggplot

Correlation R distribution histogram

```R
cr_str = c(0,0,0,15,15,15,30,30,30,45,45,45)
cor_r = apply(probe[,2:13],1,function(x) {
    cor(x,cr_str,method="pearson")
})
cor_r = data.frame(Corr=cor_r)
print(dim(cor_r)); summary(cor_r)

library(ggplot2)
ggplot(cor_r, aes(Corr))+theme_bw()+
    geom_histogram(breaks=seq(-1,1,by=0.05),
                   col="black",fill="gray",
                   show.legend=F)
ggsave(paste(date,"corr_hist.png"),width=10,height=7,units="cm")
```





# 170412 WED - ggplot

## Input AffyEC_report

```R
chp_EC = read.delim("AffyEC_Report_Probe signal.tsv")
group = c("CD","CR85","CR55","CR70",
          "CD","CR55","CR70","CR85",
          "CD","CR55","CR70","CR85")
group = factor(group, levels=c("CD","CR85","CR70","CR55"))
colnames(chp_EC) = c("AffyID",as.vector(group))
```

## Histogram

```R
library(reshape)
tmp1 = subset(chp_EC,AffyID%in%Affyid_anti)
tmp2 = subset(chp_EC,AffyID%in%Affyid_intron)
tmp3 = subset(chp_EC,AffyID%in%Affyid_exon)
tmp = rbind(tmp1,tmp2,tmp3)
tmp$ann = c(rep("Antigenomic",length(Affyid_anti)), rep("Intron",length(Affyid_intron)), rep("Exon",length(Affyid_exon)))

tmp_melt = melt(tmp,id.vars=c("AffyID","ann"))
tmp_melt_CD = subset(tmp_melt,variable=="CD")

library(ggplot2)
ggplot(tmp_melt_CD,aes(value,fill=ann,colour=ann))+
  geom_histogram(aes(y=..density..), binwidth=1, alpha=.7)+
  ggtitle("Histogram of antigenomic signal")

rm(tmp1,tmp2,tmp3,tmp,tmp_melt,tmp_melt_CD)
```
