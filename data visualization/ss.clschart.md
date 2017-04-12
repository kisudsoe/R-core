---
title: "Yeast KO strain CLS"
author: "KimSS"
date: "2017년 3월 9일"
output: html_document
---

```{r setup, include=FALSE, id:"j1emgpby"}
knitr::opts_chunk$set(echo = TRUE)
```

```{r package attachment, id:"j1emgpbz"}
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape))
```

# 170309 THU

## 1. Load data

```{r id:"j1emgpc0"}
cls = read.csv("170111_yeast cls.csv")
```

## 2. Plot strains

```{r id:"j1emgpc3"}
clsplot = function(strn) {
  date_n=subset(str_date,cls.strain%in%strn)$cls.date
  x = c(0,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,19,20,21,23,26)

  #data filtering---------
  flt = subset(cls,strain%in%c("BY4741",strn)&date%in%date_n&sample!="'no stain")
  flt_info = flt[,c(1:4,26)]
  flt_cls = flt[,5:25]
  print(subset(flt_info,note!="NA"))
  #print(flt_cls)

  #stat calculation---------
  flt_mn = aggregate(flt_cls,list(flt_info$strain,flt_info$date),mean)
  colnames(flt_mn) = c("strain","date",x)
  flt_sd = aggregate(flt_cls,list(flt_info$strain,flt_info$date),sd)
  st.err = function(x) sd(x)/sqrt(length(x))
  flt_se = aggregate(flt_cls,list(flt_info$strain,flt_info$date),st.err)
  colnames(flt_se) = c("strain","date",x)
  #print(flt_mn)

  #data rearrangement---------
  input = data.frame(melt(flt_mn,id=c("strain","date")),
                     melt(flt_se,id=c("strain","date"))$value)
  colnames(input) = c("strain","date","day","viability","se")
  input = transform(input,day=as.numeric(day)) # day as numeric
  input = subset(input,viability!="NA")
  #print(input)

  #graph visualization---------
  gp = ggplot(input,aes(x=day,y=viability,color=strain,group=strain))+
    geom_errorbar(aes(ymin=viability-se,ymax=viability+se),width=.5)+
    geom_jitter(aes(shape=strain),size=3,width=0)+
    geom_line(aes(linetype=strain),size=1)+
    #scale_color_manual(values=c("gray","black"))+
    theme_bw()+facet_grid(.~date)
  print(gp)
}
clsplot(c("FPK1","SNF4"))
```
