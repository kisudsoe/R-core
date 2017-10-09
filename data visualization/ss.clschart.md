---
title: "Yeast KO strain CLS plots"
author: "KimSS"
date: "2017-03-09"
---

# Version info

* v1.0  - 170309 THU for Mouse muscle aging (SY's cls data)
* v1.0a - 170728 for Mouse CR liver (DE's cls data)
* v1.0b - 170908, Add title (Mouse CR liver)

```r
clsplot = function(cls, x, strn, f_name, width, height) {
  date_n=subset(cls,strain%in%strn)$date

  # data filtering---------
  flt = subset(cls, strain%in%c("WT-CR","WT-NR",strn) & date%in%date_n)
  flt_info = flt[,c(1:5,25)]
  flt_cls = flt[,6:24]
  print(subset(flt_info,Note!="NA"))
  #print(flt_cls)

  # stat calculation---------
  flt_mn = aggregate(flt_cls,list(flt_info$strain,flt_info$date),mean)
  colnames(flt_mn) = c("strain","date",x)
  flt_sd = aggregate(flt_cls,list(flt_info$strain,flt_info$date),sd)
  st.err = function(x) sd(x)/sqrt(length(x))
  flt_se = aggregate(flt_cls,list(flt_info$strain,flt_info$date),st.err)
  colnames(flt_se) = c("strain","date",x)
  #rint(flt_mn)

  # data rearrangement---------
  library(reshape2)
  input = data.frame(melt(flt_mn,id=c("strain","date")),
                     melt(flt_se,id=c("strain","date"))$value)
  colnames(input) = c("strain","date","day","viability","se")
  input = transform(input,day=as.numeric(day)) # day as numeric
  input = subset(input,viability!="NA")
  #print(input)

  # graph visualization---------
  library(ggplot2)
  g = ggplot(input,aes(x=day,y=viability,color=strain,group=strain))+
    scale_y_continuous(limits=c(0,100))+
    geom_errorbar(aes(ymin=viability-se,ymax=viability+se),width=.5)+
    geom_jitter(aes(shape=strain),size=3,width=0)+
    geom_line(aes(linetype=strain),size=1)+
    ggtitle(f_name)+
    #scale_color_manual(values=c("gray","black"))+
    theme_bw()+facet_grid(.~date)
  print(g)
  if(!is.null(f_name)) {
    ggsave(paste0(f_name,".png"),plot=g,width=width,height=height)
    dev.off()
    cat("Draw cls plot done.\n")
  }
}
```
```r
# Data load
cls = read.csv("yeast od_170908.csv"); print(dim(cls))
colnames(cls)[1:5] = c("date","id","media","gene","strain") # by plot by date and strain
x = c(0,3,6,9,12,14,15,16,18,19,22,24,25,27,32,34,39,41,42)
print(colnames(cls)); head(cls)
strains = levels(cls$strain); print(strains)
# Function execution
d = 170908; strain = strains[3:4]; print(strain)
f_name = paste0(d,"_",strain)[1]; print(f_name)
clsplot(cls,x,strain,f_name, width=6,height=5)
```
