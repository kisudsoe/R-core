---
title:  "ss.circos.md"
by:     "Kim SS"
date:   "2017-06-27 TUE"
note:   "This file was first written for 2017 Spring Aging Consorcium oral presentation"
---

# TOC
[TOC]


## 2017-07-15 SAT

* Chord diagram by circlize
* See: http://zuguang.de/circlize_book/book

```r
#install.packages('RDocumentation', repos="http://cran.us.r-project.org")
library(RDocumentation)
install_package("circlize",1)
```
```r
library(circlize)
# later..
```


## 2017-06-27 TUE

* This code was first written for 2017 Spring Aging Consorcium oral presentation
* Chord diagram by OmicCircos package in Bioconductor

```r
source("https://bioconductor.org/biocLite.R")
biocLite("OmicCircos")
options(stringsAsFactors=F)
library(OmicCircos)
```

* Draw OmicCircos for gene categories

```r
f = read.csv("seg.frame.csv"); print(dim(f))
a = read.csv("seg.arc.csv")[,1:5]; print(dim(a))
l = read.csv("seg.link.csv")[1:36,]; print(dim(l))
l.pg = read.csv("seg.link.pg.csv")[1:36,]; print(dim(l.pg))

n = levels(as.factor(f[,1])); print(length(n))
db = segAnglePo(seg.dat=f, seg=n)
colors.ss = rainbow(15)
colors.hi = rainbow(25)
colors = c(colors.ss,"white",colors.hi,"white")
```
```r
pdf("OmicCircos.pdf",8,8)
plot(c(1,800),c(1,800), type="n", axes=F, xlab="", ylab="", main="")
circos(R=300, cir=db, type="chr", col=colors, print.chr.lab=T, W=20, scale=F)
#circos(R=280, cir=db, W=40, mapping=a, type="arc2", B=F, col=colors, lwd=10, cutoff=0) # 안된다?!
circos(R=290, cir=db, W=40, mapping=l, type="link", lwd=2, col=c(rep("red",8),rep("blue",28))) # link.pg가 구현이 안된다?!
dev.off()
```
