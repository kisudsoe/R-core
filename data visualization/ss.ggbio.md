---
title: "ss.ggbio"
author: "KimSS"
date: "2016-09-23 FRI"
output: html_document
note: Codes in this document are imcomplete yet.
---

```{r id:"j1emfhpc"}
source("https://bioconductor.org/biocLite.R") # run only once
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggbio)
library(dplyr)
```

```{r id:"j1emfhpe"}
data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
model <- exonsBy(txdb, by = "tx")
model17 <- subsetByOverlaps(model, genesymbol["RBM17"])
exons <- exons(txdb)
exon17 <- subsetByOverlaps(exons, genesymbol["RBM17"])
## reduce to make sure there is no overlap
## just for example
exon.new <- reduce(exon17)
## suppose
values(exon.new)$sample1 <- rnorm(length(exon.new), 10, 3)
values(exon.new)$sample2 <- rnorm(length(exon.new), 10, 10)
values(exon.new)$score <- rnorm(length(exon.new))
values(exon.new)$significant <- sample(c(TRUE,FALSE), size = length(exon.new),replace = TRUE)
## data ready
exon.new
exon.new2 = exon.new[,1:7]
exon.df = tbl_df(exon.new)
```

```{r id:"j1emfhpe"}
p17 <- autoplot(txdb, genesymbol["RBM17"])
plotRangesLinkedToData(exon.new, stat.y = c("sample1", "sample2"), annotation = list(p17))
```
