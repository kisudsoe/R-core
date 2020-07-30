---
title: "Gene set analysis using Piano library"
author: "KimSS"
date: "2017-12-21 THU"
note: Original file at `2016.04 Yeast HD.LD.Rho0\Gene set analysis.ipynb`
---

* https://bioconductor.org/packages/release/bioc/html/piano.html

```r
source("https://bioconductor.org/biocLite.R")
biocLite("piano")
library(piano)
?piano
?runGSA
```

## 1. Load category data

```r
cat = read.csv("(Archive) R process/171214_GSA_Categories.csv")
cat$Cat3 = paste(cat$Tukey3,cat$Cat0, sep="_")
cat$Cat4 = paste(cat$Tukey3,cat$Cat, sep="_")
print(dim(cat)); cat[sample(nrow(cat),5),]
data.frame(table(cat$Cat0))
```

## 2. Permutation test using lmPerm and coin packages

* http://rcompanion.org/handbook/K_01.html

### 2-1. Binary table generation

* See below `ss.permutation` function
```r
install.packages('lmPerm', repos="http://cran.us.r-project.org")
install.packages('coin', repos="http://cran.us.r-project.org")
install.packages('tictoc', repos="http://cran.us.r-project.org")
```
* Try permutation tests
```r
# parameters
n = 1; m = 19 # n= tukey; m= category
val = cat_val_df[,n]; nm = factor(cat_cat_df[,m])
table(val)
table(nm)

library(lmPerm)
summary(aovp(val~nm))
cat("===========================================\n")
library(coin)
independence_test(val~nm) # test as independent samples
#symmetry_test(val~nm) # test as paired data
```
* Permutation operation function
```r
ss.permutation = function(cat=NULL,cat_option=NULL) {
    library(tictoc)
    tic("tukey binary table generation")
    if(cat_option=="Cat0") idx = cat$Cat0
    else if(cat_option=="Cat") idx = cat$Cat

    # tukey binary table
    cats = levels(idx); tuk = levels(cat$Tukey3)
    print(paste0("n=",length(idx),", levels=",length(cats)))
    cat_tuk = NULL
    for(i in 1:length(tuk)) {
        id = unlist(tuk[i])
        cat_tuk = cbind(cat_tuk,id[match(cat$Tukey3,id)])
    }
    colnames(cat_tuk) = tuk
    cat_tuk_df = (cat_tuk!="") # If value exist, then TRUE
    cat_tuk_df[!is.na(cat_tuk)] = 1 # If value exist, then 1
    cat_tuk_df[is.na(cat_tuk)] = 0 # If value is null, then 0
    rownames(cat_tuk_df) = idx #cat$Cat0, cat$Tukey3
    print(paste(c("nrow=","ncol="),dim(cat_tuk_df)))
    toc()

    # category binary table
    tic("category binary table generation")
    cat_cat = NULL
    for(i in 1:length(cats)) {
        cat_id = unlist(cats[i])
        cat_cat = cbind(cat_cat,cat_id[match(idx,cat_id)])
    }
    colnames(cat_cat) = cats
    cat_cat_df = cat_cat
    cat_cat_df[is.na(cat_cat)] = "Others"
    rownames(cat_cat_df) = idx
    print(paste(c("nrow=","ncol="),dim(cat_cat_df)))
    toc()

    # purmutation test
    tic("permutation test")
    library(coin)
    out = NULL
    for(i in 1:length(tuk)) {
        val = cat_tuk_df[,i]
        pval = apply(cat_cat_df,2,function(col) {
            nm = factor(col)
            return(pvalue(independence_test(val~nm,
                                            #tetstat="maximum",
                                            distribution=approximate(B=1000))))
        })
        out = cbind(out,pval)
    }
    colnames(out) = tuk
    rownames(out) = cats
    print(paste(c("nrow=","ncol="),dim(out)))
    toc()
    return(out)
}
```
* Perform the permutation test
```r
cat_pval = ss.permutation(cat=cat,cat_option="Cat0")
cat_pval2 = ss.permutation(cat=cat,cat_option="Cat")
```
* Cat0 result arrangement and save
```r
cat_pval_ = cat_pval[c(19,11,9,1,10, 2,5,12,6, 3,14,16,15,20, 17,13,18,4,8, 7,21),
                     c(12,4,6,3,11,5, 2,7,9,1,10,8)]
print(paste(c("nrow=","ncol="),dim(cat_pval_)))
head(cat_pval_,10)
write.csv(cat_pval_,"171218 cat_permut_pval.csv")
```

## 3. Permutation test using piano package

* Using **runGSA** function
* http://www.bioconductor.org/packages/3.7/bioc/vignettes/piano/inst/doc/piano-vignette.pdf
```r
library(piano)
cat_ = cat$Cat4 # parameter

genes2genesets = data.frame(gene=cat$AFFYid,cat=cat_)
#head(genes2genesets)
gsumm = data.frame(table(genes2genesets$cat))
#gsumm_ = subset(gsumm,Freq>4) # filter the data
head(gsumm); print(paste("nrow=",nrow(gsumm)))

myGsc = loadGSC(genes2genesets); #myGsc
```
### GSA options
* https://www.rdocumentation.org/packages/piano/versions/1.12.0/topics/runGSA
* If geneSetStat is set to `"fisher"`, `"stouffer"`, `"reporter"` or "`tailStrength" `only p-values are allowed as geneLevelStats.
* If geneSetStat is set to `"maxmean"`, `"gsea"` or `"page"` only t-like `geneLevelStats` are allowed (e.g. t-values, fold-changes).
* For geneSetStat set to `"fisher"`, `"stouffer"`, `"reporter"`, `"wilcoxon"` or `"page"`, the gene set p-values can be calculated from a theoretical null-distribution, in this case, set `signifMethod="nullDist"`.
* For geneSetStat set to `"fisher"`, `"stouffer"`, `"reporter"`, `"wilcoxon"` or `"page"`, the gene set p-values can be calculated from a theoretical null-distribution, in this case, set `signifMethod="nullDist"`.
* For all methods `signifMethod="geneSampling"` or `signifMethod="samplePermutation"` can be used.
```r
myStats = cat$HD.LD
myStats2= cat$LD.Rho
names(myStats) = cat$AFFYid
names(myStats2)= cat$AFFYid

library(tictoc)
tic("runGSA - HD/LD")
gsaRes_hdld = runGSA(myStats,geneSetStat="gsea",
                     signifMethod="geneSampling",
                     gsc=myGsc,
                     #gsSizeLim=c(5,Inf),
                     nPerm=10000)
toc()
tic("runGSA - LD/Rho")
gsaRes_ldrho= runGSA(myStats2,geneSetStat="gsea",
                     signifMethod="geneSampling",
                     gsc=myGsc,
                     #gsSizeLim=c(5,Inf),
                     nPerm=10000)
toc()
```
* print runGSA performance
```r
gsaRes_hdld
#print(names(gsaRes))
```
* Extract pvalues from the results
```r
print(head(cat$Cat4,4))
#gs = geneSetSummary(gsaRes,"a.a.b_Protein catabolism"); gs
result_hdld = GSAsummaryTable(gsaRes_hdld)
result_ldrho= GSAsummaryTable(gsaRes_ldrho)
#GSAsummaryTable(gsaRes, save=T, file="171218 piano.tsv")
```
* Gathering the results
```r
out = NULL
out[[1]] = paste(c("nrow=","ncol="),dim(result_hdld))
out[[2]] = result_hdld
head(out[[2]],12)
out[[3]] = paste(c("nrow=","ncol="),dim(result_ldrho))
out[[4]] = result_ldrho
head(out[[4]],12)
```
* Extract pval to summary from raw results
```r
#install.packages('reshape', repos="http://cran.us.r-project.org")
result_tb = NULL
for(i in 1:nrow(result_hdld)) {
    row_hdld = result_hdld[i,]
    row_ldrho= result_ldrho[i,]
    name_hdld = unlist(strsplit(row_hdld$Name,"_"))
    name_ldrho= unlist(strsplit(row_ldrho$Name,"_"))
    if(name_hdld[1]%in%c("a.b.c","a.b.b","a.c.b")) pval = row_hdld[7]
    else if(name_hdld[1]%in%c("c.b.a","b.a.a","c.a.b")) pval = row_hdld[5]
    else if(name_ldrho[1]%in%c("b.a.b","a.a.b","b.a.c")) pval = row_ldrho[7]
    else if(name_ldrho[1]%in%c("a.b.a","b.b.a","b.c.a")) pval = row_ldrho[5]

    row = c(unlist(name_hdld[2]),unlist(name_hdld[1]),unlist(pval))
    result_tb = rbind(result_tb,row)
}
colnames(result_tb) = c("cat","tuk","pval")
rownames(result_tb) = c(1:nrow(result_tb))
result_tb = data.frame(result_tb, stringsAsFactors = FALSE)
cat("result_tb >> ")
print(paste(c("nrow=","ncol="),dim(result_tb)))
head(result_tb)

library(reshape)
result_df = cast(result_tb,cat~tuk,fun.aggregate=min,value="pval") #<- default "length"를 피하는 법
cat("result_df >> ")
print(paste(c("nrow=","ncol="),dim(result_df)))
head(result_df)
```
* Edit and Rearrange the summary table
```r
result_df2 = result_df[,-1]
rownames(result_df2) = result_df[,1]
print(colnames(result_df2)); print(head(rownames(result_df2)),12)
print(paste(c("nrow=","ncol="),dim(result_df2)))
## option 1
#result_df_ = result_df2[c(19,11,9,1,10, 2,5,12,6, 3,14,16,15,20, 17,13,18,4,8, 7,21),c(12,4,6,3,11,5, 2,7,9,1,10,8)]
#rownames(result_df_) = rownames(result_df2)[c(19,11,9,1,10, 2,5,12,6, 3,14,16,15,20, 17,13,18,4,8, 7,21)]
## option 2
result_df_ = result_df2[,c(12,4,6,3,11,5, 2,7,9,1,10,8)]
rownames(result_df_) = rownames(result_df2)

head(result_df_)
out[[5]] = result_df_
print(head(out[[2]][,1],12))
```
* Save as csv file
```r
lapply(out,function(x)
    write.table(data.frame(x),"171220 piano_gsea_subcat_fdr.csv",append=T,sep=","))
```
