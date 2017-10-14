---
Title: read_write.md
Date: 2017-07-15 SAT
Version: v3.1
Note : First written of this file was at 2013
---

* v3.0  - 2017-07-15 SAT, Update install pakcages
* v3.1  - 2017-10-14 SAT, Update keyboard input

## 1. Data input

### R update to v3.3.1

* R update in R command

```r
## R session
install.packages("installr")
library(installr)
updateR()

## Jupyter R kernel (in CMD)
conda update --all
```

### Install packages

* Install packages in Jupyter kernel (in Atom)

```r
# CRAN
install.packages('multicompView', repos="http://cran.us.r-project.org")

# Bioconductor
source("https://bioconductor.org/biocLite.R"); biocLite("limma")

# Using RDocumentation packages
install.packages('RDocumentation', repos="http://cran.us.r-project.org")
library(RDocumentation)
install_package("pbapply", 1) # type 1: CRAN; type 2: BioConductor; type 3: GitHub; type 4: default part of R
```

### Read from keyboard input
```r
a=readLines(con=file("stdin"))
```

### Read file
```r
read.delim('filename.txt', header=TRUE) # Read txt file delimited by tab.
read.csv('filename.csv',header=TRUE) 	# Read csv file.
```

### Write file
```r
write.table(fdr,"fdr.csv",sep=",",row.names=TRUE) # save into your workspace
lapply(pca.AST,function(x)
       write.table(data.frame(x),"pca.AST.csv",append=T,sep=",")) # save list as csv file
```

### winProgressBar
* Updated at 2017

```r {.lineNo}
time1=Sys.time(); n=length(iter)
pb=winProgressBar(title="Progress",label="description",min=0,max=n,width=500)
for(i in 1:300) {
	## Progress time ##
	d=difftime(Sys.time(),time1,unit="sec")
	if(d<60) d=paste(round(d,2),"sec")
	else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
	else if(d>=3600) d=paste(round(d/3600,1),"hr")
	setWinProgressBar(pb,i,title=paste0(round(i/n*100,1)," % (",i,"/",n,") done for ",d))
	###################
}
print(difftime(Sys.time(),time1,units="auto")); close(pb)
```
