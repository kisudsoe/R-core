---
Title: read_write.md
Date: 2017-07-15 SAT
Version: v3.1
Note : First written of this file was at 2013
---

* `170715` v3.0  - Update install packages
* `171014` v3.1  - Update keyboard input



# R update

## R update to latest version

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
```



# Data read/write

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



# ProgressBar

### winProgressBar

* Updated at 2017

```r {.lineNo}
pdtime = function(time) {
    t=Sys.time()
    d=difftime(t,time,unit="sec")
    if(d<60) d=paste(round(d,1),"sec")
    else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
    else if(d>=3600) d=paste(round(d/3600,1),"hr")
    out = paste("Job done:",t,"for",d)
    return(out)
}

t=Sys.time(); n=length(iter)
pb=winProgressBar(title="Progress",label="description",min=0,max=n,width=500)
for(i in 1:n) {
	## Progress time ##
	setWinProgressBar(pb,i,title=paste0(round(i/n*100,1)," % (",i,"/",n,") done for ",pdtime(t)))
	###################
}
close(pb); print(pdtime(t))
```
