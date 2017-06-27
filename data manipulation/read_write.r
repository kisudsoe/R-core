## 1. Data input ##
# First written 2013
# Last update 2017-03-10 FRI
# Version info 2.1
# Note : In case of "clipboard", you make sure :
#    1) no space between words
#    2) 'NA' for empty values.

## R update to v3.3.1
install.packages("installr")
library(installr)
updateR()

## Read file
read.delim('filename.txt', header=TRUE) # Read txt file delimited by tab.
read.csv('filename.csv',header=TRUE) 	# Read csv file.

## Write file
write.table(fdr,"fdr.csv",sep=",",row.names=TRUE) # save into your workspace
lapply(pca.AST,function(x) write.table(data.frame(x),"pca.AST.csv",append=T,sep=",")) # save list as csv file

## winProgressBar
time1=Sys.time()
pb=winProgressBar(title="Progress",label="description",min=0,max=3000,width=500)
for(i in 1:300) {
	setWinProgressBar(pb,i,title=paste(round(i/3000*100,1),"% done for",round(Sys.time()-time1,1),"sec"))
}
