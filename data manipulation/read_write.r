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
