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

## Input function
input = function(i) {
	if(i == 0) {
		data = read.table("clipboard") # Data only. No colnames and rownames
	} else if(i == 1) {
		data = read.table("clipboard", header=T) # Set colnames
	} else if(i == 2) {
		data = read.table("clipboard", row.names=1) # No colnames. Set rownames
	} else if(i == 3) {
		data = read.table("clipboard", header=T, row.names=1) # Set colnames and rownames
	}
	return(data)
}

## Read file
read.delim('filename.txt', header=TRUE) # Read txt file delimited by tab.
read.csv('filename.csv',header=TRUE) 	# Read csv file.

## Write file
write.table(fdr,"fdr.csv",sep=",",row.names=TRUE) # save into your workspace
lapply(pca.AST,function(x) write.table(data.frame(x),"pca.AST.csv",append=T,sep=",")) # save list as csv file
