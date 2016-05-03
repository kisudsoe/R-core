## 1. Data input ##
# First written 2013
# Last update 2016-04-29 FRI
# Version info 2.0a
# Note : In case of "clipboard", you make sure :
#    1) no space between words
#    2) 'NA' for empty values.

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

read.delim('filename.txt', header=TRUE) # Read txt file delimited by tab.
read.csv('filename.csv',header=TRUE) 	# Read csv file.

write.table(fdr,"fdr.csv",sep=",",row.names=TRUE) # save into your workspace