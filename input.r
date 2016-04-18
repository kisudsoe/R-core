## 1. Data input ##
# First written 2013
# Last update 2014-10-09 THR
# Version info 2.0a

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