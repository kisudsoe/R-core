# 2016-01-04 | somgrp function created
somgrp = function(data.som,i,j) {
  x = data.som$visual$x
  y = data.som$visual$y
  xax = LETTERS[1:i]
  yax = c(1:j)
  xlen = length(x)
 
  for(l in 1:xlen) {
    grpx = xax[x[l]+1]
    grpy = yax[y[l]+1]
    pas = paste(grpx,grpy,sep="") 
      # sep="" for no space-interuption between characters
    
    if(l==1)
      out = pas
    else
      out = c(out,pas)
  }
  return(out)
}
dadta209.somGroup = somgrp(data209.som,6,6)
 
 
#result.som$visual 
  # Provides the assignment of rows items to the SOM clusters. 
  # Remember that the coordinate counting starts here at zero!

# Export data option
write.table(result.som2$visual,"som.visual.csv",sep=",",row.names=TRUE)
write.table(result.som2$data,"som.data.csv",sep=",",row.names=TRUE)
write.table(result.som$code.sum,"som.sum.csv",sep=",",row.names=TRUE)