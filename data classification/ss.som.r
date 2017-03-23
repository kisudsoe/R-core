## Self-Organizing Map (SOM) ##
install.packages("som")
library(som)
data.s = t(scale(t(data))) # row-wise scale : center to zero
data.som = som(data.s, xdim=4, ydim=4, topol="rect", neigh="gaussian")
	# Performs SOM clustering
	# topol = "hexa", "rect"

plot(data.som)
#data.som.group = data.som$visual$x*10+data.som$visual$y+1

# v1.0  - 160104, somgrp function created
# v1.1  - 161115, change x-symbol to number for over 24 dimension support
ss.somgrp = function(data.som,i,j) { # 161115 ver1.1
  x = data.som$visual$x
  y = data.som$visual$y
  xax = c(1:i)
  yax = c(1:j)
  xlen = length(x)

  out = NULL
  for(l in 1:xlen) {
    grpx = xax[x[l]+1]
    grpy = yax[y[l]+1]
    if(grpx>=10&&grpy<10) pas=paste0(grpx,'x0',grpy)
    else if(grpx<10&&grpy>=10) pas=paste0('0',grpx,'x',grpy)
    else if(grpx>=10&&grpy>=10) pas=paste0(grpx,'x',grpy)
    else pas = paste0(grpx,grpy)
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

#-----------------
# 2016-11-02 WED
# Ver 1.0
ss.som.sdbox = function(data,range=c(3:6),topol="rect", neigh="gaussian") {
  data.s = t(scale(t(data)))
  n = length(range)
  for(k in 1:n) {
    i = as.numeric(range[k])
    data.som = som(data.s,xdim=i,ydim=i,topol=topol,neigh=neigh)
    x = data.som$visual$x
    y = data.som$visual$y
    xax = LETTERS[1:i]
    yax = c(1:i)

    som.group = NULL
    for(l in 1:length(x)) {
		  if(yax[y[l]+1]>=10) pas = paste0(xax[x[l]+1],yax[y[l]+1])
		  else pas = paste0(xax[x[l]+1],'0',yax[y[l]+1])
		som.group = c(som.group,pas)
    }
    som.group = as.factor(som.group)
    group = levels(som.group)
    #cat("i=",i,', ',"SOM groups =",length(group),'\n')

    num = sapply(group,function(el) {
      id = which(som.group==el)
      length(id)
    })
    sd = sapply(group,function(el) {
      id = which(som.group==el)
      if(length(id)>3) a = sd(unlist(data[id,]))
      else a = NA
      return(a)
    })
    if(k==1) {
      som.sd = data.frame(Dimension=factor(rep(length(group),length(sd))),SD=sd)
      som.num = data.frame(Dimension=factor(rep(length(group),length(num))),Number=num)
    } else {
      sd = data.frame(Dimension=factor(rep(length(group),length(sd))),SD=sd)
      num = data.frame(Dimension=factor(rep(length(group),length(num))),Number=num)
      som.sd = rbind(som.sd,sd)
      som.num = rbind(som.num,num)
    }

    #######################
    # Create progress bar #
    #######################
    if(k==1) { # set progress bar
      cat('\nProcess iteration =',n,'\n')
      pb = txtProgressBar(min=0,max=100,width=30,initial=0,style=3)
      time1 = Sys.time()
    } else { setTxtProgressBar(pb, (k/n)*100) } # show progress bar
    if(k==n) { # Duration time check
      close(pb)
      time2 = Sys.time()
      cat('\n')
      print(time2-time1)
    }
    #######################
  }
  out = list(som.sd,som.num)
  return(out)
}
# this function take long time.
som_dtm = ss.som.sdbox(data=cerev.sgl[-1],range=c(2:32))
