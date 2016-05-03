# 2015-01-06 WED
ss.somline = function(data,data.som,i=4,j=4) {
 
	# 1-1.Make SOM group list as factor
	x = data.som$visual$x
	y = data.som$visual$y
	len = length(x)
	xax = LETTERS[1:i]
	yax = c(1:j)
	for(l in 1:len) {
		pas = paste(xax[x[l]+1],yax[y[l]+1],sep="")
		# sep="" for no space-interuption between characters
		if(l==1)
			som.group = pas
		else
			som.group = c(som.group,pas)
	}
	som.group = as.factor(som.group)
	print("Level of SOM groups =")
	print(levels(som.group))
 
	# 1-2.Make total SOM group
	for(m in 1:i) {
		for(n in 1:j) {
			pas2 = paste(xax[m],yax[n],sep="")
			if(m==1&&n==1)
				group = pas2
			else
				group = c(group,pas2)
		}
	}
	group = as.factor(group)
 
	# 2.Draw plot using SOM group list
	#dev.new(with=10,height=10)
	par(mfrow=c(i,j)) # c(row,column)
	n = length(levels(group))
	print(paste("plot number=",n))
 
	for(k in 1:n) {
		sel = levels(group)[k]
		id = which(som.group==sel)
		title=paste(sel," (n=",length(id),")",sep="")
 
		if(k%%(j)==1) {
			par(mar=c(0.5,2,1.5,0))
			plot(as.numeric(data[id[1],]),
				main=title,
				col="gray50",
				type="l",
				xaxt="n",
				ylim=c(-3,3))
		} else {
			par(mar=c(0.5,1,1.5,0))
			plot(as.numeric(data[id[1],]),
				 main=title,
				 col="gray50",
				 type="l",
				 xaxt="n",
				 yaxt="n",
				 ylim=c(-3,3))
		}
 
		if(length(id)!=0) {
			rect(3.5,-3.1,6.5,3.1,col="gray90",lty=0)
			rect(9.5,-3.1,12.5,3.1,col="gray90",lty=0)
			rect(15.5,-3.1,18.5,3.1,col="gray90",lty=0)
			rect(21.5,-3.1,24.5,3.1,col="gray90",lty=0)
			rect(27.5,-3.1,30.5,3.1,col="gray90",lty=0)
			text(2,-2.8,"14")
			text(5,-2.8,"29")
			text(8,-2.8,"44")
			text(11,-2.8,"53")
			text(14,-2.8,"73")
			text(17,-2.8,"83")
			text(20,-2.8,"95")
			text(23,-2.8,"109")
			text(26,-2.8,"120")
			text(29,-2.8,"135")
		} else if(length(id)==0) {
			text(15,0,"None",cex=2)
		}
 
		for(l in 1:length(id)) {
			lines(as.numeric(data[id[l],]),col="gray50")
		}
 
		if(length(id)>3) {
			av = apply(data[id,],2,mean)
			lines(as.numeric(av),col="red",lwd=2)
		}
 
		#boxplot(data[id,],xaxt='n',yaxt='n',add=TRUE)
    
		#legend("bottomleft",title)
		#text(3,-2.5,title,font=5)
	}
}
 
ss.somline(data209.s,data209.som,6,6)