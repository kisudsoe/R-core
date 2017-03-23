# 2016-11-02 WED for Yeast dynamic netowrk project
ss.somline2 = function(data=NULL,data.som=NULL,group=NULL,
                      i=4,j=4,pdf=F,clstr=NULL) {

	# 1-1.Make SOM group list as factor
	x = data.som$visual$x
	y = data.som$visual$y
	len = length(x)
	xax = LETTERS[1:i]
	yax = c(1:j)
	som.group = NULL
	for(l in 1:len) {
		if(yax[y[l]+1]>=10) pas = paste0(xax[x[l]+1],yax[y[l]+1])
		else pas = paste0(xax[x[l]+1],'0',yax[y[l]+1])
		som.group = c(som.group,pas)
	}
	som.group = as.factor(som.group)
	cat("Level of SOM groups =",length(levels(som.group)),'\n')

	# 1-2.Make total SOM group
	for(m in 1:i) {
		for(n in 1:j) {
			if(yax[n]>=10) pas2 = paste0(xax[m],yax[n])
			else pas2 = paste0(xax[m],'0',yax[n])
			if(m==1&&n==1)
				group = pas2
			else
				group = c(group,pas2)
		}
	}
	group = as.factor(group)

	# 2.Draw plot using SOM group list
	#dev.new(with=10,height=10)
	if(length(clstr)==0) {
	  if(pdf) pdf("plot.pdf", width=i*2, height=j*2)
	  par(mfrow=c(i,j)) # c(row,column)
	  n = length(levels(group))
	  cat("plot number=",n,'\n')

	  for(k in 1:n) {
		  sel = levels(group)[k]
		  id = which(som.group==sel)
		  #cat('k=',k,', id length=',length(id),'\n')
		  title=paste(sel," (n=",length(id),")",sep="")
		  n_id = which(group.info$glucose[-1]==2)
		  cr1_id = which(group.info$glucose[-1]==1)
		  cr2_id = which(group.info$glucose[-1]==0.5)
		  cr2_id = c(cr2_id[16:18],cr2_id[1:15])
		  cr3_id = which(group.info$glucose[-1]==0.25)
		  rapa_id = c(n_id[16:18])
		  rho_id = c(n_id[22:26])
		  HD_id = c(n_id[27:31])
		  LD_id = c(n_id[32:36])
		  nr_id = c(n_id[19:21],n_id[1:15],n_id[27:31])
		  dat = as.numeric(data[id[1],nr_id])

		  if(k%%(j)==1) {
			  par(mar=c(0.5,2,1.5,0))
			  plot(dat,
				  main=title,
				  col="white",
				  type="l",
				  xaxt="n",
				  ylim=c(-3.5,3.5))
		  } else {
			  par(mar=c(0.5,1,1.5,0))
			  plot(dat,
				  main=title,
				  col="white",
				  type="l",
				  xaxt="n",
				  yaxt="n",
				  ylim=c(-3.5,3.5))
		  }

		  if(length(id)!=0) {
		    rect(3.5,-3.7,6.5,3.7,col="gray90",lty=0)
		    rect(9.5,-3.7,12.5,3.7,col="gray90",lty=0)
		    rect(15.5,-3.7,18.5,3.7,col="gray90",lty=0)

		    for(l in 1:length(id)) {
	  		  lines(as.numeric(data[id[l],nr_id[1:18]]),col="gray60")
		      lines(as.numeric(data[id[l],cr1_id]),col="gray70")
		      lines(as.numeric(data[id[l],cr2_id]),col="gray70")
		      lines(as.numeric(data[id[l],cr3_id]),col="gray70")
		      lines(x=c(16:18),
		            y=as.numeric(data[id[l],rapa_id]),col="gray50")
		      lines(x=c(19:23),
		            y=as.numeric(data[id[l],HD_id]),col="gray60")
		      lines(x=c(19:23),
		            y=as.numeric(data[id[l],LD_id]),col="gray70")
		      lines(x=c(19:23),
		            y=as.numeric(data[id[l],rho_id]),col="gray50")
		    }
	  	}

		  colr = rainbow(7)
		  if(length(id)>3) {
			  lines(as.numeric(apply(data[id,nr_id[1:18]],2,mean)),
		          col=colr[1],lwd=2)
		    lines(as.numeric(apply(data[id,cr1_id],2,mean)),
			        col=colr[2],lwd=2)
		    lines(as.numeric(apply(data[id,cr2_id],2,mean)),
			        col=colr[3],lwd=2)
		    lines(as.numeric(apply(data[id,cr3_id],2,mean)),
			        col=colr[4],lwd=2)
		    lines(x=c(16:18),
		          y=as.numeric(apply(data[id,rapa_id],2,mean)),
			        col=colr[5],lwd=2)
		    lines(x=c(19:23),
		          y=as.numeric(apply(data[id,HD_id],2,mean)),
			        col=colr[5],lwd=2)
		    lines(x=c(19:23),
		          y=as.numeric(apply(data[id,LD_id],2,mean)),
			        col=colr[6],lwd=2)
		    lines(x=c(19:23),
		          y=as.numeric(apply(data[id,rho_id],2,mean)),
			        col=colr[7],lwd=2)
		  }

		  if(length(id)!=0) {
			  text(2,-3.2,"4,5")
			  text(5,-3.2,"12")
			  text(8,-3.2,"18")
			  text(11,-3.2,"24")
			  text(14,-3.2,"48")
			  text(17,-3.2,"120")
			  text(21,-3.2,"72")
		  } else if(length(id)==0) {
			  text(12.5,0,"None",cex=2)
		  }
	  }
	  if(pdf) dev.off()
	} else {
	  par(mfrow=c(1,1),mar=c(0.5,2,1.5,7),xpd=TRUE)
	  sel = levels(group)[clstr] # clstr=integer
	  cat('Draw cluster=',sel,'\n')
		id = which(som.group==sel)
		title=paste(sel," (n=",length(id),")",sep="")

		n_id = which(group.info$glucose[-1]==2)
		cr1_id = which(group.info$glucose[-1]==1)
		cr2_id = which(group.info$glucose[-1]==0.5)
		cr2_id = c(cr2_id[16:18],cr2_id[1:15])
		cr3_id = which(group.info$glucose[-1]==0.25)
		rapa_id = c(n_id[16:18])
		rho_id = c(n_id[22:26])
		HD_id = c(n_id[27:31])
		LD_id = c(n_id[32:36])
		nr_id = c(n_id[19:21],n_id[1:15],n_id[27:31])
		dat = as.numeric(data[id[1],nr_id])

		plot(dat, main=title,
				 col="white",
				 type="l",
				 xaxt="n",
				 ylim=c(-3.5,3.5))

		if(length(id)!=0) {
		  rect(3.5,-3.7,6.5,3.7,col="gray90",lty=0)
		  rect(9.5,-3.7,12.5,3.7,col="gray90",lty=0)
		  rect(15.5,-3.7,18.5,3.7,col="gray90",lty=0)

		  for(l in 1:length(id)) {
	  	  lines(as.numeric(data[id[l],nr_id[1:18]]),col="gray60")
		    lines(as.numeric(data[id[l],cr1_id]),col="gray70")
		    lines(as.numeric(data[id[l],cr2_id]),col="gray70")
		    lines(as.numeric(data[id[l],cr3_id]),col="gray70")
		    lines(x=c(16:18),
		          y=as.numeric(data[id[l],rapa_id]),col="gray50")
		    lines(x=c(19:23),
		          y=as.numeric(data[id[l],HD_id]),col="gray60")
		    lines(x=c(19:23),
		          y=as.numeric(data[id[l],LD_id]),col="gray70")
		    lines(x=c(19:23),
		          y=as.numeric(data[id[l],rho_id]),col="gray50")
  	  }
		}

		colr = rainbow(8)
		if(length(id)>3) {
			lines(as.numeric(apply(data[id,nr_id[1:18]],2,mean)),
		        col=colr[1],lwd=2)
		  lines(as.numeric(apply(data[id,cr1_id],2,mean)),
			       col=colr[2],lwd=2)
		  lines(as.numeric(apply(data[id,cr2_id],2,mean)),
			       col=colr[3],lwd=2)
		  lines(as.numeric(apply(data[id,cr3_id],2,mean)),
			       col=colr[4],lwd=2)
		  lines(x=c(16:18),
		        y=as.numeric(apply(data[id,rapa_id],2,mean)),
			       col=colr[5],lwd=2)
		  lines(x=c(19:23),
		        y=as.numeric(apply(data[id,HD_id],2,mean)),
			       col=colr[6],lwd=2)
		  lines(x=c(19:23),
		        y=as.numeric(apply(data[id,LD_id],2,mean)),
			       col=colr[7],lwd=2)
		  lines(x=c(19:23),
		        y=as.numeric(apply(data[id,rho_id],2,mean)),
			       col=colr[8],lwd=2)
		}

		if(length(id)!=0) {
			 text(2,-3.2,"4,5 hr")
			 text(5,-3.2,"12 hr")
			 text(8,-3.2,"18 hr")
			 text(11,-3.2,"24 hr")
			 text(14,-3.2,"48 hr")
			 text(17,-3.2,"120 hr")
			 text(21,-3.2,"72 hr")
			 legend(24,3.7,
		       c("NR","CR0.25","CR0.5","CR1.0",
		         "Rapa","HD","LD","Rho0"),
		       col=colr,lty=1,lwd=2)
		} else if(length(id)==0) {
			 text(12.5,0,"None",cex=2)
		}
	}
}

## someline matrix
ss.somline(data = data.s,
           data.som = yeast.cerev_som,
           group = group.info,
           a,a,pdf=T)
## someline single plot
ss.somline(data = data.s,
           data.som = yeast.cerev_som,
           group = group.info,
           a,a,clstr=1)


# 2015-01-06 WED
# 2016-11-01 TUE for Yeast dynamic network project
ss.somline = function(data,data.som,i=4,j=4) {

	# 1-1.Make SOM group list as factor
	x = data.som$visual$x
	y = data.som$visual$y
	len = length(x)
	xax = LETTERS[1:i]
	yax = c(1:j)
	som.group = NULL
	for(l in 1:len) {
		if(yax[y[l]+1]>=10) pas = paste0(xax[x[l]+1],yax[y[l]+1])
		else pas = paste0(xax[x[l]+1],'0',yax[y[l]+1])
		som.group = c(som.group,pas)
	}
	som.group = as.factor(som.group)
	cat("Level of SOM groups =",levels(som.group),'\n')

	# 1-2.Make total SOM group
	for(m in 1:i) {
		for(n in 1:j) {
			if(yax[n]>=10) pas2 = paste0(xax[m],yax[n])
			else pas2 = paste0(xax[m],'0',yax[n])
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
	cat("plot number=",n,'\n')

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
			rect(12.5,-3.1,45.5,3.1,col="gray90",lty=0)
			#rect(9.5,-3.1,12.5,3.1,col="gray90",lty=0)
			#text(2,-2.8,"14")
			#text(5,-2.8,"29")
		} else if(length(id)==0) {
			text(15,0,"None",cex=2)
		}

		for(l in 1:length(id)) {
			lines(as.numeric(data[id[l],]),col="gray50")
		}
		#boxplot(data[id,],xaxt='n',yaxt='n',add=TRUE)

		if(length(id)>3) {
			av = apply(data[id,],2,mean)
			lines(as.numeric(av),col="red",lwd=2)
		}
		#legend("bottomleft",title)
		#text(3,-2.5,title,font=5)
	}
}

ss.somline(data.s[,c(34:45,1:33,46:60)],yeast.cerev_som,14,14)
colnames(data.s[,c(34:45,1:33,46:60)])
