# 2016-05-01 SUN----
## Tukey for post hoc test
ss.tukey = function(data,groups) { # function 'ss.tukey' start
	library(multcompView) # package multcompView install needed
	result = 0
	row = rownames(data)
	n = length(row)
	groups = as.vector(groups)
	out = NULL
  
	for(i in 1:n) { # 'for' statement start
		values = as.numeric(t(data)[,i])
    
		# post hoc tests - tukey #
		a1 = aov(values~groups)
		posthoc = TukeyHSD(x=a1, conf.level=0.95) # 95% confidence interval
		#plot(posthoc) # generate graph
		#text(0,10.7,cex=1.2,labels=row[i],xpd=TRUE) # ID on the graph
    
		# using multcompView package #
		ex.pval = extract_p(posthoc$groups)
		sigroup = multcompLetters(ex.pval)
    
		# stack iterative results
		out = rbind(out,t(data.frame(sigroup$Letters)))
    
		#######################
		# Create progress bar #
		#######################
		if (i==1) { # set progress bar
			cat('\nProcess iteration =',n,'\n')
			pb = txtProgressBar(min=0,max=100,width=30,initial=0,style=3)
			time1 = Sys.time()
		} else { setTxtProgressBar(pb, (i/n)*100) } # show progress bar
		if(i==n) { # Duration time check
			close(pb)
			time2 = Sys.time()
			cat('\nProcess done.\n')
			print(time2-time1)
		} 
		#######################
	} # 'for' statement end
	colnames(out) = rownames(data.frame(sigroup$Letters))
	rownames(out) = row
	return(out)
} # function 'ss.tukey' end

sgnl.cerev_tukey = ss.tukey(sgnl.cerev,group)