# 2016-05-01 SUN
## ver 1.0	- 160501, Yeast HD LD Rho project
## ver 1.0a	- 160504, Bugfix: as.factor(group) at line 5

ss.fc = function(data,group) {
	library(pbapply) # package for progressbar, install needed
	group.lv = levels(as.factor(group)) # bugfix 160504
	n = length(group.lv)
	cat('Level(s) of group =',group.lv,'\n\n')
	avs = NULL

	# Step 1. Calculate average of each group
	for(i in 1:n) {
		av = pbapply(data[,which(group==group.lv[i])],1,mean) # apply: 1-rows, 2-columns
		avs = cbind(avs,av)
	}
	colnames(avs) = group.lv # set colnames of avs df
	print(head(avs))

	# Step 2. Calculate fold-changes by every combination
	comb = combn(n,2) # compute combination
	m = length(comb[1,])
	fcs = NULL
	fc.col = NULL
	for(i in 1:m) {
		fc.col = c(fc.col,paste(group.lv[comb[1,i]],'/',
								group.lv[comb[2,i]],sep=""))
		cat('Process =',i,'/',m,'\n')
		fc = avs[,comb[1,i]]-avs[,comb[2,i]]
		fcs = cbind(fcs, fc)
	}
	colnames(fcs) = fc.col
	print(head(fcs))
	fcs2 = apply(fcs,2,function(x){ifelse(x<0,-2^(-x),2^x)}) # apply: 1-rows, 2-columns
	print(head(fcs2))

	# Step 3. Collect results
	#out = as.data.frame(cbind(avs,fcs,fcs2))
	return(as.data.frame(fcs2))
}

fc.cerev = as.data.frame(ss.fc(sgnl.cerev,group))
