# 2016-05-01 SUN
## ver 1.0	- 160501, Yeast HD LD Rho project
## ver 1.0a	- 160504, Bugfix: as.factor(group)
## ver 1.1  - 161114, add option: could set base group

ss.fc = function(data,group,base=NULL) {
	library(pbapply) # package for progressbar
	group.lv = levels(as.factor(group)) # bugfix 160504
	n = length(group.lv)
	cat('Level(s) of group =',group.lv,'\n\n')
	avs = NULL

	# Step 1. Calculate average of each group
	for(i in 1:n) {
		av = pbapply(data[,which(group==group.lv[i])],1,mean)
		avs = cbind(avs,av)
	}
	colnames(avs) = group.lv # set colnames of avs df
	print(head(avs))

	# Step 2. Calculate fold-changes by every combination
	if(length(base)==0) comb = combn(n,2) # compute combination
	else if(length(base)>0) { # add code 161114
	  n1 = which(group.lv==base)
	  comb = t(matrix(c(setdiff(c(1:n),n1),rep(n1,n-1)),ncol=2))
	  print(comb)
	}
	m = length(comb[1,])
	cat('Process iteration = ',m,'\n')
	fcs = NULL
	fc.col = NULL
	for(i in 1:m) {
		fc.col = c(fc.col,paste(group.lv[comb[1,i]],'/',
								group.lv[comb[2,i]],sep=""))
		#cat('Process =',i,'/',m,'\n')
		fc = avs[,comb[1,i]]-avs[,comb[2,i]]
		fcs = cbind(fcs, fc)
	}
	colnames(fcs) = fc.col
	print(tbl_df(fcs))
	fcs2 = pbapply(fcs,2,function(x){ifelse(x<0,-2^(-x),2^x)})
	print(tbl_df(fcs2))

	# Step 3. Collect results
	return(as.data.frame(fcs2))
}
cerev.fc = ss.fc(data=cerev.sgl[,-1],
                 group=droplevels(group.info$groups[-1]),
                 base="NR_005h")
