ss.gene.selec = function(data.fc,fc=1.5,data.fdr,fcfdr) {
	n= length(colnames(data.fc))
  
	for(i in 1:n) {
		fdr1 = which(data.fdr<=0.01 & abs(data.fc[,i])>fc)
		f1 = length(fdr1)
		fdr5 = which(data.fdr<=0.05 & abs(data.fc[,i])>fc)
		f5 = length(fdr5)
		fdr10 = which(data.fdr<=0.1 & abs(data.fc[,i])>fc)
		f10 = length(fdr10)
		num = c(f1,f5,f10)
    
		if(i==1) {
			result = num
			fdr1.id = fdr1
			fdr5.id = fdr5
			fdr10.id = fdr10
		} else {
			result = cbind(result,num)
			fdr1.id = union(fdr1.id,fdr1)
			fdr5.id = union(fdr5.id,fdr5)
			fdr10.id = union(fdr10.id,fdr10)
		}
	}
  
	num = c(length(fdr1.id),length(fdr5.id),length(fdr10.id))
	result = cbind(result,num)
	colnames(result) = c(colnames(data.fc),"Union")
	#rownames(result) = c("fc2&fdr1","fc2&fdr5","fc2&fdr10")
	rownames(result) = fcfdr
	print(fcfdr)

	result.summ = list(result=result,
                       fdr1.id=fdr1.id,
                       fdr5.id=fdr5.id,
                       fdr10.id=fdr10.id)
  
	return(result.summ)
}

fcfdr2 = c("fc2&fdr1","fc2&fdr5","fc2&fdr10")
gene.num2 = ss.gene.selec(fc.cerev[,1:3],fc=2,
                          pval_cerev$stat.Bonferroni,fcfdr2) # 2015-09-03 fc.fdr.selec
View(gene.num2$result)