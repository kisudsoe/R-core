# 2016-09-22 THU ver.
## Ver 1.0	- 150903, original version of this function
## ver 1.1	- 160429, using for Yeast HD LD Rho project
## ver 1.2  - 160922, bugfix, using for X-ALD project

ss.gene.selec = function(data.fc,fc=1.5,data.fdr,row.names) {
	n = length(colnames(data.fc)) # Column 1 should have AffyID
	result = NULL; fdr1.id = NULL; fdr5.id = NULL; fdr10.id = NULL;
	for(i in 2:n) {
		fdr1 = which(data.fdr<=0.01 & abs(data.fc[,i])>fc)
		f1 = length(fdr1)
		fdr5 = which(data.fdr<=0.05 & abs(data.fc[,i])>fc)
		f5 = length(fdr5)
		fdr10 = which(data.fdr<=0.1 & abs(data.fc[,i])>fc)
		f10 = length(fdr10)
		num1 = c(f1,f5,f10)

		result = cbind(result,num1)
		fdr1.id = union(fdr1.id,fdr1)
		fdr5.id = union(fdr5.id,fdr5)
		fdr10.id = union(fdr10.id,fdr10)
	}

	num2 = c(length(fdr1.id),length(fdr5.id),length(fdr10.id))
	result = cbind(result,num2)
	colnames(result) = c(colnames(data.fc)[2:n],"Union")
	rownames(result) = row.names
	print(row.names)

	result.summ = list(result=result,
                     fdr1.id=fdr1.id,
                     fdr5.id=fdr5.id,
                     fdr10.id=fdr10.id)

	return(result.summ)
}

fcfdr2 = c("fc2&fdr1","fc2&fdr5","fc2&fdr10")
gene.num2 = ss.gene.selec(fc.cerev[,1:3],fc=2,
                          pval_cerev$stat.Bonferroni,fcfdr2)
View(gene.num2$result)
