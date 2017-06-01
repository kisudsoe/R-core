# 2016-09-22 THU ver.
## Ver 1.0	- 150903, original version of this function
## ver 1.1	- 160429, using for Yeast HD LD Rho project
## ver 1.2  - 160922, bugfix, using for X-ALD project
## ver 1.2a - 160930, minor edit
## ver 1.3  - 170531, add uniononly toggle and ttest for Yeast HD/LD/Rho project

ss.gene.selec = function(data.fc,fc=c(2,3), data.pval, uniononly=F,test=NULL) {
  n= length(fc)
  m= length(colnames(data.fc))

  result=NULL; uninum=NULL
  for(i in 1:n) {
    fc.result=NULL;fdr1.id=NULL; fdr5.id=NULL; fdr10.id=NULL
    for(j in 1:m) {
      if(test=="anova") {
        fdr1 = which(data.pval<=0.01 & abs(data.fc[,j])>fc[i])
        fdr5 = which(data.pval<=0.05 & abs(data.fc[,j])>fc[i])
        fdr10 = which(data.pval<=0.1 & abs(data.fc[,j])>fc[i])
      } else if (test=="ttest") {
        fdr1 = which(data.pval[,j]<=0.01 & abs(data.fc[,j])>fc[i])
        fdr5 = which(data.pval[,j]<=0.05 & abs(data.fc[,j])>fc[i])
        fdr10 = which(data.pval[,j]<=0.1 & abs(data.fc[,j])>fc[i])
      } else { stop("please input --> test='anova' or 'ttest'"); }
      f1 = length(fdr1)
      f5 = length(fdr5)
      f10 = length(fdr10)
      col = c(f1,f5,f10)

      fc.result = cbind(fc.result,col)
      fdr1.id = union(fdr1.id,fdr1)
      fdr5.id = union(fdr5.id,fdr5)
      fdr10.id = union(fdr10.id,fdr10)
    }
    uni = c(length(fdr1.id),length(fdr5.id),length(fdr10.id))

    if(uniononly==F) {
      result = rbind(result,fc.result)
      uninum = c(uninum,uni)
    } else if(uniononly==T) {
      result = cbind(result,uni)
    }
  }

  if(uniononly==F) {
    result = cbind(result,uninum)
    colnames(result) = c(colnames(data.fc),"Union")

    thrsh_fc = c(paste0("fc",fc))
    thrsh_pval = c(paste0(names(data.pval)[1],c('_1','_5','10')))
    thrsh = as.vector(outer(thrsh_pval,thrsh_fc,paste,sep="."))
    print("Significance thresholds are fixed at 1%, 5%, 10%.")
    rownames(result) = thrsh
  } else if(uniononly==T) {
    colnames(result) = c(paste0('FC',fc))
    rownames(result) = paste0(names(data.pval)[1],c('_1','_5','10'))
    print("Significance thresholds are fixed at 1%, 5%, 10%.")
  }

  write.csv(result,"ss.gene.selec_result.csv")
  return(result)
}

fc = data.frame(rma_df$LD.HD,rma_df$Rho0.HD,rma_df$Rho0.LD)
bonf = rma_df$stat.Bonferroni
names(bonf) = "Bonf"
fc_thrsh = c(1.5,2,3,5,10)

ss.gene.selec(data.fc=fc,fc=fc_thrsh, data.pval=bonf, uniononly=T,test='ttest')
# test='anova' or 'ttest'
