# 2016-09-22 THU ver.
## Ver 1.0	- 150903, original version of this function
## ver 1.1	- 160429, using for Yeast HD LD Rho project
## ver 1.2  - 160922, bugfix, using for X-ALD project
## ver 1.2a - 160930, minor edit
## ver 1.3  - 170531, add uniononly toggle and ttest for Yeast HD/LD/Rho project
## ver 1.4  - 170907, add fc criteria option

ss.gene.selec = function(data.fc, fc=c(2,3), data.pval, pval=c(0.01,0.05,0.1), uniononly=F,test=NULL) {
  n= length(fc); m= ncol(data.fc); l= length(pval)
  result=NULL; uninum=NULL
  for(i in 1:n) { # iteration for fc
    fc.result=NULL;fdr.id=NULL#;fdr1.id=NULL; fdr5.id=NULL; fdr10.id=NULL
    for(j in 1:m) {
      fdr=NULL; col = NULL
      for(k in 1:l) {
        if(test=="anova") {
          fdr[[k]] = which(data.pval<=pval[k] & abs(data.fc[,j])>fc[i])
        } else if(test=="ttest") {
          fdr[[k]] = which(data.pval[,j]<=pval[k] & abs(data.fc[,j])>fc[i])
        } else { stop("please input --> test='anova' or 'ttest'"); }
        col = c(col,length(fdr[[k]]))
        if(j==1) { fdr.id[k]=list(NULL) }
        fdr.id[[k]] = union(fdr.id[[k]],fdr[[k]])
      }
      fc.result = cbind(fc.result,col)
    }

    uni = NULL
    for(k in 1:l) { uni = c(uni,length(fdr.id[[k]])) }

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

    thrsh_fc = c(paste0("FC",fc))
    thrsh_pval = c(paste0(names(data.pval)[1],pval))
    thrsh = as.vector(outer(thrsh_pval,thrsh_fc,paste,sep="."))
    rownames(result) = thrsh
  } else if(uniononly==T) {
    colnames(result) = c(paste0('FC',fc))
    rownames(result) = paste0(names(data.pval)[1],pval)
  }

  write.csv(result,"ss.gene.selec_result.csv")
  print("Result saved as 'ss.gene.selec_result.csv'.")
  return(result)
}

metab_fc = data.frame(metab_$HD.LD,metab_$LD.Rho0)
fc_thrsh = c(1.5,2,3,5)
metab_ttest = data.frame(metab_$t.test1,metab_$t.test2); names(metab_ttest) = "ttest"
metab_fdr = metab_$FDR; names(metab_fdr) = "FDR"

ss.gene.selec(data.fc=metab_fc,fc=fc_thrsh, data.pval=metab_ttest,uniononly=F,test='ttest')
# test='anova' or 'ttest'
