## Input group information & t-test
group = c("CC-ALD","CC-ALD","CC-ALD","AMN","AMN","AMN","AMN","CC-ALD")
ss.ttest.fdr = function(data, group) {
  n = length(rownames(data))
  group = as.factor(group)
  print(summary(group))
  grp1_id = which(group==levels(group)[1])
  grp2_id = which(group==levels(group)[2])
  
  raw.p = NULL
  rown = rownames(data)
  print(head(rown))
  
  # STEP 1. Calculate t-test p-value
  for(i in 1:n) {
    pval = t.test(data[i,grp1_id],data[i,grp2_id])$p.value
    raw.p = c(raw.p,pval)
    
    #######################
    # Create progress bar #
    #######################
    if (i==1) { # set progress bar
      cat('\nProcess iteration =',n,'\n')
      pb = txtProgressBar(min=0,max=100,initial=0,style=3)
      time1 = Sys.time()
    } else { setTxtProgressBar(pb, (i/n)*100) } # show progress bar
    if(i==n) { # Duration time check
      close(pb)
      time2 = Sys.time()
      cat('\n')
      print(time2-time1)
    } 
    #######################
  }
  cat('----------------------------------')
  
  # STEP 2. Calculate adjust p-values
  raw.p = as.vector(raw.p)
  library(fdrtool)
  #pval.estimate.eta0(p1,method="bootstrap") # plot 그리기
  q = fdrtool(raw.p,statistic="pvalue")
  fdrtool.qval = q$qval
  cat('(1/9) fdrtool.tfdr\t-> done\n')
  fdrtool.lfdr = q$lfdr
  cat('(2/9) fdrtool.dfdr\t-> done\n')
  
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("qvalue")
  library(qvalue)
  q2 = qvalue(raw.p)
  qvalue.fdr = q2$qvalues
  cat('(3/9) qvalue.fdr\t-> done\n')
  
  stat.BH.fdr = p.adjust(raw.p,"BH")
  cat('(4/9) stat.fdr.BH\t-> done\n')
  #stat.Hommel = p.adjust(raw.p,"hommel") # take time too long
  cat('(5/9) stat.Hommel\t-> skip, this takes time too long.\n')
  stat.BY = p.adjust(raw.p,"BY")
  cat('(6/9) stat.BY\t\t-> done\n')
  stat.Hochberg = p.adjust(raw.p,"hochberg")
  cat('(7/9) stat.Hochberg\t-> done\n')
  stat.Holm = p.adjust(raw.p,"holm")
  cat('(8/9) stat.Holm\t\t-> done\n')
  stat.Bonferroni = p.adjust(raw.p,"bonferroni")
  cat('(9/9) stat.Bonferroni\t-> done\n')
  
  result = data.frame(raw.p,
                      fdrtool.qval,
                      fdrtool.lfdr,
                      qvalue.fdr,
                      stat.BH.fdr,
                      #stat.Hommel,
                      stat.BY,
                      stat.Hochberg,
                      stat.Holm,
                      stat.Bonferroni)
  rownames(result) = rown # set row names
  
  ## Plot draw start ##
  result.r = result[order(raw.p),] # re-ordered by raw.p
  
  matplot(result.r$raw.p, result.r, xlim=c(0,1),ylim=c(0,1),
          main="Various multiple test correction methods",
          xlab="Raw p-value",
          ylab="Adjusted p-value",
          type="l",
          col=c("black",rainbow(9)), # 색지정
          lty=1,
          lwd=2)
  legend('bottomright',
         legend=c("raw.p",
                  "fdrtool.tfdr",
                  "fdrtool.dlfdr",
                  "qvalue.fdr",
                  "stat.BH.fdr",
                  #"stat.Hommel",
                  "stat.BY",
                  "stat.Hochberg",
                  "stat.Holm",
                  "stat.Bonfferoni"),
         col=c("black",rainbow(9)), # 색지정
         cex=1,
         pch=16)
  
  x11()
  par(mfrow=c(2,2)) # 2행 2열로 plot 파티션 나누기
  hist(raw.p, xlim=c(0,1))
  hist(fdrtool.qval, xlim=c(0,1))
  hist(stat.BH.fdr, xlim=c(0,1))
  hist(stat.Bonferroni, xlim=c(0,1))
  
  cat('\nPlot draw\t\t-> done\n')
  cat('----------------------------------')
  ## plot draw done ##
  
  #######################
  # Duration time check #
  #######################
  time2 = Sys.time()
  print(time2-time1)
  #######################
  
  return(result)
}