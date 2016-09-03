args = commandArgs(TRUE)
if (length(args) < 2) {
  cat("\n less argument \n")
  cat(" Usage : Rscript mk_generation.r input_file_name output_file_name \n")
  q()
}
df = read.csv(file=args[1],header=T)

corr.pval = function(df.input) {
  df1 = df.input[,7:13]
  df2 = df.input[,14:32]
  df1_nm = colnames(df.input)[7:13]
  df2_nm = colnames(df.input)[14:32]
  n1 = length(df1[1,])
  n2 = length(df2[1,])
  corr = NULL
  pval = NULL
  for(i in 1:n1) {
    corr_col = NULL
    pval_col = NULL
    for(j in 1:n2) {
      corr.t = cor.test(df1[,i],df2[,j])
      corr_col = c(corr_col,corr.t$estimate)
      pval_col = c(pval_col,corr.t$p.value)
    }
    corr = cbind(corr,corr_col)
    pval = cbind(pval,pval_col)
  }
  colnames(corr) = paste0('Corr_',df1_nm)
  rownames(corr) = df2_nm
  colnames(pval) = paste0('pval_',df1_nm)
  out = cbind(corr,pval)
  #print(out)
  write.csv(out,paste0(args[2],"_result.csv"))
  return()
}

corr.pval(df)
