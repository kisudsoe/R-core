# 2016-05-06 ver
## ver 2.0	- 160503, Advanced process speed than 'ss.anova.fdr' function
## ver 2.1	- 160506, set rownames to result for Yeast HDLDRho project
## ver 2.1a - 160811, eliminate qvalue function (unavailable package)
## ver 2.2  - 160927, using ggplot2
## ver 2.3  - 170907, image save as png for Yeast HD/LD/Rho project

ss.anova.fdr2 = function(data,group) {
	library(pbapply) # package for progressbar, install needed
	n = length(rownames(data))
	out = NULL
	raw.p = NULL
	time1 = Sys.time()

	# Step 1. ANOVA iteration
	raw.p = pbapply(data,1,function(row) {
								if(!any(is.na(row))) {
									fit = lm(row~group,na.action=na.exclude)
									anovaP = anova(fit)$'Pr(>F)'[1]
								} else { anovaP = NA } })
	time2 = Sys.time()
	cat('Process iteration =',n,'\n')
	print(time2-time1)
	cat('----------------------------------\n\n')

	# Step 2. Adjusted p-values
	raw.p = as.vector(raw.p)
	length(raw.p)

	library(fdrtool)
	#pval.estimate.eta0(p1,method="bootstrap") # Draw plot
	q = fdrtool(raw.p,statistic="pvalue")
	fdrtool.qval = q$qval
	cat('(1/9) fdrtool.tfdr\t-> done\n')
	fdrtool.lfdr = q$lfdr
	cat('(2/9) fdrtool.dfdr\t-> done\n')

	#library(qvalue)
	#q2 = qvalue(raw.p)
	#qvalue.fdr = q2$qvalues
	cat('(3/9) qvalue.fdr\t-> skip, not available package (for R version 3.3.1)\n')

	stat.BH.fdr = p.adjust(raw.p,"BH")
	cat('(4/9) stat.fdr.BH\t-> done\n')
	#stat.Hommel = p.adjust(raw.p,"hommel") # take time too long
	cat('(5/9) stat.Hommel\t-> skip, Hommel takes time too long.\n')
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
                      #qvalue.fdr,
                      stat.BH.fdr,
                      #stat.Hommel,
                      stat.BY,
                      stat.Hochberg,
                      stat.Holm,
                      stat.Bonferroni)
	rownames(result) = rownames(data) # set row names

	## Plot draw start ##
	result.r = result[order(raw.p),] # re-ordered by raw.p

	matplot(result.r$raw.p, result.r, xlim=c(0,1),ylim=c(0,1),
            main="Various multiple test correction methods",
            xlab="Raw p-value",
            ylab="Adjusted p-value",
            type="l",
            col=c("black",rainbow(9)), # Set color
            lty=1,
            lwd=2)
	legend('bottomright',
           legend=c("raw.p",
                    "fdrtool.tfdr",
										"fdrtool.dlfdr",
                    #"qvalue.fdr",
                    "stat.BH.fdr",
                    #"stat.Hommel",
                    "stat.BY",
                    "stat.Hochberg",
                    "stat.Holm",
                    "stat.Bonfferoni"),
		   col=c("black",rainbow(9)), # Set color
           cex=1,
           pch=16)
  dev.copy(png,filename="ss.anova.fdr2_matplot.png",width=7,height=7,units="in",res=100)
  dev.off()

	pvals = stack(result)
	library(ggplot2)
  p = ggplot(pvals,aes(ind,values))+theme_bw()+
	  geom_point(aes(colour=ind),alpha=0.2,position="jitter")+
    geom_boxplot(alpha=0,colour="black")+
    labs(title="Distribution of adjusted p-values",
         x="Algorithms", y="p-values")+
    theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(filename="ss.anova.fdr2_pvaldist.png",plot=p,width=7,height=7,units="in")
  #print(p)

	#x11()
	par(mfrow=c(2,2)) # 2x2 plot partition
	hist(raw.p, xlim=c(0,1))
	hist(stat.BH.fdr, xlim=c(0,1))
	hist(stat.Holm, xlim=c(0,1))
	hist(stat.Bonferroni, xlim=c(0,1))

	cat('\nPlot draw\t\t-> done\n')
	cat('----------------------------------\n')

	## plot draw done ##

	#######################
	# Duration time check #
	#######################
	time2 = Sys.time()
	print(time2-time1)
	#######################

	return(result)
}

pval_gene = ss.anova.fdr2(gene_sgnl,group)
