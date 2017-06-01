### ROC curve ###
# Script by Kim, Seung-Soo
# ver 1.0 written 2014-10-27 TUE

### Bibliography ###
# http://publicifsv.sund.ku.dk/~tag/Teaching/Barcelona/lecturenotes/How-to-Roc-Barcelona.html
# http://rpubs.com/cardiomoon/64987

id.neg_control=which(data.ann$mRna...Description=="neg_control ")
id.pos_control=which(data.ann$mRna...Description=="pos_control ")
tmp_neg=data.frame(data[id.neg_control,],group="neg_control")
tmp_pos=data.frame(data[id.pos_control,],group="pos_control")
data_control=rbind(tmp_neg,tmp_pos)
 
# Histogram generation
library(ggplot2)
ggplot(data_control, aes(x=W14M1, fill=group, colour=group))+
	geom_histogram(aes(y=..density..), binwidth=1, alpha=.7,
				   position="identity")+
	ggtitle("Histogram of controls")


# ROC curve generation
library(Epi)
data_control.ROC=ROC(form=group~W14M1,data=data_control,plot="ROC")
 
optimal_lr.eta=function(x){
	no=which.max(x$res$sens+x$res$spec)[1]
	result=x$res$lr.eta[no]
	result
}
optimal_lr.eta(data_control.ROC)
 
optimal_cutpoint=function(x){
	y=optimal_lr.eta(x)
	b0=unname(x$lr$coeff[1])
	b1=unname(x$lr$coeff[2])
	result=(-log(1/y-1)-b0)/b1
	result
}
optimal_cutpoint(data_control.ROC)
 
# optimal_cutpoint vector
ssROC = function(target) {
	library(Epi)
	leng=length(colnames(target))-1
	coln = colnames(target)
	cat(paste(" data column length =",leng,"\n",
              "data length =",length(target[,1]),"\n"))
 
	for(i in 1:leng) {
		observ = target[,i]
		result.ROC = ROC(test=observ,stat=target$group,plot="ROC")
 
		#Get optimal_lr.eta?
		no=which.max(result.ROC$res$sens+result.ROC$res$spec)[1]
		y = result.ROC$res$observ[no]
 
		# Get optimal_cutpoint
		#b0=unname(result.ROC$lr$coeff[1])
		#b1=unname(result.ROC$lr$coeff[2])
		#out.r = (-log(1/y-1)-b0)/b1
		out.ROC = y
 
		if(i==1)
			out = out.ROC
		else
			out = rbind(out,out.ROC)
	}
	#colnames(out) = c("lr.eta","cutpoint")
	#rownames(out) = coln
  
	return(out)
}
tmp=ssROC(data_control)