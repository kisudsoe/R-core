# 2016-06-29 WED, 06-30 THU----
## R update to v3.3.1
install.packages("installr")
library(installr)
updateR()

## Multi factor analysis for dog characters
## Using FactoMineR library
data = read.table("clipboard",header=T,row.names=1)

install.packages("FactoMineR")
install.packages("missMDA")
library(FactoMineR)
library(missMDA) # handle missing value for FactoMineR

res1 = PCA(data_res) # 11 categories

res_impute.comp = imputePCA(data_res)
res_impute = PCA(res_impute.comp$completeObs)

res_prcomp = prcomp(data_res,scale=TRUE) # Should be no 'NA' value in data
biplot(res_prcomp)
fa1 = factanal(data_res, factor=3, rotation="varimax") # Should be no 'NA' value in data
fa1


# 2016-07-01 FRI----
## PCAplot function for displaying PCA result
## chull function for draw borderline of dots in scatter plot
## ver 1.0	- 160701, Dog project
PCAplot = function(PCAout,criteria,range=NULL,
                   crite2=NULL,range2=NULL,hulls=TRUE,labels=TRUE) {
  #par(mar=c(5,5,5,5), xpd=TRUE) # (Bottom, Left, Top, Right)
  par(xpd=FALSE)
  plot(PCAout$ind$coord[,1:2], pch=19, col="black",
       xlab=paste("Dim 1(",round(PCAout$eig[1,2],2),"%)"),
       ylab=paste("Dim 2(",round(PCAout$eig[2,2],2),"%)"),
       main="Individuals factor map (PCA)")
  abline(h=0,lty=2)
  abline(v=0,lty=2)
  
  criteria.title = names(criteria)[1]
  if(is.null(range)){
    criteria = as.factor(criteria)
    criteria.lv = levels(criteria)
  } else {
    criteria.lv = range # For criteria range
  }
  cols = rainbow(length(criteria.lv))
  colstrans = rainbow(length(criteria.lv), alpha=0.3)
  
  ## This section for crite2 border hull polygons
  if(!is.null(crite2)){
    crite2.title = names(crite2)[1]
    if(is.null(range2))
      crite2.lv = levels(as.factor(crite2))
    else
      crite2.lv = range2
    cols2 = gray.colors(length(crite2.lv))
    colstrans2 = gray.colors(length(crite2.lv),alpha=0.3)
    
    for(i in 1:length(crite2.lv)){
      # getting the convex hull of each unique point set
      if(is.null(range2)){
        crite2Ids = which(crite2==crite2.lv[i])
      } else {
        if(i==1){
          crite2Ids = which(0<crite2 & crite2<=crite2.lv[i])
        }else{
          crite2Ids = which(crite2.lv[i-1]<crite2 & crite2<=crite2.lv[i])
        }
      }
      Coord = PCAout$ind$coord[crite2Ids,1:2]
      hullcoor = Coord[chull(Coord[,1],Coord[,2]),]
      #polygon(hullcoor,col=colstrans2[i], border=cols2[i])
      polygon(hullcoor,col=colstrans2[i], border="black")
    }
    legend("topright",title=crite2.title,legend=crite2.lv,pch=15,col=cols2)
  }
  
  for(i in 1:length(criteria.lv)){
    if(is.null(range)) {
      criteriaIds = which(criteria==criteria.lv[i])
    } else {
      if(i==1){
        criteriaIds = which(0<criteria & criteria<=criteria.lv[i])
      }else{
        criteriaIds = which(criteria.lv[i-1]<criteria & criteria<=criteria.lv[i])
      }
    }
    Coord = PCAout$ind$coord[criteriaIds,1:2]
    points(Coord, pch=19, col=cols[i])
    
    # getting the convex hull of each unique point set
    if(hulls==TRUE) {
      hullcoor = Coord[chull(Coord[,1],Coord[,2]),]
      polygon(hullcoor,col=colstrans[i], border=colstrans[i]) # border=cols[i]
    }
  }
  
  par(xpd=TRUE)
  legend("topleft",title=criteria.title,legend=criteria.lv,pch=19,col=cols)
  if(labels==TRUE) {
    text(PCAout$ind$coord,labels=rownames(res1$call$X),pos=3)
  }
}

category = setNames(res1$call$X$Age,"Age")
category2 = setNames(res1$call$X$Species,"Species")
PCAplot(PCAout=res1,criteria=category,
        #range=c(0.65,0.7,0.75,0.85), # AT ratio
        #range=c(20,25,30,35,45), # Weight
        range=c(4,8,14), # Age
        crite2=category2,
        #range2=c(0.65,0.7,0.75,0.85), # AT ratio
        hulls=FALSE,labels=FALSE)

summary(res1$call$X)
summary(as.factor(category))
hist(category)
table(cut(category,breaks=c(0,4,8,14)))