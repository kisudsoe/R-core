---
title:  "GO analysis"
author: "Kim SS"
date:   "2017-06-12 MON"
description: "This file is written for Yeast HD/LD/Rho project with DAVID results"
---

```r
#install.packages("irr",repos="http://cran.us.r-project.org")
ss.mdsMatrix = function(go) { # Matrix for Multidimensional scaling
  library(irr)
  terms = go$Term
  agr = NULL; kappa=NULL; prop=NULL
  time1=Sys.time(); n=length(terms); i=1
  pb=winProgressBar(title="ss.goGroup",label="ss.goGroup function progress", min=0,max=n,width=500)
  for(t in terms) {
    t.df = subset(go,Term==t)
    t.genes = strsplit(as.character(t.df$Genes),", ")[[1]]

    # union gene list
    union.genes.list = sapply(go$Genes,function(x) {
      strsplit(as.character(x),", ")[[1]]
    })
    union.genes = Reduce(union, union.genes.list)

    # Agreement proportion
    t.agr = apply(go[,c(3,7)],1,function(row) { # Select column= Term, Genes
      row.genes = strsplit(as.character(row[2]),", ")[[1]]
      int = intersect(t.genes,row.genes)
      uni.diff = setdiff(union.genes,union(t.genes,row.genes))
      agree = (length(int)+length(uni.diff))/length(union.genes)
      return(agree)
    })
    agr = cbind(agr,t.agr)

    # Kappa score
    t.ka = apply(go[,c(3,7)],1,function(row) {
      row.genes = strsplit(as.character(row[2]),", ")[[1]]
      g1 = t.genes[match(union.genes,t.genes)]
      g2 = row.genes[match(union.genes,row.genes)]
      g3 = cbind(g1,g2)
      g = g3!="" # Value exist -> TRUE
      g[is.na(g)]=F # NA -> FALSE
      return(kappa2(g)$value)
      #return(agree(g)$value)
    })
    kappa = cbind(kappa,t.ka)

    # Proportion
    t.pr = apply(go[,c(3,7)],1,function(row) {
      row.genes = strsplit(as.character(row[2]),", ")[[1]]
      int = intersect(t.genes,row.genes)
      pr = length(int)/length(t.genes)
      return(pr)
    })
    prop = cbind(prop,t.pr)

    ## Progress time ##
    d=difftime(Sys.time(),time1,unit="sec")
    if(d<60) d=paste(round(d,1),"sec")
    else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
    else if(d>=3600) d=paste(round(d/3600,1),"hr")
    setWinProgressBar(pb,i,title=paste0(round(i/n*100,1)," % (",i,"/",n,") done for ",d))
    ###################
    i=i+1
  }
  colnames(agr)=terms; rownames(agr)=terms
  colnames(kappa)=terms; rownames(kappa)=terms
  colnames(prop)=terms; rownames(prop)=terms
  print(difftime(Sys.time(),time1,units="auto")); close(pb)
  return(out=list(agreement=agr,kappa=kappa,proportion=prop))
}

#install.packages("geomnet",repos="http://cran.us.r-project.org")
ss.gonet = function(mdsMat,goinfo,group=NULL,threshold=0.5,fname="ss.gonet",title="", width=10,height=8) {
  n = length(colnames(mdsMat))
  net = NULL
  #for(i in 1:n) {
  #  c = mdsMat[,i]; c.name = colnames(mdsMat)[i]; m = length(c)
  #  net.sub = cbind(rep(c.name,m),names(c),c)
  #  net = rbind(net,net.sub)
  #}
  net = data.frame(as.table(mdsMat))
  colnames(net) = c("from_id","to_id","value")
  rownames(net) = NULL
  net = data.frame(net)
  net$value = as.numeric(as.character(net$value))
  net = subset(net, value >threshold)

  net.ann = merge(x=net,y=goinfo, by_x="from_id",by_y="from_id")
  title = paste0(title," (",names(mdsMat)[1]," <",threshold,")")
  #write.csv(net.ann,paste0(fname,".csv"),row.names=F)

  library(geomnet)
  if(group=="Category") Category_ = factor(net.ann$Category)
  else if(group=="Geneset") Category_ = factor(net.ann$Geneset)
  g1=ggplot(data=net.ann,aes(from_id=from_id,to_id=to_id))+
    geom_net(aes(size=GeneN,shape=Category_,colour=Category_),
             layout.alg="kamadakawai",
             fontsize=2.5,directed=F,labelon=T,vjust=0,ecolour="grey70",
             labelcolour="black",linewidth=1.5,selfloops=F,repel=T)+
    ggtitle(title)+
    scale_shape_manual(values=1:nlevels(Category_))+
    theme(plot.title=element_text(size=8,face="bold"))+
    theme_net()
  #print(g)
  g2=ggplot(data=net.ann,aes(from_id=from_id,to_id=to_id))+
    geom_net(aes(size=GeneN,shape=Category_,colour=Category_),
             layout.alg="kamadakawai",
             fontsize=2.5,directed=F,labelon=F,vjust=0,ecolour="grey70",
             labelcolour="black",linewidth=1.5,selfloops=F,repel=T)+
    ggtitle(title)+
    scale_shape_manual(values=1:nlevels(Category_))+
    theme(plot.title=element_text(size=8,face="bold"))+
    theme_net()

  ggsave(paste0("gonet_",fname,"_",threshold,".png"),g1,width=width,height=height)
  ggsave(paste0("gonet_",fname,"_",threshold,"_n.png"),g2,width=width,height=height)
  return(net.ann)
}

#install.packages("ggrepel",repos="http://cran.us.r-project.org")
ss.gomds = function(mdsmat, goinfo, threshold=0.5, title="", fname="ss.gomds",width=8,height=7) {
  mds = cmdscale(dist(mdsmat),eig=T,k=2)
  mds.coor = data.frame(mds$points)
  mds.coor$from_id = rownames(mds.coor)
  colnames(mds.coor) = c("x","y","from_id")
  mds.coor.ann = merge(x=mds.coor,y=goinfo, by_x="from_id",by_y="from_id")

  library(ggplot2); library(ggrepel)
  title = paste0(gs," (by ",names(mdsmat)[1],")")
  g1=ggplot(data=mds.coor.ann,aes(x=x,y=y))+theme_bw()+
     geom_point(aes(size=GeneN,shape=Category,colour=Category))+
     geom_text_repel(aes(label=from_id),size=4)+
     scale_shape_manual(values=1:nlevels(Category_))+
     ggtitle(title)
  g2=ggplot(data=mds.coor.ann,aes(x=x,y=y))+theme_bw()+
     geom_point(aes(size=GeneN,shape=Category,colour=Category))+
     scale_shape_manual(values=1:nlevels(Category_))+
     ggtitle(title)
  ggsave(paste0("gomds_",fname,".png"),g1,width=width,height=height)
  ggsave(paste0("gomds_",fname,"_n.png"),g2,width=width,height=height)
}

## To generage GO groups by kappa scores
ss.goGroup = function(net) {
  n=nrow(net); group=NULL; time1=Sys.time()
  pb=winProgressBar(title="ss.goGroup",label="ss.goGroup function progress", min=0,max=n,width=500)
  ann = unique(net[,c(-2,-3)]); net_ = unique(net[,1:2])
  colnames(ann)=c("Goid","Geneset","Category","GeneN","TermN")
  for(i in 1:n) {
    go = c(as.character(net_[i,1]),as.character(net_[i,2]))
    if(i==1) group[[1]] = go
    else {
      m=length(group)
      for(j in 1:m) {
        int = length(intersect(group[[j]],go))
        if(int>0) { group[[j]] = union(group[[j]],go); break }
        else if(int==0 && j==m) group[[j+1]] = go
    } }
    ## Progress time ##
    d=difftime(Sys.time(),time1,unit="sec")
    if(d<60) d=paste(round(d,1),"sec")
    else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
    else if(d>=3600) d=paste(round(d/3600,1),"hr")
    setWinProgressBar(pb,i,title=paste0(round(i/n*100,1)," % (",i,"/",n,") done for ",d))
    ###################
  }
  group2=NULL
  for(i in 1:length(group)) {
    if(i==1) group2[[1]] = group[[i]]
    for(j in 1:length(group2)) {
      int = length(intersect(group[[i]],group2[[j]]))
      if(int>0) { group2[[j]] = union(group[[i]],group2[[j]]); break }
      else if(int==0 && j==length(group2)) group2[[j+1]] = group[[i]]
    }
  }
  names(group2) = c(1:length(group2))
  out = stack(group2); colnames(out) = c("Goid","GO_groups")
  out_ann = unique(merge(x=out,y=ann, by_x="Goid",by_y="Goid"))
  print(difftime(Sys.time(),time1,units="auto")); close(pb)
  return(out_ann)
}
```

## 1. Load data

* If single geneset analysis: goto **section 2**
* If multiple geneset analysis: goto **section 3**

```r
go = read.csv("170718_DAVID_1738degs_tukey2.csv")
print(dim(go))
colnames(go) = c("Geneset","Category","Term","GeneN","Percentile","PValue","Genes", "DEGn","TotalN","TotalGene","Enrichment","Bonferroni","Benjamini","FDR")
str(go)

go_5 = subset(go,FDR<0.05)
print(dim(go_5))
data.frame(table(go_5$Geneset))
data.frame(table(go_5$Category))
```

## 2. GeneSet Enrichment Analysis (GSEA)

1. Cohen's kappa score (irr package) `<-` comparing agreement proportion later
2. Draw GO netowrk (geomnet package)
3. Draw Multidimensional scaling (MDS) plot (cmdscale function)

```r
## 170718
genesets=c("aab","aba","abb","abc","acb",
           "baa","bab","bac","bba","bca","cab","cba")
fnames = genesets
print(length(genesets))

## 1. Filter geneset & Set filename ##
d = "170718"; indx = 1
#gs=genesets[indx]; fname = paste0(d,"_",fnames[indx])
gs = genesets; fname = paste0(d,"_tukey2")
cat(paste(">> Geneset name=",paste(gs,collapse=","),"\n"))
go_ = subset(go_5,Geneset %in% gs)
cat(paste("1/6. Filtered data. row=", dim(go_)[1], "col=",dim(go_)[2],"\n"))

## 2. Cohen's kappa score calculation using irr package ##
## ss.mdsMatrix: [[1]]=agreement; [[2]]=kappa; [[3]]=proportion
mdsmat = ss.mdsMatrix(go=go_)
cat(paste("2/6. MDS matrix process done. row=", dim(mdsmat[[2]])[1], "col=",dim(mdsmat[[2]])[2],"\n"))

## 3. Merge GO information for GO network ##
goinfo = unique(data.frame(from_id=go_$Term, Geneset=go_$Geneset, Category=go_$Category,
                           GeneN=go_$GeneN, TotalN=go_$TotalN))
cat("3/6. GO information merge done.\n")

## 4. Draw GO network using Geomnet package ##
names(mdsmat[[1]])="agreement"; names(mdsmat[[2]])="kappa"; names(mdsmat[[3]])="proportion"
net = ss.gonet(mdsMat=mdsmat[[2]], goinfo=goinfo, group="Geneset", threshold=0.6,
               title=gs, fname=fname, width=7, height=6)
cat(paste("4/6. Saved GO net result.\n>> File name=",fname,"\n"))

## 5. Draw Multidimensional scaling (MDS) plot ## <- 별로다. 구분이 잘 안됨
ss.gomds(mdsmat=mdsmat[[2]],goinfo=goinfo, title=gs, fname=fname,
         width=9,height=7)
cat(paste("5/6. Saved GO MDS plot.\n>> File name=",fname,"\n"))

## 6. Export GO table with kappa group
goGroup = ss.goGroup(net)
write.csv(goGroup,paste0(fname,"_GOgroups.csv"),row.names=F)
cat(paste("6/6. Saved GO groups by kappa threshold.\n>> File name=",fname,"_GOgroups.csv\n"))
```

## 3. GSEA for Multi-genesets

1. Cohen's kappa score (irr package) `<-` comparing agreement proportion later
2. Draw GO netowrk (geomnet package)

```r
## 170629
genesetlist=list(c('mmm','mmm_aab','mmm_baa','mmm_cab'),
              c('mpm','mpm_baa','mpm_bab','mpm_cab'),
              c('mpp','mpp_aab','mpp_bab','mpp_cab'),
              c('pmm','pmm_aab','pmm_cab'),
              c('pmp','pmp_baa','pmp_bab','pmp_cab'),
              c('ppp','ppp_aab','ppp_baa','ppp_cab'))
fnamelist = c('mmm','mpm','mpp','pmm','pmp','ppp')

## 1. Filter geneset & Set filename ##
d = "170629"; indx_ = 6
gs=genesetlist[[indx_]]; fname = paste0(d,"_",fnamelist[indx_])
cat(paste("Geneset names=",paste(gs,collapse=", "),"\n"))
go_ = subset(go_5,Geneset%in%gs)
cat(paste("1/4. Filtered data. row=", dim(go_)[1], "col=",dim(go_)[2],"\n"))

## 2. Cohen's kappa score calculation using irr package ##
## ss.mdsMatrix: [[1]]=agreement; [[2]]=kappa; [[3]]=proportion
mdsmat = ss.mdsMatrix(go=go_)
cat(paste("2/4. MDS matrix process done. row=", dim(mdsmat[[2]])[1], "col=",dim(mdsmat[[2]])[2],"\n"))

## 3. Merge GO information for GO network ##
goinfo = unique(data.frame(from_id=go_$Term, Geneset=go_$Geneset, Category=go_$Category,
                           GeneN=go_$GeneN, TotalN=go_$TotalN))
cat("3/4. GO information merge done.\n")

## 4. Draw GO network using Geomnet package ##
names(mdsmat[[1]])="agreement"; names(mdsmat[[2]])="kappa"; names(mdsmat[[3]])="proportion"
net = ss.gonet(mdsMat=mdsmat[[2]], goinfo=goinfo, group="Geneset", threshold=0.6,
               title=gs, fname=fname, width=9, height=8)
cat(paste("4/4. Saved GO net result.\n  >> File name=",fname,".png\n"))

## 5. Export GO table with kappa group
goGroup = ss.goGroup(net)
write.csv(goGroup,paste0(fname,"_GOgroups.csv"),row.names=F)
cat(paste("6/6. Saved GO groups by kappa threshold.\n>> File name=",fname,"_GOgroups.csv\n"))
```
