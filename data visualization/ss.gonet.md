---
title:    "ss.gonet.md"
author:   "Kim SS"
date:     "2017-06-19 MON"
note:     "This file was written for Yeast HD/LD/Rho project"
---

* Version info
* 1.0  - 170619, first release

```r {.lineNo}
ss.mdsMatrix = function(go) { # Matrix for Multidimensional scaling
  library(irr)
  terms = go$Term
  agr = NULL; kappa=NULL; prop=NULL
  time1=Sys.time()
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
  }
  colnames(agr)=terms; rownames(agr)=terms
  colnames(kappa)=terms; rownames(kappa)=terms
  colnames(prop)=terms; rownames(prop)=terms
  print(paste("Process time=",round(Sys.time()-time1,1),"sec"))
  return(out=list(agreement=agr,kappa=kappa,proportion=prop))
}

#install.packages("geomnet",repos="http://cran.us.r-project.org")
ss.gonet = function(mdsMat,goinfo,threshold=0.5,fname="ss.gonet",title="", width=10,height=8) {
  n = length(colnames(mdsMat))
  net = NULL
  for(i in 1:n) {
    c = mdsMat[,i]
    c.name = colnames(mdsMat)[i]
    m = length(c)
    net.sub = cbind(rep(c.name,m),names(c),c)
    net = rbind(net,net.sub)
  }
  colnames(net) = c("from_id","to_id","value")
  rownames(net) = NULL
  net = data.frame(net)
  net$value = as.numeric(as.character(net$value))
  net = subset(net, value >threshold)

  net.ann = merge(x=net,y=goinfo, by_x="from_id",by_y="from_id")
  title = paste0(title," (",names(mdsMat)[1]," threshold=",threshold,")")

  library(geomnet)
  g1=ggplot(data=net.ann,aes(from_id=from_id,to_id=to_id))+
    geom_net(aes(size=GeneN,shape=Category,colour=Category),
             layout.alg="kamadakawai",
             fontsize=4,directed=F,labelon=T,vjust=0,ecolour="grey70",
             labelcolour="black",linewidth=1.5,selfloops=F,repel=T)+
    ggtitle(title)+
    theme(plot.title=element_text(size=8,face="bold"))+
    theme_net()
  #print(g)
  g2=ggplot(data=net.ann,aes(from_id=from_id,to_id=to_id))+
    geom_net(aes(size=GeneN,shape=Category,colour=Category),
             layout.alg="kamadakawai",
             fontsize=4,directed=F,labelon=F,vjust=0,ecolour="grey70",
             labelcolour="black",linewidth=1.5,selfloops=F,repel=T)+
    ggtitle(title)+
    theme(plot.title=element_text(size=8,face="bold"))+
    theme_net()
  ggsave(paste0("gonet_",fname,".png"),g1,width=width,height=height)
  ggsave(paste0("gonet_",fname,"_n.png"),g2,width=width,height=height)
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
     ggtitle(title)
  g2=ggplot(data=mds.coor.ann,aes(x=x,y=y))+theme_bw()+
     geom_point(aes(size=GeneN,shape=Category,colour=Category))+
     ggtitle(title)
  ggsave(paste0("gomds_",fname,".png"),g1,width=width,height=height)
  ggsave(paste0("gomds_",fname,"_n.png"),g2,width=width,height=height)
}
```

1. Load data
------------

```r {.lineNo}
go = read.csv("170611_DAVID results.csv")
print(dim(go))
colnames(go) = c("Geneset","Category","Term","GeneN","Percentile","PValue","Genes", "DEGn","TotalN","TotalGene","Enrichment","Bonferroni","Benjamini","FDR")
str(go)

go_5 = subset(go,FDR<0.05)
print(dim(go_5))
data.frame(table(go_5$Geneset))
```

2. GeneSet Enrichment Analysis (GSEA)
-------------------------------------

1. Cohen's kappa score (irr package) `<-` comparing agreement proportion later
2. Draw GO netowrk (geomnet package)
3. Draw Multidimensional scaling (MDS) plot (cmdscale function)

```r {.lineNo}
genesets=c('0a_1789_total','0b_1204_total','0c_835_total',
          '1a_1789_HD+/LD','1b_1204_HD+/LD','1c_835_HD+/LD',
          '2a_1789_HD-/LD','2b_1204_HD-/LD',
          '3a_1789_LD+/Rho','3b_1204_LD+/Rho','3c_835_LD+/Rho',
          '4a_1789_LD-/Rho','4b_1204_LD-/Rho','4c_835_LD-/Rho',
          '5a_1789_HD+/Rho','5b_1204_HD+/Rho','5c_835_HD+/Rho',
          '6a_1789_HD-/Rho','6b_1204_HD-/Rho','6c_835_HD-/Rho')
fnames=c('0a_1789_total','0b_1204_total','0c_835_total',
         '1a_1789_HD.p.LD','1b_1204_HD.p.LD','1c_835_HD.p.LD',
         '2a_1789_HD.m.LD','2b_1204_HD.m.LD',
         '3a_1789_LD.p.Rho','3b_1204_LD.p.Rho','3c_835_LD.p.Rho',
         '4a_1789_LD.m.Rho','4b_1204_LD.m.Rho','4c_835_LD.m.Rho',
         '5a_1789_HD.p.Rho','5b_1204_HD.p.Rho','5c_835_HD.p.Rho',
         '6a_1789_HD.m.Rho','6b_1204_HD.m.Rho','6c_835_HD.m.Rho')

## 1. Filter geneset & Set filename ##
d = "170614"; indx = 20
gs=genesets[indx]; fname = paste0(d,"_",fnames[indx])
cat(paste(">> Geneset name=",gs,"\n"))
go_ = subset(go_5,Geneset==gs)
cat(paste("1/5. Filtered data. row=", dim(go_)[1], "col=",dim(go_)[2],"\n"))

## 2. Cohen's kappa score calculation using irr package ##
## ss.mdsMatrix: [[1]]=agreement; [[2]]=kappa; [[3]]=proportion
mdsmat = ss.mdsMatrix(go=go_)
cat(paste("2/5. MDS matrix process done. row=", dim(mdsmat[[2]])[1], "col=",dim(mdsmat[[2]])[2],"\n"))

## 3. Merge GO information for GO network ##
goinfo = unique(data.frame(from_id=go_$Term, Category=go_$Category,
                           GeneN=go_$GeneN, TotalN=go_$TotalN))
cat("3/5. GO information merge done.\n")

## 4. Draw GO network using Geomnet package ##
names(mdsmat[[1]])="agreement"; names(mdsmat[[2]])="kappa"; names(mdsmat[[3]])="proportion"
net = ss.gonet(mdsMat=mdsmat[[2]], goinfo=goinfo, threshold=0.6,
               title=gs, fname=fname, width=10, height=8)
cat(paste("4/5. Saved GO net result.\n>> File name=",fname,"\n"))

## 5. Draw Multidimensional scaling (MDS) plot ## <- 별로다. 구분이 잘 안됨
ss.gomds(mdsmat=mdsmat[[2]],goinfo=goinfo, title=gs,fname=fname,
         width=9,height=7)
cat(paste("5/5. Saved GO MDS plot.\n>> File name=",fname,"\n"))
```
