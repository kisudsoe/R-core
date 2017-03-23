# 2016-11-04 FRI for Yeast dynamic network
ss.heatmap3 = function(mat,group.info,clstr=NULL,hc.order=NULL,yeast.cerev_somgrp=NULL) {
  mat = as.matrix(yeast.cerev.sgl[-1])
  if(length(clstr)>0) {
    mat_id = which(yeast.cerev_somgrp[,3] %in% clstr)
    mat.s = t(scale(t(mat)))[mat_id,]
  } else mat.s = t(scale(t(mat)))

  col.rg = c(-3,0,3)
  colGroup = data.frame(#Group=group.info$groups[-1],
                      #Group2=group.info$groups2[-1],
                      Time=group.info$time[-1],
                      #Glucose=group.info$glucose[-1],
                      Treat=group.info$treat[-1])
  colGroup2 = data.frame(Group2=group.info$groups2[-1])
  ha = HeatmapAnnotation(df=colGroup)
  ha2 = HeatmapAnnotation(df=colGroup2)

  hc = hclust(as.dist(1-cor(mat, method="spearman")),method="complete")
  hc.dd = as.dendrogram(hc)
  #order.dendrogram(hc.dd)
  hc.dd = rotate(hc.dd,hc.order)
  hr = hclust(as.dist(1-cor(t(mat),method="spearman")),method="complete")
  hr.dd = as.dendrogram(hr)

  ## filter somgrp
  if(length(clstr)>0) yeast.cerev_somgrp = yeast.cerev_somgrp[mat_id,]

  ## sort data.s by somgrp
  mat.s = as.data.frame(mat.s)
  mat.s_som = data.frame(cbind(mat.s,som=yeast.cerev_somgrp[,3]))
  mat.s_som2 = mat.s_som[order(mat.s_som$som),]
  mat.s = apply(mat.s_som2[,1:60],2,function(x) as.numeric(as.character(x)))

  ## set rowAnnotation
  l = NULL
  rainb = list()
  yeast.cerev_somgrp = yeast.cerev_somgrp[order(yeast.cerev_somgrp[,3]),]
  for(i in 1:length(colnames(yeast.cerev_somgrp))) {
    l1 = length(levels(factor(yeast.cerev_somgrp[,i])))
    l = c(l,l1)
    rainb[[i]] = rainbow(l1)
    names(rainb[[i]]) = levels(factor(yeast.cerev_somgrp[,i]))
  }
  nms = colnames(yeast.cerev_somgrp)
  rowA = rowAnnotation(df=yeast.cerev_somgrp,
                       col=list(som9=rainb[[1]],
                                som25=rainb[[2]],
                                som100=rainb[[3]]))

  cell.cols = colorRamp2(col.rg,c("Cyan","black","Yellow"))
  Heatmap(mat.s, col=cell.cols,
			    name = "Probeset\nintensity\n(z-score)",
			    cluster_rows=F,
			    cluster_columns=hc.dd,
			    top_annotation=ha,
			    bottom_annotation =ha2,
			    column_dend_height=unit(3,"cm"),
			    column_dend_reorder=FALSE,
			    row_dend_width=unit(3,"cm"),
			    row_names_side="left",
			    row_names_gp=gpar(fontsize=7))+rowA
}
ss.heatmap2(mat=yeast.cerev.sgl[1],group.info=group.info,
           clstr=c("A01","A02","A03","A04","A05",
                   "A06","A07","A08","A09","A10"),
           hc.order=c(1:12,31:51,53:60,52,28:30,13:27),
           yeast.cerev_somgrp=yeast.cerev_somgrp)


# v1.1  160720 - Color edit and label added
install.packages("heatmap.plus")
ss.heatmap = function(data,hr=TRUE,row.col=NULL,hc=TRUE,hc.order=NULL) {
	library("gplots"); library("devtools")
	library("heatmap.plus")
	x = as.matrix(data)
	cols = rainbow(5)
	col.cols = c(rep(cols[1],3),rep(cols[2],3),rep(cols[3],3),rep(cols[4],3),rep(cols[5],3))
	print(col.cols)
	if(length(row.col)>0) {
		#row.fac = c("grey80","yellow")
		row.fac = rainbow(length(levels(factor(row.col))))
		row.cols = row.fac[as.vector(row.col)]
	} else {
		row.cols = NULL
	}

	# Hierarchical clustering #
	## Correlation methods = "Pearson", "Kendall", "Spearman"
	## hclust methods = "complete", "average", "median", "centroid", etc.
	data.mat = as.matrix(data)
	if(hr==TRUE) {
		hr = hclust(as.dist(1-cor(t(data.mat),method="pearson")),method="complete")
		hr.dd = as.dendrogram(hr)
	} else {
		hr.dd = FALSE
	}
	if(hc==TRUE) {
		hc = hclust(as.dist(1-cor(data.mat, method="spearman")),method="complete")
		hc.dd = as.dendrogram(hc)
		print(order.dendrogram(hc.dd))
		if(length(hc.order)>0) {
			hc.dd = reorder(hc.dd,hc.order,mean)
			print(order.dendrogram(hc.dd))
		}
	}

	#source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
	# Color function to generate blue-yellow heat maps
	#png(filename="heatmap.png", width=400, height=1080) # Set session address before use
	my.colorFct <- function(n = 50, low.col = 0.5, high.col=0.17, saturation = 1) {
		if (n < 2) stop("n must be greater than 2")
		n1 <- n%/%2
		n2 <- n - n1
		c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2)))
	}

	# load heatmap.3 from url
	#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R", sha1 = "e0a2ae0c0a4c2ddaefb2b3e9bb4551c430535a98")
	out = heatmap.2(x,
			col=my.colorFct,
			## rows
			Rowv=hr.dd,
			RowSideColors=row.cols,
			## columns
			Colv=hc.dd,
			ColSideColors=col.cols,
			scale="row", # scale 기준방향 설정
			trace="non",
			## color key + density info
			key.title="Cell Colors",
			keysize=1,
			density.info="histogram")
	#dev.off()
	return()
}

tmp = ss.heatmap(data.arr.heat,
                 row.col=data.arr$Marker,
                 hr=FALSE, hc=TRUE,
                 hc.order=c(1,2,3,5,6,4,7,8,9,11,12,10,13,14,15))

write.table(heatmap,"heatmap.csv",sep=",",row.names=TRUE)

## Bibliography ##
# http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/heatmap/
# ?heatmap; http://127.0.0.1:23779/library/stats/html/heatmap.html
# ?heatmap.2; http://127.0.0.1:23779/library/gplots/html/heatmap.2.html


## ss.heatmap version 2
## v2.0  160721 - Fisrt written
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
install.packages("dendsort")
install.packages("dendextendRcpp")
install.packages("circlize")

ss.heatmap2 = function(data,colrange=NULL,hr=TRUE,row.col=NULL,row.col2=NULL,hc=TRUE,hc.order=NULL) {
	library(ComplexHeatmap); library(circlize); library(dendsort); library(dendextend)

	mat = as.matrix(data)
	mat.s = t(scale(t(mat))) # scale by row-lines
	mat.s_max = max(mat.s)
	mat.s_min = min(mat.s)
	col.rg = c(mat.s_min,0,mat.s_max)
	if(length(colrange)>0) {col.rg = colrange}

	colGroup = data.frame(Sample=c(rep("WL",3),rep("Og",3),rep("Bl",3),rep("Hw",3),rep("Ya",3)))
	ha = HeatmapAnnotation(df=colGroup) #,
				#col=list(Samples=c("WL"="FF0000","Og"="#CCFF00",
				#"Bl"="#00FF66","HW"="#0066FF","Ya"="#CC00FF")))
	#ra = rowAnnotation(row.col,
		#    col=list(Group=c("W"=col.fa[1],"O"=col.fa[2],"B"=col.fa[3],"H"=col.fa[4],"Y"=col.fa[5]),
		#    Marker=c("-"="grey80","mk"="yellow")), width=unit(1,"cm"))
	if(hc==TRUE) {
		hc = hclust(as.dist(1-cor(mat, method="spearman")),method="complete")
		hc.dd = as.dendrogram(hc)
	} else {hc.dd=FALSE}
	if(hr==TRUE) {
		hr = hclust(as.dist(1-cor(t(mat),method="pearson")),method="complete")
		hr.dd = as.dendrogram(hr)
	} else {hr.dd=FALSE}
	if(length(hc.order)>0) {
		print(order.dendrogram(hc.dd))
		hc.dd = rotate(hc.dd,hc.order) #rotate dendrogram leafs
		print(order.dendrogram(hc.dd))
	}
	cell.cols = colorRamp2(col.rg,c("Cyan","black","Yellow"))

	Heatmap(mat.s, col=cell.cols,
			name = "Probeset\nintensity\n(Scaled)",
			cluster_rows=hr.dd,cluster_columns=hc.dd,
			top_annotation=ha,
			column_dend_height=unit(3,"cm"),
			column_dend_reorder=FALSE,
			row_dend_width=unit(3,"cm"),
			row_names_side = "left",
			row_names_gp = gpar(fontsize=7)) +
	Heatmap(row.col2,name="Group",col=c("-"="black","+"="#990000"),width=unit(3,"cm"))+
	Heatmap(row.col,name="Marker",col=c("-"="grey80","marker"="yellow"),
			show_row_names=TRUE,
			row_names_gp = gpar(fontsize=7),
			width=unit(0.5,"cm"))
	#draw(ra + ht)
}
ss.heatmap2(data.arr.heat, hr=FALSE,hc=TRUE,
			row.col=Marker,
			row.col2=Group,
			hc.order=c("WL.1","WL.2","WL.3","Og.2","Og.3","Og.1",
			"Bl.1","Bl.2","Bl.3","Hw.2","Hw.3","Hw.1","Ya.1","Ya.2","Ya.3"))

ss.heatmap2(data.arr.heat,colrange=c(-2,0,2), hr=FALSE,hc=TRUE)

-----
## ss.heatmap version 3
## 161025 - For X-ALD project
## 161101 - Adapted for Yeast dynamic network project
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dendsort))
suppressPackageStartupMessages(library(dendextend))

sgl = read.csv("161025_cluster_heapmap_input.csv")
sgl.n = sgl[,-1]

mat = as.matrix(sgl.n)
mat.s = t(scale(t(mat)))
mat.s_max = max(mat.s); mat.s_min = min(mat.s)
col.rg = c(mat.s_min,0,mat.s_max)

colGroup = data.frame(Group=c(rep("AMN",4),rep("CCALD",4),rep("Normal",3)))
ha = HeatmapAnnotation(df=colGroup)

# hclust method="ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
# cor method="pearson", "kendall", "spearman"
# as.dist method="euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hc = hclust(as.dist(1-cor(mat, method="spearman")),method="complete")
hc.dd = as.dendrogram(hc)
hr = hclust(as.dist(1-cor(t(mat),method="spearman")),method="complete")
hr.dd = as.dendrogram(hr)

cell.cols = colorRamp2(col.rg,c("Cyan","black","Yellow"))
Heatmap(mat.s, col=cell.cols,
			name = "Probeset\nintensity\n(Scaled)",
			cluster_rows=hr.dd,cluster_columns=hc.dd,
			top_annotation=ha,
			column_dend_height=unit(3,"cm"),
			column_dend_reorder=FALSE,
			row_dend_width=unit(3,"cm"),
			row_names_side = "left",
			row_names_gp = gpar(fontsize=7))
