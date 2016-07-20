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


## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")

ss.heatmap2 = function(data,colrange=NULL,hr=TRUE,row.col=NULL,hc=TRUE,hc.order=NULL) {
	library(ComplexHeatmap); library(circlize)
	mat = as.matrix(data)
	mat.s = t(scale(t(mat))) # scale by row-lines
	mat.s_max = max(mat.s)
	mat.s_min = min(mat.s)
	col.rg = c(mat.s_min,0,mat.s_max)
	if(length(colrange)>0) {col.rg = colrange}

	colGroup = data.frame(Sample=c(rep("WL",3),rep("Og",3),rep("Bl",3),rep("Hw",3),rep("Ya",3)))
	col.fa = rainbow(5)
	print("col.fa=rainbow(5)")
	ha = HeatmapAnnotation(df=colGroup,
				col=list(Samples=c("WL"=col.fa[1],"Og"=col.fa[2],
				"Bl"=col.fa[3],"HW"=col.fa[4],"Ya"=col.fa[5])))
	
	if(length(row.col)>0) {
		ra = rowAnnotation(row.col,
		    col=list(Group=c("W"=col.fa[1],"O"=col.fa[2],"B"=col.fa[3],"H"=col.fa[4],"Y"=col.fa[5]), Marker=c("-"="grey80","mk"="yellow")), width=unit(1,"cm"))
	}

	if(hc==TRUE) {hc = hclust(as.dist(1-cor(mat, method="spearman")),method="complete")}
	if(hr==TRUE) {hclust(as.dist(1-cor(t(mat),method="pearson")),method="complete")}
	column_dend_reorder = hc.order #<-??
	cell.cols = colorRamp2(col.rg,c("Cyan","black","Yellow"))

	ht = Heatmap(mat.s, col=cell.cols,
			name = "Probeset\nintensity\n(Scaled)",
			cluster_rows=hr,cluster_columns=hc,
			top_annotation=ha,
			column_dend_height=unit(3,"cm"),
			column_dend_reorder=column_dend_reorder,
			row_dend_width=unit(3,"cm"))
	draw(ra + ht)
}
ss.heatmap2(data.arr.heat, hr=FALSE,hc=TRUE,
			row.col=data.frame(Marker,Group),
			hc.order=c(1,2,3,5,6,4,7,8,9,11,12,10,13,14,15))
			#hc.order=c(1:15))
ss.heatmap2(data.arr.heat,colrange=c(-2,0,2), hr=FALSE,hc=TRUE)