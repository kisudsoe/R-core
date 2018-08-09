### [R] Clustering.R ###
# v1.0  - 150821 FRI,First written
# v1.0a - 160502 MON, update
# v1.1  - 170606, Update for Atom
# to find 'k' how about to use NbClust?

## Bibliography ##
# http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#R_clustering

## K-means Clustering ##

##################################################
# data =
# ID	CD	CD	CD	CR85	CR85	CR85	CR70	CR70	CR70	CR55	CR55	CR55
# 10413492	60.09154042	46.09345481	48.55203969	78.32549082	30.5870698	55.6334242	87.25706243	84.94929092	63.63462164	618.9117923	368.4589277	200.0532145
# 10605666	76.37454632	39.64590628	55.2353513	20.51944777	19.68916467	32.4590754	3.990280042	3.755330806	11.47196005	3.563932468	5.039008885	3.592305385
# 10537062	4654.788816	3379.513823	3827.535385	815.5728843	2723.601656	2901.533761	392.9195114	275.7429136	751.6201844	58.35158746	116.8885658	71.89053824
# 10419568	146.9885533	121.0265701	136.6896313	16.03664006	48.78277761	35.32677143	18.50149569	25.46811008	14.80025369	18.43774098	15.73385186	19.84632529
# 10467744	419.7190108	315.4570974	446.5630275	75.72774164	325.424712	409.2021741	45.49323114	49.07274265	60.52922608	41.11434671	29.76242139	42.24415982
##################################################

ss.kmeans = function(data,k,hr=TRUE,hc=TRUE) {
	km = kmeans(t(scale(t(data))),k)
		# K-menas clustering
		# row-wise scaling returns a mean close to zero and
		# a standard deviation of one.
	#km$cluster
	result.km = as.matrix(km$cluster)

	# Designate color to each k-cluster #
	mycol <- sample(rainbow(256))
	mycol <- mycol[as.vector(km$cluster)]

	# Hierarchical clustering #
	## Correlation methods = "Pearson", "Kendall", "Spearman"
	## hclust methods = "complete", "average", "median", "centroid", etc.
	data.mat = as.matrix(data)
	if(hr==TRUE) {
		hr = hclust(as.dist(1-cor(t(data.mat),method="pearson")),method="complete")
		hr.dd = as.dendrogram(hr)
	} else { hr.dd=FALSE }
	if(hc==TRUE) {
		hc = hclust(as.dist(1-cor(data.mat, method="spearman")),method="complete")
		hc.dd = as.dendrogram(hc)
	} else { hc.dd=FALSE }

	#source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
	# Color function to generate blue-yellow heat maps
	my.colorFct <- function(n=50, low.col=0.5, high.col=0.17, saturation=1) {
		if (n<2) stop("n must be greater than 2")
		n1 <- n%/%2
		n2 <- n - n1
		c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2)))
	}

	# Draw heatmap #
	heatmap(data.mat,
			    Rowv=hr.dd,
			    Colv=hc.dd,
			    col=my.colorFct(),
			    scale="row",
			    #ColSideColors=heat.colors(length(hc$labels)),
			    RowSideColors=mycol)
	# This heatmap shows how the obtained K-means clusters can be
	# compared with the hierarchical clustering results by highlighting
	# them in the heatmap color bar.
  dev.copy(png,"ss.kmeans_fig.png",width=8,height=8,units="in",res=100)
  dev.off()

	return(result.km)
}

result.km = ss.kmeans(data,2)

write.table(result.km,"kmeans.csv",sep=",",row.names=TRUE)
