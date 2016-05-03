### [R] Clustering.R ###
# First written 2015-08-21 FRI
# Last update 2016-05-02 MON
# Version info 1.0a

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

data = input(3)
 
ss.kmeans = function(data,i) {
	km = kmeans(t(scale(t(data))),i) 
		# K-menas clustering
		# row-wise scaling returns a mean close to zero and 
		# a standard deviation of one.
	#km$cluster
	result.km = as.matrix(km$cluster)
 
	# Designate color to each k-cluster #
	mycol <- sample(rainbow(256))
	mycol <- mycol[as.vector(km$cluster)]
 
	# Hierarchical clustering #
	data.mat = as.matrix(data)
	hr = hclust(as.dist(1-cor(t(data.mat),method="pearson")),method="complete")
	hc = hclust(as.dist(1-cor(data.mat, method="spearman")),method="complete")
		# Correlation methods = "Pearson", "Kendall", "Spearman"
		# hclust methods = "complete", "average", "median", "centroid", etc.
	source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
	#heatmap(data.mat, 
	#        Rowv=as.dendrogram(hr), 
	#        Colv=as.dendrogram(hc), 
	#        col=my.colorFct(), 
	#        scale="row")
  
	# Draw heatmap #
	heatmap(data.mat, 
			Rowv=as.dendrogram(hr), 
			Colv=as.dendrogram(hc), 
			col=my.colorFct(), 
			scale="row", 
			#ColSideColors=heat.colors(length(hc$labels)), 
			RowSideColors=mycol) 
		# This heatmap shows how the obtained K-means clusters can be 
		# compared with the hierarchical clustering results by highlighting 
		# them in the heatmap color bar.
  
	return(result.km)
}
 
result.km = ss.kmeans(data,2)
 
write.table(result.km,"kmeans.csv",sep=",",row.names=TRUE)