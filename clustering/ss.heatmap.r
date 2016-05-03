## Daw a Heat Map ##

# heatmap(x, Rowv = NULL, Colv = if(symm)"Rowv" else NULL,
#         distfun = dist, hclustfun = hclust,
#         reorderfun = function(d, w) reorder(d, w),
#         add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
#         scale = c("row", "column", "none"), na.rm = TRUE,
#         margins = c(5, 5), ColSideColors, RowSideColors,
#         cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
#         labRow = NULL, labCol = NULL, main = NULL,
#         xlab = NULL, ylab = NULL,
#         keep.dendro = FALSE, verbose = getOption("verbose"), ...)

# rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
# heat.colors(n, alpha = 1)
# terrain.colors(n, alpha = 1)
# topo.colors(n, alpha = 1)
# cm.colors(n, alpha = 1)

ssheatmap = function(data) {
	library("gplots")
 
	x = as.matrix(data)
	color.group = c(rep("#FF0000",3),rep("#FFFF00",3),rep("#FF00FF",3),rep("#00FFFF",3),rep("#888800",3))
 
	#png(filename="heatmap.png", width=400, height=1080) # 실행전에 작업공간 설정할 것
	out = heatmap.2(x,
					col=topo.colors(100),
					#ColSideColors=color.group,
					scale = "row", # scale 기준방향 설정
					trace="non",
					density.info="none"
					)
	#dev.off()

	return(t(out$carpet))
}
 
heatmap.fc = ssheatmap(data.fc)
 

## Bibliography ##
# http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/heatmap/
# ?heatmap; http://127.0.0.1:23779/library/stats/html/heatmap.html
# ?heatmap.2; http://127.0.0.1:23779/library/gplots/html/heatmap.2.html

write.table(heatmap,"heatmap.csv",sep=",",row.names=TRUE)