# Function start (draw 2-D boxplot)-----
#id = which(factor[,1] == "0-0")
#result[id,]

data=input(3) # scale data
group=input(1) # group by SOM or k-means

groupbox = function(data,group,i,j) {
	group = as.factor(group)
 
	#dev.new(with=10,height=10)
	par(mfrow=c(i,j)) # c(row,column)
	n = length(levels(group))
	print(paste("plot number=",n))
 
	for(k in 1:n) {
		sel = levels(group)[k]
		id = which(group==sel)
		title=paste(sel," (n=",length(id),")",sep="")
 
		if(k%%j==1) {
			par(mar=c(0,2,2,0)) # (bottom,left,up,right)
			boxplot(data[id,],
					main=title,
					xaxt='n',
					ylim=c(-3,3))
		} else {
			par(mar=c(0,1,2,1)) # (bottom,left,up,right)
			boxplot(data[id,],
					main=title,
					xaxt='n',
					yaxt='n', 
					ylim=c(-3,3))
		}
		rect(3.5,-3.5,6.5,3.5,col="grey90",lty=0)
		rect(9.5,-3.5,12.5,3.5,col="grey90",lty=0)
		rect(15.5,-3.5,18.5,3.5,col="grey90",lty=0)
		rect(21.5,-3.5,24.5,3.5,col="grey90",lty=0)
		rect(27.5,-3.5,30.5,3.5,col="grey90",lty=0)
		text(2,-2.8,"14")
		text(5,-2.8,"29")
		text(8,-2.8,"44")
		text(11,-2.8,"53")
		text(14,-2.8,"73")
		text(17,-2.8,"83")
		text(20,-2.8,"95")
		text(23,-2.8,"109")
		text(26,-2.8,"120")
		text(29,-2.8,"135")
 
		boxplot(data[id,],
				xaxt='n',
				yaxt='n',
				add=TRUE)
 
		#legend("bottomleft",title)
		#text(3,-2.5,title,font=5)
	}
 
}
 
groupbox(result,factor$som.4,3,2)
# Function end-----