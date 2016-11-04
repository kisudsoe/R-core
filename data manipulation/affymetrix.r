## Boxplot for selected GEO samples ##

library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE60596", GSEMatrix =TRUE, getGPL=FALSE)

# group names for all samples in a series
gsms <- "000333222111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("CD","CR85","CR70","CR55")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf","#f2cb98","#dcdaa5", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE60596", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
