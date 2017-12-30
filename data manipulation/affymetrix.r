# CEL file process in R

## install package
source("https://bioconductor.org/biocLite.R")
biocLite("affy")

## load packages and CEL files
library(affy)
setwd("C:/Users/KimSS-Work/Desktop/GSE60596_RAW")
c = ReadAffy()

image(c) # draw image plot

## RMA
c_rma = rma(c) # RMA by rma
c_ex = expresso(c, bgcorrect.method="rma", # RMA by expresso
                normalize.method="quantiles",
                pmcorrect.method="pmonly",
                summary.method="medianpolish")

## Result visualization
# (1) extract expression table
cl = log2(exprs(c)) # before normalization
ces = exprs(c_ex) # after normalization

# (2) boxplot
par(mfrow=c(1,2),
    mar=c(17,2,3,1)) # bottom, left, top, right
boxplot(cl,main="Before",ylim=c(0,15),las=2)
boxplot(ces,main="After",ylim=c(0,15),las=2)

# (3) MA plot
mva.pairs(cl[,1:3])
mva.pairs(ces[,1:3])


# Data download from GEO DB using R code
# load series and platform data from GEO
## GSE60596: Data of mouse CR fat paper
library(Biobase)
library(GEOquery)

gset <- getGEO("GSE60596", GSEMatrix =TRUE, getGPL=FALSE)

# group names for all samples in a series
gsms <- "000333222111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

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

# 1. Group to colnames
group = c(rep("CD",3),rep("CR85",3),rep("CR70",3),rep("CR55",3))
group = factor(group, levels=c("CD","CR85","CR70","CR55"))

View(g.ex)

# 2. calculate ANOVA p-value
pval = apply(g.ex, 1, function(row) {
  fit = lm(row~group)
  anovaP = anova(fit)$'Pr(>F)'[1]
})

# 3. calculate Adjusted p-values
pval.adj = p.adjust(pval,"BH") # adjust by FDR

library(ggplot2)
p.df = stack(data.frame(pval,pval.adj))
ggplot(p.df,aes(x=ind,y=values))+theme_bw()+
  geom_boxplot(aes(fill=ind))+
  labs(title="Distribution of pval and p.adj (FDR)",
       x="Algorithms", y="p-values")

# 4. calculate fold change values
## average signal by groups
g.av = NULL
for(i in levels(group)) {
  a.av = apply(g.ex[,which(group==i)],1,mean)
  g.av = cbind(g.av,a.av)
}
colnames(g.av) = levels(group)

## calculate fold change
g.fc = NULL
for(i in 1:length(colnames(g.av)[-1])) {
  a.fc.log = g.av[,i+1]-g.av[,1] # log2 form
  a.fc = sapply(a.fc.log, function(x) ifelse(x>0,2^x,-(2^(-x)))
  g.fc = cbind(g.fc, a.fc)
}
colnames(g.fc) = levels(group)[-1]

# 5. Filter differentially expressed genes
## Statistical criteria by FDR <0.05 & |FC| >1.5
g.fc.max = apply(abs(g.fc),1,max)
id_deg = which(pval.adj<0.05 & g.fc.max>1.5)
length(id_deg)
head(id_deg)

# 6. Correlation histogram
cr = c(rep(0,3),rep(15,3),rep(30,3),rep(45,3))
g.deg = g.ex[id_deg,] # filter signals of DEGs
g.corr = apply(g.deg,1,function(row) {
  a.corr = cor(cr, row)
})
head(g.corr)
hist(g.corr, breaks=seq(-1,1,by=0.1),
     main="Histogram of correlation", xlab="Correlation r value")

# 7. SOM clustering
library(som)
g.deg.s = t(scale(t(g.deg)))
g.som = som(g.deg.s, xdim=4,ydim=4, topol="rect", neigh="gaussian")

plot(g.som, main="SOM clustering")

## box plot of a cluster
id_x = which(g.som$visual$x==3)
id_y = which(g.som$visual$y==0)
id_xy = intersect(id_x,id_y)
clst = g.deg.s[id_xy,]

boxplot(clst, ylim=c(-3,3), las=2,
        main=paste0("Pattern of cluster"))

# 8. Results export
p.adj.deg = pval.adj[id_deg]
g.som.xy = data.frame(g.som$visual$x, g.som$visual$y)

out = cbind(g.deg, p.adj.deg, g.fc.deg, g.corr, g.som.xy)
write.csv(out, "analysis.csv")
