# ss.heatmap



## 180129 Yeast HD-LD project 2

update several options, modify to display only a subcategory and so on..

```R
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap") # 왜 안되지?

ss.heatmap4 = function(mat,title="", type=NULL, annot=NULL, location=F, hr.order=NULL,
                       split=NULL, split_order=NULL, ra_width=10, width=7, height=8, img_name=NULL) {
  library(ComplexHeatmap)
  tukey_order = c("cab","baa","cba","abc","abb","acb", "aba","bab", "bca","bba","aab","bac")
  #if(type=="metab") tukey_order = c("b.b.a","b.a.a","c.b.a","c.a.b","b.a.c", "a.b.a","a.c.b","a.b.c","a.a.b")
  # Sort matrix by tukey and HD/LD FC
  if(type=="metab" && nrow(mat)>1) {
    mat.s1 = t(scale(t(mat[,1:24])))
    mat.s1.HD = apply(mat.s1[,1:8],1,mean)
    mat = mat[order(factor(mat$Tukey,levels=tukey_order),-mat.s1.HD),]
    mat_name = apply(data.frame(mat$Tukey2,mat$Name2),1,function(row) {
      #if(as.numeric(row[1])<2 & as.numeric(row[1])>-2) paste0(row[2],"*")
      if(row[1]%in%c("b.b.a","a.a.b")) paste0(row[2],"*")
      else row[2]
    })
  } else if(nrow(mat)>1) {
    mat.s1 = t(scale(t(mat[,2:16])))
    mat.s1.HD = apply(mat.s1[,1:5],1,mean); #mat.s1.LD = apply(mat.s1[,6:10],1,mean)
    mat = mat[order(factor(mat$Tukey,levels=tukey_order),-mat.s1.HD),] #-mat$`HD/LD`, -mat.s1.LD
  }

  # Concatenate Row names
  if(type=="none") {
    mat_ = mat[,2:16]; tukey_ = data.frame(Tukey=mat$Tukey2)
    annot_ = paste0(mat$AFFYid," (",mat$Symbol1,") ",mat$Tukey2)
    row.names(mat_) = annot_
    row.names(tukey_) = annot_
  } else if(type%in%c("some","summ")) {
    mat_ = mat[,2:16]; tukey_ = data.frame(Tukey=mat$Tukey2)
    if(location) location_ = mat[,24:ncol(mat)]
    ann1 = paste0(mat$Tukey2,"-",mat$AFFYid," (",mat$Symbol1,") ")
    if(annot=="Category") annot_ = paste0(ann1,mat$Category)
    else if(annot=="Cat0") annot_ = paste0(ann1,mat$Cat0)
    else if(annot=="Cat1") annot_ = paste0(ann1,mat$Cat1)
    else if(annot=="Cat2") annot_ = paste0(ann1,mat$Cat2)
    else if(annot=="Cat3") annot_ = paste0(ann1,mat$Cat3)
    else if(annot=="Cat") annot_ = paste0(ann1,mat$Cat)
    else if(annot=="Name") annot_ = ann1
    else annot_ = paste0(mat$AffyID," ",mat$SGDid)
    row.names(mat_) = annot_
    row.names(tukey_) = annot_#mat$SGDid
    if(location) row.names(location_) = annot_
  } else if(type=="metab") {
    mat_ = mat[,1:24]; tukey_ = data.frame(Tukey=mat$Tukey)
    ann2 = paste0(mat$Tukey,"-",mat$KEGG," (",mat$Name.y,")")
    if(annot=="Category") annot_ = paste0(mat$Category)
    else if(annot=="Cat0") annot_ = paste0(mat$Cat0)
    else if(annot=="Cat1") annot_ = paste0(mat$Cat1)
    else if(annot=="Cat2") annot_ = paste0(mat$Cat2)
    else if(annot=="Cat3") annot_ = paste0(mat$Cat3)
    else if(annot=="Name") annot_ = paste0(mat_name)
    else annot_ = paste0(mat$Category,", ",mat$Name)
    annot_ = paste0(ann2," ",annot_)
    row.names(mat_) = annot_
    row.names(tukey_) = annot_
  } else {
    mat_ = mat
  }
  
  # heatmap split by annotation
  if(split=="Tukey") {split = factor(mat$Tukey2); print(">> Split by tukey <<")}
  else if(split=="Location") {split = factor(mat$location); print(">> Split by Location <<")}
  else if(split=="Cat0") {split = factor(mat$Cat0); print(">> Split by Cat0 <<")}
  else if(split=="Cat1") {split = factor(mat$Cat1); print(">> Split by Cat1 <<")}
  else if(split=="Cat2") {split = factor(mat$Cat2); print(">> Split by Cat2 <<")}
  if(!is.null(split_order)) {split = factor(split,levels=split_order)}
  mat.s = t(scale(t(mat_)))

  #install.packages('circlize',repos="http://cran.us.r-project.org")
  library(circlize)
  #mat.s_max = max(mat.s); mat.s_min = min(mat.s)
  #col.rg = c(mat.s_min,0,mat.s_max)
  col.rg = c(-2.2,0,2.2)
  if(type=="metab") cell.cols = colorRamp2(col.rg,c("Green","Black","Red"))
  else cell.cols = colorRamp2(col.rg,c("Cyan","black","Yellow"))

  if(type!="metab") {
    hc = hclust(as.dist(1-cor(mat_,method="pearson")),method="complete")
    hc.dd = as.dendrogram(hc)
    if(!is.null(hr.order)) hc.dd = reorder(hc.dd,hr.order,mean)
  }
  
  if(type%in%c("none","summ")) {
    hm = Heatmap(matrix              = mat.s,
                 column_title        = title,
                 name                = "Z-score",
                 col                 = cell.cols,
                 cluster_columns     = F, #hc.dd,
                 cluster_rows        = F,
                 column_dend_reorder = F,
                 column_dend_height  = unit(2,"cm"),
                 column_names_side   = "top",
                 show_row_names      = F,
                 split               = split,
                 row_title_gp        = gpar(cex=1),
                 row_names_gp        = gpar(cex=1),
                 row_names_max_width = unit(ra_width,"cm"),
                 show_heatmap_legend = T)
  } else if(type=="some") {
    hm = Heatmap(matrix              = mat.s,
                 column_title        = title,
                 name                = "z-score",
                 col                 = cell.cols,
                 cluster_columns     = F, #hc.dd,
                 cluster_rows        = F,
                 column_dend_reorder = F,
                 column_dend_height  = unit(1,"cm"),
                 column_names_side   = "top",
                 row_dend_width      = unit(1,"cm"),
                 show_row_names      = F,
                 split               = split,
                 row_title_gp        = gpar(cex=1),
                 row_names_gp        = gpar(cex=1),
                 row_names_max_width = unit(ra_width,"cm"),
                 show_heatmap_legend = T)
    #hm = hm1+ra
    #print(draw(hm, heatmap_legend_side = "left"))
  } else if(type=="metab") {
    hm = Heatmap(matrix              = mat.s,
                 column_title        = title,
                 name                = "z-score",
                 col                 = cell.cols,
                 cluster_columns     = F, #hc.dd,
                 cluster_rows        = F,
                 column_dend_reorder = F,
                 column_dend_height  = unit(1,"cm"),
                 column_names_side   = "top",
                 show_row_names      = F,
                 split               = split,
                 row_title_gp        = gpar(cex=1),
                 row_names_gp        = gpar(cex=1),
                 row_names_max_width = unit(ra_width,"cm"),
                 show_heatmap_legend = T)
  } else {
    hm = Heatmap(matrix              = mat.s,
                 column_title        = title,
                 name                = "Probeset\nintensity\n(Scaled)",
                 col                 = cell.cols,
                 cluster_columns     = hc.dd,
                 column_dend_reorder = F,
                 column_dend_height  = unit(2,"cm"),
                 column_names_side   = "top",
                 row_dend_width      = unit(2,"cm"),
                 show_row_names      = F,
                 split               = split,
                 row_title_gp        = gpar(cex=0.8))
    #print(hm)
  }
  
  # 2nd heatmap for tukey pattren
  # Tukey pattern heatmap
  col2 = c("bba"="goldenrod1","aab"="dodgerblue", # Rho-specific genes
           "bca"="goldenrod1","bac"="dodgerblue",  # Buffered Rho genes
           "aba"="goldenrod1","bab"="dodgerblue", # LD-specific genes
           "cab"="lightpink", "acb"="lightblue1",  # Buffered HD genes
           "baa"="red",       "abb"="blue",  # HD-specific genes
           "cba"="red",       "abc"="blue")
  if(type%in%c("none","summ")) show_row_names = F
  else if(location) show_row_names = F
  else show_row_names = T
  h2=Heatmap(tukey_, name="Tukey",
             cluster_columns=F,
             cluster_rows=F,
             show_row_names=show_row_names,
             column_names_side="top",
             row_names_max_width = unit(ra_width,"cm"),
             col=col2)
  # 3rd heatmap for location table
  if(location) {
    col3 = c("YES"="forestgreen","NO"="grey80")
    if(type%in%c("none","summ")) show_row_names = F
    else show_row_names = T
    h3=Heatmap(location_, name="Location",
               cluster_columns=F,
               cluster_rows=F,
               column_names_side="top",
               show_row_names=show_row_names,
               row_names_max_width = unit(ra_width,"cm"),
               col=col3)
  }
  
  # Draw heatmap
  if(type=="summ") {
    subset = sample(11811,500)
    labels = mat$Cat0[subset]
    ra = rowAnnotation(link=row_anno_link(at=subset,labels=labels), width=unit(1,"cm")+max_text_width(labels))
    if(location) print(hm+h2+h3+ra)
    else print(hm+h2+ra)
  } else {
    
    if(location) print(hm+h2+h3)
    else print(hm+h2)
  }
  
  if(!is.null(img_name)) {
    dev.copy(png,img_name,width=width,height=height,units="cm",res=100)
    dev.off()
    cat("Draw heatmap done.\n")
  }
}

# Function for Location table by Affyid
ss.affyLocation = function(annot_) {
    # strsplit for location of genes
    locations = annot_$Location
    locations_ = NULL
    for(i in 1:length(locations)) {
        if(length(grep(", ",locations[i]))>0) {
            loc_split = strsplit(as.character(locations[i]),split=", ")
            rows = data.frame(id=rep(annot_$AFFYid[i],length(loc_split)),loc=loc_split)
            colnames(rows) = c("id","location")
            locations_ = rbind(locations_,rows)
        } else {
            row = data.frame(id=annot_$AFFYid[i],location=locations[i])
            locations_ = rbind(locations_,row)
        }
    }
    locations_ = unique(locations_)
    loc.df = data.frame(table(factor(as.character(locations_$loc))))
    #loc.df = loc.df[order(loc.df$Freq,decreasing=T),] # Custom location order
    print(loc.df)
    
    # affyid list of locations
    affyid_list=NULL; loc_names=NULL;
    loc = loc.df$Var1
    for(i in 1:length(loc)) {
        loc_sub = subset(locations_,location %in% loc[i])
        loc_ids = as.character(unique(loc_sub$id))
        if(i==1&length(loc_ids)==1) affyid_list = list(loc_ids) # bugfix_171114
        else affyid_list[[i]] = loc_ids
        loc_names = c(loc_names,paste0(loc[i]," (",length(loc_ids),")"))
    }
    names(affyid_list) = loc_names
    
    # Union table generation
    affyid_uni = Reduce(union,affyid_list)
    affyid_tbl = NULL
    for(i in 1:length(loc_names)) {
        ids = affyid_list[[i]]
        affyid_tbl = cbind(affyid_tbl, ids[match(affyid_uni,ids)])
    }
    affyid_df = (affyid_tbl!="") # If value exist, then TRUE
    affyid_df[!is.na(affyid_tbl)] = "YES"
    affyid_df[is.na(affyid_df)] = "NO" # If value is null, then FALSE
    rownames(affyid_df) = affyid_uni; colnames(affyid_df) = loc_names
    return(affyid_df)
}
```





## 170609 Yeast HD-LD project

```r
ss.heatmap4 = function(mat,title="", hr.order=NULL, split=NULL, img_name=NULL) {
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("ComplexHeatmap")
  library(ComplexHeatmap)
  mat = mat
  mat.s = t(scale(t(mat)))

  #n = levels(factor(dat$tukey))
  #colr = sample(rainbow(256),size=length(n))
  #colr.grp = colr[as.numeric(deg_1738$tukey)]

  #install.packages('circlize',repos="http://cran.us.r-project.org")
  library(circlize)
  mat.s_max = max(mat.s); mat.s_min = min(mat.s)
  col.rg = c(mat.s_min,0,mat.s_max)
  #col.rg = c(-1.5,0,1.5)
  cell.cols = colorRamp2(col.rg,c("Cyan","black","Yellow"))

  hc = hclust(as.dist(1-cor(mat,method="pearson")),method="complete")
  hc.dd = as.dendrogram(hc)
  if(!is.null(hr.order)) hc.dd = reorder(hc.dd,hr.order,mean)

  hm = Heatmap(matrix=mat.s,
               column_title=title,
               name="Probeset\nintensity\n(Scaled)",
               col=cell.cols,
               cluster_columns=hc.dd,
               column_dend_reorder=F,
               column_dend_height=unit(2,"cm"),
               row_dend_width=unit(2,"cm"),
               show_row_names=F,
               split=split,
               row_title_gp=gpar(cex=0.8))
  print(hm)
  if(!is.null(img_name)) {
    dev.copy(png,img_name,width=7,height=8,units="in",res=100)
    dev.off()
    cat("Draw heatmap done.\n")
  }
}
```

```R
ss.heatmap4(deg_1738[,2:16],
            title="Heatmap of DEG 1,738 by Tukey",
            hr.order=c(1:5,6:10,11:15),
            split=tukey,
            img_name="170609_heatmap_deg1738_tukey.png")
```



## 161104 Yeast dynamic network

```R
ss.heatmap3 = function(mat,group.info,clstr=NULL,hc.order=NULL,yeast.cerev_somgrp=NULL) {
  library(ComplexHeatmap); library(circlize)
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
```

```R
ss.heatmap2(mat=yeast.cerev.sgl[1],group.info=group.info,
           clstr=c("A01","A02","A03","A04","A05",
                   "A06","A07","A08","A09","A10"),
           hc.order=c(1:12,31:51,53:60,52,28:30,13:27),
           yeast.cerev_somgrp=yeast.cerev_somgrp)
```



## 161025 X-ALD project

* ss.heatmap version 3
* 161101 Adapted for Yeast dynamic network project

```R
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
```



## 160721 Chicken project

- ss.heatmap version 2

```R
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
```

```R
ss.heatmap2(data.arr.heat, hr=FALSE,hc=TRUE,
			row.col=Marker,
			row.col2=Group,
			hc.order=c("WL.1","WL.2","WL.3","Og.2","Og.3","Og.1", "Bl.1","Bl.2","Bl.3","Hw.2","Hw.3","Hw.1","Ya.1","Ya.2","Ya.3"))
ss.heatmap2(data.arr.heat,colrange=c(-2,0,2), hr=FALSE,hc=TRUE)
```



## 160720 version 1

```R
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
```

```R
tmp = ss.heatmap(data.arr.heat,
                 row.col=data.arr$Marker,
                 hr=FALSE, hc=TRUE,
                 hc.order=c(1,2,3,5,6,4,7,8,9,11,12,10,13,14,15))
write.table(heatmap,"heatmap.csv",sep=",",row.names=TRUE)
```

