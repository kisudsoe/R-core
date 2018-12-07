# 2017-07-04 TUE, ss.venn
- ver1.0 - 170704, Generate for Dog project
```r
ss.venn = function(grouplist, main="") { # 170704 ver - id for multi-groups
	# grouplist = list(group1,group2,group3)
	library(limma)
	# Generate union group from input three groups
	unionlist = Reduce(union, grouplist)
	# Get names from input vectors
	g.title = names(grouplist)
	
	# Togethering the vectors to one dataFrame list
	unionPr=NULL; title=NULL
 	for(i in 1:length(grouplist)) {
  		unionPr = cbind(unionPr,grouplist[[i]][match(unionlist,grouplist[[i]])])
   		title = c(title,paste(g.title[i],'\n',length(grouplist[[i]])))
  	}
  	rownames(unionPr) = unionlist

	# Make TRUE/FALSE table to match with ID list
	union = (unionPr != "") 			# Transform values to TRUE, if ID exists
	union[is.na(union)] = FALSE			# Transform NA to FALSE value
	union = as.data.frame(union) 		# Make 'union' to data.frame form
	colnames(union) = title 			# Names attach to venn diagram
	out = list(list=union, vennCounts=vennCounts(union))

  	## Generate Venn Diagram
	v = vennDiagram(union,main=paste0(main," (",nrow(union)," genes)"), circle.col=rainbow(length(g.title)))
  	print(v)
  	dev.copy(png,"ss.venn.png",width=8,height=8,units="in",res=100)
  	graphics.off()
	return(out)
}

ortho_egid_canis = list(ortho_ensem_egid_canis,ortho_inpara_egid_canis,ortho_homo_egid_canis)
names(ortho_egid_canis) = c("Ensembl_dog","Inparanoid_dog","HomoloGene_dog")
ortho_union_egid_canis = ss.venn(ortho_egid_canis, main="Dog-Yeast orthologs (Dog) from DBs")
#--------------------
affyid1 = unlist(read.delim("clipboard")); print(length(affyid1)) # Protein_catabolism
affyid2 = unlist(read.delim("clipboard")); print(length(affyid2)) # Mitochondrion
affyids = list(affyid1,affyid2); names(affyids) = c("Protein_catabolism","Mitochondrion")
out = ss.venn(affyids,main="Genes: Protein catabolism - Mitochondrion")
```



# 2016-05-04 WED, ss.tukey3Venn

- ver 1.0	- 160504, First version using for Human X-ALD project, Yeast project
- Ver 1.0a - 160505, Add function for print process information out
```r
ss.tukey3Venn = function(tukey) { # 160504 row[2] is centered.
	library(pbapply)
	time1 = Sys.time()
	# Make T/F table
	print(tukey[1:10,])
	out = NULL
	n = length(rownames(tukey))
	cat('\n------------------------------------------------------\n')
	cat('Process iteration =',n,'\n')
	out = pbapply(tukey,1,function(row){
		cp1 = ifelse(length(grep(row[2],row[1]))==0,TRUE,FALSE)
		cp2 = ifelse(length(grep(row[2],row[3]))==0,TRUE,FALSE)
		cp3 = ifelse(length(grep(row[1],row[3]))==0,TRUE,FALSE)
		output = c(cp1, cp2, cp3)
		return(output)
	})
	cat('------------------------------------------------------\n\n')
	out = t(out)
	coln = colnames(tukey)
	colnames(out) = c(paste(coln[2],"!=",coln[1],sep=""),
					  paste(coln[2],"!=",coln[3],sep=""),
					  paste(coln[1],"!=",coln[3],sep=""))
	print(out[1:10,])
	
	# Make Venn Diagram
	library(limma)
	vennDiagram(out) # Generate Venn Diagram
	dev.copy(png,"ss.tukey3Venn_fig.png",width=8,height=7,units="in",res=100)
  dev.off()

	# Make return table
	out = data.frame(tukey,out)
	time2 = Sys.time()
	print(time2-time1)
	return(out)
}

sgnl.cerev_tukey_venn = ss.tukey3Venn(sgnl.cerev_tukey)
```



# 2016-07-04, ss.venn3

- Install limma package through biocLite

```r
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
```

- Coded for general
- v1.0 2015 - original version
- v1.1 160704 - Add names of venn from vector
- v1.2 170606 - Bug fix for Atom
```r
ss.venn3 = function(group1, group2, group3, main="") { # 160704 ver - id for each 3-groups

	# Generate union group from input three groups
	unionlist = Reduce(union, list(group1, group2, group3))
	#unionlist = union(group1,group2)
	#unionlist = union(unionlist,group3)
	
	# Get names from input vectors
	title = c(names(group1)[1], names(group2)[1], names(group3)[1])
  	print(title)

	# Togethering the vectors to one dataFrame list
	unionPr = data.frame(list=unionlist,
						 g1=character(length(unionlist)),
						 g2=character(length(unionlist)),
						 g3=character(length(unionlist)))
	
	unionPr$g1 = group1[match(unionPr$list,group1)]
	unionPr$g2 = group2[match(unionPr$list,group2)]
	unionPr$g3 = group3[match(unionPr$list,group3)]
	
	rownames(unionPr) = unionPr[,1] 	# column 1 list as rownames
	#print(head(rownames(unionPr)))
	unionPr[1] = NULL 					# delete column 1
	
	# Make TRUE/FALSE table to match with ID list
	union = (unionPr != "") 			# Transform values to TRUE, if ID exits
	union[is.na(union)] = FALSE			# Transform NA to FALSE value
	
	union = as.data.frame(union) 		# Make 'union' to data.frame form
	title = c(paste(title[1],'\n',length(group1)),
			      paste(title[2],'\n',length(group2)),
			      paste(title[3],'\n',length(group3)))
	colnames(union) = title 			# Names attach to venn diagram
	
	library(limma)
	union = list(list=union, vennCounts=vennCounts(union))
	vennDiagram(union$list,main=main) 			# Generate Venn Diagram
  	dev.copy(png,"ss.venn3_fig.png",width=8,height=7,units="in",res=100)
  	dev.off()

	return(union)
}

sgd.Affy_venn = ss.venn2(sgd.id$gene.primaryIdentifier,
                         sgd.Affyid$SGD.accession.number)
```



# 2016-05-16 MON, ss.Venn4

- ver 1.0 - 160516, Coded for general
```r
ss.venn4 = function(group1,group2,group3,group4) {
	unionlist = union(group1,union(group2,union(group3,group4))) # Venn for 4-groups
	unionPr = data.frame(list=unionlist,
	                     g1=character(length(unionlist)),
						 g2=character(length(unionlist)),
						 g3=character(length(unionlist)),
	                     g4=character(length(unionlist)))
	unionPr$g1 = group1[match(unionPr$list,group1)]
	unionPr$g2 = group2[match(unionPr$list,group2)]
	unionPr$g3 = group3[match(unionPr$list,group3)]
	unionPr$g4 = group4[match(unionPr$list,group4)]
	rownames(unionPr) = unionPr[,1]
	unionPr[1] = NULL
	union = (unionPr != "") # TRUE, if id exist.
	union[is.na(union)] = FALSE # FALSE, if id null.
	print(head(union))
	
	library(limma)
	colnames(union) = c(length(group1),length(group2),length(group3),length(group4))
	union = list(list=union, vennCounts=vennCounts(union))
	vennDiagram(union$list)
	
	return(union)
}
```