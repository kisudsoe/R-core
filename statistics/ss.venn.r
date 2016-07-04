# 2016-05-04 WED
## ver 1.0	- 160504, First version using for Human X-ALD project, Yeast project
## Ver 1.0a - 160505, Add function for print process information out
ss.tukey3Venn = function(tukey,group) { # 160504 row[2] is centered.
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
	
	# Make return table
	out = data.frame(tukey,out)
	time2 = Sys.time()
	print(time2-time1)
	return(out)
}

sgnl.cerev_tukey_venn = ss.tukey3Venn(sgnl.cerev_tukey,group)


# Venn Diagram (3 circles)
## Install limma package through biocLite
source("https://bioconductor.org/biocLite.R")
biocLite("limma")

## Coded for general
## v1.0 2015 - original version
## v1.1 160704 - Add names of venn from vector
ss.venn3 = function(group1, group2, group3,main="") { # 160704 ver - id for each 3-groups
  
	# Generate union group from input three groups
	unionlist = union(group1,group2)
	unionlist = union(unionlist,group3)

	# Get names from input vectors
	title = c(names(group1)[1], names(group2)[1], names(group3)[1])
  
	# Togethering the vectors to one dataFrame list
	unionPr = data.frame(list=unionlist,
						 g1=character(length(unionlist)),
						 g2=character(length(unionlist)),
						 g3=character(length(unionlist)))
  
	unionPr$g1 = group1[match(unionPr$list,group1)]
	unionPr$g2 = group2[match(unionPr$list,group2)]
	unionPr$g3 = group3[match(unionPr$list,group3)]
  
	rownames(unionPr) = unionPr[,1] 	# column 1 list as rownames
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
  
	return(union)
}

fc_union_id = ss.venn3(fc_Rho.HD_id,fc_Rho.LD_id,fc_HD.LD_id)


# 2016-05-09 MON
## ver 1.0 - 160509, Coded for general
ss.venn2 = function(group1,group2) { #Venn for 2-groups
	unionlist = union(group1,group2)
	unionPr = data.frame(list=unionlist,
                         g1=character(length(unionlist)),
                         g2=character(length(unionlist)))
	unionPr$g1 = group1[match(unionPr$list,group1)]
	unionPr$g2 = group2[match(unionPr$list,group2)]
	rownames(unionPr) = unionPr[,1]
	unionPr[1] = NULL
	union = (unionPr != "") # TRUE, if id exist.
	union[is.na(union)] = FALSE # FALSE, if id null.
	print(head(union))
  
	library(limma)
	colnames(union) = c(length(group1),length(group2))
	union = list(list=union, vennCounts=vennCounts(union))
	vennDiagram(union$list)
  
	return(union)
}

sgd.Affy_venn = ss.venn2(sgd.id$gene.primaryIdentifier,
                         sgd.Affyid$SGD.accession.number)
						 
# 2016-05-16 MON
## ver 1.0 - 160516, Coded for general
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