

ss.venn3 = function(group1, group2, group3) { # data with 3-group columns
  
	# Generate union group from input three groups
	unionlist = union(group1,group2)
	unionlist = union(unionlist,group3)
  
	# Togethering the vectors to one dataFrame list
	unionPr = data.frame(list=unionlist,
						 g1=character(length(unionlist)),
						 g2=character(length(unionlist)),
						 g3=character(length(unionlist)))
  
	unionPr$g1 = group1[match(unionPr$list,group1)]
	unionPr$g2 = group2[match(unionPr$list,group2)]
	unionPr$g3 = group3[match(unionPr$list,group3)]
  
	rownames(unionPr) = unionPr[,1] # column 1 list as rownames
	unionPr[1] = NULL # delete column 1
  
	# Make TRUE/FALSE table to match with ID list
	union = (unionPr != "") # ID 있으면 모두 TRUE로 데이터 형 변환
	union[is.na(union)] = FALSE # NA는 모두 FALSE로 변환
  
	union = as.data.frame(union) # union을 data.frame 형식으로 변환
	colnames(union) = c(length(group1),length(group2),length(group3))
  
	#coln = colnames(data.frame(group1[1],group2[1],group3[1])) # Names of groups
	#colnames(union) = coln # Names attach to venn diagram
  
	library(limma)
	union = list(list=union, vennCounts=vennCounts(union))
	vennDiagram(union$list) # Generate Venn Diagram
  
	return(union)
}

fc_union_id = ss.venn3(fc_Rho.HD_id,fc_Rho.LD_id,fc_HD.LD_id)


ss.tukey3Venn2 = function(tukey,group) {
	# Make T/F table
	print(tukey[1:10,])
	out = NULL
	n = length(rownames(tukey))
	out = apply(tukey,1,function(row){
				cp1 = ifelse(length(grep(row[2],row[1]))==0,TRUE,FALSE)
				cp2 = ifelse(length(grep(row[2],row[3]))==0,TRUE,FALSE)
				cp3 = ifelse(length(grep(row[1],row[3]))==0,TRUE,FALSE)
				output = c(cp1, cp2, cp3)
				return(output)
			})
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
	return(out)
}