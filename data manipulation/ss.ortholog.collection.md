---
title:  "Yeast orthologs"
author: "KimSS"
date:   "2017-07-03 MON"
note:   "This code was written for ortholog collections from Ensembl Compara, Inparanoid and NCBI HomoloGene"
---

# Ortholog Collection

Code logs
---------

v1.1  - 161125, Yeast network project
v1.2  - 170417, Mice CR liver project
v1.3  - 170703, Dog project

* Get R packages for useMart
* https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html

```R
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("hom.Sc.inp.db")
#biocLite("hom.Hs.inp.db")
biocLite("annotationTools")
biocLite("org.Hs.eg.db")
biocLite("org.Sc.sgd.db")
```

* Custom functions

```R
## Ensembl search of 'Ortho_inpara' to get EntrezGene ID
ens.getBM = function(query_ids,attri=NULL,filter,marts,labels){
  if(length(attri)==0) attributes = c('ensembl_peptide_id','ensembl_gene_id',
                 'external_gene_name','entrezgene')
  else attributes = attri
  i = 1; n = length(marts)
  out = data.frame()
  for(mart in marts) {
    rst = getBM(attributes=attributes, filters=filter,
                values=query_ids, mart=mart)
    speci = rep(labels[i],length(rst[,1]))
    rst = data.frame(rst,Species=speci)
    out = rbind(out,rst)
    cat(paste(i,'=',labels[i],'>>',length(rst[,1]),'\n'))
    i = i+1
  }
  return(out)
}

## Compiling EntrezGene ID annotation to Inparanoid orthologs
ss.annot = function(origin,annot,colRange=c(2:5)){ # First columns of both entries should be IDs
  out = data.frame()
  i=1; n=nrow(origin); time1=Sys.time()
  pb=winProgressBar(title="Progress",label="description",min=0,max=n,width=500)
  print(paste0('iteration= ',n))
  for(id in origin[,1]) {
    ann_row = which(annot[,1]==id)
    annNum = length(ann_row)
    # 160717_bugfix: fill blanks
    if(length(apply(origin[i,],2, function(x) which(x=="")))>0){
      for(j in 1:length(origin[i,])) {
        if(origin[i,j]==""){ origin[i,j]="NA" }
    } } # end script #
    if(annNum==0) {
      annot_row = data.frame(origin[i,],ensembl_gene_id=NA,
                             external_gene_name=NA,entrezgene=NA,
                             Species=NA,annNum=annNum)
    } else {
      annot_row = data.frame(origin[i,],annot[ann_row,colRange],annNum=annNum)
    }
    out = rbind(out,annot_row)
    ## Progress time ##
  	d=difftime(Sys.time(),time1,unit="sec")
  	if(d<60) d=paste(round(d,2),"sec")
  	else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
  	else if(d>=3600) d=paste(round(d/3600,1),"hr")
  	setWinProgressBar(pb,i,title=paste0(round(i/n*100,1)," % (",i,"/",n,") done for ",d))
  	###################
    i = i+1
  }
  print(difftime(Sys.time(),time1,units="auto")); close(pb)
  return(out)
}

## Filter HomoloGene DB by taxonomy
homo.filter = function(homo_db, tax1="9606", tax2="4932") {
  time1=Sys.time()
  rows1 = which(homo_db$TaxId==tax1)
  rows2 = which(homo_db$TaxId==tax2)
  rows_union = union(rows1,rows2)

  hid1 = homo_db$HID[rows1]
  hid2 = homo_db$HID[rows2]
  hid = intersect(hid1, hid2)

  n = length(hid)
  pb=winProgressBar(title="Progress",label="homo.filter function",min=0,max=n,width=500)
  out = data.frame()
  for(i in 1:n){
    hid_rows = which(homo_db$HID==hid[i])
    rows = intersect(rows_union,hid_rows)
    out = rbind(out,homo_db[rows,])

    ## Progress time ##
  	d=difftime(Sys.time(),time1,unit="sec")
  	if(d<60) d=paste(round(d,2),"sec")
  	else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
  	else if(d>=3600) d=paste(round(d/3600,1),"hr")
  	setWinProgressBar(pb,i,title=paste0(round(i/n*100,1)," % (",i,"/",n,") done for ",d))
  	###################
  }
  print(difftime(Sys.time(),time1,units="auto")); close(pb)
  return(out)
}

# 170704 ver - grouplist = list(group1,group2,group3) <- id list
ss.venn = function(grouplist, main="") {
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
	union = (unionPr != "") 			# Transform values to TRUE, if ID exits
	union[is.na(union)] = FALSE			# Transform NA to FALSE value
	union = as.data.frame(union) 		# Make 'union' to data.frame form
	colnames(union) = title 			# Names attach to venn diagram
	out = list(list=union, vennCounts=vennCounts(union))

  ## Generate Venn Diagram
	v = vennDiagram(union$list,main=paste0(main," (",nrow(union)," genes)"))
  print(v)
  dev.copy(png,"ss.venn3.png",width=8,height=8,units="in",res=100)
  graphics.off()
	return(out)
}
```

# 1. Orthologs from Ensembl Compara

```R
library(biomaRt)
human_ds = "hsapiens_gene_ensembl"; human_fl = "with_homolog_scer"
mouse_ds = "mmusculus_gene_ensembl"; mouse_fl = "with_scerevisiae_homolog"
canis_ds = "cfamiliaris_gene_ensembl"; canis_fl = "with_scerevisiae_homolog"
yeast_ds = "scerevisiae_gene_ensembl"

ensem1 = useMart("ensembl", dataset=yeast_ds)
ensem2 = useMart("ensembl", dataset=canis_ds) # 170704

## To see list of ensembl data ##
ensembl=useMart("ensembl")
listDatasets(ensembl)             # To find database list
print(listAttributes(ensem2)$name) # To find attribute list
listFilters(ensem2)                # To find filter list
##-----------------------------##

#attributes = c('ensembl_gene_id','hgnc_symbol','entrezgene')
attributes = c('ensembl_gene_id','external_gene_name','entrezgene')
attributesL = c('ensembl_gene_id','external_gene_name','entrezgene')
ortho_ensem = getLDS(attributes,filters=canis_fl,
                     values=TRUE,mart=ensem2, # values=TRUE : all exist values
                     attributesL=attributesL,martL=ensem1)
write.csv(ortho_ensem,"ss.ortholog/ortho_ensem.csv",row.names=F)
ortho_ensem = read.csv("ss.ortholog/ortho_ensem.csv")
head(ortho_ensem)

# Ensembl orthologs
#ortho_ensem_ensid_human = na.omit(unique(ortho_ensem$Gene.stable.ID)) # bugfix 170417
ortho_ensem_ensid_canis = na.omit(unique(ortho_ensem$Gene.stable.ID)) # 170703
ortho_ensem_ensid_yeast = na.omit(unique(ortho_ensem$Gene.stable.ID.1))
print(length(ortho_ensem_ensid_canis))
print(length(ortho_ensem_ensid_yeast))

# Orthologs Entrez ID
#ortho_ensem_egid_human = na.omit(unique(ortho_ensem$NCBI.gene.ID)) # bugfix 170417
ortho_ensem_egid_canis = na.omit(unique(ortho_ensem$NCBI.gene.ID)) # 170703
ortho_ensem_egid_yeast = na.omit(unique(ortho_ensem$NCBI.gene.ID.1))
print(length(ortho_ensem_egid_canis))
print(length(ortho_ensem_egid_yeast))

# > head(ortho_ensem)
#   Ensembl.Gene.ID HGNC.symbol EntrezGene.ID Ensembl.Gene.ID.1 Associated.Gene.Name EntrezGene.ID.1
# 1 ENSG00000155366        RHOC           389           YPR165W                 RHO1          856294
# 2 ENSG00000119688       ABCD4          5826           YPL147W                 PXA1          855956
# 3 ENSG00000090238       YPEL3         83719           YBL049W                 MOH1          852231
# 4 ENSG00000184985      SORCS2         57537           YCR101C                               850465
# 5 ENSG00000152556        PFKM          5213           YGR240C                 PFK1          853155
# 6 ENSG00000099365       STX1B        112755           YMR183C                 SSO2          855221
```

# 2. Orthologs from Inparanoid yeast

- https://bioconductor.org/packages/release/data/annotation/html/hom.Hs.inp.db.html

* https://bioconductor.org/packages/release/data/annotation/html/hom.Sc.inp.db.html

```R
library(DBI) # To utilize sql query

library(hom.Hs.inp.db) # Human
library(hom.Mm.inp.db) # Mouse
library(hom.Sc.inp.db) # Yeast for Dog project, 170704

## Check data ##
print(ls("package:hom.Sc.inp.db"))
print(dbListTables(hom.Sc.inp_dbconn())) # List species
dbListFields(hom.Sc.inp_dbconn(),"Canis_familiaris") # Search species

tmp_inpara = dbGetQuery(hom.Sc.inp_dbconn(), "SELECT * FROM Canis_familiaris;")
print(dim(tmp_inpara))
head(ortho_inpara2)
#as.list(hom.Hs.inpSACCE)[1:4] # What data in this?
#print(head(as.list(hom.Sc.inpCANFA)))
#mget("ENSCAFP00000007145",hom.Sc.inpCANFA) # Example for catch a data, 외않되?
##------------##

## sql query grammers ##
# tmpsql = "SELECT * FROM Canis_familiaris WHERE inp_id='ENSCAFP00000007145';"
# tmpsql = "SELECT * FROM Saccharomyces_cerevisiae WHERE inp_id='YDR288W';"
# tmpsql = "SELECT * FROM Saccharomyces_cerevisiae WHERE clust_id='2143';"
##--------------------##

ortho_inpara = dbGetQuery(hom.Sc.inp_dbconn(), "SELECT * FROM Canis_familiaris;")

## get orthologs from Inparanoid yeast
library(hom.Sc.inp.db)

speci1="SACCE"; label1="S_cerevisiae"
speci2="CANFA"; label2="CL_familiaris"
## Ensembl search of 'Ortho_inpara' to get EntrezGene ID
attributesI = c('ensembl_peptide_id','ensembl_gene_id','external_gene_name','entrezgene')
ortho_inpara_names = ens.getBM(query_ids=ortho_inpara$inp_id,
                               filter="ensembl_peptide_id",
                               marts=c(yeast,canis),
                               labels=c('S_cerevisiae','CL_familiaris'))
write.csv(ortho_inpara_names,"ss.ortholog/ortho_inpara_names.csv",row.names=F)
ortho_inpara_names = read.csv("ss.ortholog/ortho_inpara_names.csv")
print(dim(ortho_inpara_names))
ortho_inpara_names[sample(nrow(ortho_inpara_names),10),]
data.frame(summary(factor(ortho_inpara_names$Species)))
data.frame(table(ortho_inpara$species))

ortho_inpara_yeast = which(ortho_inpara$species==speci1)
ortho_inpara_canis = which(ortho_inpara$species==speci2)
ortho_inpara_inpid_yeast = na.omit(unique(ortho_inpara$inp_id[ortho_inpara_yeast]))
ortho_inpara_inpid_canis = na.omit(unique(ortho_inpara$inp_id[ortho_inpara_canis]))
print(length(ortho_inpara_yeast))
print(length(ortho_inpara_canis))
print(length(ortho_inpara_inpid_yeast))
print(length(ortho_inpara_inpid_canis))

# InParanoid orthologs
#ortho_inpara_human_id = which(ortho_inpara_names$Species=="H_sapiens")
ortho_inpara_yeast_id = which(ortho_inpara_names$Species==label1)
ortho_inpara_canis_id = which(ortho_inpara_names$Species==label2) # 170704
#ortho_inpara_egid_human = na.omit(unique(ortho_inpara_names$entrezgene[ortho_inpara_human_id]))
ortho_inpara_egid_yeast = na.omit(unique(ortho_inpara_names$entrezgene[ortho_inpara_yeast_id]))
ortho_inpara_egid_canis = na.omit(unique(ortho_inpara_names$entrezgene[ortho_inpara_canis_id])) # 170704


## Compiling EntrezGene ID anntation to Inparanoid orthologs
ortho_inpara_ann = ss.annot(ortho_inpara, ortho_inpara_names, colRange=c(2:5))
head(ortho_inpara)
head(ortho_inpara_names)
write.csv(ortho_inpara_ann, "ss.ortholog/ortho_inpara_ann.csv", row.names=F) # ignore warning message
ortho_inpara_ann = read.csv("ss.ortholog/ortho_inpara_ann.csv")
head(ortho_inpara_ann)
data.frame(table(ortho_inpara_ann$annNum))
```

# 3. Get orthologs from HomoloGene

```R
#library(annotationTools)
#browseVignettes("annotationTools")

homologene = read.delim("ss.ortholog/homologene_build68.data.tsv", header=FALSE)
colnames(homologene) = c('HID','TaxId','GeneId','GeneSymbol','Proteingi','Proteinacc')
print(paste("rows=",dim(homologene)[1],", columns=",dim(homologene)[2])); head(homologene)
taxids = read.delim("ss.ortholog/HomoloGene_taxid_taxname.txt", header=FALSE)
taxids

### Homologene - TaxId ###
# 9606	Homo sapiens
# 4932	Saccharomyces cerevisiae
# 10090	Mus musculus
# 9615	Canis lupus familiaris
############################

tax1="4932"; label1='S_cerevisiae' # Saccharomyces cerevisiae
tax2="9615"; label2='CL_familiaris' # Canis lupus familiaris
ortho_homo = homo.filter(homologene, tax1=tax1,tax2=tax2) # Yeast & Dog - 9.5 min
write.csv(ortho_homo, "ss.ortholog/ortho_homo.csv",row.names=F)
ortho_homo = read.csv("ss.ortholog/ortho_homo.csv"); print(dim(ortho_homo))
ortho_homo[sample(nrow(ortho_homo),10),]
data.frame(table(ortho_homo$TaxId))

#ortho_homo_human_id = which(ortho_homo$TaxId=="9606")
ortho_homo_yeast_id = which(ortho_homo$TaxId==tax1)
ortho_homo_canis_id = which(ortho_homo$TaxId==tax2)
#ortho_homo_egid_human = na.omit(unique(ortho_homo$GeneId[ortho_homo_human_id]))
ortho_homo_egid_yeast = na.omit(unique(ortho_homo$GeneId[ortho_homo_yeast_id]))
ortho_homo_egid_canis = na.omit(unique(ortho_homo$GeneId[ortho_homo_canis_id]))
print(length(ortho_homo_yeast_id))
print(length(ortho_homo_canis_id))
print(length(ortho_homo_egid_yeast))
print(length(ortho_homo_egid_canis))

## Get ensembl_gene_id from entrez Id of homologene orthologs
attributesH = c('entrezgene','ensembl_gene_id','external_gene_name')
ortho_homo_names = ens.getBM(query_ids=ortho_homo$GeneId,
                             attri=attributesH,
                             filter="entrezgene",
                             marts=c(yeast,canis),
                             labels=c(label1,label2))
ortho_homo_names[sample(nrow(ortho_homo_names),5),]
head(ortho_homo)
ortho_homo_ann = ss.annot(ortho_homo[,c(3,4,1,2)], ortho_homo_names, colRange=c(1:4)); # 11.6s
print(dim(ortho_homo_ann)) # ignore warning message
ortho_homo_ann[sample(nrow(ortho_homo_ann),10),]
data.frame(table(ortho_homo_ann$annNum))
```

# 4. Venn analysis of orthologs

* using 'ss.venn3' function - 170704 renewal ver

```R
ortho_egid_canis = list(ortho_ensem_egid_canis,ortho_inpara_egid_canis,ortho_homo_egid_canis)
names(ortho_egid_canis) = c("Ensembl_dog","Inparanoid_dog","HomoloGene_dog")

#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
ortho_union_egid_canis = ss.venn3(ortho_egid_canis, main="Dog-Yeast orthologs (Dog) from DBs")
print(dim(ortho_union_egid_canis$list))

ortho_egid_yeast = list(ortho_ensem_egid_yeast,ortho_inpara_egid_yeast,ortho_homo_egid_yeast)
names(ortho_egid_yeast) = c("Ensembl_yeast","Inparanoid_yeast","HomoloGene_yeast")

ortho_union_egid_yeast = ss.venn3(ortho_egid_yeast, main="Dog-Yeast orthologs (Yeast) from DBs")
print(dim(ortho_union_egid_yeast$list))
```

# 5. Search target orthologs

```R
ortho.genes = function(targets=NULL, dblist=c("ensembl","inparanoid","homologene"), Tbform=1){ # Tbform=1(cluster table)/2(1:1 match table)
  time1 = Sys.time()
  ensem = data.frame(); inpara = data.frame()
  homol = data.frame(); out = data.frame()
  n=length(targets); i=1
  pb=winProgressBar(title="Progress",label="ortho.genes functoon",min=0,max=n,width=500)
  for(target in targets){
    if(!is.numeric(target)){
      target = toupper(target) # Capital charater
      #target = paste0('^',toupper(target),'$')
    } else if(i==1){ print("target is numeric") }

    ## Search ensembl DB
    if("ensembl" %in% dblist){
      if(is.numeric(target)){
        ensem_etz = data.frame(ortho_ensem$NCBI.gene.ID,ortho_ensem$NCBI.gene.ID.1) # bugfix 170417
        ensem_rows = which(apply(ensem_etz,1,function(x) any(which(x==target))))
      } else {
        ensem_rows = which(apply(ortho_ensem,1,function(x) any(which(x==target))))
      }
      if(Tbform==1 & length(ensem_rows>0)) {
        ensem_tmp = data.frame(ortho_ensem[ensem_rows,],Group=paste0("q",i))
        names(ensem_tmp) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID",
                             "Ensembl.Gene.ID.1","Gene.Symbol.1","EntrezGene.ID.1",
                             "Group")
        ensem = rbind(ensem,ensem_tmp)
      } else if(Tbform==2) {
        ensem_hu = unique(ortho_ensem[ensem_rows,1:3])
        if(length(ensem_hu[,1])!=0){
          ensem_hu = data.frame(ensem_hu,Species=label2,Group=paste0("q",i))
          names(ensem_hu) = c("Ensembl.Gene.ID","Gene.Symbol",
                              "EntrezGene.ID","Species","Group")
        } else { ensem_hu = NULL }
        ensem_ye = unique(ortho_ensem[ensem_rows,4:6])
        if(length(ensem_ye[,1])!=0){
          ensem_ye = data.frame(ensem_ye,Species=label1,Group=paste0("q",i))
          names(ensem_ye) = c("Ensembl.Gene.ID","Gene.Symbol",
                              "EntrezGene.ID","Species","Group")
        } else { ensem_ye = NULL }
        ensem = rbind(ensem,ensem_hu,ensem_ye)
    } }

    ## Search inparanoid DB with annotataion
    if("inparanoid" %in% dblist){
      if(is.numeric(target)){
        inp_id_row = which(ortho_inpara_ann$entrezgene==target)
      } else {
        inp_id_row = which(apply(ortho_inpara_ann,1,function(x) any(which(x==target))))
      }
      inp_cid = unique(ortho_inpara_ann[inp_id_row,2])
      inp_cid_rows = which(ortho_inpara_ann$clust_id==inp_cid)
      if(Tbform==1) {
        #inpara = rbind(inpara,ortho_inpara_ann[inp_cid_rows,])
        inpara_tmp = unique(ortho_inpara_ann[inp_cid_rows,])
        hu_id = which(inpara_tmp$species=="CANFA") # "MUSMU", "HOMSA"
        ye_id = which(inpara_tmp$species=="SACCE")
        for(id in hu_id) {
          if(length(ye_id)!=0) {
            inpara_hu.ye = cbind(inpara_tmp[id,6:8],inpara_tmp[ye_id,6:8],
                                 inpara_tmp[ye_id,2])
            names(inpara_hu.ye) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID",
                                    "Ensembl.Gene.ID.1","Gene.Symbol.1","EntrezGene.ID.1","cid")
          } else inpara_hu.ye = NULL
          inpara = rbind(inpara,inpara_hu.ye)
        }
      } else if(Tbform==2) {
        inpara_tmp = unique(ortho_inpara_ann[inp_cid_rows,c(6:9,2)])
        names(inpara_tmp) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID","Species","cid")
        inpara = rbind(inpara,inpara_tmp)
    } }

    ## Search homologene DB with annotation
    if("homologene" %in% dblist){
      if(is.numeric(target)){
        homo_id_row = which(ortho_homo_ann$GeneId==target)
      } else {
        homo_id_row = which(apply(ortho_homo_ann,1,function(x) any(which(x==target))))
      }
      homo_hid = unique(ortho_homo_ann[homo_id_row,3])
      homo_hid_rows = which(ortho_homo_ann$HID==homo_hid)
      if(Tbform==1) {
        #homol = rbind(homol,ortho_homo_ann[homo_hid_rows,])
        homol_tmp = unique(ortho_homo_ann[homo_hid_rows,])
        hu_id = which(homol_tmp$Species==label2)
        ye_id = which(homol_tmp$Species==label1)
        for(id in hu_id) {
          if(length(ye_id)!=0) {
            homol_hu.ye = cbind(homol_tmp[id,c(6,7,5)],homol_tmp[ye_id,c(6,7,5)],
                                homol_tmp[ye_id,3])
              names(homol_hu.ye) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID",
                                     "Ensembl.Gene.ID.1","Gene.Symbol.1","EntrezGene.ID.1","HID")
          } else homol_hu.ye = NULL
          homol = rbind(homol,homol_hu.ye)
        }
      } else if(Tbform==2) {
        homol_tmp = unique(ortho_homo_ann[homo_hid_rows,c(6,2,1,8,3)])
        names(homol_tmp) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID","Species","HID")
        homol = rbind(homol,homol_tmp)
    } }

    ## Progress time ##
    d=difftime(Sys.time(),time1,unit="sec")
    if(d<60) d=paste(round(d,2),"sec")
    else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
    else if(d>=3600) d=paste(round(d/3600,1),"hr")
    setWinProgressBar(pb,i,title=paste0(round(i/n*100,1)," % (",i,"/",n,") done for ",d))
    ###################
    i = i+1
  }

  if(Tbform==1) { # Compiling union table as 1:1 match
    ## bugfix_160718
    Tb = data.frame(Ensembl.Gene.ID=character(0),
                    Gene.Symbol=character(0),
                    EntrezGene.ID=character(0),
                    Ensembl.Gene.ID.1=character(0),
                    Gene.Symbol.1=character(0),
                    EntrezGene.ID.1=character(0))
    if(length(ensem)==0) ensem = data.frame(Tb, Group=character(0))
    if(length(inpara)==0) inpara = data.frame(Tb, cid=character(0))
    if(length(homol)==0) homol = data.frame(Tb, HID=character(0))
    cola=7 # ortholog group id
  } else if(Tbform==2) { cola=5 } # Compiling union table as group result
  #print(paste("ensem, length=",length(rownames(ensem))))
  #print(ensem)
  #print(paste("inpara, length=",length(rownames(inpara))))
  #print(inpara)
  #print(paste("homol, length=",length(rownames(homol))))
  #print(homol)
  unionlist = rbind(ensem[,-cola],inpara[,-cola],homol[,-cola]) # change 170417
  unionlist = unique(unionlist)
  unionTb = data.frame(unionlist,
                       Ensembl=character(length(unionlist[,1])),
                       Inparanoid=character(length(unionlist[,1])),
                       Homologene=character(length(unionlist[,1])))
  unionTb$Ensembl    = ensem[match(unionTb[,1],ensem[,1]),cola]
  unionTb$Inparanoid = inpara[match(unionTb[,1],inpara[,1]),cola]
  unionTb$Homologene = homol[match(unionTb[,1],homol[,1]),cola]
  out = unionTb

  print(difftime(Sys.time(),time1,units="auto")); close(pb)
  return(out)
}
```
```R
#ortho.genes(target=c("PXA1","PXA2"),Tbform=2)
#ortho.genes(target="ABCD1")
#ortho.genes(target=c(6520,215,5594))
#ortho.genes(target=c("SLC3A2","ABCD1","MAPK1"))
egid_canis = Reduce(union, ortho_egid_canis); print(length(egid_canis))
egid_yeast = Reduce(union, ortho_egid_yeast); print(length(egid_yeast))
ortho_canis = ortho.genes(target=egid_canis)
print(dim(ortho_canis)) # 7,167
write.csv(ortho_canis,"170704_ortholog_canis.csv",row.names=F)
ortho_yeast = ortho.genes(target=egid_yeast)
print(dim(ortho_yeast)) # 8,950
write.csv(ortho_yeast,"170704_ortholog_yeast.csv",row.names=F)

# why did different result show? Yeast results contains all the canis results.
print(length(setdiff(ortho_yeast$Ensembl.Gene.ID,ortho_canis$Ensembl.Gene.ID))) # 1,102
print(length(setdiff(ortho_yeast$Ensembl.Gene.ID.1,ortho_canis$Ensembl.Gene.ID1))) # 2,817
```
