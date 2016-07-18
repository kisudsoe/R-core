# Get R packages for useMart
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("hom.Sc.inp.db")
biocLite("hom.Hs.inp.db")
biocLite("annotationTools")
biocLite("org.Hs.eg.db")
biocLite("org.Sc.sgd.db")

# 1. Orthologs from Ensembl Compara
library(biomaRt)

human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
yeast = useMart("ensembl", dataset="scerevisiae_gene_ensembl")

attributes = c('ensembl_gene_id','hgnc_symbol','entrezgene')
attributesL = c('ensembl_gene_id','external_gene_name','entrezgene')
ortho_ensem = getLDS(attributes,filters="with_homolog_scer",
                     values=TRUE,mart=human, # values=TRUE : all exist values
                     attributesL=attributesL,martL=yeast)

ortho_ensem_egid_human = na.omit(unique(ortho_ensem$EntrezGene.ID))
ortho_ensem_egid_yeast = na.omit(unique(ortho_ensem$EntrezGene.ID.1))

# > head(ortho_ensem)
#   Ensembl.Gene.ID HGNC.symbol EntrezGene.ID Ensembl.Gene.ID.1 Associated.Gene.Name EntrezGene.ID.1
# 1 ENSG00000155366        RHOC           389           YPR165W                 RHO1          856294
# 2 ENSG00000119688       ABCD4          5826           YPL147W                 PXA1          855956
# 3 ENSG00000090238       YPEL3         83719           YBL049W                 MOH1          852231
# 4 ENSG00000184985      SORCS2         57537           YCR101C                               850465
# 5 ENSG00000152556        PFKM          5213           YGR240C                 PFK1          853155
# 6 ENSG00000099365       STX1B        112755           YMR183C                 SSO2          855221

### To see list of ensembl data ###
ensembl=useMart("ensembl")
View(listDatasets(ensembl)) # To find database list

View(listAttributes(human)) # To find attribute list
View(listFilters(human))    # To find filter list


# 2. Orthologs from Inparanoid human
## https://bioconductor.org/packages/release/data/annotation/html/hom.Hs.inp.db.html
library(hom.Hs.inp.db)
ls("package:hom.Hs.inp.db")
as.list(hom.Hs.inpSACCE)[1:4] # What data in this?
mget("ENSP00000348965",hom.Hs.inpSACCE) # Example for catch a data

library(DBI) # To utilize sql query
dbListTables(hom.Hs.inp_dbconn()) # What species for orthologs?
dbListFields(hom.Hs.inp_dbconn(),"Saccharomyces_cerevisiae") # What columns in sql table

tmpsql = "SELECT * FROM Saccharomyces_cerevisiae WHERE inp_id='ENSP00000348965';"
tmpsql = "SELECT * FROM Saccharomyces_cerevisiae WHERE inp_id='YDR288W';"
tmpsql = "SELECT * FROM Saccharomyces_cerevisiae WHERE clust_id='2143';"

ortho_inpara = dbGetQuery(hom.Hs.inp_dbconn(),
                          "SELECT * FROM Saccharomyces_cerevisiae;")

## get orthologs from Inparanoid yeast
# library(hom.Sc.inp.db)
# ls("package:hom.Sc.inp.db")
# dbListTables(hom.Sc.inp_dbconn())
# dbListFields(hom.Sc.inp_dbconn(),"Homo_sapiens")
# ortho_inpara2 = dbGetQuery(hom.Sc.inp_dbconn(),
#                            "SELECT * FROM Homo_sapiens;") # Same as data from human

## Ensembl search of 'Ortho_inpara' to get EntrezGene ID
ens.getBM = function(query_ids,filter,marts,labels){
  attributes = c('ensembl_peptide_id','ensembl_gene_id',
                 'external_gene_name','entrezgene')
  i = 1
  n = length(marts)
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

ortho_inpara_names = ens.getBM(query_ids=ortho_inpara$inp_id,
                               filter="ensembl_peptide_id",
                               marts=c(human,yeast),
                               labels=c('H_sapiens','S_cerevisiae'))
head(ortho_inpara_names)
dim(ortho_inpara_names)
summary(factor(ortho_inpara_names$Species))

ortho_inpara_human_id = which(ortho_inpara_names$Species=="H_sapiens")
ortho_inpara_yeast_id = which(ortho_inpara_names$Species=="S_cerevisiae")
ortho_inpara_egid_human = na.omit(unique(ortho_inpara_names$entrezgene[ortho_inpara_human_id]))
ortho_inpara_egid_yeast = na.omit(unique(ortho_inpara_names$entrezgene[ortho_inpara_yeast_id]))

## Compiling EntrezGene ID annotation to Inparanoid orthologs
ss.annot = function(origin,annot){ # First columns of both entries should be IDs
  out = data.frame()
  i = 1
  n = length(origin[,1])
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
      annot_row = data.frame(origin[i,],annot[ann_row,2:5],annNum=annNum)
    }
    out = rbind(out,annot_row)
    #######################
    # Create progress bar #
    #######################
    if(i==1) { # set progress bar
      cat('\nProcess iteration =',n,'\n')
      pb = txtProgressBar(min=0,max=100,width=30,initial=0,style=3)
      time1 = Sys.time()
    } else { setTxtProgressBar(pb, (i/n)*100) } # show progress bar
    if(i==n) { # Duration time check
      close(pb)
      time2 = Sys.time()
      cat('\n')
      print(time2-time1)
    } 
    #######################
    i = i+1
  }
  return(out)
}
ortho_inpara_ann = ss.annot(ortho_inpara, ortho_inpara_names) # ignore warning message
summary(factor(ortho_inpara_ann$annNum))


## 3. get orthologs from HomoloGene
library(annotationTools)
browseVignettes("annotationTools")

homologene = read.delim("D:/KimSS-NAS/LFG/Works/2016.04 Human X-ALD/(Archive) Orthologs/Hs-Sc/homologene_build68.data.tsv", header=FALSE)
colnames(homologene) = c('HID','TaxId','GeneId','GeneSymbol','Proteingi','Proteinacc')

homo.filter = function(homo_db) {
  human_rows = which(homo_db$TaxId=="9606")
  yeast_rows = which(homo_db$TaxId=="4932")
  hu.ye_rows_union = union(human_rows,yeast_rows)
  
  human_hid = homo_db$HID[human_rows]
  yeast_hid = homo_db$HID[yeast_rows]
  hu.ye_hid = intersect(human_hid,yeast_hid)
  
  n = length(hu.ye_hid)
  out = data.frame()
  for(i in 1:n){
    hu.ye_hid_rows = which(homo_db$HID==hu.ye_hid[i])
    hu.ye_rows = intersect(hu.ye_rows_union,hu.ye_hid_rows)
    out = rbind(out,homo_db[hu.ye_rows,])
    #######################
    # Create progress bar #
    #######################
    if(i==1) { # set progress bar
      cat('\nProcess iteration =',n,'\n')
      pb = txtProgressBar(min=0,max=100,width=30,initial=0,style=3)
      time1 = Sys.time()
    } else { setTxtProgressBar(pb, (i/n)*100) } # show progress bar
    if(i==n) { # Duration time check
      close(pb)
      time2 = Sys.time()
      cat('\n')
      print(time2-time1)
    } 
    #######################
  }

  return(out)
}

ortho_homo = homo.filter(homologene)
summary(factor(ortho_homo$TaxId))

ortho_homo_human_id = which(ortho_homo$TaxId=="9606")
ortho_homo_yeast_id = which(ortho_homo$TaxId=="4932")
ortho_homo_egid_human = na.omit(unique(ortho_homo$GeneId[ortho_homo_human_id]))
ortho_homo_egid_yeast = na.omit(unique(ortho_homo$GeneId[ortho_homo_yeast_id]))

### Homologene - TaxId ###
# 9606	Homo sapiens
# 4932	Saccharomyces cerevisiae
# ---   ---
# 10090	Mus musculus
# 10116	Rattus norvegicus
# 28985	Kluyveromyces lactis
# 318829	Magnaporthe oryzae
# 33169	Eremothecium gossypii
# 3702	Arabidopsis thaliana
# 4530	Oryza sativa
# 4896	Schizosaccharomyces pombe
# 5141	Neurospora crassa
# 6239	Caenorhabditis elegans
# 7165	Anopheles gambiae
# 7227	Drosophila melanogaster
# 7955	Danio rerio
# 8364	Xenopus (Silurana) tropicalis
# 9031	Gallus gallus
# 9544	Macaca mulatta
# 9598	Pan troglodytes
# 9615	Canis lupus familiaris
# 9913	Bos taurus

## Get ensembl_gene_id from entrez Id of homologene orthologs
attributesH = c('entrezgene','ensembl_gene_id','external_gene_name')
ortho_homo_names = ens.getBM(query_ids=ortho_homo$GeneId,
                             attri=attributesH,
                             filter="entrezgene",
                             marts=c(human,yeast),
                             labels=c('H_sapiens','S_cerevisiae'))
ortho_homo_ann = ss.annot(ortho_homo, ortho_homo_names) # ignore warning message
summary(factor(ortho_homo_ann$annNum))


## 4. Venn analysis of orthologs
### using 'ss.venn3' function - 160704 ver
ortho_ensem_egid_human = setNames(ortho_ensem_egid_human,"Ensembl_human")
ortho_inpara_egid_human = setNames(ortho_inpara_egid_human,"Inparanoid_human")
ortho_homo_egid_human = setNames(ortho_homo_egid_human,"HomoloGene_human")

source("https://bioconductor.org/biocLite.R")
biocLite("limma")
ortho_union_egid_human = ss.venn3(ortho_ensem_egid_human,ortho_inpara_egid_human,
                    ortho_homo_egid_human, main="Human-Yeast orthologs (Human) from DBs")
length(ortho_union_egid_human$list[,1])

ortho_ensem_egid_yeast = setNames(ortho_ensem_egid_yeast,"Ensembl_yeast")
ortho_inpara_egid_yeast = setNames(ortho_inpara_egid_yeast,"Inparanoid_yeast")
ortho_homo_egid_yeast = setNames(ortho_homo_egid_yeast,"HomoloGene_yeast")

ortho_union_egid_yeast = ss.venn3(ortho_ensem_egid_yeast,ortho_inpara_egid_yeast,
                            ortho_homo_egid_yeast, main="Human-Yeast orthologs (Yeast) from DBs")
length(ortho_union_egid_yeast$list[,1])


## 5. Search target orthologs
ortho.genes = function(targets,dblist=c("ensembl","inparanoid","homologene"),
                       Tbform=1){ # Tbform=1(cluster table)/2(1:1 match table)
  ensem = data.frame(); inpara = data.frame()
  homol = data.frame(); out = data.frame()
  n=length(targets); i=1
  for(target in targets){
    if(!is.numeric(target)){
      target = toupper(target) # Capital charater
      #target = paste0('^',toupper(target),'$')
    } else if(i==1){ print("target is numeric") }
    
    ## Search ensembl DB
    if("ensembl" %in% dblist){
      if(is.numeric(target)){
        ensem_etz = data.frame(ortho_ensem$EntrezGene.ID,ortho_ensem$EntrezGene.ID.1)
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
          ensem_hu = data.frame(ensem_hu,Species="H_sapiens",Group=paste0("q",i))
          names(ensem_hu) = c("Ensembl.Gene.ID","Gene.Symbol",
                              "EntrezGene.ID","Species","Group")
        } else { ensem_hu = NULL }
        ensem_ye = unique(ortho_ensem[ensem_rows,4:6])
        if(length(ensem_ye[,1])!=0){
          ensem_ye = data.frame(ensem_ye,Species="S_cerevisiae",Group=paste0("q",i))
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
        hu_id = which(inpara_tmp$species=="HOMSA")
        ye_id = which(inpara_tmp$species=="SACCE")
        for(id in hu_id) {
          inpara_hu.ye = cbind(inpara_tmp[id,6:8],inpara_tmp[ye_id,6:8],
                               inpara_tmp[ye_id,2])
          names(inpara_hu.ye) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID",
                                  "Ensembl.Gene.ID.1","Gene.Symbol.1","EntrezGene.ID.1","cid")
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
        hu_id = which(homol_tmp$Species=="H_sapiens")
        ye_id = which(homol_tmp$Species=="S_cerevisiae")
        for(id in hu_id) {
          homol_hu.ye = cbind(homol_tmp[id,c(6,7,5)],homol_tmp[ye_id,c(6,7,5)],
                              homol_tmp[ye_id,3])
          names(homol_hu.ye) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID",
                                 "Ensembl.Gene.ID.1","Gene.Symbol.1","EntrezGene.ID.1","HID")
          homol = rbind(homol,homol_hu.ye)
        }
      } else if(Tbform==2) {
        homol_tmp = unique(ortho_homo_ann[homo_hid_rows,c(6,2,1,8,3)])
        names(homol_tmp) = c("Ensembl.Gene.ID","Gene.Symbol","EntrezGene.ID","Species","HID")
        homol = rbind(homol,homol_tmp)
      } }
    
    #######################
    # Create progress bar #
    #######################
    if(i==1) { # set progress bar
      cat('\nProcess iteration =',n,'\n')
      pb = txtProgressBar(min=0,max=100,width=30,initial=0,style=3)
      time1 = Sys.time()
    } else { setTxtProgressBar(pb, (i/n)*100) } # show progress bar
    if(i==n) { # Duration time check
      close(pb)
      time2 = Sys.time()
      cat('\n')
      print(time2-time1)
    } 
    #######################
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
    if(length(ensem)==0) {
      ensem = data.frame(Tb, Group=character(0)) }
    if(length(inpara)==0) {
      inpara = data.frame(Tb, cid=character(0)) }
    if(length(homol)==0) {
      homol = data.frame(Tb, HID=character(0)) }
    cola=6; colb=7
  } else if(Tbform==2) { cola=4; colb=5 } # Compiling union table as group result
  
  unionlist = rbind(ensem[,1:cola],inpara[,1:cola],homol[,1:cola])
  unionlist = unique(unionlist)
  unionTb = data.frame(unionlist,
                       Ensembl=character(length(unionlist[,1])),
                       Inparanoid=character(length(unionlist[,1])),
                       Homologene=character(length(unionlist[,1])))
  unionTb$Ensembl    = ensem[match(unionTb[,1],ensem[,1]),colb]
  unionTb$Inparanoid = inpara[match(unionTb[,1],inpara[,1]),colb]
  unionTb$Homologene = homol[match(unionTb[,1],homol[,1]),colb]
  out = unionTb
  
  return(out)
}

ortho.genes(target=c("PXA1","PXA2"))
ortho.genes(target="ABCD1")
ortho.genes(target=c(6520,215,5594))
ortho.genes(target=c("SLC3A2","ABCD1","MAPK1"))