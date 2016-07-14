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

ortho_ensem = getLDS(attributes,filters="with_homolog_scer",values=TRUE,mart=human,
                     attributesL=attributesL,martL=yeast) # values=TRUE : all exist values

# > head(ortho_ensem)
#   Ensembl.Gene.ID HGNC.symbol EntrezGene.ID Ensembl.Gene.ID.1 Associated.Gene.Name EntrezGene.ID.1
# 1 ENSG00000155366        RHOC           389           YPR165W                 RHO1          856294
# 2 ENSG00000119688       ABCD4          5826           YPL147W                 PXA1          855956
# 3 ENSG00000090238       YPEL3         83719           YBL049W                 MOH1          852231
# 4 ENSG00000184985      SORCS2         57537           YCR101C                               850465
# 5 ENSG00000152556        PFKM          5213           YGR240C                 PFK1          853155
# 6 ENSG00000099365       STX1B        112755           YMR183C                 SSO2          855221

## To see list of ensembl data ##
ensembl=useMart("ensembl")
View(listDatasets(ensembl)) # To find database list

View(listAttributes(human)) # To find attribute list
View(listFilters(human))    # To find filter list


# 2. Orthologs from Inparanoid human
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
ens.getBM = function(inp_id,species){
  attributes = c('ensembl_gene_id','external_gene_name','entrezgene')
  if(species=="HOMSA"){
    BM = getBM(attributes=attributes, filters="ensembl_peptide_id",
                   values=inp_id, mart=human)
    nrow = round(length(unlist(BM))/3,0)
    out = data.frame(matrix(unlist(BM),nrow=nrow,byrow=F))
  } else if(species=="SACCE"){
    BM = getBM(attributes=attributes, filters="ensembl_gene_id",
                   values=inp_id, mart=yeast)
    nrow = round(length(unlist(BM))/3,0)
    out = data.frame(matrix(unlist(BM),nrow=nrow,byrow=F))
  }
  if(length(out)==0){
    nrow = 0
    out = data.frame(X1=inp_id,X2='NA',X3='NA',X4='NA',X5=nrow)
  } else { out = cbind(rep(inp_id,nrow),out,rep(nrow)) }
  names(out) = c('inp_id',names(BM),'nrow')
  return(out)
}

ens.multi = function(ortho_inpara){
  n = length(ortho_inpara[,1])
  nms = data.frame()
  for(i in 1:n){
    inp_id = ortho_inpara[i,1]
    species = ortho_inpara[i,3]
    nms.tmp = ens.getBM(inp_id,species)
    nms = rbind(nms,nms.tmp)
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
  return(nms)
}

ortho_inpara_ensembl = ens.multi(ortho_inpara)
summary(factor(ortho_inpara_names$nrow))