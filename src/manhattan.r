snp = readRDS('Manhattan plot_input.rds')
snp$CHR = as.numeric(gsub('X',23,snp$CHR))

library(qqman)
name_pdf = '66,419_MitoSetPS3[1471].pdf'
name_png = '66,419_MitoSetPS3[1471].png'
#pdf(name_pdf,width=14,height=6,paper='special')
png(name_png,width=23,height=8,units="in",res=100)
manhattan(snp,chr="CHR",bp="POS_START",snp="SNP_rsID",p="P",
          chrlabs = c(1:22,'X'),
          suggestiveline=F,
          genomewideline=-log10(1e-04),
          annotatePval=0.008,
          annotateTop=F)
dev.off