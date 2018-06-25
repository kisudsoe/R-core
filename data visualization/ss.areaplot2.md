---
title: "ss.areaplot2"
author: "KimSS"
date: "2016-09-26 MON"
output: html_document
---

- 2016-09-26 ver for X-ALD project/splicing analysis

```{r id:"j1emhtea"}
ss.areaPlot2 = function(gen.sgl,id,grp) {
  n = length(gen.sgl)
  gen.sgl.trans = cbind(stack(gen.sgl),sample=rep(colnames(gen[4:14]),n),
                        Group=rep(grp,n)) %>% tbl_df
  #print(head(gen.sgl.trans,10))
  gen_sum = summarise(group_by(gen.sgl.trans,ind,Group),
                      Max=max(values),Min=min(values)) %>% tbl_df
  gen_sum = arrange(gen_sum,Group) # input data preparation is done.

  #library(ggplot2)
  p = ggplot()+theme_bw()+
    geom_ribbon(data=gen_sum,
                aes(x=ind,ymax=Max,ymin=Min,group=Group,fill=Group),
                alpha=0.4)+
    stat_summary(data=gen.sgl.trans,
                 aes(ind,values,group=Group,colour=Group),
                 fun.y=mean,geom="line",size=1)+
    geom_point(data=gen.sgl.trans,
               aes(ind,values,colour=Group,shape=Group),
               size=2)+
    labs(title=paste0("Splicing analysis of AffyID: ",id),
         x="AffyExonID",y="Probeset intensity")+
    theme(axis.text.x=element_text(angle=45, hjust=1))
  return(p)
}

id = 16957951
gen = filter(union_sgnl,AffyGeneID==id) %>% tbl_df
gen.sgl = gen[,4:14] %>% t %>% data.frame %>% tbl_df
colnames(gen.sgl) = as.vector(gen$AffyExonID)
gen.sgl

ss.areaPlot2(gen.sgl,id,group.info$Group)
'''
```