---
title: "ss.areaplot2"
author: "KimSS"
date: "2016-09-26 MON"
output: html_document
---



## 160926 X-ALD project/splicing analysis

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
<<<<<<< HEAD
```


```
## 160324 Mouse muscle aging project

- Line area plot for fibrosis network (2)
- This code edited from line plot "liness (2016-03-21)"
- Minor edit 01 (2016-03-31 THU) -
- Minor edit 02 (2016-04-01 FRI) - Add 'data_total'

​```R
data = input(1) # Input signal data
data_total = input(1) # input total data
group = input(0) # input group info : 14	14	14	29	29	29	44	44	44...

ss.areaPlot = function(data,srch='',type=1,group) {
  # Search string in srch for row id
  row = which(apply(data[,1:6],1, function(x) any(grepl(srch,x))))[1]

  # Get data from the row
  datrow = data[row,]
  print(datrow)
  title=paste(datrow[1],' ',unlist(datrow[2]),' ',
              unlist(datrow[6]),' (',unlist(datrow[3]),')')
  title2=paste(unlist(datrow[6]),'\n',unlist(datrow[2]),
               ' (',unlist(datrow[3]),')') # FC
  title3=paste(unlist(datrow[6]),'\n',unlist(datrow[2]),
               ' (',unlist(datrow[4]),')') # FD
  title4=paste(unlist(datrow[6]),'\n',unlist(datrow[2]),
               ' (',unlist(datrow[5]),')') # Corr
  datr1 = length(datrow)
  datr2 = datr1-29
  datlin = as.numeric(datrow[datr2:datr1])
  ymax = max(datlin)
  ymin = min(datlin)
  if(ymax-ymin<1.5) {
    yrng = c(ymin-0.2,ymax+0.2)
  } else
    yrng = c(ymin-0.5,ymax+0.5)

  group.lv = levels(as.factor(unlist(group)))
  n = length(group.lv)
  print(group.lv)

  for(l in 1:n) {
    id = which(group==group.lv[l])
    if(l==1) {
      datlin.max = max(datlin[id])
      datlin.min = min(datlin[id])
      datlin.mid = median(datlin[id])
    } else {
      datlin.max = c(datlin.max,max(datlin[id]))
      datlin.min = c(datlin.min,min(datlin[id]))
      datlin.mid = c(datlin.mid,median(datlin[id]))
    }
  }
  if(type==1) {
    par(mar=c(0.5,2,0.5,0.5))
    plot(datlin.max,
         main='',
         col="gray50",
         type="l",
         xaxt="n",
         ylim=yrng)
  } else if(type==2) {
    par(mar=c(0.5,2.5,8,0.5))
    plot(datlin.max,
        main=title3,
        cex.main=3.5,
        col="gray50",
        type="l",
        xaxt="n",
        ylim=yrng)
  } else if(type==3) {
    par(mar=c(0.5,2.5,8,0.5))
    plot(datlin.max,
         main=title4,
         cex.main=3.5,
         col="gray50",
         type="l",
         xaxt="n",
         ylim=yrng)
  }

  rect(1.5,yrng[1],2.5,yrng[2],col="gray90",lty=0)
  rect(3.5,yrng[1],4.5,yrng[2],col="gray90",lty=0)
  rect(5.5,yrng[1],6.5,yrng[2],col="gray90",lty=0)
  rect(7.5,yrng[1],8.5,yrng[2],col="gray90",lty=0)
  rect(9.5,yrng[1],10.3,yrng[2],col="gray90",lty=0)
  text(1,yrng[1]+0.1,"14")
  text(2,yrng[1]+0.1,"29")
  text(3,yrng[1]+0.1,"44")
  text(4,yrng[1]+0.1,"53")
  text(5,yrng[1]+0.1,"73")
  text(6,yrng[1]+0.1,"83")
  text(7,yrng[1]+0.1,"95")
  text(8,yrng[1]+0.1,"109")
  text(9,yrng[1]+0.1,"120")
  text(10,yrng[1]+0.1,"135")
  if(type==1) {
    text(5.5,yrng[2]-0.1,title,cex=1.2,font=2)
  }

  lines(datlin.max,col="gray50",lwd=1)
  lines(datlin.min,col="gray50",lwd=1)

  datlin.minr = rev(datlin.min)
  graytrans = rgb(50,50,50,100,maxColorValue=255)
  polygon(c(1:10,10:1),c(datlin.max,datlin.minr),col=graytrans,border=NA)

  lines(datlin.mid,col="black",lwd=3)
}

ss.areaPlot(data,'Cry2',type=2,group) # type=1/2=fd/3=corr
ss.areaPlot(data_total,'bmal1',type=2,group_total)
```

>>>>>>> 6df2c2f77f43481f13ee7ada7852930408bf2f44
