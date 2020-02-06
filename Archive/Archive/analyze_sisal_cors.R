### Correlations

library(nest)
library(corrplot)

outdir<-"~/R/spielwiese/SISAL_analyses/output/"
TAB<-read.csv("~/R/spielwiese/SISAL_analyses/data/past_millenium_SISALv1b.csv")

ids<-unique(TAB$entity_id)
getts<-function(TAB,i){ts1<-(TAB[TAB$entity_id==ids[i],c(3,4)]);ret<-zoo(ts1[,2],order.by=ts1[,1])}

TS<-list()
for (i in 1:length(ids)){
  TS[[i]]<-getts(TAB,i)
}
names(TS)<-ids


TSpm<-lapply(TS,window,end=1100,start=-40)

sapply(TSpm,length)


# I don't have the SISAL metadata here, so I'm making lats/lons up.
metafake<-data.frame(rep(0,length(TSpm)),rep(0,length(TSpm)),names(TSpm))
colnames(metafake)<-c("Lat","Lon","Name")

# now call the function that checks the sampling.
OUT<-quality_check(TSpm,metafake,T0=-40,T1=1100,maxhiat = 200,min.res=100)
#
# one could set different start/end periods, then multiple time series splits would be returned, here we just get one:
plot(OUT$WBT.out$`188`$tsplit[[1]])
OUT$ind.ok # these fulfill the abovementioned criteria

TSfilt<-lapply(OUT$WBT.out,function(x){x$tsplit[[1]]})

# remove those with less than 50 samples
TSfilt<-TSfilt[-which(sapply(TSfilt,length)<50)]


C<-matrix(NA,nrow=length(TSfilt),ncol=length(TSfilt))
colnames(C)<-rownames(C)<-names(TSfilt)
C.detr<-C.detr.ci<-N.detr<-N<-P<-P.detr<-C


################################################################################
################################################################################
## Correlations for the unvariate time series, both raw and detrended with a 500-year trend removed
ts.detr<-list()
for (i in 1:(length(TSfilt)-1)){
  # effective filtering to 500a
  x.detr<-mean(TSfilt[[i]])+TSfilt[[i]]-gaussbandpass(TSfilt[[i]],1,500)$trend

  ts.detr[[i]]<-x.detr
  for (j in (i+1):length(TSfilt)){
    y.detr<-mean(TSfilt[[j]])+TSfilt[[j]]-gaussbandpass(TSfilt[[j]],1,500)$trend

    temp<-nexcf_ci(x.detr,y.detr,conflevel=0.1)
    C.detr[i,j]<-temp$rxy
    P.detr[i,j]<-temp$pval
    if (temp$rxy<=0) {C.detr.ci[i,j]<-temp$ci[1]} else {C.detr.ci[i,j]<-temp$ci[2]}
    N.detr[i,j]<-temp$neff
    rm(temp)

    temp<-nexcf_ci(TSfilt[[i]],TSfilt[[j]],conflevel=0.1)

    C[i,j]<-temp$rxy
    P[i,j]<-P[j,i]<-temp$pval
    N[i,j]<-temp$neff
    C[j,i]=C[i,j]
    C.detr[j,i]=C.detr[i,j]
    rm(temp)
  }
  }



## Plot the time series just to check
plot(TSfilt[[1]])
lines(ts.detr[[1]],col="red")


## find the highest correlation pairs
srtC<-sort(C)
ij.min<-which((sort(C))[1]==C,arr.ind = TRUE)
(unique(names(TSfilt)[ij.min]))
ij.max<-which(C==max(C,na.rm = TRUE),arr.ind = TRUE)
(unique(names(TSfilt)[ij.max]))

# Die Entities 305, 390 und 286 scheinen sehr interessant zu sein!





## Visualize correlation matrix
p.mat<-P
M<-C

corrplot(C, type="upper", order="hclust",
         p.mat = P, sig.level = 0.1)
mtext("Correlations, sig. level 0.1",line=2,side=1)

corrplot(C.detr,type="upper",order="hclust",p.mat=P.detr,sig.level = 0.1)
mtext("Detrended correlations",line=2,side=1)


C.comb<-C
C.comb[lower.tri(C.comb)]<-C.detr[lower.tri((C.detr))]


P.comb<-P
P.comb[lower.tri(P.comb)]<-t(P.detr[upper.tri((P.detr))])
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")))

col3<-colorRampPalette(c("blue","white","red"))

pdf(file = paste(outdir,"combined_matrix_colr.pdf",sep=""))
corrplot(C.comb,type="full",p.mat=P.comb,sig.level = 0.1,diag=FALSE,is.corr=TRUE,col=col2(200))
dev.off()
pdf(file = paste(outdir,"combined_matrix_number.pdf",sep=""))
corrplot(C.comb,type="full",p.mat=P.comb,sig.level = 0.1,diag=FALSE,is.corr=TRUE,col=col2(200),method="number")
dev.off()
