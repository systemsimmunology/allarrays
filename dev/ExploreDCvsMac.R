
##########################################
## BMDM vs DC signatures 
#########################################

pdata.dir <- file.path(Sys.getenv("AA"),"processed_data/20091015")
util.dir <- file.path(Sys.getenv("AA"),"utils")

load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))

source(paste(util.dir,"utilitiesPlot.R",sep="/"))
source(paste(util.dir,"utilitiesMetaData.R",sep="/"))

load(paste(pdata.dir,"dm.3prime.RData",sep="/"))
load(paste(pdata.dir,"dm.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.3prime.RData",sep="/"))
load(paste(pdata.dir,"mm.RData",sep="/"))

all.conds <- colnames(mm)

v1 <- unlist(lapply(CSSs.tc.exon,"[[","Cell Type"))
v2 <- unlist(lapply(CSSs.tc.3prime,"[[","Cell Type"))
v <- c(v1,v2)
# which is the same as this v
v <- unlist(lapply(CSSs.tc,"[[","Cell Type"))
bmdm.conds <- names(which(v=="BMDM"))
dc.conds <- names(which(v=="Dendritic Cell"))
n.bmdm <- length(bmdm.conds)
n.dc <- length(dc.conds)
n.conds <- length(v)

## bmdm set on exon arrays: timecourses and data matrices
bex <- intersect(bmdm.conds,names(CSSs.tc.exon))
CSSs.tc.bex <- CSSs.tc.exon[bex]
dm.exon.bmdm <- dm.exon[,unlist(lapply(CSSs.tc.bex,"[[","DM Column"))]

## dc set on exon arrays: timecourses and data matrices
dcex <- intersect(dc.conds,names(CSSs.tc.exon))
CSSs.tc.dcex <- CSSs.tc.exon[dcex]
dm.exon.dc <- dm.exon[,unlist(lapply(CSSs.tc.dcex,"[[","DM Column"))]

## bmdm set on 3 prime arrays: timecourses and data matrices
b3p <- intersect(bmdm.conds,names(CSSs.tc.3prime))
CSSs.tc.b3p <- CSSs.tc.3prime[b3p]
dm.3prime.bmdm <- dm.3prime[,unlist(lapply(CSSs.tc.b3p,"[[","DM Column"))]

## Attempt simple probabilistic classifier
## scores for "DC-ness"
w <- (length(bmdm.conds) - apply(mm[,bmdm.conds],1,sum)) + apply(mm[,dc.conds],1,sum)
##w <- w/(length(bmdm.conds)+length(dc.conds))
we <- names(sort(w,decreasing=TRUE)[1:25])

## Ranking by "BMDM-ness"
we <- names(sort(w,decreasing=FALSE)[1:25])

vdc <- apply(mm[,dc.conds],1,sum)/n.dc
vmac <- apply(mm[,bmdm.conds],1,sum)/n.bmdm
plot(vmac,vdc,xlab="DC frequency",ylab="Mac frequency")
text(vmac,vdc,labels=gene.symbol[rownames(mm)],pos=3)


eid <- we[2]

sum(mm[eid,dc.conds])/length(dc.conds)
sum(mm[eid,bmdm.conds])/length(bmdm.conds)
paste(mm[eid,dc.conds],collapse="")
paste(mm[eid,bmdm.conds],collapse="")

we <- names(sort(w,decreasing=TRUE)[1:25])

x11()
main <- paste(gene.symbol[eid],": DCs, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.dcex, data.matrix = dm.exon, main=main,labvec=mm[eid,], ymax=max(dm.exon.dc[eid,]))

x11()
main <- paste(gene.symbol[eid],": BMDMs, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.bex, data.matrix = dm.exon, main=main, labvec=mm[eid,], ymax=max(dm.exon.bmdm[eid,]))

main <- paste(gene.symbol[eid],": BMDMs, 3 Prime Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.b3p, data.matrix = dm.3prime,main=main, labvec=mm[eid,], ymax=max(dm.3prime.bmdm[eid,]))

### Which DC stims were also done for BMDM

cc.bex <- character() ## conditions in common BMDM conds
cc.dcex <- character() ## conditions in common DC conds
for ( cond in dc.conds ){
  query <-  CSSs.tc[[cond]][c("Strain","Stimulus 1","Stimulus 2")]
  v <- inListSoft(query, CSSs.tc.exon[bex])
  if ( length(v) > 0 ) {
    cat("***",cond, "***\n")
    cat(names(CSSs.tc.exon[bex[v]]),"\n")
    cc.bex <- c(cc.bex,bex[v])
    cc.dcex <- c(cc.dcex,cond)
  }
}
## Note that the two are different in length as we
## have both Male and Female for BMDM in some cases




eids <- rownames(mm.cdag)
eids <- intersect(cdag,rownames(dm.exon))

##
## Magnitudes
##

dcm <- meanAbsMultipleCSSTs(eids,CSSs.tc.exon[dcex],data.matrix=dm.exon)
bm <- meanAbsMultipleCSSTs(eids,CSSs.tc.exon[bex],data.matrix=dm.exon)

dcmm <- apply(dcm,1,mean)
bmm <- apply(bm,1,mean)

plot(dcmm,bmm,xlab="DC arrays, mean expression",ylab="BMDM arrays, mean expression",main="Expression of CD Antigens")
text(dcmm,bmm,label=gene.symbol[eids],pos=4)
abline(0,1)

plot(log2(dcmm),log2(bmm),xlab="DC arrays, log2(mean expression)",ylab="BMDM arrays, log2(mean expression)",main="Expression of CD Antigens")
text(log2(dcmm),log2(bmm),label=gene.symbol[eids],pos=4)
abline(0,1)

##
## Induction Ratios
##

maxRatioMultipleCSSTs 

dcm <- maxRatioMultipleCSSTs(eids,CSSs.tc.exon[cc.dcex],data.matrix=dm.exon)
bm <- maxRatioMultipleCSSTs(eids,CSSs.tc.exon[cc.bex],data.matrix=dm.exon)

dcmm <- apply(dcm,1,mean) ## Mean induction ratio for DCs, if that makes any sense
bmm <- apply(bm,1,mean) ## Mean induction ratio for BMDMs, if that makes any sense

x11()
plot(dcmm,bmm,xlab="DC arrays, mean induction",ylab="BMDM arrays, mean induction",main="Expression of CD Antigens", xlim=c(0,10),ylim=c(0,10))
text(dcmm,bmm,label=gene.symbol[eids],pos=4)
abline(0,1)

eid <- gene.eid["Tlr3"]
x11()
main <- paste(gene.symbol[eid],": DCs, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.dcex[cc.dcex], data.matrix = dm.exon, main=main,labvec=mm[eid,], ymax=max(dm.exon.dc[eid,]))
x11()
main <- paste(gene.symbol[eid],": BMDMs, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.bex[cc.bex], data.matrix = dm.exon, main=main, labvec=mm[eid,], ymax=max(dm.exon.bmdm[eid,]))



##
## Random trials
##
gene.symbol[eids],pos=4)
abline(0,1)
dc.conds.rand <- all.conds[sample(1:n.conds,n.dc)]
bmdm.conds.rand <- setdiff(all.conds,dc.conds.rand)
w.rand <- (n.bmdm - apply(mm[,bmdm.conds.rand],1,sum)) + apply(mm[,dc.conds.rand],1,sum)
##w.rand <- w.rand/n.conds
hist(w.rand, breaks=50)
sort(w.rand,decreasing=TRUE)[1:10]

x11()
par(mfrow=c(2,1))
hist(w, breaks=1:n.conds, xlim=c(100,n.conds),ylim = c(0,100))
hist(w.rand, breaks=1:n.conds, xlim=c(100,n.conds),ylim = c(0,100))


x11()
par(mfrow=c(2,1))
hist(w, breaks=1:n.conds, xlim=c(90,100))
hist(w.rand, breaks=1:n.conds, xlim=c(90,100))

