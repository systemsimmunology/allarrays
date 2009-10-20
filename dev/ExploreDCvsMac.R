

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

## Attempt simple probabilistic classifier
w <- (length(bmdm.conds) - apply(mm[,bmdm.conds],1,sum)) + apply(mm[,dc.conds],1,sum)
##w <- w/(length(bmdm.conds)+length(dc.conds))
names(sort(w,decreasing=TRUE)[1:10])

sum(mm[eid,dc.conds])/length(dc.conds)
sum(mm[eid,bmdm.conds])/length(bmdm.conds)

## bmdm set on exon arrays: timecourses and data matrices
bex <- intersect(bmdm.conds,names(CSSs.tc.exon))
CSSs.tc.bex <- CSSs.tc.exon[bex]
dm.exon.bmdm <- dm.exon[,unlist(lapply(CSSs.tc.bex,"[[","DM Column"))]

## bmdm set on 3 prime arrays: timecourses and data matrices
b3p <- intersect(bmdm.conds,names(CSSs.tc.3prime))
CSSs.tc.b3p <- CSSs.tc.3prime[b3p]
dm.3prime.bmdm <- dm.3prime[,unlist(lapply(CSSs.tc.b3p,"[[","DM Column"))]

## bmdm set on exon arrays: timecourses and data matrices
dcex <- intersect(dc.conds,names(CSSs.tc.exon))
CSSs.tc.dcex <- CSSs.tc.exon[dcex]
dm.exon.dc <- dm.exon[,unlist(lapply(CSSs.tc.dcex,"[[","DM Column"))]

we <- names(sort(w,decreasing=TRUE)[1:25])
eid <- we[2]


x11()
main <- paste(gene.symbol[eid],": DCs, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.dcex, data.matrix = dm.exon, main=main,labvec=mm[eid,], ymax=max(dm.exon.dc[eid,]))

x11()
main <- paste(gene.symbol[eid],": BMDMs, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.bex, data.matrix = dm.exon, main=main, labvec=mm[eid,], ymax=max(dm.exon.bmdm[eid,]))

main <- paste(gene.symbol[eid],": BMDMs, 3 Prime Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.b3p, data.matrix = dm.3prime,main=main, labvec=mm[eid,], ymax=max(dm.3prime.bmdm[eid,]))

### Which DC stims were also done for BMDM
for ( cond in dc.conds ){
  query <-  CSSs.tc[[cond]][c("Strain","Stimulus 1","Stimulus 2")]
  v <- inListSoft(query, CSSs.tc.exon[bex])
  if ( length(v) > 0 ) {
    cat("***",cond, "***\n")
    cat(names(CSSs.tc.exon[bex[v]]),"\n")
  }
}

for ( cond in dc.conds ){
  query <-  CSSs.tc[[cond]][c("Strain","Stimulus 1","Stimulus 2")]
  v <- inListSoft(query, CSSs.tc.3prime[b3p])
  if ( length(v) > 0 ) {
    cat("***",cond, "***\n")
    cat(names(CSSs.tc.3prime[b3p[v]]),"\n")
  }
}

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


