

pdata.dir <- file.path(Sys.getenv("AA"),"processed_data/20091015")
util.dir <- file.path(Sys.getenv("AA"),"utils")

load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.eid.RData",sep="/"))

source(paste(util.dir,"utilitiesPlot.R",sep="/"))
source(paste(util.dir,"utilitiesMetaData.R",sep="/"))

load(paste(pdata.dir,"dm.3prime.RData",sep="/"))
load(paste(pdata.dir,"dm.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.3prime.RData",sep="/"))
load(paste(pdata.dir,"mm.RData",sep="/"))

all.conds <- colnames(mm)

rt <- read.table("/Users/thorsson/allarrays/dev/StimulusClasses.tsv",as.is=TRUE,sep='\t',header=TRUE)
gclass <- rt[["General.Class"]]
names(gclass) <- rt[["Stimulus"]]

v1 <- unlist(lapply(CSSs.tc.exon,"[[","Stimulus 1"))
v2 <- unlist(lapply(CSSs.tc.3prime,"[[","Stimulus 1"))
v <- c(v1,v2)
# which (should be?) the same as this v
v <- unlist(lapply(CSSs.tc,"[[","Stimulus 1"))

v.conds <- names(v[which(gclass[v]=="Viral")])
b.conds <- names(v[which(gclass[v]=="Bacterial")])

nb <- length(b.conds)
nv <- length(v.conds)
nall <- nb+nv

## b set on exon arrays: timecourses and data matrices
bex <- intersect(b.conds,names(CSSs.tc.exon))
CSSs.tc.bex <- CSSs.tc.exon[bex]
dm.exon.b <- dm.exon[,unlist(lapply(CSSs.tc.bex,"[[","DM Column"))]

## b set on 3 prime arrays: timecourses and data matrices
b3p <- intersect(b.conds,names(CSSs.tc.3prime))
CSSs.tc.b3p <- CSSs.tc.3prime[b3p]
dm.3prime.b <- dm.3prime[,unlist(lapply(CSSs.tc.b3p,"[[","DM Column"))]

## v set on exon arrays: timecourses and data matrices
vex <- intersect(v.conds,names(CSSs.tc.exon))
CSSs.tc.vex <- CSSs.tc.exon[vex]
dm.exon.v <- dm.exon[,unlist(lapply(CSSs.tc.vex,"[[","DM Column"))]

## v set on 3 prime arrays: timecourses and data matrices
v3p <- intersect(v.conds,names(CSSs.tc.3prime))
CSSs.tc.v3p <- CSSs.tc.3prime[v3p]
dm.3prime.v <- dm.3prime[,unlist(lapply(CSSs.tc.v3p,"[[","DM Column"))]


## Attempt simple probabilistic classifier
w <- (length(b.conds) - apply(mm[,b.conds],1,sum)) + apply(mm[,v.conds],1,sum)
w <- w/(length(b.conds)+length(v.conds))

vv <- apply(mm[,v.conds],1,sum)/nv
bb <- apply(mm[,b.conds],1,sum)/nb
plot(bb,vv,xlab="Bacterial frequency",ylab="Viral frequency")
text(bb,vv,labels=gene.symbol[rownames(mm)],pos=3)


## Condition weight??
cw <- apply(mm,2,sum)/nrow(mm)
## Normalize to a "probability"
cw <- cw/sum(cw)
mmw <- t(t(mm)*cw)
vvw <- apply(mmw[,v.conds],1,sum)/nv
bbw <- apply(mmw[,b.conds],1,sum)/nb
plot(bbw,vvw,xlab="Weighted Bacterial frequency",ylab="Weighted Viral frequency")
text(bbw,vvw,labels=gene.symbol[rownames(mm)],pos=3)


we <- names(sort(w,decreasing=TRUE)[1:25])


we <- names(sort(w,increasing=TRUE)[1:25])


eid <- we[2]
sum(mm[eid,v.conds])/length(v.conds)
sum(mm[eid,b.conds])/length(b.conds)
paste(mm[eid,v.conds],collapse="")
paste(mm[eid,b.conds],collapse="")


x11()
main <- paste(gene.symbol[eid],": Viral activation, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.vex, data.matrix = dm.exon, main=main,labvec=mm[eid,], ymax=max(dm.exon.v[eid,]))

x11()
main <- paste(gene.symbol[eid],": Bacterial activation, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.bex, data.matrix = dm.exon, main=main, labvec=mm[eid,], ymax=max(dm.exon.b[eid,]))

x11()
main <- paste(gene.symbol[eid],": Bacterial activation, 3 Prime Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.b3p, data.matrix = dm.3prime,main=main, labvec=mm[eid,], ymax=max(dm.3prime.b[eid,]))

x11()
main <- paste(gene.symbol[eid],": Viral activation, 3 Prime Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.v3p, data.matrix = dm.3prime,main=main, labvec=mm[eid,], ymax=max(dm.3prime.v[eid,]))

n.conds <- length(v)

