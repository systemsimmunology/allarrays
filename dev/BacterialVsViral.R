

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

## Attempt simple probabilistic classifier
w <- (length(b.conds) - apply(mm[,b.conds],1,sum)) + apply(mm[,v.conds],1,sum)
##w <- w/(length(bmdm.conds)+length(dc.conds))
names(sort(w,decreasing=TRUE)[1:10])

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

we <- names(sort(w,decreasing=TRUE)[1:25])
eid <- we[1]

x11()
main <- paste(gene.symbol[eid],": Viral, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.vex, data.matrix = dm.exon, main=main,labvec=mm[eid,], ymax=max(dm.exon.v[eid,]))

x11()
main <- paste(gene.symbol[eid],": Bacterial, Exon Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.bex, data.matrix = dm.exon, main=main, labvec=mm[eid,], ymax=max(dm.exon.bmdm[eid,]))

x11()
main <- paste(gene.symbol[eid],": Bacterial, 3 Prime Arrays", collapse=" ")
gridPlotCSS(eid, CSSs.tc.b3p, data.matrix = dm.3prime,main=main, labvec=mm[eid,], ymax=max(dm.3prime.bmdm[eid,]))


n.conds <- length(v)

