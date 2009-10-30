
util.dir <- file.path(Sys.getenv("AA"),"utils")
pdata.dir <- file.path(Sys.getenv("AA"),"processed_data/20091015") 
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))
source(paste(util.dir,"utilitiesSigTest.R",sep="/"))
load(paste(pdata.dir,"dm.3prime.RData",sep="/"))
load(paste(pdata.dir,"dm.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.3prime.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.RData",sep="/"))

eids.on.exon.array <- rownames(dm.exon)
eids.on.3prime.array <- rownames(dm.3prime)
eids.on.both <- intersect(eids.on.exon.array,eids.on.3prime.array)
eids <- eids.on.both

### assumption is that these will have have no colnames in common
mm.exon <- binarizeCSSTs(eids.on.both,CSSs.tc.exon,data.matrix=dm.exon,ratio.cutoff=2.5, abs.cutoff=200)
mm.3prime <- binarizeCSSTs(eids.on.both,CSSs.tc.3prime,data.matrix=dm.3prime,ratio.cutoff=2.5, abs.cutoff=200)
##mm  <- cbind(mm.3prime[,threeprime.set],mm.exon[,exon.set])
mm <-cbind(mm.3prime,mm.exon)

save(mm,file=paste(pdata.dir,"mm.RData",sep="/"))

## Per gene, how many conditions is it expressed in 
n.expressed.conds <- apply(mm,1,sum)
genes.not.expressed <- names(which( n.expressed.conds == 0 ))
genes.always.expressed <- names(which( n.expressed.conds == length(CSSs.tc) ))
variable.genes <- setdiff(rownames(mm),union(genes.not.expressed,genes.always.expressed))
