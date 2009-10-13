
util.dir <- file.path(Sys.getenv("AA"),"utils")
pdata.dir <- file.path(Sys.getenv("AA"),"processed_data/20090507") 
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

###
### Functional groups
###
go.dir <- file.path(Sys.getenv("DATA_DIR"),"GeneOntology")
ca <- as.character(read.table(paste(go.dir,"CytokineActivity.tsv",sep="/"),as.is=TRUE)$V1)
cb <- as.character(read.table(paste(go.dir,"CytokineBinding.tsv",sep="/"),as.is=TRUE)$V1)
extcell <- as.character(read.table(paste(go.dir,"ExtracellularRegion.tsv",sep="/"),as.is=TRUE)$V1)
tfa <-  as.character(read.table(paste(go.dir,"TranscriptionFactorActivity.tsv",sep="/"),as.is=TRUE)$V1)
##m2a <- as.character(read.table("/Users/thorsson/data/MacPolarization/WoundHealingM2a.tsv",as.is=TRUE,sep='\t',header=TRUE)[,"Gene.ID"])

## Intersecting extcell with ca seems to have little effect now
## Just skip it
eids <- intersect(ca,variable.genes)

## requires mm
ncbi.dir <- file.path(Sys.getenv("DATA_DIR"),"ncbi")
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))

ofile <- paste(pdata.dir,"CAboole.tsv",sep="/")
write.table(file=ofile,mm[intersect(ca,eids.on.both),],quote=FALSE,sep="\t")

## For MATLAB write separate files with values and labels
set <- intersect(ca,eids.on.both)
ofile <- paste(pdata.dir,"CAbooleVals.tsv",sep="/")
write.table(file=ofile,mm[set,],quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
ofile <- paste(pdata.dir,"CAbooleGeneIDs.tsv",sep="/")
write.table(file=ofile,cbind(set,gene.symbol[set]),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
