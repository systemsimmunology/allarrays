
util.dir <- file.path(Sys.getenv("AA"),"utils")
pdata.dir <- file.path(Sys.getenv("AA"),"processed_data/20091015") 
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))

load(paste(pdata.dir,"mm.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.RData",sep="/"))

##
## All genes
## 
## Per gene, how many conditions is it expressed in 
n.expressed.conds <- apply(mm,1,sum)
genes.not.expressed <- names(which( n.expressed.conds == 0 ))
genes.always.expressed <- names(which( n.expressed.conds == length(CSSs.tc) ))
variable.genes <- setdiff(rownames(mm),union(genes.not.expressed,genes.always.expressed))

## Per condition, how many genes are expressed
n.expressed.genes <- apply(mm,2,sum)
blank.conditions <- names(which(n.expressed.genes==0))

##
## Genes in Subset of interest
##
go.dir <- file.path(Sys.getenv("DATA_DIR"),"GeneOntology")
ca <- as.character(read.table(paste(go.dir,"CytokineActivity.tsv",sep="/"),as.is=TRUE)$V1)

eids <- intersect(ca,variable.genes)
blank.conds <- names(which(apply(mm[eids,],2,sum)==0))
var.conds <- names(which(apply(mm[eids,],2,sum)!=0))

#
# 
#

genes.for.plot <- eids
heatmap(t(mm[genes.for.plot,]),scale="none",margins=c(15,15),cexCol=0.9,cexRow=1.2,labCol=gsym[genes.for.plot])

maxes <- matrix(nrow=length(eids),ncol=length(CSSs.tc))
rownames(maxes) <- eids
colnames(maxes) <- names(CSSs.tc)
## This may be meaningless for two array types ???

for ( csst in CSSs.tc.3prime ){
  relcols <- csst[["DM Column"]]
  maxvec <- apply(dm.3prime[eids,relcols],1,max)
  maxes[eids,csst[["name"]]] <- maxvec
}

for ( csst in CSSs.tc.exon ){
  relcols <- csst[["DM Column"]]
  maxvec <- apply(dm.exon[eids,relcols],1,max)
  maxes[eids,csst[["name"]]] <- maxvec
}


genes.for.plot <- setdiff(eids, genes.not.expressed )
heatmap(t(log(maxes[genes.for.plot,])),scale="none",margins=c(15,15),cexCol=0.9,cexRow=1.2,labCol=gsym[genes.for.plot])


## Explore classes
table(unlist(lapply(CSSs.tc,"[[","Strain")))

##
## Explore results of influence graph (with Timo, Ilya)
## 

minf <- as.matrix(read.table("/Volumes/ILYA LAB/Cytokine\ network/Matlab/G.tsv",sep="\t",as.is=TRUE))      
ca <- rownames(mm.ca)
colnames(minf) <- ca
rownames(minf) <- ca

## top influences


winners <- matrixExtrema(minf,10,decreasing=TRUE)$pairs

cl <- integer()
for ( i in 1:nrow(winners) ){
  one <- winners[i,1]
  two <- winners[i,2]
  ds <- dist.genes.mat[one,two]
  cat(gs[one],"->",gs[two],"\t",minf[one,two],"\t",ds,"\n")
  cl <- c(cl,ds)
}
