
eids.on.exon.array <- rownames(dm.exon)
eids.on.3prime.array <- rownames(dm.3prime)
eids.on.both <- intersect(eids.on.exon.array,eids.on.3prime.array)

mm.exon <- binarizeCSSTs(eids.on.both,CSSs.tc.exon,data.matrix=dm.exon,ratio.cutoff=2.5, abs.cutoff=200)
mm.3prime <- binarizeCSSTs(eids.on.both,CSSs.tc.3prime,data.matrix=dm.3prime,ratio.cutoff=2.5, abs.cutoff=200)
mm  <- cbind(mm.3prime[,threeprime.set],mm.exon[,exon.set])
n.expressed.conds <- apply(mm,1,sum)
genes.not.expressed <- eids[which( n.expressed.conds == 0 )]
genes.always.expressed <- eids[which( n.expressed.conds == length(CSSs.tc) )]
n.expressed.genes <- apply(mm,2,sum)

###
### Functional groups
###
ca <- as.character(read.table("/Users/thorsson/data/GeneOntology/CytokineActivity.tsv",as.is=TRUE)$V1)
cb <- as.character(read.table("/Users/thorsson/data/GeneOntology/CytokineBinding.tsv",as.is=TRUE)$V1)
extcell <- as.character(read.table("/Users/thorsson/data/GeneOntology/ExtracellularRegion.tsv",as.is=TRUE)$V1)
tfa <- as.character(read.table("/Users/thorsson/data/GeneOntology/TranscriptionFactorActivity.tsv",as.is=TRUE)$V1)
m2a <- as.character(read.table("/Users/thorsson/data/MacPolarization/WoundHealingM2a.tsv",as.is=TRUE,sep='\t',header=TRUE)[,"Gene.ID"])
## Intersecting extcell with ca seems to have little effect now
## Just skip it
eids <- intersect(ca,eids.on.both)


## requires mm 
write.table(file="CAboole.tsv",mm[intersect(ca,eids.on.both),],quote=FALSE,sep="\t")
set <- intersect(ca,eids.on.both)
write.table(file="CAbooleVals.tsv",mm[set,],quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(file="CAbooleGeneIDs.tsv",cbind(set,gene.symbol[set]),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
## For MATLAB write separate files with values and labels
set <- intersect(ca,eids.on.both)
write.table(file="CAbooleVals.tsv",mm[set,],quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(file="CAbooleGeneIDs.tsv",cbind(set,gene.symbol[set]),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

##########################################
## M2-ness scores
##########################################
v <- apply(mm,2,sum)

###################
##### Older code
####################
 
mm <- matrix(nrow=length(eids),ncol=length(CSSs.tc))
rownames(mm) <- eids
colnames(mm) <- names(CSSs.tc)             
for ( eid in eids ){
  bvec.exon <- binarizeCSSTs(eid,CSSs.tc.exon,data.matrix=dm.exon,cutoff=5)
  bvec.3prime <- binarizeCSSTs(eid,CSSs.tc.3prime,data.matrix=dm.3prime,cutoff=5)
  mm[eid,names(bvec.exon)] <- bvec.exon
  mm[eid,names(bvec.3prime)] <- bvec.3prime
}

n.expressed.conds <- apply(mm,1,sum)
genes.not.expressed <- eids[which( n.expressed.conds == 0 )]
genes.always.expressed <- eids[which( n.expressed.conds == length(CSSs.tc) )]



genes.for.plot <- setdiff(eids,union(genes.not.expressed,genes.always.expressed))
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
