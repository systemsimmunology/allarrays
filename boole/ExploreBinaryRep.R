

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
