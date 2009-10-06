
dMat <- 1-cor(t(lps.ratios.w.zero[targSet,]))
dr <- dist(dMat)

hc <- hclust( dr , "complete" )
x11()
plot(hc, main="complete")


## Genes expressed in more than some number of  conditions

eids <- names(which(apply(mm,1,sum)>10))

dMat.exon <- 1-cor(dm.exon[eids,])
dr.exon <- dist(dMat.exon)
hc.exon <- hclust( dr.exon , "complete" )
x11()
plot(hc.exon, main="complete")



dMat.3prime <- 1-cor(dm.3prime[eids,])
dr.3prime <- dist(dMat.3prime)
hc.3prime <- hclust( dr.3prime , "complete" )
x11()
plot(hc.3prime, main="complete")


## Hamming distance and dendrogram of conditions
hd.conds <- hamming.distance(t(mm[eids,])) ## requires library(e1071)

dist.conds <- as.dist(hd.conds)
#### equivalently, use manhattan distance
dist.conds <- dist(t(mm[eids,]),method='manhattan')
### 'euclidian' also equivalent, after squaring the result ! 
###  'binary' seems to be #XOR/#OR so don't use that
## the e1071 function seems slow
## system.time(a <- as.dist(hamming.distance(t(mm[eids,]))))
## loses to
## system.time(a<- dist(t(mm[eids,]),method='manhattan'))
hc.conds <- hclust( dist.conds, "complete" )
x11()
hp <- plot(hc.conds, main="complete")
dendrorder.conds <- hc.conds$labels[hc.conds$order] ## Conditions in dendrogram order

# Hamming distance and dendrogram of genes
dist.genes <- dist(mm[eids,var.conds],method='manhattan')
hc.genes <- hclust( dist.genes, "complete" )
dendrorder.genes <- hc.genes$labels[hc.genes$order] ## Genes in dendrogram order
dg <- dendrorder.genes


x11()
hp <- plot(hc.genes, main="complete")

k=4
hc.conds <- hclust( dist.conds, "complete" )
cluster.members.conds <- cutree(hc.conds,k)
cmoc <- cluster.members.conds[hc.conds$order]

## White off Black on
collars <- c("white","black")



c1 <- names(which(cmoc==1))
image(mm[dendrorder.genes,c1],col=collars,main="Cluster 1")



c2 <- names(which(cmoc==2))
image(mm[dendrorder.genes,c2],col=collars,main="Cluster 2")


c6 <- names(which(cmoc==6))
image(mm[dendrorder.genes,c6],col=collars,main="Cluster 6")



c8 <- names(which(cmoc==8))
image(mm[dendrorder.genes,c8],col=collars,main="Cluster 8")



##c1 should be split into three!
dc1 <- dist(t(mm[eids,c1]),method='manhattan')
hc1 <- hclust(dc1,"complete")
cluster.members.conds.1 <- cutree(hc1,3)
cmoc.1 <- cluster.members.conds.1[hc1$order]
c11 <- names(which(cmoc.1==1))
image(mm[dendrorder.genes,c11],col=collars,main="Cluster 1.1")
c13 <- names(which(cmoc.1==3))
image(mm[dendrorder.genes,c13],col=collars,main="Cluster 1.3")




## Careful, this can be  too big for X11
Rowv  <- as.dendrogram(hc.mm)
hv <- heatmap(t(mm[eids,]),Rowv=Rowv, scale="none",margins=c(15,15),cexCol=0.9,cexRow=1.2,labCol=gsym[eids])
so <- colnames(mm)[hv$rowInd][129:1]


orderFromLabel <- hc.mm$order
names(orderFromLabel) <- hc.mm$labels



## this will get cluster members in plot order
cmo <- cluster.members[hc.mm$order]
c3 <- names(which(cmo==3))
image(mm[eids,c3])

c6 <- colnames(mm)[which(cluster.members==6)]
image(t(mm[eids,c6]))
    
## Temporary code
for ( i in 1:68 ){
  cat( i, hamming.distance(mm[dg[i],],mm[dg[i+1],]), "\n")
}

