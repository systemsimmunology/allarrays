
### Clustering of binary vectors

## Hamming distance and dendrogram of conditions
hd.conds <- hamming.distance(t(mm[eids,var.conds])) ## requires library(e1071)

dist.conds <- as.dist(hd.conds)
#### equivalently, use manhattan distance
dist.conds <- dist(t(mm[eids,var.conds]),method='manhattan')
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


x11()
par(mfrow=c(2,2))

c1 <- names(which(cmoc==1))
image(mm[dendrorder.genes,c1],col=collars,main="Cluster 1")

c2 <- names(which(cmoc==2))
image(mm[dendrorder.genes,c2],col=collars,main="Cluster 2")

c3 <- names(which(cmoc==3))
image(mm[dendrorder.genes,c3],col=collars,main="Cluster 3")

c4 <- names(which(cmoc==4))
image(mm[dendrorder.genes,c4],col=collars,main="Cluster 4")


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

