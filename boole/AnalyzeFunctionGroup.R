
##Set mmg to be the binary matrix for group of interest
mmg <- mm.ca
eids <- rownames(mm.ca)

dist.genes <- dist(mmg,method='manhattan')
hc.genes <- hclust( dist.genes, "complete" )
hc.genes$labels <- gene.symbol[hc.genes$labels]
dendrorder.genes <- hc.genes$labels[hc.genes$order] ## Genes in dendrogram order
dg <- dendrorder.genes
dist.genes.mat <- as.matrix(dist.genes)

## Horizontal plot
plot(as.dendrogram(hc.genes),horiz=T)

##
## Multidimensional scaling
## http://www.statmethods.net/advstats/mds.html
fit <- cmdscale(dist.genes,eig=TRUE,k=2)


# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type="n")
text(x, y, labels=gene.symbol[eids],cex=0.9) 

# Nonmetric MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

library(MASS)
 
svec <- vector(length=10)
svec <- 0

for ( k in 1:10 ){
  fit <- isoMDS(dist.genes, k=k) # k is the number of dim
  svec[k] <- fit[["stress"]] 
}

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
x11()
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Nonmetric MDS", type="n")
text(x, y, labels=gene.symbol[eids],cex=0.9)

##
## Heatmap
##
genes.for.plot <- eids
heatmap(t(mmg[genes.for.plot,]),scale="none",margins=c(15,15),cexCol=0.9,cexRow=1.2,labCol=gsym[genes.for.plot])

