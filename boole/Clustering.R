

blank.conds <- names(which(apply(mm[eids,],2,sum)==0))
var.conds <- names(which(apply(mm[eids,],2,sum)!=0))
## We will only use the var.conds from now on
var.conds <- setdiff(var.conds, "BMDM_Bl6_salmonella.flgBpEm62__Male" ) ## male and female identical?
var.conds <- setdiff(var.conds, "BMDM_Bl6_Estradiol_Ifn gamma_Female" ) ## Identical to BMDM_Bl6_Ifn gamma__Male
var.conds <- setdiff(var.conds, "mTEC_Bl6_X31 influenza__Female" ) ## Identical to mTEC_Bl6_PR8 influenza__Female

### Clustering of binary vectors

## Hamming distance and dendrogram of conditions

##
## Hamming distance
##
## Can be computed using the package e1071
### hd.conds <- hamming.distance(t(mm[eids,var.conds])) ## requires library(e1071)
### dist.conds <- as.dist(hd.conds)
### 
### equivalently, use manhattan distance
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
cell.type <- as.character(lapply(CSSs.tc[dendrorder.conds],"[[","Cell Type"))
strain <- as.character(lapply(CSSs.tc[dendrorder.conds],"[[","Strain"))
stimulus <- as.character(lapply(CSSs.tc[dendrorder.conds],"[[","Stimulus 1"))
condlabs <- paste(paste(cell.type,strain,sep="_"),stimulus,sep="_")

# Hamming distance and dendrogram of genes
dist.genes <- dist(mm[eids,var.conds],method='manhattan')
hc.genes <- hclust( dist.genes, "complete" )
hc.genes$labels <- gene.symbol[hc.genes$labels]
dendrorder.genes <- hc.genes$labels[hc.genes$order] ## Genes in dendrogram order
dg <- dendrorder.genes
dist.genes.mat <- as.matrix(dist.genes)

x11()
hp <- plot(hc.genes, main="complete")

## Horizontal plot

plot(as.dendrogram(hc.genes),horiz=T)

## Clustering 

##
## Clusters, genes
##  

k=3
cluster.members.genes <- cutree(hc.genes,k)
cmoc <- cluster.members.genes[hc.genes$order]

## White off Black on
collars <- c("white","black")

x11()
par(mfrow=c(2,2))

c1 <- gene.eids[names(which(cmoc==1))]
image(t(mm[c1,dendrorder.conds]),col=collars,main="Cluster 1")

c2 <- gene.eids[names(which(cmoc==2))]
image(t(mm[c2,dendrorder.conds]),col=collars,main="Cluster 2")

c3 <- gene.eids[names(which(cmoc==3))]
image(t(mm[c3,dendrorder.conds]),col=collars,main="Cluster 3")

##
## Multidimensional scaling
## http://www.statmethods.net/advstats/mds.html
fit <- cmdscale(dist.genes,eig=TRUE,k=2)
fit 

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
## Clusters, conditions
##

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


fit <- isoMDS(dist.conds, k=2) # k is the number of dim

# plot solution
x11()
x <- fit$points[,1]
y <- fit$points[,2]

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Nonmetric MDS", type="n")

##plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Nonmetric MDS")

text(x, y, labels=var.conds,cex=0.9)


## Careful, this can be  too big for X11
Rowv  <- as.dendrogram(hc.mm)
hv <- heatmap(t(mm[eids,]),Rowv=Rowv, scale="none",margins=c(15,15),cexCol=0.9,cexRow=1.2,labCol=gsym[eids])
so <- colnames(mm)[hv$rowInd][129:1]

