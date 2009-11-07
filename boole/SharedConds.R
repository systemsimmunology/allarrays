
## for a group of genes
## thought to be induced under
## subset of conditions

mmat <- mm.tfa
conds <- colnames(mmat)
genes <- rownames(mmat)

## what is/are the 'shared meta data feature(s)'
dist.genes <- dist(mmat,method='manhattan')
dist.genes.mat <- as.matrix(dist.genes)

es <- gene.eid[c("Stat2","Irf7","Stat1","Znfx1","Batf2")]

reduceByMetaData <- function ( mmat, mstring="Cell Type"){

  conds <- colnames(mmat)
  genes <- rownames(mmat)

  mdatapercond <- unlist(lapply(CSSs.tc[conds],"[[",mstring))
  names(mdatapercond) <- conds
  
  cb <- table(mdatapercond)
  cts <- names(cb)
  vab <- as.vector(cb)
  names(vab) <- cts ## tally conditions available for each metadata value
  
  bb <- matrix(nrow=length(genes),ncol=length(cts))
  rownames(bb) <- genes
  colnames(bb) <- cts
  for ( ct in cts ){
    hh <- mmat[,which(mdatapercond==ct)]
    if ( !is.vector(hh)) {
      bb[,ct] <- apply(hh,1,sum)
    } else {
      bb[,ct] <- hh
    }  
  }
  olist <- list()
  olist[["mmeta"]] <- bb
  olist[["tmeta"]] <- vab
  olist
}

scoreMetaRep <- function ( eids, mmeta, tmeta ){
  bbb <- mmeta[eids,]
  apply(t(t(bbb)/tmeta),2,mean)
}
                          

res <- reduceByMetaData(mmat,"Cell Type")
sort(scoreMetaRep(es,res[["mmeta"]],res[["tmeta"]]),decreasing=TRUE)

res <- reduceByMetaData(mmat,"Strain")
sort(scoreMetaRep(es,res[["mmeta"]],res[["tmeta"]]),decreasing=TRUE)

res <- reduceByMetaData(mmat,"Stimulus 1")
sort(scoreMetaRep(es,res[["mmeta"]],res[["tmeta"]]),decreasing=TRUE)


bbb <- bb[es,]

## aggregate measure

## normalize counts by total possible counts
## average over genes set

apply(t(t(bbb)/vab),2,mean)


###
### Attic
###


mv <- mmat[es,]

e1 <- names(which(mv[1,]==1))
e2 <- names(which(mv[2,]==1))
e3 <- names(which(mv[3,]==1))
e4 <- names(which(mv[4,]==1))

w1 <- table(strain[e1])/length(e1)
w2 <- table(strain[e2])/length(e2)
w3 <- table(strain[e3])/length(e3)
w4 <- table(strain[e4])/length(e4)

wb <- table(strain)/137

w1/wb[names(w1)]
w2/wb[names(w2)]
w3/wb[names(w3)]
w4/wb[names(w4)]


c1 <- table(celltype[e1])
c2 <- table(celltype[e2])
c3 <- table(celltype[e3])
c4 <- table(celltype[e4])


c1/cb[names(c1)]
c2/cb[names(c2)]
c3/cb[names(c3)]
c4/cb[names(c4)]
