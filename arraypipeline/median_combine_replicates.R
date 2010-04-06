## Invoke using
## R --vanilla --slave --args "Expression matrix with replicates" "Group File" "Output File" < median_combine_replicates.R

ca <- commandArgs()
expfile <- ca[5]
groupfile <- ca[6]
outfile <- ca[7]

## Load Group Files
rt <- read.table(groupfile,as.is=TRUE)
v <- rt$V2
names(v) <- rt$V1
uvs <- unique(v)

## Load Expression Matrix before replicate combination
load("~/allarrays/tmp/emat.RData")

rt <- as.matrix(read.table("~/allarrays/tmp/emat.tsv",as.is=TRUE,colClasses = "character"))

## 
h <- numeric()
for ( uv in uvs ){
  cols <- names(which(v==uv))
  g <- apply(emat[,cols],1,median)
  h <- cbind(h,g)
}
colnames(h) <- uvs

function( matrix,topLeftString="" ){
  rownames <- rownames(matrix)
  colnames <- colnames(matrix)
  outMat <- cbind(rownames,matrix)
  outMat <- rbind(c(topLeftString,colnames),outMat)
  rownames(outMat) <- NULL
  colnames(outMat) <- NULL
  return(outMat)
}

eemat <- h

write.table(matrixPrintFormat(h),file="eemat.tsv",sep="\t",quote=FALSE,row.names = FALSE, col.names = FALSE)
save(eemat,file="eemat.RData")

