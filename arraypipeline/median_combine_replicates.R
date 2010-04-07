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

emat <- as.matrix(read.table(expfile,as.is=TRUE))
## More faithful version of columns headers takes care of
## R incorrectly inserts X before colheaders that begin with a number
## and converting "-" to "."
colnames(emat) <- strsplit(readLines(expfile,n=1),'\t')[[1]][-1]

## 
h <- numeric()
for ( uv in uvs ){
  cols <- names(which(v==uv))
  g <- apply(emat[,cols],1,median)
  h <- cbind(h,g)
}
colnames(h) <- uvs

matrixPrintFormat <- 
function( matrix,topLeftString="" ){
  rownames <- rownames(matrix)
  colnames <- colnames(matrix)
  outMat <- cbind(rownames,matrix)
  outMat <- rbind(c(topLeftString,colnames),outMat)
  rownames(outMat) <- NULL
  colnames(outMat) <- NULL
  return(outMat)
}

write.table(matrixPrintFormat(h),file=outfile,sep="\t",quote=FALSE,row.names = FALSE, col.names = FALSE)
##save(eemat,file="eemat.RData")

