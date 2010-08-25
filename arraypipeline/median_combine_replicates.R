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


emat <- read.matrix(expfile)

h <- numeric()
for ( uv in uvs ){
  cols <- names(which(v==uv))
  g <- apply(emat[,cols],1,median)
  h <- cbind(h,g)
}
colnames(h) <- uvs

write.table(matrixPrintFormat(h),file=outfile,sep="\t",quote=FALSE,row.names = FALSE, col.names = FALSE)
eemat <- h
save(eemat,file="eemat.RData")

