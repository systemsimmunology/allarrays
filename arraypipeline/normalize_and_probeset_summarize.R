## Invoke using
## R --vanilla --slave --args " "CEL list, no path""CEL file folder" < normalize_and_probeset_summarize.R

library(affy)
library(mouse4302mmentrezgcdf) ## R CMD INSTALL mouse4302mmentrezgcdf_12.1.0.tar.gz

library(mogene10stv1mmentrezgcdf)

ca <- commandArgs()

filelist <- ca[5]
celfile.path <- ca[6]
  
## Arguments for ReadAffy
## Filename to specify files and filename order
filenames <- as.vector(read.table(filelist,as.is=TRUE)$V1)
## celfile.path

## Read probevel data into AffyBatch object
rawdata<-ReadAffy(filenames=filenames,celfile.path=celfile.path )

## Change to custom CDF 
rawdata@cdfName <- "Mouse4302_Mm_ENTREZG"
rawdata@cdfName <- "MoGene10stv1_Mm_ENTREZG"

## Background correct, Normalize, Probeset Summarize
eset <- rma(rawdata) ## rma or gcrma
emat <- exprs(eset)


## Insert the leading cell of matrix for printing
matrixPrintFormat <- function( matrix,topLeftString="" ){
  rownames <- rownames(matrix)
  colnames <- colnames(matrix)
  outMat <- cbind(rownames,matrix)
  outMat <- rbind(c(topLeftString,colnames),outMat)
  rownames(outMat) <- NULL
  colnames(outMat) <- NULL
  return(outMat)
}

## Save
save(emat,file="emat.RData")
write.table(matrixPrintFormat(emat),file="emat.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
