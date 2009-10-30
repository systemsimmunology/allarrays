##
## sliceAndDiverify
##
### combines sliceByAvailableRow and diversify
sliceAndDiversify <- function (rowlabs,inMat,col.on.min=NULL,col.on.max=NULL, row.on.min=NULL,row.on.max=NULL ){
  inMat.f1 <- sliceByAvailableRows(rowlabs,inMat)
  outMat <- diversify( inMat.f1, col.on.min,col.on.max,row.on.min,row.on.max)
}

##
## sliceByAvailableRows
##
## For given query of row labels, return a sliced matrix if rows are availabe 
sliceByAvailableRows <- function ( rowlabs, inMat ) {
  rowlabs.f1 <- intersect(rowlabs,rownames(inMat))
  inMat.f1 <- inMat[rowlabs.f1,]
  cat("Of ",length(rowlabs), "rows, ",length(rowlabs.f1)," were found in input matrix\n") 
  inMat.f1
}

##
## diversify
##
## Input boolean matrix and return boolean matrix satifying the following diversity conditions
##
## col.on.min: For any row, the *fraction* of 1s must exceed this value
## col.on.max: For any row, the *fraction* of 1s must not exceed this value

## NB: Denominators for fractions are from original matrix 


## If no thresholds are given, return matrix in which all rows and columns are variable ( not all zero or all one )

diversify <- function( inMat, col.on.min=NULL,col.on.max=NULL, row.on.min=NULL,row.on.max=NULL ){

  tiny <- 1e-10
  ## tiny threshold amounts to >0 and <1 
  if ( is.null(col.on.min) ){ col.on.min <- tiny }
  if ( is.null(col.on.max) ){ col.on.max <- 1-tiny }
  if ( is.null(row.on.min) ){ row.on.min <- tiny }
  if ( is.null(row.on.max) ){ row.on.max <- 1-tiny }
  
  ## Per row, how many cols are one
  n.on.cols <- apply(inMat,1,sum)
  ## The fraction of rows in which as col is 1 
  frac.on.cols <- n.on.cols/ncol(inMat) 

  row.keepers <- names(which( (frac.on.cols>col.on.min ) & (frac.on.cols<col.on.max ) ))
  
  ## Per column, how many rows are 1 
  n.rows.on <- apply(inMat,2,sum)
  ## The fraction rows on in each col
  frac.rows.on <- n.rows.on/nrow(inMat) 

  col.keepers <- names(which( (frac.rows.on>row.on.min ) & (frac.rows.on<col.on.max ) ))
                    
  outMat <- inMat[row.keepers,col.keepers]

  cat(length(row.keepers), "rows, and ",length(col.keepers)," cols pass threshold.\n")

  outMat

}
  

writeBooleMat <- function ( inMat, prefix = "", outdir ){
  ofile <- paste(outdir,"/",prefix,"booleMat.tsv",sep="")
  mm.ca.mpf <- matrixPrintFormat(mm.ca,topLeftString="GeneID")
  write.table(file=ofile,mm.ca.mpf,sep="\t",row.names=FALSE,col.names=FALSE)
  set <- rownames(mm.ca)
  ## For MATLAB write separate files with values and labels
  ofile <- paste(outdir,"CAbooleVals.tsv",sep="/")
  write.table(file=ofile,mm.ca,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  ofile <- paste(outdir,"CAbooleGeneIDs.tsv",sep="/")
  write.table(file=ofile,cbind(set,gene.symbol[set]),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  ofile <- paste(outdir,"CAbooleMatWithGeneSymbols.tsv",sep="/")
  mm.ca.mpf[,1] <- c("GeneSymbol",gene.symbol[set])
  write.table(file=ofile,mm.ca.mpf,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

}
