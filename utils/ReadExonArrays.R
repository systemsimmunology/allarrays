
## Read expression sets
## Usage
## 

##R --slave --vanilla --arch x86_64 < ReadExonArrays.R 


##
## Read Data_Matrix to dm
##  
dm.file <- "Data_Matrix.tsv"
dmd <- as.data.frame(read.table(dm.file,sep="\t",as.is=TRUE,row.names=2,header=TRUE,strip.white=TRUE))
## The "row.names=2,header=TRUE" seems to mangle strings containing a dash (converts to ".")
## Correct this as follows
dm.header <- as.character(unlist(read.table(dm.file,sep="\t",as.is=TRUE,nrows=1)))
conds <- dm.header[9:length(dm.header)]
eids.onarray <- rownames(dmd)
##This is commented out as we use most recent NCBI mappings
##genesymbol <- dm[,2] 
##names(genesymbol) <- eids.onarray
dm <- data.matrix(dmd[,-(1:7)])
colnames(dm) <- conds 
dm <- 2^dm ## Natural units
save(dm,file="dm.RData")

##
## Read expression_set to es
##  
es.file <- "expression_set.tsv"
esd <- as.data.frame(read.table(es.file,sep="\t",as.is=TRUE,row.names=1,header=TRUE,strip.white=TRUE))
es <- as.matrix(esd)

## The "row.names=2,header=TRUE" seems to mangle strings containing a dash (converts to ".")
## Correct this as follows
es.header <- as.character(unlist(read.table(es.file,sep="\t",as.is=TRUE,nrows=1)))
samps <- es.header
colnames(es) <- samps
# check if you like
identical(rownames(es),as.character(sapply(eids.onarray,paste,"_at",sep="")))
rownames(es) <- eids.onarray
es <- 2^es ## Natural units
save(es,file="es.RData")

## Input sample to meta mapping file is called Mouse_Exon.txt
im.file <- "Mouse_Exon.txt"
imd <- read.table(im.file,sep="\t",as.is=TRUE,header=FALSE,strip.white=TRUE)
## strip path in front of cel files
im.samps <- sapply(imd$V1,substr,47,10000) ## FRAGILE !!
im.grps <- imd$V2
ugs <- unique(im.grps)
im <- list()
for ( ug in ugs ){
  im[[ug]] <- as.character(im.samps[which(im.grps==ug)])
}
save(im,file="im.RData")
