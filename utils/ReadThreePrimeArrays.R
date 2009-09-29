## Read expression sets
## Usage
## 

##R --vanilla --slave --arch x86_64 < ReadThreePrimeArrays.R 

ag.file <- "AllGenes.txt"
agd <- as.data.frame(read.table(ag.file,sep="\t",as.is=TRUE,row.names=2,header=TRUE,strip.white=TRUE))
ag.header <- as.character(unlist(read.table(ag.file,sep="\t",as.is=TRUE,nrows=1)))
conds <- ag.header[10:length(ag.header)]
ag <- as.matrix(agd[9:length(agd)])
colnames(ag) <- conds
ag <- 2^ag ## Natural units
dm <- ag
save(dm,file='dm.RData')

