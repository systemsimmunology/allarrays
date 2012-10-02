##
## Retrieve list required to rerun normalization of allarrays
## for BrainArray annotations
##

library(rjson)
source("~/allarrays/utils/httpget.R")
source("~/allarrays/utils/utilitiesMetaData.R")
source("~/allarrays/utils/utilitiesSampleGroup.R")

## Location of output files from automated pipeline
repo.dir <- file.path(Sys.getenv("ILYA_LAB"),"Vesteinn/data/ImmunoRepository/sampleData/microarray/runs/Aderem Three Prime Arrays_2009-05-07_at_00.06.21/output/Mouse 430 2.0")

load(paste(c(repo.dir),"sglist.RData",sep="/")) ## 228 sgs, same length as dm.columns
load(paste(c(repo.dir),"dm.columns.RData",sep="/")) 

## Step through list by index
sglist.names <- names(sglist)
n.sglist <- length(sglist.names)
sg.info <- character() ## for sample-group specifics
sample.info <- character() ## for sample specifics
for ( i in 1:n.sglist ){
  repository.name <- sglist.names[i]
  dmcol.name <- as.character(dm.columns[repository.name])
  gno <- sglist[[i]] ## sg object
  sgname <- nullToBlank("sample_group_name",gno)
  vecy <- c(gno[["Cell Type"]],gno[["Strain"]],
            nullToBlank("Stimulus 1",gno),nullToBlank("Time 1",gno),
            gno[["Sex"]])
  vecyy <- c(gno[["Cell Type"]],gno[["Strain"]],
            nullToBlank("Stimulus 1",gno),nullToBlank("Time 1",gno),
            nullToBlank("Stimulus 2",gno),nullToBlank("Time 2",gno),
            gno[["Sex"]])
  csst.name <- paste(vecyy,collapse="_")
  vecy <- c(c(csst.name,repository.name,dmcol.name,sgname),vecy)
  if ( length(vecy) != 9 ){stop()}
  sg.info <- rbind(sg.info,vecy)
  for(uuid in gno$sampleUUIDs){
    rdp <- getNodeObject(uuid)[["raw_data_path"]]
    celfile <- basename(rdp)
    cinfo <- c(celfile,dmcol.name,rdp)
    sample.info <- rbind(sample.info,cinfo)
  }
}
rownames(sg.info) <- NULL
rownames(sample.info) <- NULL
colnames(sg.info) <- c("CSST","sglist Name","Data_Matrix Column","sample_group_name","Cell Type","Strain","Stimulus 1","Time 1","Sex")

##
## Write sample group information into .TSV
## 

write.table(sg.info,file="AllThreePrimeArraysSampleGroupInfo.tsv",quote=FALSE,sep="\t",row.name=F)

write.table(sample.info,file="AllThreePrimeArraysSampleInfo.tsv",quote=FALSE,row.names=F,col.names=F,sep="\t")

## Get unique rows after the fact using unix 
