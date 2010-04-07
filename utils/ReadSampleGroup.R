
##
## Read sample groups from file. Arrange into time-courses
##
## Invoke using 
## R --arch x86_64 --slave --args "Sample Group MetaData (tab-delimited)" < ReadSampleGroups.R 

ca <- commandArgs()
infile  <- ca[4]

source("/Users/thorsson/allarrays/utils/utilitiesSampleGroup.R")
source("/Users/thorsson/allarrays/utils/utilitiesMetaData.R")

sglist <- readSGfromfile(infile)

CSSs <- getUniqueCSSs(sglist)

CSSs.tc <- findTimeCourse(CSSs, sglist)

zeroTimes <- findTimeZero( CSSs, sglist)

### Add zero Times to time course objects
CSSs.tc <- includeZeroTime( CSSs.tc, zeroTimes )

##
## sample group to expression matrix column headers
## Commented out for now as there is a simple exact match
##
if(FALSE){
if( file.access(paste(data.path,"eemat.tsv",sep="/")) ==0 ){
  dm.file <- paste(data.path,"eemat.tsv",sep="/")
  dm.header <- as.character(unlist(read.table(dm.file,sep="\t",as.is=TRUE,nrows=1)))
  conds <- dm.header[2:length(dm.header)]
  dm.columns <- conds
  names(dm.columns) <- conds
  ##  dm.columns <- DMcolsFromSG( sglist, conds, exon.form=FALSE )
}
}

dm.columns <- names(sglist)
names(dm.columns) <- names(sglist)
### It will be convenient to carry these mappings with SCCs.tc
for ( tc in CSSs.tc ){
  dmc <- as.character(dm.columns[tc[["Sample Group"]]])
  CSSs.tc[[tc$name]][["DM Column"]] <- dmc
}


save(sglist,file="sglist.RData")
save(CSSs,file="CSSs.RData")
save(CSSs.tc,file="CSSs.tc.RData")

