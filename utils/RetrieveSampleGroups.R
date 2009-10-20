
##
## Read sample groups. Arrange into time-courses
##
## Invoke using 
## R --arch x86_64 --slave --args "Full JCR Path To Sample Group" < RetrieveSampleGroups.R 

## E.g
## R --arch x86_64 --slave --args "/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - BMDM Arrays_2009-03-29_at_09.20.12/output/Mouse Exon/sampleGlossary" < /Volumes/thorsson/MetaAnalysis2009/RShared/RetrieveSampleGroups.R
ca <- commandArgs()
sg.path <- ca[4]
data.path <- "."

##test values

data.path <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-10-15_at_16.51.44"
sg.path <- "/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-10-15_at_16.51.44/output/Mouse Exon/sampleGlossary"

##data.path <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-08-22_at_09.20.09"
##sg.path <- "/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-08-22_at_09.20.09/output/Mouse Exon/sampleGlossary"

##data.path <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-05-02_at_09.20.05"
##sg.path <- "/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-05-02_at_09.20.05/output/Mouse Exon/sampleGlossary"

##data.path <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/ExpressionSets/GeneLevel/Genomics Expression Public Dataset"
##sg.path <- "/ExpressionSets/GeneLevel/Genomics Expression Public Dataset/sampleGlossary"

##test values
##data.path <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/sampleData/microarray/runs/Aderem Three Prime Arrays_2009-05-07_at_00.06.21/output/Mouse 430 2.0"
##sg.path <- "/sampleData/microarray/runs/Aderem Three Prime Arrays_2009-05-07_at_00.06.21/output/Mouse 430 2.0/sampleGlossary"

library(rjson)
util.dir <- file.path(Sys.getenv("AA"),"utils")
sources <- c("httpget.R","utilitiesMetaData.R","utilitiesSampleGroup.R")
for ( src in sources ){source(paste(util.dir,src,sep="/"))}

## Get Sample Group object
sg.obj <- getNodeObjectByPath(sg.path)
## parse out child nodes of this object
## to yield the actual sample groups
sglist.new <- parseSampleGroup( sg.obj )

sglist.1stim <- removeUnequalTime( sglist.new )
  
sglist.nn <- repairSGgender( sglist.1stim )

bsgname <- paste(data.path,"BadSampleGroups.txt",sep="/")
if(file.access(bsgname)==0) { ## if file was found 
  baddies <- read.table(file=bsgname,as.is=TRUE,sep='\t')$V1
  sglist.nn <- sglist.nn[setdiff(names(sglist.nn),baddies)]
}
 
sglist <- sglist.nn

save(sglist,file="sglist.RData")

## Get unique CSS objects
## Stimulated
CSSs <- getUniqueCSSs( sglist )

bcname <- paste(data.path,"BadCSS.txt",sep="/")
if(file.access(bcname)==0) {
  baddies <- read.table(file=bcname,as.is=TRUE,sep='\t')$V1
  CSSs <- CSSs[setdiff(names(CSSs),baddies)]
}

save(CSSs,file="CSSs.RData")

## Unstimulated Stimulus 1 = ""
CSSs.unstim <- getUnstimmedCSSs( sglist )
## Never used?

### Create time course objects
CSSs.tc <- findTimeCourse(CSSs, sglist )

## Zero time sglabel (sample group name)for every stimmed CSS
ztname <- paste(data.path,"ZeroTimeChoices.tsv",sep="/")
if(file.access(ztname)==0) {
  zeroTimes <- findTimeZero( CSSs, sglist, zfile = ztname )
} else {
  zeroTimes <- findTimeZero( CSSs, sglist, zfile = NULL )
}

### Add zero Times to time course objects
CSSs.tc <- includeZeroTime( CSSs.tc )

save(CSSs.tc,file="CSSs.tc.RData")

##
##  This relies on the existence of an expression matrix file

if(file.access(paste(data.path,"Data_Matrix.tsv",sep="/"))==0 ) { ## Exon Arrays
  dm.file <- paste(data.path,"Data_Matrix.tsv",sep="/")
  dm.header <- as.character(unlist(read.table(dm.file,sep="\t",as.is=TRUE,nrows=1)))
  conds <- dm.header[9:length(dm.header)]
  dm.columns <- DMcolsFromSG( sglist, conds )
}
if( file.access(paste(data.path,"AllGenes.txt",sep="/")) ==0 ){
  dm.file <- paste(data.path,"AllGenes.txt",sep="/")
  dm.header <- as.character(unlist(read.table(dm.file,sep="\t",as.is=TRUE,nrows=1)))
  conds <- dm.header[10:length(dm.header)]
  dm.columns <- DMcolsFromSG( sglist, conds, exon.form=FALSE )
}

save(dm.columns,file="dm.columns.RData")

### It will be convenient to carry these mappings with SCCs.tc
for ( tc in CSSs.tc ){
  dmc <- as.character(dm.columns[tc[["Sample Group"]]])
  CSSs.tc[[tc$name]][["DM Column"]] <- dmc
}


## Before writing up, check for length consistencies and
## absence of nulls

for ( csstc in CSSs.tc ){
  if ( is.null(csstc$Sex) ){
    sto("Missing Sex attribute in",csstc$name)
  }
  s <- csstc[["Sample Group"]]
  t <- csstc[["Time 1"]]
  r <- csstc[["Replicate Count"]]
  if ( length(s) != length(t) ){
    stop("Sample Group and Time 1 of unequal length in ",csstc$name)
  }  
  if ( length(s) != length(t) ){
    stop("Sample Group and Replicate Count of unequal length in ",csstc$name)
  }
  if (length(which(is.null(s))>0)){
    stop("NULL Sample Group in ",csstc$name)
  }
  if (length(which(is.null(t))>0)){
    stop("NULL Time 1 in ",csstc$name)
  }
  if (length(which(is.null(r))>0)){
    stop("NULL Replicate Count in ",csstc$name)
  }  
}

save(CSSs.tc,file="CSSs.tc.RData")

  


