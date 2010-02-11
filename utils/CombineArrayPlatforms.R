##
## Load and combine 3'prime arrays and exon arrays into single summary matrix and CSSs
##
##

indir.3prime <- "20090507.3prime"
indir.exon <- "20091015.exon"
outdir <- "20091015"

util.dir <- file.path(Sys.getenv("AA"),"utils")
ncbi.dir <- file.path(Sys.getenv("DATA_DIR"),"ncbi")
data.dir <- file.path(Sys.getenv("AA"),"data") # local repository for "raw" data
pdata.dir <- file.path(Sys.getenv("AA"),"processed_data") 
## Augment the above to full path
indir.exon <- paste(data.dir,indir.exon,sep="/")
indir.3prime <- paste(data.dir,indir.3prime,sep="/")
outdir <- paste(pdata.dir,outdir,sep="/")

load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.eid.RData",sep="/"))

source(paste(util.dir,"utilitiesPlot.R",sep="/"))
source(paste(util.dir,"utilitiesSigTest.R",sep="/"))
source(paste(util.dir,"utilitiesMetaData.R",sep="/"))

data.path.exon <- indir.exon
data.path.3prime <- indir.3prime

##Use this to read off network directory
## Have currently moved all input files to local
##data.path.exon <- paste(Sys.getenv("DATA_DIR"),"ImmunoRepository/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-05-02_at_09.20.05",sep="/")
##data.path.3prime <- paste(Sys.getenv("DATA_DIR"),"ImmunoRepository/sampleData/microarray/runs/Aderem Three Prime Arrays_2009-05-07_at_00.06.21/output/Mouse 430 2.0",sep="/")

## Year 1 release of portal
## data.path.3prime <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/ExpressionSets/GeneLevel/Genomics Expression Public Dataset"

### 
### Load expression data (saved locally as network transfer is slow )
###
 
load(paste(indir.3prime,"dm.3prime.allarrays.RData",sep="/"))
dm.3prime <- dm; rm(dm) 

##load("data/dm.3prime.released.RData")
##dm.3prime <- ag; rm(ag) # called ag in older files

load(paste(indir.exon,"dm.RData",sep="/"))
dm.exon <- dm

### 
### Load expression column mappings
###
load(paste(data.path.3prime,"dm.columns.RData",sep="/"))
dm.columns.3prime <- dm.columns; rm(dm.columns)
load(paste(data.path.exon,"dm.columns.RData",sep="/"))
dm.columns.exon <- dm.columns; rm(dm.columns)

##
## Load exon timecourses
##

load(paste(data.path.exon,"CSSs.tc.RData",sep="/"))
## 5-2-09 Exon Set: 82 timecourses
## 8-22-09 Exon Set: 98 timecourses
## 10-15-09 Exon Set: 101 timecoures
## 1-25-10. That last number now appears to be 98, but differs from 8-22-09
##lung.conds <- names(which(unlist(lapply(CSSs.tc,"[[","Cell Type"))=="Lung"))
##CSSs.tc <- CSSs.tc[setdiff(names(CSSs.tc),lung.conds)]
## removed 3
##mtec.conds <- names(which(unlist(lapply(CSSs.tc,"[[","Cell Type"))=="mTEC"))
##CSSs.tc <- CSSs.tc[setdiff(names(CSSs.tc),mtec.conds)]
## removed 3
so <- sort(names(CSSs.tc),index.return=TRUE)$ix
CSSs.tc <- CSSs.tc[so]
CSSs.tc.exon <- CSSs.tc ; rm(CSSs.tc)
 
##
## Load 3 prime timecourses
##
load(paste(data.path.3prime,"CSSs.tc.RData",sep="/"))
so <- sort(names(CSSs.tc),index.return=TRUE)$ix
CSSs.tc <- CSSs.tc[so]
CSSs.tc.3prime <- CSSs.tc ; rm(CSSs.tc)
## Have 60 with 2009-05-07
## Some time courses are really short?
maxtimes <- unlist(lapply(lapply(CSSs.tc.3prime,"[[","Time 1"),max))
CSSs.tc.3prime <- CSSs.tc.3prime[setdiff(names(CSSs.tc.3prime),names(which(maxtimes < 120)))]
## 54 tcs remain (2009-05-07)
## No exon arrays have maxtimes < 120 

## Can be convenient to get rid of "unneeded" columns in data matrices i.e. those that are not in the current set of timecourses
needed.cols.3prime <- unique(unlist(lapply(CSSs.tc.3prime,"[[","DM Column")))
baddies <- setdiff(needed.cols.3prime,colnames(dm.3prime))
if ( length(baddies) > 0 ){
  cat("Error. Bad Time points: ",baddies,"\n")
} else {
  dm.3prime <- dm.3prime[,needed.cols.3prime]
}

needed.cols.exon <- unique(unlist(lapply(CSSs.tc.exon,"[[","DM Column")))
baddies <- setdiff(needed.cols.exon,colnames(dm.exon))
if ( length(baddies) > 0 ){
  cat("Error. Bad Time points: ",baddies,"\n")
} else {
  dm.exon <- dm.exon[,needed.cols.exon]
}

## Unique
exon.only <- setdiff(names(CSSs.tc.exon),names(CSSs.tc.3prime))
threeprime.only <- setdiff(names(CSSs.tc.3prime),names(CSSs.tc.exon))
CSSs.tc <- c(CSSs.tc.exon[exon.only],CSSs.tc.3prime[threeprime.only]) ## (First definition of this )
### c *can* be unpredictable for lists


## In common
inc <- intersect(names(CSSs.tc.exon),names(CSSs.tc.3prime))
### Make choice[s for which arrays to get repsentative values
### Attempt to combine all data into a single set
## There are eleven CSSs for both array types, five using the public data
## Choose, based e.g. on sample number
### See below for diagnostic tools
arraychoice <- character()
arraychoice["BMDM_Bl6_AC-LDL__Female"] <- "3prime" ## Both at 24 hrs, but most LDLs were done for 3' (?)
arraychoice["BMDM_Bl6_Ifn beta__Female"] <- "3prime" ## more time points. Tmax is 24hrs for exon, though.
arraychoice["BMDM_Bl6_Ifn gamma__Female"] <- "exon" ## Tough call, see below
arraychoice["BMDM_Bl6_LPS__Female"] <- "3prime" ## Longer coverage
arraychoice["BMDM_Bl6_PAM3_Poly IC_Female"] <- "3prime" ## more coverage
arraychoice["BMDM_Bl6_PAM3__Female"] <- "3prime" ## Longer
arraychoice["BMDM_Bl6_Poly IC__Female"] <- "3prime" ## Longer
arraychoice["BMDM_Bl6_salmonella.wt__Female"] <- "exon" ## longer TC
arraychoice["BMDM_Cebp/d null_LPS__Female"] <- "exon" ## Longer tc
arraychoice["BMDM_Cebp/d null_PAM3__Female"] <- "exon" ## Longer tc
arraychoice["BMDM_Trif null_PAM3_Poly IC_Female"] <- "exon" ## Not sure here (both are at 12hrs only)
for ( co in inc ) {
  if ( arraychoice[co] == "exon" ){
    CSSs.tc[co] <- CSSs.tc.exon[co]
  } else if (  arraychoice[co] == "3prime" )
    CSSs.tc[co] <- CSSs.tc.3prime[co]
}
so <- sort(names(CSSs.tc),index.return=TRUE)$ix
CSSs.tc <- CSSs.tc[so]
## Length 119 in June 2009
## Length 131 in September 2009
## Length 141 in September 2009

## threeprime.set: Those conditions defined by the three prime set
threeprime.conds <- sort(c(threeprime.only,names(which(arraychoice=="3prime"))))
## exon.set: Those conditions defined by the exon set
exon.conds <- sort(c(exon.only,names(which(arraychoice=="exon"))))

## REDEFINE, REMOVING REDUNDANCIES
## Remember to re-read objects from disk if you don't like this
CSSs.tc.3prime <- CSSs.tc.3prime[threeprime.conds]
CSSs.tc.exon <- CSSs.tc.exon[exon.conds]

save(CSSs.tc,file=paste(outdir,"CSSs.tc.RData",sep="/"))
save(CSSs.tc.3prime,file=paste(outdir,"CSSs.tc.3prime.RData",sep="/"))
save(CSSs.tc.exon,file=paste(outdir,"CSSs.tc.exon.RData",sep="/"))
save(dm.exon,file=paste(outdir,"dm.exon.RData",sep="/"))
save(dm.3prime,file=paste(outdir,"dm.3prime.RData",sep="/"))

