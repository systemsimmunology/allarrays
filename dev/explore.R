   
load("/Volumes/ILYA LAB/Vesteinn/data/ncbi/gene.symbol.RData")
load("/Volumes/ILYA LAB/Vesteinn/data/ncbi/gene.eid.RData")

source("Rshared/utilitiesPlot.R")
source("Rshared/utilitiesSigTest.R")
source("Rshared/utilitiesMetaData.R")

data.path.exon <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/sampleData/microarray/runs/AutomaticTasks - Gene-Level Exon Pipeline - All Exon Arrays_2009-05-02_at_09.20.05"

## Recent 3 prime runs
data.path.3prime <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/sampleData/microarray/runs/Aderem Three Prime Arrays_2009-05-07_at_00.06.21/output/Mouse 430 2.0"

## Year 1 release of portal
## data.path.3prime <- "/Volumes/ILYA LAB/Vesteinn/data/ImmunoRepository/ExpressionSets/GeneLevel/Genomics Expression Public Dataset"

### 
### Load expression data (saved locally as network transfer is slow )
###

load("data/dm.3prime.allarrays.RData")
dm.3prime <- dm; rm(dm) 

##load("data/dm.3prime.released.RData")
##dm.3prime <- ag; rm(ag) # called ag in older files

load("data/dm.exon.RData")
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

lung.conds <- names(which(unlist(lapply(CSSs.tc,"[[","Cell Type"))=="Lung"))
CSSs.tc <- CSSs.tc[setdiff(names(CSSs.tc),lung.conds)]
## removed 3
mtec.conds <- names(which(unlist(lapply(CSSs.tc,"[[","Cell Type"))=="mTEC"))
CSSs.tc <- CSSs.tc[setdiff(names(CSSs.tc),mtec.conds)]
## removed 3
so <- sort(names(CSSs.tc),index.return=TRUE)$ix
CSSs.tc <- CSSs.tc[so]
## Have 76
CSSs.tc.exon <- CSSs.tc ; rm(CSSs.tc)

##
## Load 3 prime timecourses
##
load(paste(data.path.3prime,"CSSs.tc.RData",sep="/"))
so <- sort(names(CSSs.tc),index.return=TRUE)$ix
CSSs.tc <- CSSs.tc[so]
CSSs.tc.3prime <- CSSs.tc ; rm(CSSs.tc)
## Have 27 with released set, 70 with 2009-05-07
## Some time courses are really short?
maxtimes <- unlist(lapply(lapply(CSSs.tc.3prime,"[[","Time 1"),max))
CSSs.tc.3prime <- CSSs.tc.3prime[setdiff(names(CSSs.tc.3prime),names(which(maxtimes < 120)))]
## 64 tcs remain (2009-05-07)
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
## Length 129 in June 2009

## threeprime.set: Those conditions defined by the three prime set
threeprime.set <- sort(c(threeprime.only,names(which(arraychoice=="3prime"))))
## exon.set: Those conditions defined by the exon set
exon.set <- sort(c(exon.only,names(which(arraychoice=="exon"))))

## REDEFINE, REMOVING REDUNDANCIES
## Remember to re-read objects from disk if you don't like this
CSSs.tc.3prime <- CSSs.tc.3prime[threeprime.set]
CSSs.tc.exon <- CSSs.tc.exon[exon.set]

###########
########## Test plots for the above
##########

  
gridPlotCSS(gene.eid["Il12b"],CSSs.tc.3prime[inc],data.matrix=dm.3prime,nx=3,ny=4)
x11()
gridPlotCSS(gene.eid["Il12b"],CSSs.tc.exon[inc],data.matrix=dm.exon,nx=3,ny=4)


for ( co in inc ){
  cat("==========", co,"============\n")
  cat("3 prime:",CSSs.tc.3prime[[co]][["Time 1"]],"\n")
  cat("Exon:",CSSs.tc.exon[[co]][["Time 1"]],"\n")  
}

  
eid <- gene.eid["Arg1"]

v <- binarizeCSSTs(eid,CSSs.tc.exon,data.matrix=dm.exon) 
gridPlotCSS(eid,CSSs.tc.exon,data.matrix=dm.exon,nx=10,ny=8,labvec=v)

eid <- gene.eid["Chi3l3"]

v <- binarizeCSSTs(eid,CSSs.tc.3prime,data.matrix=dm.3prime) 
 
gridPlotCSS(eid,CSSs.tc.3prime[so],data.matrix=dm.3prime,nx=6,ny=5,labvec)

####

both <- intersect(names(CSSs.tc.3prime),names(CSSs.tc.exon))
gridPlotCSS(eid,CSSs.tc.exon[both],data.matrix=dm.exon,nx=3,ny=2)
x11()
gridPlotCSS(eid,CSSs.tc.3prime[both],data.matrix=dm.3prime,nx=3,ny=2)







### Make choices for which arrays to get repsentative values
CSSs.tc["BMDM_Bl6_AC-LDL__Female"] <- CSSs.tc.3prime["BMDM_Bl6_AC-LDL__Female"] ## Both at 24 hrs, but most LDLs were done for 3' (?)
CSSs.tc["BMDM_Bl6_Ifn beta__Female"] <- CSSs.tc.3prime["BMDM_Bl6_Ifn beta__Female"] ## more time points. Tmax is 24hrs for exon, though.
CSSs.tc["BMDM_Bl6_Ifn gamma__Female"] <- CSSs.tc.exon["BMDM_Bl6_Ifn gamma__Female"] ## Tough call, see below
## CSSs.tc.3prime[["BMDM_Bl6_Ifn gamma__Female"]][["Time 1"]] = 0   60  120  240  480 1440
## CSSs.tc.exon[["BMDM_Bl6_Ifn gamma__Female"]][["Time 1"]] == 0  240  720 1020 1440
CSSs.tc["BMDM_Bl6_LPS__Female"] <- CSSs.tc.3prime["BMDM_Bl6_LPS__Female"] ## Longer coverage
CSSs.tc["BMDM_Bl6_PAM3_Poly IC_Female"] <- CSSs.tc.3prime["BMDM_Bl6_PAM3_Poly IC_Female"] ## more coverage
CSSs.tc["BMDM_Bl6_PAM3__Female"] <- CSSs.tc.3prime["BMDM_Bl6_PAM3__Female"] ## Longer
CSSs.tc["BMDM_Bl6_Poly IC__Female"] <- CSSs.tc.3prime["BMDM_Bl6_Poly IC__Female"] ## Longer
CSSs.tc["BMDM_Bl6_salmonella.wt__Female"] <- CSSs.tc.exon["BMDM_Bl6_salmonella.wt__Female"] ## longer TC
CSSs.tc["BMDM_Cebp/d null_LPS__Female"] <- CSSs.tc.exon["BMDM_Cebp/d null_LPS__Female"] ## Longer tc
CSSs.tc["BMDM_Cebp/d null_PAM3__Female"] <- CSSs.tc.exon["BMDM_Cebp/d null_PAM3__Female"] ## Longer tc
CSSs.tc["BMDM_Trif null_PAM3_Poly IC_Female"] <- CSSs.tc.exon["BMDM_Trif null_PAM3_Poly IC_Female"] ## Not sure here (both are at 12hrs only)
