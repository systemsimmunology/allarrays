library(rjson)
library(httpRequest)
library(chron)

util.dir <- file.path(Sys.getenv("AA"),"utils")
source(file.path(Sys.getenv("HOME"),"bin/R/MatrixPrintFormat.R"))
source(paste(util.dir,"httpget.R",sep="/"))
source(paste(util.dir,"utilitiesMetaData.R",sep="/"))
source(paste(util.dir,"utilitiesSampleGroup.R",sep="/"))

timefilter <- FALSE 
debug <- TRUE

if ( timefilter ){
  ##the period from March 1, 2009 - August 30, 2009 (Not sure about this!)
  ## format m/d/y
  startDate <- "3/1/09"
  endDate <- "8/30/09"
}

term.array="/sampleData/microarray/chips/Mouse Exon"

##term.array="/sampleData/microarray/chips/Mouse 430 2.0"

##term.array="/sampleData/microarray/chips/Isbimm"

##term.array="/sampleData/microarray/chips/Mouse Promoter 1.0R"

### First version: Agnostic to Project
### Get all arrays ?

### Find all unique combinations of
### Cell Type

### Strain

### Ignore (in this round) Sex, Stim 2, 

## 1. Find all unique combinations

##term.ct="Dendritic Cell"
##term.ct="BMDM"
##names(term.ct)="Cell Type"
##term.str="Bl6"
##names(term.str)="Strain"
##term.stim="R848"
##names(term.stim)="Stimulus 1"
##term.t="\"0120\""
##names(term.t)="Time 1" 


names(term.array) = "chip"
termlist <- list()
termlist[[1]] <- term.array

obj.1 <- searchByNameValue(termListHttpForm(termlist))
##length(obj.1) is 538 on 5-14-09 for exon arrays
##length(obj.1) is 600 on 7-24-09 for exon arrays
##length(obj.1) is 648 on 8-30-09 for exon arrays
##length(obj.1) is 672 on 9-30-09 for exon arrays

if ( timefilter ){
  ## Filter by time if desired
  obj <- filterArrayListByDate(obj.1,startDate,endDate)
  cat(strsplit(term.array,split="/")[[1]][5],"\n")
  cat("Start date:",startDate,"End date:",endDate,"\n")
  cat("Before time range filtering:",length(obj.1),"\n")
  cat("After time range filtering:",length(obj),"\n")  
} else {  
  ## or not
  obj <- obj.1
  cat("Total arrays by specified array type:",length(obj),"\n")  
}

all.cell.types <- unique(unlist(lapply(obj,"[[","Cell Type")))
cat("All Cell Types:",paste(all.cell.types,collapse=", "),"\n")  
all.times <- sort(unique(unlist(lapply(obj,"[[","Time 1"))))
all.times <- all.times[sort(as.numeric(all.times),index.return=TRUE)$ix]
## last one as the string sort can fail


## Cell type for all samples
## 62 samples with no cell type
all.cell.type <- as.character(lapply(obj,"[[","Cell Type"))
cell.type.null <- which(as.character(lapply(obj,"[[","Cell Type"))=="NULL")
## unlist(lapply(obj[cell.type.null],"[[","name")) ## names of the arrays with no cell type
cat("Arrays with no Cell Type:",length(cell.type.null),"\n")
## Lack ofo metadata, including Investigator
## Seems Thunder data is here!
## The others might be VL's titration experiments

## Time 1 for all samples
## 84 experiments have no Time 1 
all.time1 <- as.character(lapply(obj,"[[","Time 1"))
time1.null <- which(as.character(lapply(obj,"[[","Time 1"))=="NULL")

## Stime 1 for all samples 
all.stim1 <- as.character(lapply(obj,"[[","Stimulus 1"))

baddies <- which(all.cell.type =="NULL" | all.stim1 == "NULL" | all.cell.type == "NULL" )
## 105 are missing one of the above attributes
cat("Arrays missing either Cell Type, Strain, Stimulus 1:",length(baddies),"\n")

## 43 with cell type but not Time1, Stim1
## length(which(all.cell.type!="NULL" & all.stim1=="NULL"))
## These all seem like "unstim", for a variety of strains
## They do not have "Stim 1==Unstim" 

cat("Arrays with Cell Type, Strain, Stimulus 1:",length(which(all.cell.type !="NULL" & all.stim1 != "NULL" & all.cell.type != "NULL")),"\n")

## 567 have Strain,Stimulus 1, Time 1
## 515, excluding double stims
## 52 doubles by this count
inds.dbls <- which(as.character(lapply(obj,"[[","Stimulus 2")) !=  "NULL" )
 
doobles <- cbind(
                 unlist(lapply(obj[inds.dbls],"[[","Cell Type")),
                 unlist(lapply(obj[inds.dbls],"[[","Strain")),
                 unlist(lapply(obj[inds.dbls],"[[","Stimulus 1")),
                 unlist(lapply(obj[inds.dbls],"[[","Time 1")),
                 unlist(lapply(obj[inds.dbls],"[[","Stimulus 2")),
                 unlist(lapply(obj[inds.dbls],"[[","Time 2"))
                 )
## Various ways to organize these
## Need to account for Stim1/Stim2 labeling arbitrarity
## In most cases, time1=time2 (sampling time)
## In some case time != time2 (administration time )


## Or April 2009, three prime arrys
##all.cell.types <- "BMDM"
##all.times <- c("0000","0020","0040","0060","0080","0100","0120","0180","0240","0300","0360","0420","0480","0720","1080","1440","2880")

## Or, May 2009, exon arrays
##all.cell.types <- c("BMDM","Dendritic Cell","Lung","MycBMM","mTEC","Splenic Dendritic Cell","Lung Dendritic Cell")
##all.times <- c("0000","0060","0120","0180","0240","0360","0480","0540","0720","1020","1080","1440")

metadataMat <- character()
timetallyMat <- numeric()

for ( cell.type in all.cell.types ){

  term.ct <-  cell.type
  names(term.ct) <- "Cell Type"
  termlist[[2]] <- term.ct

  objj <- searchByNameValue(termListHttpForm(termlist))
  ## Can include samples with No Stimulus 1 , Time 1
  
  ## Filter by time if desired
  if ( timefilter ){
    objj <- filterArrayListByDate(objj,startDate,endDate)  
  }
  
  ## Hopefully this filter can be ditched some day
  ## Addama is return matches by substring for Cell Type
  baddies <- which(unlist(lapply(objj,"[[","Cell Type")) != cell.type)
  if ( length(baddies) > 0 ){
    objj <- objj[-baddies]
  }
  
  if ( debug ){
    cat("***",cell.type,"***\n")
    b2 <- unlist(lapply(obj.1[which(all.cell.type ==cell.type & all.stim1 != "NULL" & all.cell.type != "NULL")],"[[","name"))
    cat("length from big query, filtered on cell type",length(b2),"\n")  
    cat("length before filtering",length(objj),"\n")
  }
  
  ## For clean indexing operations below, remove the nulls
  time1.null <- which(as.character(lapply(objj,"[[","Time 1"))=="NULL")
  if ( length(time1.null) > 0){
    objj <- objj[-time1.null]
  }
  stim1.null <- which(as.character(lapply(objj,"[[","Stimulus 1"))=="NULL")
  if ( length(stim1.null) > 0){
    objj <- objj[-stim1.null]
  }
  strain.null <- which(as.character(lapply(objj,"[[","Strain"))=="NULL")
  if ( length(strain.null) > 0){
    objj <- objj[-strain.null]
  }
  if ( debug ) { cat("after removing nulls",length(objj),"\n") }

  ## Ditch samples with Stim2 for now
  objj <- removeStim2(objj)
  if ( debug ) { cat("after removing doubles",length(objj),"\n") }

  stim1s <- unique(unlist(lapply(objj,"[[","Stimulus 1")))
  for ( stim1 in stim1s ){
    inds <- which(unlist(lapply(objj,"[[","Stimulus 1"))==stim1)
    objjj <- objj[inds] ## now has unique cell type and stim
    strains <- unique(unlist(lapply(objjj,"[[","Strain")))
    for ( strain in strains ){
      inds <- which(unlist(lapply(objjj,"[[","Strain"))==strain)
      timetally.1 <- table(unlist(lapply(objjj[inds],"[[","Time 1")))
      timetally.2 <- as.numeric(timetally.1)
      names(timetally.2) <- names(timetally.1)
      timetally.3 <- vector(length=length(all.times))
      names(timetally.3) <- all.times
      timetally.3[names(timetally.2)] <- as.numeric(timetally.2)
      if ( !is.null(ncol(timetallyMat))){
        if ( ncol(timetallyMat) != length(timetally.3) ){
          stop("Mismatch\n")
        }
      }
      timetallyMat <- rbind(timetallyMat,timetally.3)
      metadataMat <- rbind(metadataMat,c(cell.type,stim1,strain))
    }
  }
}
rownames(timetallyMat) <- NULL
 
tt2 <- cbind(timetallyMat,apply(timetallyMat,1,sum))
tt3 <- rbind(tt2,apply(tt2,2,sum))

mm3 <- rbind(metadataMat,c("Column Total","",""))

outMat <- cbind(mm3,tt3)
outMat <- rbind(c("Cell Type","Stimulus 1","Strain",all.times,"Row Total"),outMat)
rownames(outMat) <- NULL
colnames(outMat) <- NULL
write.table(outMat,file="AllMetaData.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

cat("Have Cell Type,Strain, Stimulus 1:",sum(timetallyMat),"\n")

outMat <- rbind(c("Cell Type","Strain","Stimulus 1","Time 1","Stimulus 2","Time 2"),doobles)
write.table(outMat,file="Doubles.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

##write.table(outMat,file="AllThreePrimeMetaData.tsv",sep="\t",quote=FALSE)


