
#
# Identify files and parameters for curated 3 prime
#

library(rjson)
library(httpRequest)
source("/Users/thorsson/allarrays/utils/utilitiesSampleGroup.R")
source("/Users/thorsson/allarrays/utils/utilitiesMetaData.R")
source("/Users/thorsson/allarrays/utils/httpget.R")

basename <- function(path){
  toks <- split1(path,splitchar="/")
  basename <- toks[length(toks)]
  basename
}

## utility function for creating cel file matrix )
celMat <- function( stim, times, zerotconds, sex="Female" ){
  tc <- list()
  tc[["Cell Type"]] <- "BMDM"
  tc[["Strain"]] <- "Bl6"
  tc[["Stimulus 1"]] <- stim
  tc[["Sex"]] <- sex
  tc[["Time 1"]] <- times
  res <- celFilesFromCSSTCs(tc)
  if ( stim == "Poly IC" ){
    stim <- "Poly-IC"
  }
  zstring <- paste(c("BMDM_Bl6_",stim,"_0000___",sex),collapse="")
  res <- rbind(cbind(zstring,zerotconds),res)
  colnames(res) <- NULL
  rb <- as.character(sapply(res[,2],basename))
  collect <- cbind(stim,res[,1],rb,res[,2])
  colnames(collect) <- NULL
  collect
}


collezion <- character()
##ls /net/arrays/Affymetrix/core/probe_data/200410/*CEL | grep R848 | grep "\-0"
zerot <-c(
          "/net/arrays/Affymetrix/core/probe_data/200410/20041011_01_R848-A-0.CEL",
          "/net/arrays/Affymetrix/core/probe_data/200410/20041012_01_R848-B-0.CEL"
          )
res <- celMat("R848",c(0,20,40,60,80,120),zerot)
collezion <- rbind(collezion,res)

##ls /net/arrays/Affymetrix/core/probe_data/200407/*CEL | grep PAM2 | grep "\-0"
zerot <-c(
          "/net/arrays/Affymetrix/core/probe_data/200407/20040706_01_PAM2A-0.CEL",
          "/net/arrays/Affymetrix/core/probe_data/200407/20040707_01_PAM2B-0.CEL",
          "/net/arrays/Affymetrix/core/probe_data/200407/20040708_01_PAM2C-0.CEL"
          )
res <- celMat("PAM2",c(0,20,40,60,80,120),zerot,sex="Male")
collezion <- rbind(collezion,res)

##ls /net/arrays/Affymetrix/core/probe_data/200407/*CEL | grep PAM3 | grep "\-0"
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200407/20040712_01_PAM3A-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200407/20040713_01_PAM3B-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200407/20040715_01_PAM3C-0.CEL"
           )
res <- celMat("PAM3",c(0,20,40,60,80,120,240,480,720),zerot,sex="Male")
collezion <- rbind(collezion,res)

## ls /net/arrays/Affymetrix/core/probe_data/200410/*CEL | grep PIC | grep "\-0"
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200410/20041004_01_PIC-A-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200410/20041005_01_PIC-B-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200410/20041005_05_PIC-C-0.CEL"
           )
res <- celMat("Poly IC",c(0,20,40,60,80,120,240,480,720),zerot)
collezion <- rbind(collezion,res)

##ls /net/arrays/Affymetrix/core/probe_data/200503/*CEL | grep CPG | grep "\-0"
#ls /net/arrays/Affymetrix/core/probe_data/200502/*CEL | grep CPG | grep "\-0"
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200503/20050317_05_CPG_C-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200502/20050224_07_CPG_A-0.CEL"
           )
res <- celMat("CpG",c(0,20,40,60,80,120),zerot)
collezion <- rbind(collezion,res)

write.table(collezion,file="/Users/thorsson/allarrays/auxfiles/try.tsv",quote=FALSE,col.name=FALSE,row.name=FALSE,sep="\t")

if (FALSE) {

tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "Myd88 null"
tc[["Stimulus 1"]] <- "LPS"
tc[["Time 1"]] <- c(0,20,60,120)
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/MyD88LPSthreePrimeCelFiles_unix")
tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "Myd88 null"
tc[["Stimulus 1"]] <- "PAM3"
tc[["Time 1"]] <- c(0,20,60,120)
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/MyD88PAM3threePrimeCelFiles_unix")

tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "Myd88 null"
tc[["Stimulus 1"]] <- "Poly IC"
tc[["Time 1"]] <- c(0,20,60,120)
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/MyD88PolyICthreePrimeCelFiles_unix")
##ls /net/arrays/Affymetrix/core/probe_data/200502/*CEL | grep MyD88 | grep "\-0"

## TRIF KO
tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "TRIF null (LPS2)"
tc[["Stimulus 1"]] <- "LPS"
tc[["Time 1"]] <- c(0,20,60,120)
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/TRIFLPSthreePrimeCelFiles_unix")
tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "TRIF null (LPS2)"
tc[["Stimulus 1"]] <- "PAM2"
tc[["Time 1"]] <- c(0,20,60,120)
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/TRIFPAM2threePrimeCelFiles_unix")

}
