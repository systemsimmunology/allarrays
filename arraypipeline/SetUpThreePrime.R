
#
# Identify files and parameters for curated 3 prime
#

library(rjson)
library(httpRequest)

util.dir <- file.path(Sys.getenv("AA"),"utils")
aux.dir <- file.path(Sys.getenv("AA"),"utils")

source(paste(util.dir,"httpget.R",sep="/"))
source(paste(util.dir,"utilitiesSampleGroup.R",sep="/"))
source(paste(util.dir,"utilitiesMetaData.R",sep="/"))

## This is the table we will build up
collezion <- character()
##
## LPS
##
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200406/20040603_03_LPS1-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200406/20040621_01_LPS2-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200406/20040622_01_LPS3-0.CEL"
           )
## leave 12hrs out for a bit
res <- celMat("LPS",c(0,20,40,60,80,120,240,360,480,1080,1440),zerot)
##/net/arrays/Affymetrix/core/probe_data/200410/20041020_04_LPS_10ng-ul_24h.CEL needs to be removed from output file
collezion <- rbind(collezion,res)

##
## R848
## 
##ls /net/arrays/Affymetrix/core/probe_data/200410/*CEL | grep R848 | grep "\-0"
zerot <-c(
          "/net/arrays/Affymetrix/core/probe_data/200410/20041011_01_R848-A-0.CEL",
          "/net/arrays/Affymetrix/core/probe_data/200410/20041012_01_R848-B-0.CEL"
          )
res <- celMat("R848",c(0,20,40,60,80,120),zerot)
collezion <- rbind(collezion,res)

##
## PAM2
## 
##ls /net/arrays/Affymetrix/core/probe_data/200407/*CEL | grep PAM2 | grep "\-0"
zerot <-c(
          "/net/arrays/Affymetrix/core/probe_data/200407/20040706_01_PAM2A-0.CEL",
          "/net/arrays/Affymetrix/core/probe_data/200407/20040707_01_PAM2B-0.CEL",
          "/net/arrays/Affymetrix/core/probe_data/200407/20040708_01_PAM2C-0.CEL"
          )
res <- celMat("PAM2",c(0,20,40,60,80,120),zerot,sex="Male")
collezion <- rbind(collezion,res)

##
## PAM3
## 
##ls /net/arrays/Affymetrix/core/probe_data/200407/*CEL | grep PAM3 | grep "\-0"
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200407/20040712_01_PAM3A-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200407/20040713_01_PAM3B-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200407/20040715_01_PAM3C-0.CEL"
           )
res <- celMat("PAM3",c(0,20,40,60,80,120,240,480,720),zerot,sex="Male")
collezion <- rbind(collezion,res)

##
## Poly IC
##
## ls /net/arrays/Affymetrix/core/probe_data/200410/*CEL | grep PIC | grep "\-0"
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200410/20041004_01_PIC-A-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200410/20041005_01_PIC-B-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200410/20041005_05_PIC-C-0.CEL"
           )
res <- celMat("Poly IC",c(0,20,40,60,80,120,240,480,720),zerot)
collezion <- rbind(collezion,res)

##
## CpG
##
##ls /net/arrays/Affymetrix/core/probe_data/200503/*CEL | grep CPG | grep "\-0"
#ls /net/arrays/Affymetrix/core/probe_data/200502/*CEL | grep CPG | grep "\-0"
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200503/20050317_05_CPG_C-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200502/20050224_07_CPG_A-0.CEL"
           )
res <- celMat("CpG",c(0,20,40,60,80,120),zerot)
collezion <- rbind(collezion,res)

#
# MyD88
#
##ls /net/arrays/Affymetrix/core/probe_data/200502/*.CEL | grep MyD88 | grep "\-0"
zerot <- c(
           "/net/arrays/Affymetrix/core/probe_data/200502/20050208_01_MyD88_LPS_M2-0.CEL",
           "/net/arrays/Affymetrix/core/probe_data/200502/20050224_08_MyD88_PIC_C-0.CEL"
           )
res <- celMat(stim="LPS",c(0,20,60,120),zerot,strain="Myd88 null",label="MyD88-LPS")
collezion <- rbind(collezion,res)
res <- celMat(stim="Poly IC",c(0,20,60,120),zerot,strain="Myd88 null",label="MyD88-PolyIC")
collezion <- rbind(collezion,res)
res <- celMat(stim="PAM3",c(0,20,60,120),zerot[1],strain="Myd88 null",label="MyD88-PAM3")
collezion <- rbind(collezion,res)


##
## TRIF
##
res <- celMat(stim="LPS",c(0,20,60,120),zerotconds=NULL,strain="TRIF null (LPS2)",label="TRIF-LPS")
collezion <- rbind(collezion,res)
## These need weeding by hand
res <- celMat(stim="PAM2",c(0,20,60,120),zerotconds=NULL,strain="TRIF null (LPS2)",label="TRIF-PAM2")
collezion <- rbind(collezion,res)


##
## All bad cases covered by looking for stim amount in cel name
##
baddies <- grep("ng",collezion[,3])
collezion <- collezion[-baddies,]

#
# Output
#
aux.dir <- file.path(Sys.getenv("AA"),"utils")
ofile <- paste(aux.dir,"ThreePrimeMasterFile.tsv",sep="/"))
write.table(collezion,file=ofile,quote=FALSE,col.name=FALSE,row.name=FALSE,sep="\t")

