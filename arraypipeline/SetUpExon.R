
#
# Identify files and parameters for curated exon 
#
# Can run in a tmp folder, e.g as
#  R --no-save < $AA/arraypipeline/SetUpExon.R
 
library(rjson)
library(httpRequest)

util.dir <- file.path(Sys.getenv("AA"),"utils")
aux.dir <- file.path(Sys.getenv("AA"),"auxfiles")

source(paste(util.dir,"httpget.R",sep="/"))
source(paste(util.dir,"utilitiesSampleGroup.R",sep="/"))
source(paste(util.dir,"utilitiesMetaData.R",sep="/"))

## This is the table we will build up
collezion <- character()

##
## PAM3 
##
zerotconds <- c(
                "/net/arrays/Affymetrix/core/probe_data/200702/20070220_01_KK1.0.4hrs.1.CEL",
                "/net/arrays/Affymetrix/core/probe_data/200702/20070228_01_KK2.0.4hrs.1.CEL",
                "/net/arrays/Affymetrix/core/probe_data/200702/20070228_05_NN.0.4hrs.1.CEL"
                )
res <- celMat("PAM3",c(0,60,240,360,720),zerotconds=zerotconds,term.array="/sampleData/microarray/chips/Mouse Exon")
collezion <- rbind(collezion,res)

### res <- celMat("PAM3",c(0,60,240,360,720),zerotconds=NULL,sex="Male",term.array="/sampleData/microarray/chips/Mouse Exon")

collezion <- rbind(collezion,res)


##
## PolyIC 
##
zerotconds <- c(
                "/net/arrays/Affymetrix/core/probe_data/200702/20070220_01_KK1.0.4hrs.1.CEL",
                "/net/arrays/Affymetrix/core/probe_data/200702/20070228_01_KK2.0.4hrs.1.CEL",
                "/net/arrays/Affymetrix/core/probe_data/200702/20070228_05_NN.0.4hrs.1.CEL"
                )
res <- celMat("Poly IC",c(0,240,720,1440),zerotconds=zerotconds,term.array="/sampleData/microarray/chips/Mouse Exon")
collezion <- rbind(collezion,res)
##res <- celMat("Poly IC",c(0,240,720,1440),zerotconds=NULL,sex="Male",term.array="/sampleData/microarray/chips/Mouse Exon")
##collezion <- rbind(collezion,res)

##
## PolyIC  MyD88-TRIF 

## Not sure what the unstimulated cel file is here

zerotconds <- "/net/arrays/Affymetrix/core/probe_data/200709/20070919_08_MyDTRIF_________BMDM_M_AY_DZ_Std.CEL"
res <- celMat("Poly IC",c(0,240,720,1440),strain="MyD88-Trif null",sex="Male",zerotconds=zerotconds,term.array="/sampleData/microarray/chips/Mouse Exon")

collezion <- rbind(collezion,res)

##
## Some bad cases covered by looking for NStd in cel name
##
baddies <- grep("NStd",collezion[,3])
collezion <- collezion[-baddies,]

#
# Output
#
aux.dir <- file.path(Sys.getenv("AA"),"auxfiles")
ofile <- paste(aux.dir,"ExonArrayMasterFile.tsv",sep="/")
write.table(collezion,file=ofile,quote=FALSE,col.name=FALSE,row.name=FALSE,sep="\t")


