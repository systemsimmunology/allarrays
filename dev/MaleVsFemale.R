
BMDM_IRF1 null_Poly IC__Female

BMDM_IRF1 null_Poly IC__Male

girl <- "BMDM_IRF1-null_Poly-IC_0720___Female"

boy <- "BMDM_IRF1-null_Poly-IC_0720___Male"
 
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

isMale <- function ( inString ) {
  !identical(inString,strsplit(inString,split="Male")[[1]])
}


isFemale <- function ( inString ) {
  !identical(inString,strsplit(inString,split="Female")[[1]])
}


## grab non-gender part of column name
preMale <- function (inString) {
  substr(inString,1,nchar(inString)-4)
}


## grab non-gender part of column name
preFemale <- function (inString) {
  substr(inString,1,nchar(inString)-6)
}

exon.conds <- colnames(dm.exon)

b <- names(which(sapply(exon.conds,isMale)))
bb <- as.character(sapply(b,preMale))


c <- names(which(sapply(exon.conds,isFemale)))
cc <- as.character(sapply(c,preFemale))


dd <- intersect(bb,cc)
## wow, 13 of these

ddm <- paste(dd,"Male",sep="")

ddf <- paste(dd,"Female",sep="")


collection <- character()
for ( i in 1:13) {
  eids <- names(sort(dm.exon[,ddf[i]]-dm.exon[,ddm[i]],decreasing=TRUE)[1:10])
  cat(gene.symbol[eids],"\n")
  collection <- c(collection,eids)
  ##plot(dm.exon[,ddm[i]],dm.exon[,ddf[i]])
}
gene.symbol[names(sort(table(collection),decreasing=TRUE)[1:10])]



collection <- character()
for ( i in 1:13) {
  eids <- names(sort(dm.exon[,ddm[i]]-dm.exon[,ddf[i]],decreasing=TRUE)[1:10])
  cat(gene.symbol[eids],"\n")
  collection <- c(collection,eids)
  ##plot(dm.exon[,ddm[i]],dm.exon[,ddf[i]])
}
gene.symbol[names(sort(table(collection),decreasing=TRUE)[1:10])]



gene.symbol[names(sort(dm.exon[,ddm[2]]-dm.exon[,ddf[2]],decreasing=TRUE)[1:10])]


gene.symbol[names(sort(dm.exon[,ddm[4]]-dm.exon[,ddf[4]],decreasing=TRUE)[1:10])]


##
## Find "male" and "female"  (top ten)
##

## "Male" genes

gene.symbol[names(sort(table(collection),decreasing=TRUE)[1:10])]
   69550    56040    12257    12861    17069    20525    22121    12314 
  "Bst2"  "Rplp1"   "Tspo" "Cox6a1"   "Ly6e" "Slc2a1" "Rpl13a"  "Calm2" 
   15254   193740 
 "Hint1" "Hspa1a" 

eid <- gene.eid["Bst2"]
dm.exon[eid,ddm]-dm.exon[eid,ddf]
dm.exon[eid,ddm]/dm.exon[eid,ddf] ## Ratios, though almost always positive are below 2


## "Female genes"
  16956    56619   627695   213742    22352    11957    12874    16414 
   "Lpl" "Clec4e" "627695"   "Xist"    "Vim"  "Atp5j"    "Cpd"  "Itgb2" 
   18792    20302 
  "Plau"   "Ccl3" 

 
eid <- gene.eid["Vim"]
dm.exon[eid,ddf]-dm.exon[eid,ddm]
dm.exon[eid,ddf]/dm.exon[eid,ddm] ## Ratios, though almost always positive are below 2






### Here's another way to score

# female-ness

femininity <- apply(dm.exon[,ddf]/dm.exon[,ddm],1,mean)

save(femininity,file="femininity.RData")

g1 <- names(sort(apply(dm.exon[,ddf]/dm.exon[,ddm],1,mean),decreasing=TRUE)[1:10])
 
ind <- 3
cbind(dm.exon[g1[ind],ddf],dm.exon[g1[ind],ddm])

### Female top ten

gene.symbol[g1]
     213742      319178       14562      627853       22042      319172 
     "Xist" "Hist1h2bb"      "Gdf3"  "EG627853"      "Tfrc" "Hist1h2ab" 
     666094       94176       21973      667280 
 "EG666094"     "Dock2"     "Top2a"  "EG667280" 

##
## Male ness
## 

masculinity <- apply(dm.exon[,ddm]/dm.exon[,ddf],1,mean)

save(masculinity,file="masculinity.RData")

b1 <- names(sort(apply(dm.exon[,ddm]/dm.exon[,ddf],1,mean),decreasing=TRUE)[1:10])
 
ind <- 3
cbind(dm.exon[b1[ind],ddf],dm.exon[b1[ind],ddm])

## Male top ten

## First one only induced in PolyIC and Leismania

gene.symbol[b1]
               26908                26900                22290 
           "Eif2s3y"              "Ddx3y"                "Uty" 
              213002            100039803               677027 
            "Ifitm6"       "LOC100039803" "OTTMUSG00000010155" 
           100039720                66101            100042429 
         "100039720"               "Ppih"       "LOC100042429" 
               69902 
             "Mrto4"

