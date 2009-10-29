
########################
##  Salmonella
########################

## Oddity: Why is there both a "BMDM_Bl6_salmonella.flgBpEm62__Female" and "BMDM_Bl6_salmonella.flgBpEm62__Male"

## This is a fix that no longer seems required
## CSSs.tc.exon[["BMDM_MyD88-Trif null_Listeria__Female"]][["Sample Group"]] <- CSSs.tc.exon[["BMDM_MyD88-Trif null_Listeria__Female"]][["Sample Group"]][2:4]
## CSSs.tc.exon[["BMDM_MyD88-Trif null_Listeria__Female"]][["Time 1"]] <- CSSs.tc.exon[["BMDM_MyD88-Trif null_Listeria__Female"]][["Time 1"]][2:4]
## CSSs.tc.exon[["BMDM_MyD88-Trif null_Listeria__Female"]][["DM Column"]] <- CSSs.tc.exon[["BMDM_MyD88-Trif null_Listeria__Female"]][["DM Column"]][2:4]

## Why is there a "MyD88-Trif-null_Listeria_in-vitro_0000___BMDM_Mouse" sample group??

stim1s.exon <- unlist(lapply(CSSs.tc.exon,"[[","Stimulus 1"))
stim1s.3prime <- unlist(lapply(CSSs.tc.3prime,"[[","Stimulus 1"))

salmonella.strings <- c("salmonella.wt","FlgEProfect","FliCProfect","salmonella.flgB","salmonella.prgH","salmonella.flgBpEm62")
listeria.strings <- c("Listeria")
 
sc.exon <- names(CSSs.tc.exon)[ stim1s.exon %in% salmonella.strings ]
sc.3prime <- names(CSSs.tc.3prime)[ stim1s.3prime %in% salmonella.strings ] ## keep in mind that we removed the wt.
lc.exon <- names(CSSs.tc.exon)[ stim1s.exon %in% listeria.strings ]

## We're down to several now, so let's just choose the plot order
sc.exon <- c( "BMDM_Bl6_salmonella.wt__Female","BMDM_MyD88-Trif null_salmonella.wt__Female","BMDM_Bl6_FlgEProfect__Female" , "BMDM_Bl6_FliCProfect__Female" )
sc.3prime <- c("BMDM_Bl6_salmonella.wt__Female","BMDM_Bl6_salmonella.flgB__Female","BMDM_Bl6_salmonella.prgH__Female","BMDM_Bl6_salmonella.flgBpEm62__Female")


#######################################
## Try to identify salmonella classifiers
########################################

salmonella.conds <- c(sc.exon,sc.3prime)
non.salmonella.conds <- setdiff(colnames(mm),salmonella.conds) ## would be better just to remove IPAF from mm

w <- (length(non.salmonella.conds) - apply(mm[,non.salmonella.conds],1,sum)) + apply(mm[,salmonella.conds],1,sum)
##w <- w/(length(non.salmonella.conds)+length(salmonella.conds))
names(sort(w,decreasing=TRUE)[1:10])

sum(mm[v,salmonella.conds])/length(salmonella.conds)

sum(mm[v,non.salmonella.conds])/length(non.salmonella.conds)

apply(mm[v,salmonella.conds],1,sum)/length(salmonella.conds)

apply(mm[v,non.salmonella.conds],1,sum)/length(non.salmonella.conds)

ww <- apply(mm[v,salmonella.conds],1,sum)/length(salmonella.conds) + (1 - apply(mm[v,non.salmonella.conds],1,sum)/length(non.salmonella.conds) )
ww <- ww/2.


##################################################
##
## Wild-type + MyD88-TRIF mouse, Salmonella, Listeria, and LPS
##  
#################################################

cois <- c("BMDM_Bl6_salmonella.wt__Female","BMDM_Bl6_Listeria__Female","BMDM_Bl6_LPS__Female","BMDM_MyD88-Trif null_salmonella.wt__Female")

eids <- rownames(mm) ## include all, or replace with filtered 

## This is for the combined array mm, which can be imperfect
## as the choice made for LPS is currently 3 prime
sa.set <- names(which(mm[eids,"BMDM_Bl6_salmonella.wt__Female"]==1))
li.set <- names(which(mm[eids,"BMDM_Bl6_Listeria__Female"]==1))
lps.set <- names(which(mm[eids,"BMDM_Bl6_LPS__Female"]==1))
bug.set <- intersect(sa.set,li.set)
bugonly <- setdiff(bug.set,lps.set)


## Separation by platform
sa.3prime <- names(which(mm.3prime[eids,"BMDM_Bl6_salmonella.wt__Female"]==1))
sa.exon <- names(which(mm.exon[eids,"BMDM_Bl6_salmonella.wt__Female"]==1))
lps.3prime <- names(which(mm.3prime[eids,"BMDM_Bl6_LPS__Female"]==1))
lps.exon <- names(which(mm.exon[eids,"BMDM_Bl6_LPS__Female"]==1))

##
##
load("/Volumes/ILYA LAB/Vesteinn/data/ncbi/gene.eid.RData")

eid <- gene.eid["Pdgfb"] ##

eid <- gene.eid["Il17ra"] ##


## Handy plotting routine for both platforms 
bugPlot <- function ( eid ){
  dev.set(2)
  main <- paste(gene.symbol[eid],': Exon Array',sep='')
  gridPlotCSS(eid,CSSs.tc.exon[c(sc.exon,lc.exon)],data.matrix=dm.exon,nx=3,ny=2, main=main)
  dev.set(3)
  main <- paste(gene.symbol[eid],': 3 prime Array',sep='')
  gridPlotCSS(eid,CSSs.tc.3prime[sc.3prime],data.matrix=dm.3prime,nx=2,ny=2,main=main)
}


## Handy plotting routine for both platforms 
bugPlot <- function ( eid ){
  dev.set(2)
  main <- paste(gene.symbol[eid],': Exon Array',sep='')
  gridPlotCSS(eid,CSSs.tc.exon[c(sc.exon,lc.exon)],data.matrix=dm.exon,nx=3,ny=2, main=main)
  dev.set(3)
  main <- paste(gene.symbol[eid],': 3 prime Array',sep='')
  gridPlotCSS(eid,CSSs.tc.3prime[sc.3prime],data.matrix=dm.3prime,nx=2,ny=2,main=main)
}



## Handy plotting routine for both platforms 
bugtlrPlot <- function ( eid ){
  dev.set(2)
  main <- paste(gene.symbol[eid],': Exon Array',sep='')
  gridPlotCSS(eid,CSSs.tc.exon[intersect(cois,names(CSSs.tc.exon)) ],data.matrix=dm.exon,nx=2,ny=2, main=main)
  dev.set(3)
  main <- paste(gene.symbol[eid],': 3 prime Array',sep='')
  gridPlotCSS(eid,CSSs.tc.3prime[intersect(cois,names(CSSs.tc.3prime))],data.matrix=dm.3prime,nx=1,ny=2,main=main)
}


## A more continuous version, based on maximum ratios

source("../utils/utilitiesSigTest.R")
mr.exon <- maxRatioMultipleCSSTs(eids,CSSs.tc.exon,data.matrix=dm.exon)
mr.3prime <- maxRatioMultipleCSSTs(eids,CSSs.tc.3prime,data.matrix=dm.3prime)

## Having moved to a choice, when two arrays are possible, some of these
## are not longer available. Can reinstate if needed

## IF using full info from both arrays 

sa.3prime.vec <- mr.3prime[eids,"BMDM_Bl6_salmonella.wt__Female"]
sa.exon.vec <- mr.exon[eids,"BMDM_Bl6_salmonella.wt__Female"]
sa.vec <- (sa.3prime.vec + sa.exon.vec)/2.

lps.3prime.vec <- mr.3prime[eids,"BMDM_Bl6_LPS__Female"]
lps.exon.vec <- mr.exon[eids,"BMDM_Bl6_LPS__Female"]
lps.vec <- (lps.3prime.vec + lps.exon.vec)/2.

## ELSE, if a choice has been made
sa.exon.vec <- mr.exon[eids,"BMDM_Bl6_salmonella.wt__Female"]
sa.vec <- sa.exon.vec

lps.3prime.vec <- mr.3prime[eids,"BMDM_Bl6_LPS__Female"]
lps.vec <- lps.3prime.vec

li.exon.vec <- mr.exon[eids,"BMDM_Bl6_Listeria__Female"]

plot(sa.vec,lps.vec,xlim=c(0,10),ylim=c(0,10))
text(sa.vec,lps.vec,labels=gene.symbol[eids])

bugoriented <- names(which(sa.exon.vec>1.5 & sa.3prime.vec>1.5 & abs(lps.exon.vec)<1.5 & abs(lps.3prime.vec)<1.5 ))

bugoriented <- names(which(sa.exon.vec> 2.5 & abs(lps.exon.vec)<2.0 ))



x11()
cond1 <- "BMDM_Bl6_salmonella.wt__Female"
cond2 <- "BMDM_Bl6_salmonella.flgB__Female"

cond1 <- "BMDM_Bl6_salmonella.flgB__Female"
cond2 <- "BMDM_Bl6_salmonella.prgH__Female"

cond1 <- "BMDM_Bl6_FlgEProfect__Female"
cond2 <- "BMDM_Bl6_FliCProfect__Female"

cond1 <- "BMDM_Bl6_salmonella.wt__Female"
cond2 <- "BMDM_MyD88-Trif null_salmonella.wt__Female"
   
dev.set(4)
maxx <- 40
mr.mat <- mr.exon
plot(mr.mat[eids,cond1],mr.mat[eids,cond2],xlab=cond1,ylab=cond2,xlim=c(0,maxx),ylim=c(0,maxx))
text(mr.mat[eids,cond1],mr.mat[eids,cond2],labels=gs[eids],pos=4)
abline(0,1)

## GO testing
library(GO.db)
library(GOstats)
library(mogene10stv1mmentrezg.db)
source("/Users/thorsson/data/GeneOntology/goUtilities.R")
load("/Users/thorsson/data/GeneOntology/mogene10stv1mmentrezg.RData")

es <- intersect(eids,mogene10stv1mmentrezg.eids.goMapped.byps)

res <- goSigTableEid(es,entrezUniverse=mogene10stv1mmentrezg.eids,annotation="mogene10stv1mmentrezg.db")
OR
res <- goSignificanceEid(es,entrezUniverse=mogene10stv1mmentrezg.eids,annotation="mogene10stv1mmentrezg.db")

v <- res$goMap[[1]]

## Longer description of gene
unlist(mget(v,org.Mm.egGENENAME))

## Yet another gene symbol mapping
mget(v,org.Mm.egSYMBOL)
mget("Tlr4",org.Mm.egSYMBOL2EG)

##
##
##
## Are there genes induced in the Myd88/Trif doulbe ko macrophages by salmonella infection?
## How do these genes act in WT macrophages?
##

sa.myd88trif.exon.vec <- mr.exon[eids,"BMDM_MyD88-Trif null_salmonella.wt__Female"]
sa.nontlr <- names(which(sa.myd88trif.exon.vec > 4))
sa.nontlr <- names(sort(sa.myd88trif.exon.vec[sa.nontlr],decreasing=TRUE))


## for lps, we essentially have the 3prime data

lps.myd88.exon.vec <- mr.3prime[eids,"BMDM_Myd88 null_LPS__Female"]
lps.trif.exon.vec <- mr.3prime[eids,"BMDM_Trif null (LPS2)_LPS__Female"]

 
x11()

eid <- gene.eid["Ccl7"]

## Plots for TRIF dependence
exon.trif.conds <- c(
"BMDM_Bl6_LPS__Female",
"BMDM_Bl6_Poly IC__Male",
"BMDM_MyD88-Trif null_Poly IC__Male",
"BMDM_Bl6_salmonella.wt__Female",
"BMDM_MyD88-Trif null_salmonella.wt__Female",
"BMDM_Bl6_Listeria__Female",
"BMDM_MyD88-Trif null_Listeria__Female"
                     )

threeprime.trif.conds <- c(
"BMDM_Bl6_LPS__Female",
"BMDM_Bl6_Poly IC__Female",
"BMDM_Trif null (LPS2)_LPS__Female",
"BMDM_Trif null_LPS__Female"
)
                           

trifPlot <- function ( eid ){
  dev.set(2)
  main <- paste(gene.symbol[eid],': Exon Array',sep='')  
  gridPlotCSS(eid,CSSs.tc.exon[exon.trif.conds],data.matrix=dm.exon,nx=4,ny=2, main=main, tmax=12.)
  dev.set(3)
  main <- paste(gene.symbol[eid],': 3 prime Array',sep='')
  gridPlotCSS(eid,CSSs.tc.3prime[threeprime.trif.conds],data.matrix=dm.3prime,nx=2,ny=2,main=main,tmax=12)
}

