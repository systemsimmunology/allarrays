
##
## Look at overall signal intensity and ratios
##
## Different versions of normalized data sets

util.dir <- file.path(Sys.getenv("AA"),"utils")
#pdata.dir <- file.path(Sys.getenv("AA"),"processed_data/20121002") 
data.dir <- file.path(Sys.getenv("AA"),"data") 
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))
source(paste(util.dir,"utilitiesSigTest.R",sep="/"))

## curated dataset, includes VERA SAM statistical testing
dir3prime.cur <- paste(data.dir,"20120926.curated.3prime",sep="/")
## all 3 prime arrays
dir3prime.all <- paste(data.dir,"20121001.3prime",sep="/")
## exon lps timecourse
direxon.cur <- paste(data.dir,"20121001.curated.exon",sep="/")
## all exon arrays 
direxon.all <- paste(data.dir,"20120928.exon",sep="/")
## chipseq uses: 3prime.cur and exon.cur
## tfinf uses: 3prime.cur
## allarrays uses: 3prime.all, exon.all

##
## Load VERA SAM processed
##
selectiveLoad("lps.lambdas",paste(dir3prime.cur,"all.lambdas.objects.RData",sep="/"))
selectiveLoad("lps.mus",paste(dir3prime.cur,"all.mus.objects.RData",sep="/"))
selectiveLoad("lps.ratios",paste(dir3prime.cur,"all.ratios.objects.RData",sep="/"))
lambda.cutoff <- 26.44526

##
## Load non-stat tested
##
load(paste(dir3prime.all,"dm.RData",sep="/"))
dm.3prime <- dm
load(paste(dir3prime.all,"CSSs.tc.RData",sep="/"))
CSSs.tc.3prime <- CSSs.tc

load(paste(direxon.all,"dm.RData",sep="/"))
dm.exon <- dm
load(paste(direxon.all,"CSSs.tc.RData",sep="/"))
CSSs.tc.exon <- CSSs.tc

load(paste(direxon.cur,"dm.RData",sep="/"))
dm.exon <- dm
load(paste(direxon.cur,"CSSs.tc.RData",sep="/"))
CSSs.tc.exon <- CSSs.tc
##Warning: order dependence above

onboth.ps <- intersect(rownames(dm.exon),rownames(lps.mus))
## works for exon curated (LPS only)
rats.exon <- log10(dm.exon[,2:ncol(dm.exon)]/dm.exon[,1])

##
## First compare levels in both the curated sets
## 

##
## Compare absolute values
##
quartz()

plot(log10(lps.mus[onboth.ps,"min0"]),log10(dm.exon[onboth.ps,"BMDM_Bl6_LPS_0000___Female"]),xlab="3 prime",ylab="exon",main="Log10(Intensity), T=0")

# exon abs values 2-3 times higher than three prime
quantile(probs=seq(0.1,0.9,0.1),x=lps.mus[onboth.ps,"min0"])
quantile(probs=seq(0.1,0.9,0.1),x=dm.exon[onboth.ps,"BMDM_Bl6_LPS_0000___Female"])
quantile(probs=seq(0.1,0.9,0.1),x=dm.exon[onboth.ps,"BMDM_Bl6_LPS_0000___Female"])/quantile(probs=seq(0.1,0.9,0.1),x=lps.mus[onboth.ps,"min0"])

abline(0,1)
abline(log10(2.953436),1) # %50 quantile

# Can also try to look at distributions with kernel smoothing
library(KernSmooth)
x <- cbind(log10(lps.mus[onboth.ps,"min0"]),log10(dm.exon[onboth.ps,"BMDM_Bl6_LPS_0000___Female"]))
est <- bkde2D( x, bandwidth=c(0.2,0.2) ) # larger bandwithds <-> wider kernel
contour(est$x1,est$x2,est$fhat,xlab="3 prime",ylab="exon",main="Log10(Intensity), T=0")
abline(log10(2.953436),1)
## two peak structure is somewhat baffling. slope close to kernel overall "axis"
## hmm maybe upper peak is the "true" expressed ones we care about ?


##
##  Same kind of analysis, at T=2 hrs this time
##
quartz()
plot(lps.mus[onboth.ps,"min120"],dm.exon[onboth.ps,"BMDM_Bl6_LPS_0120___Female"])
plot(log10(lps.mus[onboth.ps,"min120"]),log10(dm.exon[onboth.ps,"BMDM_Bl6_LPS_0120___Female"]))
# exon abs values 2.6-2.8 higher than three prime
quantile(probs=seq(0.25,0.75,0.25),x=lps.mus[onboth.ps,"min120"])
quantile(probs=seq(0.25,0.75,0.25),x=dm.exon[onboth.ps,"BMDM_Bl6_LPS_0120___Female"])
quantile(probs=seq(0.25,0.75,0.25),x=dm.exon[onboth.ps,"BMDM_Bl6_LPS_0120___Female"])/quantile(probs=seq(0.25,0.75,0.25),x=lps.mus[onboth.ps,"min120"])

abline(0,1)
abline(log10(2.838369),1)
x <- cbind(log10(lps.mus[onboth.ps,"min120"]),log10(dm.exon[onboth.ps,"BMDM_Bl6_LPS_0120___Female"]))
est <- bkde2D( x, bandwidth=c(0.2,0.2) ) # larger bandwithds <-> wider kernel
contour(est$x1,est$x2,est$fhat,xlab="3 prime",ylab="exon",main="Log10(Intensity), T=120")
abline(log10(2.838369),1)
#
# Ratios
#

quartz()
plot(lps.ratios[onboth.ps,"min120"],rats.exon[onboth.ps,"BMDM_Bl6_LPS_0120___Female"])
abline(0,1)
x <- cbind(lps.ratios[onboth.ps,"min120"],rats.exon[onboth.ps,"BMDM_Bl6_LPS_0120___Female"])
est <- bkde2D( x, bandwidth=c(0.2,0.2) ) # larger bandwithds <-> wider kernel
contour(est$x1,est$x2,est$fhat,xlab="3 prime",ylab="exon",main="Log10 Ratios, T=120")
## interesting that this is mostly circular, though induced genes are on diagonal
## probably has to do with how few induced genes there are relative to whole

plot(lps.ratios[onboth.ps,"min60"],rats.exon[onboth.ps,"BMDM_Bl6_LPS_0060___Female"])
abline(0,1)

plot(lps.ratios[onboth.ps,"hr4"],rats.exon[onboth.ps,"BMDM_Bl6_LPS_0240___Female"],xlab="3 prime",ylab="exon",main="Log10 Ratios, T=240")
abline(0,1)
x <- cbind(lps.ratios[onboth.ps,"hr4"],rats.exon[onboth.ps,"BMDM_Bl6_LPS_0240___Female"])
est <- bkde2D( x, gridsize = c(1000L, 1000L),bandwidth=c(0.2,0.2) ) # larger bandwithds <-> wider kernel
contour(est$x1,est$x2,est$fhat,xlab="3 prime",ylab="exon",main="Log10 Ratios, T=240",nlevels=50)
# weird, elongation only happens when gridsize upped from default of 51L to 500L or more 
## And what is this L thing?

v3 <- quantile(x=lps.ratios[onboth.ps,"hr4"],probs=seq(0.1,0.9,0.1))
ve <- quantile(x=rats.exon[onboth.ps,"BMDM_Bl6_LPS_0240___Female"],probs=seq(0.1,0.9,0.1))
plot(v3,ve)
abline(0,1)
## Fairly consistent
## For very low ratios ve deciles are somewhat higher


## At the very highest percentiles, we do start to see some differences
v3 <- quantile(x=lps.ratios[onboth.ps,"hr4"],probs=seq(0.01,0.99,0.01))
ve <- quantile(x=rats.exon[onboth.ps,"BMDM_Bl6_LPS_0240___Female"],probs=seq(0.01,0.99,0.01))
plot(v3,ve)
abline(0,1)

## Here's the full QQ-plot
qqplot(lps.ratios[onboth.ps,"hr4"],rats.exon[onboth.ps,"BMDM_Bl6_LPS_0240___Female"])

##
## Lambda cutoff in terms of ratios
##

quartz()
plot(lps.ratios[onboth.ps,"min120"],lps.lambdas[onboth.ps,"min120"])
abline(h=lambda.cutoff)

## ratio of nearest probeset to cutoff
log10ratcut <- abs(lps.ratios[which.min(abs(lps.lambdas[,"min120"]-lambda.cutoff)),"min120"])
## 0.1738
ratcut <- 10^log10ratcut
## 1.492107
## This probeset has a mu value of 436.9811


## Not sure about the analysis below

## Summary
## Use log-ratio cutoff of 1.492107 ( or round off to 1.5? ) corresponding to lambda.cutoff 
## "mu.cutoff" should be about 
##quantile(probs=seq(0.1,0.9,0.1),x=lps.mus[onboth.ps,"min0"])

##      10%       20%       30%       40%       50%       60%       70%       80% 
#  9.77417  13.29508  18.25626  26.57390  41.94745  70.95654 122.94988 218.43454 
#      90% 
#454.80523 
# quantile(probs=seq(0.1,0.9,0.1),x=dm.exon[onboth.ps,"BMDM_Bl6_LPS_0000___Female"])
#       10%        20%        30%        40%        50%        60%        70% 
#  20.03423   32.44529   47.60614   70.63696  123.88910  223.81326  365.21303 
#       80%        90% 
# 579.60069 1104.71831 
#quantile(probs=seq(0.1,0.9,0.1),x=dm.exon[onboth.ps,"BMDM_Bl6_LPS_0000___Female"])/quantile(probs=seq(0.1,0.9,0.1),x=lps.mus[onboth.ps,"min0"])
#     10%      20%      30%      40%      50%      60%      70%      80% 
#2.049711 2.440398 2.607661 2.658133 2.953436 3.154230 2.970422 2.653430 
#     90%
#2.428992 
# Use one of the above as scaling factors. Eg. mu.cutoff of 100 for 3prime would be 297 for exon! 


#--------------------------------------------
# Let's see if this changes the correspondence between constitutive counts on the two platforms

mu.cutoff <- 200 
mu.high.cutoff <- mu.cutoff
high.expressors <- names(which(apply(lps.mus[,1:imax]>mu.high.cutoff,1,sum)==imax))
constitutive.3prime.ps <- setdiff(high.expressors,lps.6hr.ps)
constitutive.3prime.eid <- as.character(ncbiID[constitutive.3prime.ps])

## binarization, with max time six hours
abs.cutoff <- 2.953436 * mu.cutoff
rat.cutoff <- 1.492107

dm.lps.exon <- dm.exon

## nt=5 corresponds to 6 hrs
nt <- 5
maxabs <- apply(dm.lps.exon[,1:nt],1,max)
abs.logvec <- ( maxabs > abs.cutoff )
tcs <- dm.lps.exon
ratmat <- (tcs/tcs[,1])[,2:nt]
maxrats <- apply(ratmat,1,max)
minrats <- apply(ratmat,1,min)
rat.logvec <- ( maxrats > rat.cutoff ) | ( minrats < (1/rat.cutoff) )
diffexp.exon.ps <- names(which(abs.logvec & rat.logvec))

mu.exon.cutoff <- 2.970422 * mu.cutoff #
# quite a bit higher than it was
mu.high.cutoff <- mu.exon.cutoff # was 300
high.expressors <- names(which(apply(dm.lps.exon[,1:nt]>mu.high.cutoff,1,sum)==nt))
constitutive.exon.ps <- setdiff(high.expressors,diffexp.exon.ps)

#compareSets(constitutive.exon.ps,constitutive.3prime.ps)
#   a    b  a^b  a-b  b-a  a+b 
# 619  932  277  342  655 1274 


#--------------------------------------------
# Older code here below
#--------------------------------------------

load(paste(pdata.dir,"dm.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.exon.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.3prime.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.RData",sep="/"))


source("~/macrophage/AffyArray/newAffy/newFunctions.R") ## for compareSets

eids.on.exon.array <- rownames(dm.exon)
eids.on.3prime.array <- rownames(dm.3prime)
eids.on.both <- intersect(eids.on.exon.array,eids.on.3prime.array)
compareSets(eids.on.exon.array,eids.on.3prime.array)
##    a     b   a^b   a-b   b-a   a+b 
## 24124 21509 16368  7756  5141 29265 


coi <- "BMDM_Bl6_LPS__Female" 

coi <- "BMDM_Bl6_salmonella.wt__Female"

set.3prime <- names(which(mm.3prime[eids,coi]==1))
set.exon <- names(which(mm.exon[eids,coi]==1))

compareSets(set.3prime,set.exon)

v <- setdiff(set.3prime,set.exon)

v <- setdiff(set.exon,set.3prime)
  
i <- 6
par(mfrow=c(2,1))
plotCSS(v[i],CSSs.tc.exon[[coi]],data.matrix=dm.exon)
plotCSS(v[i],CSSs.tc.3prime[[coi]],data.matrix=dm.3prime)
## In 3' not exon: seems due to basal level
## Rhoq has higher intesity in 3'
## In exon, not 3'
## Steap4: Subtle abs cuttoff dependence
## Basal level
## 3' "misses"
## Given that probes are "better" in exon arrays
## Perhaps a rule that favors exons in those cases?
## Sp110 is poster case for this

##
## Main factors in this comparison
## - Subtleties of ration or abs threshold 
## - t=0 levels
## - bad probes for 3' ???




