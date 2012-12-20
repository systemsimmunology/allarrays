

source(paste(file.path(Sys.getenv("AA"),"utils"),"VERAandSAMutils.R",sep="/"))

out.dir <- Sys.getenv("PWD")
expression.dir <- out.dir
annotation.dir <- out.dir
##expression.dir <- file.path(Sys.getenv("TFINF"),"expression_data")
##annotation.dir <- file.path(Sys.getenv("TFINF"),"annotations")

##
## Annotation objects
##

rt <- read.table(file.path(Sys.getenv("AA"),"annotation/Mouse4302_Mm_ENTREZG_mapfile"),sep="\t",as.is=TRUE,header=TRUE)
probeset <- rt$ProbesetID
ncbiID <- as.character(rt$EntrezID)
names(ncbiID) <- probeset
cname.compare <- rt$GeneSymbol
names(cname.compare) <- probeset
## aliases
gene.names.all <- cname.compare
commonName <- cname.compare
gene.ids.all <- probeset

load("~/data/ncbi/np.RData")
names(np) <- paste(names(np),"_at",sep="")

annotation.objects <- c("probeset",
                        "gene.ids.all",
                        "commonName",
                        "np",
                        "ncbiID",
                        "gene.names.all",
                        "cname.compare")

# identical: commonName, gene.names.all
# identical: probeset, gene.ids.all
# keep both, as functions may rely on them

save (list=annotation.objects, file=paste(annotation.dir,"annotation.objects.RData",sep="/"))

repProbes.cname <- probeset
names(repProbes.cname) <- cname.compare
repProbes.ncbiID <- probeset
names(repProbes.ncbiID) <- ncbiID
repProbes.np <- names(np)
names(repProbes.np) <- np

save(list=c("repProbes.cname","repProbes.np","repProbes.ncbiID"),file=paste(annotation.dir,"representativeProbes.RData",sep="/"))



##
## Read and save VERA and SAM files
##

conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr6","hr8","hr18","hr24")
res <- readVERAandSAMfiles("LPS",expression.dir,niceheaders=conds)
lps.mus <- res$mus
lps.ratios <- res$ratios
lps.lambdas <- res$lambdas
lps.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
res <- readVERAandSAMfiles("PAM2",expression.dir,niceheaders=conds)
pam2.mus <- res$mus
pam2.ratios <- res$ratios
pam2.lambdas <- res$lambdas
pam2.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
#original (male - female mix)
#conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr8","hr12")
res <- readVERAandSAMfiles("PAM3",expression.dir,niceheaders=conds)
pam3.mus <- res$mus
pam3.ratios <- res$ratios
pam3.lambdas <- res$lambdas
pam3.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
#original (male - female mix)
#conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr8","hr12")
res <- readVERAandSAMfiles("R848",expression.dir,niceheaders=conds)
r848.mus <- res$mus
r848.ratios <- res$ratios
r848.lambdas <- res$lambdas
r848.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120")
res <- readVERAandSAMfiles("CpG",expression.dir,niceheaders=conds)
cpg.mus <- res$mus
cpg.ratios <- res$ratios
cpg.lambdas <- res$lambdas
cpg.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min40","min60","min80","min120","hr4","hr8","hr12")
res <- readVERAandSAMfiles("Poly-IC",expression.dir,niceheaders=conds)
polyIC.mus <- res$mus
polyIC.ratios <- res$ratios
polyIC.lambdas <- res$lambdas
polyIC.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min60","min120")
res <- readVERAandSAMfiles("MyD88-LPS",expression.dir,niceheaders=conds)
myd88KOlps.mus <- res$mus
myd88KOlps.ratios <- res$ratios
myd88KOlps.lambdas <- res$lambdas
myd88KOlps.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min60","min120")
res <- readVERAandSAMfiles("MyD88-PolyIC",expression.dir,niceheaders=conds)
myd88KOpolyIC.mus <- res$mus
myd88KOpolyIC.ratios <- res$ratios
myd88KOpolyIC.lambdas <- res$lambdas
myd88KOpolyIC.muStdErrs <- res$muStdErrs

conds <- c("min0","min20","min60","min120")
myd88KOpam3.mus <- read.matrix(paste(c(expression.dir,"/","MyD88-PAM3.mus"),collapse=""))
colnames(myd88KOpam3.mus) <- conds

conds <- c("min0","min60","min120")
trifKOpam2.mus <- read.matrix(paste(c(expression.dir,"/","TRIF-PAM2.mus"),collapse=""))
colnames(trifKOpam2.mus) <- conds

conds <- c("min0","min60","min120")
trifKOlps.mus <- read.matrix(paste(c(expression.dir,"/","TRIF-LPS.mus"),collapse=""))
colnames(trifKOlps.mus) <- conds

all.ratios.objects <- ls()[grep("ratios$",ls())]
save ( list=all.ratios.objects, file=paste(expression.dir,"all.ratios.objects.RData",sep="/"))
all.mus.objects <- ls()[grep("mus$",ls())]
save ( list=all.mus.objects, file=paste(expression.dir,"all.mus.objects.RData",sep="/"))
all.lambdas.objects <- ls()[grep("lambdas$",ls())]
save ( list=all.lambdas.objects, file=paste(expression.dir,"all.lambdas.objects.RData",sep="/"))
all.muStdErrs.objects <- ls()[grep("muStdErrs$",ls())]
save ( list=all.muStdErrs.objects, file=paste(expression.dir,"all.muStdErrs.objects.RData",sep="/"))

probeset <- rownames(lps.mus)
nps <- nrow(lps.mus)
lambdaCube <- rep(NA,6*nps*5)
dim(lambdaCube) <- c(6,nps,5)
dimnames(lambdaCube)[[3]] <- c("min20","min40","min60","min80","min120")
dimnames(lambdaCube)[[2]] <- probeset[1:nps]
dimnames(lambdaCube)[[1]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
lambdaCube["LPS",,] <- lps.lambdas[1:nps,1:5]
lambdaCube["PAM2",,] <- pam2.lambdas[1:nps,1:5]
lambdaCube["PAM3",,] <- pam3.lambdas[1:nps,1:5]
lambdaCube["PolyIC",,] <- polyIC.lambdas[1:nps,1:5]
lambdaCube["R848",,] <- r848.lambdas[1:nps,1:5]
lambdaCube["CpG",,] <- cpg.lambdas[1:nps,1:5]

ratioCube <- rep(NA,6*nps*5)
dim(ratioCube) <- c(6,nps,5)
dimnames(ratioCube)[[3]] <- c("min20","min40","min60","min80","min120")
dimnames(ratioCube)[[2]] <- probeset[1:nps]
dimnames(ratioCube)[[1]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")

ratioCube["LPS",,] <- lps.ratios[1:nps,1:5]
ratioCube["PAM2",,] <- pam2.ratios[1:nps,1:5]
ratioCube["PAM3",,] <- pam3.ratios[1:nps,1:5]
ratioCube["PolyIC",,] <- polyIC.ratios[1:nps,1:5]
ratioCube["R848",,] <- r848.ratios[1:nps,1:5]
ratioCube["CpG",,] <- cpg.ratios[1:nps,1:5]

muCube <- rep(NA,6*nps*6)
dim(muCube) <- c(6,nps,6)
dimnames(muCube)[[3]] <- c("min0","min20","min40","min60","min80","min120")
dimnames(muCube)[[2]] <- probeset[1:nps]
dimnames(muCube)[[1]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")

muCube["LPS",,] <- lps.mus[1:nps,1:6]
muCube["PAM2",,] <- pam2.mus[1:nps,1:6]
muCube["PAM3",,] <- pam3.mus[1:nps,1:6]
muCube["PolyIC",,] <- polyIC.mus[1:nps,1:6]
muCube["R848",,] <- r848.mus[1:nps,1:6]
muCube["CpG",,] <- cpg.mus[1:nps,1:6]

all.Cube.objects <- ls()[grep("Cube$",ls())]
save( list=all.Cube.objects, file=paste(expression.dir,"all.Cube.objects.RData",sep="/"))

##
## Intensities scaled to maximum
##

max.intensity <- apply(lps.mus,1,max) 
lps.mat.max1 <- lps.mus/max.intensity ## this one will be modified for nucloc

##max.intensity <- apply(pam2.mus,1,max) 
pam2.mat.max1 <- pam2.mus/max.intensity

##max.intensity <- apply(pam3.mus,1,max) 
pam3.mat.max1 <- pam3.mus/max.intensity

##max.intensity <- apply(polyIC.mus,1,max) 
polyIC.mat.max1 <- polyIC.mus/max.intensity

##max.intensity <- apply(r848.mus,1,max) 
r848.mat.max1 <- r848.mus/max.intensity

##max.intensity <- apply(cpg.mus,1,max) 
cpg.mat.max1 <- cpg.mus/max.intensity

ofile <- paste(expression.dir,"scaled.mus.objects.RData",sep="/")
savelist <- c("max.intensity","lps.mat.max1","pam2.mat.max1","pam3.mat.max1","polyIC.mat.max1","r848.mat.max1","cpg.mat.max1" )

save(list=savelist, file=ofile)


##########################################################
# Overall signal integral
##########################################################

imax <- 6
m <- lps.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
lps.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#lps.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1


imax <- 6
m <- pam2.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
pam2.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#pam2.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

imax <- 6
m <- pam3.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
pam3.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#pam3.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

imax <- 6
m <- polyIC.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
polyIC.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#polyIC.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

imax <- 6
m <- r848.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
r848.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#r848.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

imax <- 6
m <- cpg.mus[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,20,40,60,80,120)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
cpg.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
#CpG.mic <- intvec/((tvec[length(tvec)]-tvec[1])*m[,1]) - 1

all.mic <- cbind(lps.mic,pam2.mic,pam3.mic,polyIC.mic,r848.mic,cpg.mic)
colnames(all.mic) <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")

save(all.mic,file=paste(expression.dir,"all.mic.RData",sep="/"))

##
## Significance lists
##

mu.cutoff <- 100
##mu.cutoff <- 0
##lambda.cutoff <- 57.2 ## earlier cutoff - leads to 3372 genes for full time-course
#lambda.cutoff <- 150

## May 2010
##lambda.cutoff <- 26.61275 ## 0.05% cutoff - leads to 4913 genes for full time-course
##lambda.cutoff <- 66.31579 ## 0.01% cutoff - leads to 3069 genes for full time-course

## September 2012
lambda.cutoff <- 26.44526

# Up to 24hrs
imax <- 11
lps.full.ps.sig <- rownames(sigSlice(lambda.cutoff,lps.mus[,1:imax],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[,1:imax]<mu.cutoff,1,sum)==imax))
lps.full.ps.sig <- setdiff(lps.full.ps.sig,low.expressors)
cat("No. LPS genes up to 24 hrs: ",length(lps.full.ps.sig),"\n")

## Up to 2hrs


imax <- 6

stimulus <- "LPS"
lps.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,lps.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
lps.ps.sig <- setdiff(lps.ps.sig,low.expressors)
lps.cname.sig <- sort(unique(cname.compare[lps.ps.sig]))

stimulus <- "PAM2"
pam2.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,pam2.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
pam2.ps.sig <- setdiff(pam2.ps.sig,low.expressors)
pam2.cname.sig <- sort(unique(cname.compare[pam2.ps.sig]))

stimulus <- "PAM3"
pam3.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,pam3.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
pam3.ps.sig <- setdiff(pam3.ps.sig,low.expressors)
pam3.cname.sig <- sort(unique(cname.compare[pam3.ps.sig]))

stimulus <- "PolyIC"
polyIC.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,polyIC.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
polyIC.ps.sig <- setdiff(polyIC.ps.sig,low.expressors)
polyIC.cname.sig <- sort(unique(cname.compare[polyIC.ps.sig]))

stimulus <- "R848"
r848.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,r848.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
r848.ps.sig <- setdiff(r848.ps.sig,low.expressors)
r848.cname.sig <- sort(unique(cname.compare[r848.ps.sig]))

stimulus <- "CpG"
cpg.ps.sig <- rownames(sigSlice(lambda.cutoff,ratioCube[stimulus,,1:(imax-1)],lambdaCube[stimulus,,1:(imax-1)]))
low.expressors <- names(which(apply(muCube[stimulus,cpg.ps.sig,1:imax]<mu.cutoff,1,sum)==imax))
cpg.ps.sig <- setdiff(cpg.ps.sig,low.expressors)
cpg.cname.sig <- sort(unique(cname.compare[cpg.ps.sig]))

all.ps.sig <- union(lps.ps.sig,union(pam2.ps.sig,union(pam3.ps.sig,union(polyIC.ps.sig,union(r848.ps.sig,cpg.ps.sig)))))
all.cname.sig <- sort(unique(cname.compare[all.ps.sig]))
n.all.ps.sig <- length(all.ps.sig)

cat("Gene counts up to and including 2 hrs\n")
cat("No. LPS genes: ",length(lps.ps.sig),"\n")
cat("No. PAM2 genes: ",length(pam2.ps.sig),"\n")
cat("No. PAM3 genes: ",length(pam3.ps.sig),"\n")
cat("No. PolyIC genes: ",length(polyIC.ps.sig),"\n")
cat("No. R848 genes: ",length(r848.ps.sig),"\n")
cat("No. CpG genes: ",length(cpg.ps.sig),"\n")
cat("No. All genes: ",length(all.ps.sig),"\n")

sigCompCube <- rep(NA,6*6*6)
dim(sigCompCube) <- c(6,6,6)
dimnames(sigCompCube)[[3]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
dimnames(sigCompCube)[[2]] <- c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
dimnames(sigCompCube)[[1]] <- c("a","b","a^b","a-b","b-a","a+b")

sigCompCube[,"LPS","PAM2"] <- compareSets(lps.cname.sig,pam2.cname.sig)
sigCompCube[,"LPS","PAM3"] <- compareSets(lps.cname.sig,pam3.cname.sig)
sigCompCube[,"LPS","PolyIC"] <- compareSets(lps.cname.sig,polyIC.cname.sig)
sigCompCube[,"LPS","R848"] <- compareSets(lps.cname.sig,r848.cname.sig)
sigCompCube[,"LPS","CpG"] <- compareSets(lps.cname.sig,cpg.cname.sig)

sigCompCube[,"PAM2","PAM3"] <- compareSets(pam2.cname.sig,pam3.cname.sig)
sigCompCube[,"PAM2","PolyIC"] <- compareSets(pam2.cname.sig,polyIC.cname.sig)
sigCompCube[,"PAM2","R848"] <- compareSets(pam2.cname.sig,r848.cname.sig)
sigCompCube[,"PAM2","CpG"] <- compareSets(pam2.cname.sig,cpg.cname.sig)

sigCompCube[,"PAM3","PolyIC"] <- compareSets(pam3.cname.sig,polyIC.cname.sig)
sigCompCube[,"PAM3","R848"] <- compareSets(pam3.cname.sig,r848.cname.sig)
sigCompCube[,"PAM3","CpG"] <- compareSets(pam3.cname.sig,cpg.cname.sig)

sigCompCube[,"PolyIC","R848"] <- compareSets(polyIC.cname.sig,r848.cname.sig)
sigCompCube[,"PolyIC","CpG"] <- compareSets(polyIC.cname.sig,cpg.cname.sig)

sigCompCube[,"R848","CpG"] <- compareSets(r848.cname.sig,cpg.cname.sig)

sigCompCube["a^b",,]/sigCompCube["a+b",,]


outMat <- matrix(NA,nrow=6,ncol=6)
rownames(outMat) <-  c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
colnames(outMat) <-  c("LPS","PAM2","PAM3","PolyIC","R848","CpG")
m1 <- sigCompCube["a^b",,]/sigCompCube["a",,]
m2 <- t(sigCompCube["a^b",,]/sigCompCube["b",,])
outMat[upper.tri(outMat)] <- m1[upper.tri(m1)]
outMat[lower.tri(outMat)] <- m2[lower.tri(m2)]

write.table(100*outMat,file="sigCompCube.tsv",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

lps.ps.sig <- rownames(sigSlice(lambda.cutoff,lps.mus[,1:imax],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[,1:imax]<mu.cutoff,1,sum)==imax))
lps.ps.sig <- setdiff(lps.ps.sig,low.expressors)
#seems identical to computation above
lps.full.ps.sig.rep <- sort(repProbes.cname[unique(sort(cname.compare[lps.full.ps.sig]))])


all.ps.lists <- ls()[grep("ps.sig$",ls())]
all.ps.lists <- c(all.ps.lists,"lps.full.ps.sig.rep")
save( list=all.ps.lists, file=paste(annotation.dir,"all.ps.list.objects.RData",sep="/"))

