
# For 3' arrays:
# Turn VERA and SAM output files into RData
# Also non-replicated samples
# 
utils.dir <- file.path(Sys.getenv("AA"),"utils")
source(paste(utils.dir,"VERAandSAMutils.R",sep="/"))

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

rt <- read.table(file.path(Sys.getenv("AA"),"annotation/MoGene4302_Mm_ENTREZG_mapfile"),sep="\t",as.is=TRUE,header=TRUE)
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

source("countSigGenes.R")



all.ps.lists <- ls()[grep("ps.sig$",ls())]
all.ps.lists <- c(all.ps.lists,"lps.full.ps.sig.rep")
save( list=all.ps.lists, file=paste(annotation.dir,"all.ps.list.objects.RData",sep="/"))

