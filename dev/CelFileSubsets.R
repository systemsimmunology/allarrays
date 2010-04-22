
source("/Users/thorsson/allarrays/utils/utilitiesSampleGroup.R")

tc <- CSSs.tc[["BMDM_Bl6_R848__Female"]]
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/R848threePrimeCelFiles_unix")
#ls /net/arrays/Affymetrix/core/probe_data/200410/*CEL | grep R848 | grep "\-0"

tc <- CSSs.tc[["BMDM_Bl6_PAM2__Male"]]
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/PAM2threePrimeCelFiles_unix")
#ls /net/arrays/Affymetrix/core/probe_data/200407/*CEL | grep PAM2 | grep "\-0"

tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "Bl6"
tc[["Stimulus 1"]] <- "PAM3"
tc[["Time 1"]] <- c(0,20,40,60,80,240,480,720)

##tc <- CSSs.tc[["BMDM_Bl6_PAM3__Female"]]
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/PAM3threePrimeCelFiles_unix")
#ls /net/arrays/Affymetrix/core/probe_data/200407/*CEL | grep PAM3 | grep "\-0"

tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "Bl6"
tc[["Stimulus 1"]] <- "Poly IC"
tc[["Time 1"]] <- c(0,20,40,60,80,240,480,720)

res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/PolyICthreePrimeCelFiles_unix")
## ls /net/arrays/Affymetrix/core/probe_data/200410/*CEL | grep PIC | grep "\-0"

tc <- list()
tc[["Cell Type"]] <- "BMDM"
tc[["Strain"]] <- "Bl6"
tc[["Stimulus 1"]] <- "CpG"
tc[["Time 1"]] <- c(0,20,40,60,80,240,480,720)

##tc <- CSSs.tc[["BMDM_Bl6_PAM3__Female"]]
res <- celFilesFromCSSTCs(tc,file="/Users/thorsson/allarrays/auxfiles/PolyICthreePrimeCelFiles_unix")

