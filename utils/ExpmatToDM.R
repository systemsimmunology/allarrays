dm <- read.matrix("expmat.tsv")
dm <- 2^dm
save(dm,file="dm.RData")
