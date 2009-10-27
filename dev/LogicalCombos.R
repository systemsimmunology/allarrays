

util.dir <- file.path(Sys.getenv("AA"),"utils")
pdata.dir <- file.path(Sys.getenv("AA"),"processed_data/20091015") 
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))

load(paste(pdata.dir,"mm.RData",sep="/"))
load(paste(pdata.dir,"CSSs.tc.RData",sep="/"))

 
all.strain <- unlist(lapply(CSSs.tc,"[[","Strain"))
all.stim1 <- unlist(lapply(CSSs.tc,"[[","Stimulus 1"))
all.stim2 <- unlist(lapply(CSSs.tc,"[[","Stimulus 2"))

ss <- cbind(all.strain,all.stim1,all.stim2)
colnames(ss) <- c("Strain","Stimulus 1","Stimulus2")
write.table(ss,file="lc.tsv",sep='\t',quote=FALSE,col.names=TRUE,row.names=TRUE)
