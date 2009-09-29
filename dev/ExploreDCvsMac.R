

##########################################
## BMDM vs DC signatures 
#########################################
v1 <- unlist(lapply(CSSs.tc.exon,"[[","Cell Type"))
v2 <- unlist(lapply(CSSs.tc.3prime,"[[","Cell Type"))
v <- c(v1,v2)

bmdm.conds <- names(which(v=="BMDM"))
dc.conds <- names(which(v=="Dendritic Cell"))

## ?Attempt simple probabilistic classifier
w <- (length(bmdm.conds) - apply(mm[,bmdm.conds],1,sum)) + apply(mm[,dc.conds],1,sum)
w <- w/(length(bmdm.conds)+length(dc.conds))
names(sort(w,decreasing=TRUE)[1:10])

sum(mm[eid,dc.conds])/length(dc.conds)
sum(mm[eid,bmdm.conds])/length(bmdm.conds)

gridPlotCSS(eid, CSSs.tc.exon[dc.conds], data.matrix = dm.exon, nx=5, ny=5)
 
bex  <- intersect(bmdm.conds,names(CSSs.tc.exon))
x11()
gridPlotCSS(eid, CSSs.tc.exon[bex], data.matrix = dm.exon, nx=8, ny=7, ylim =c(0,max(dm.exon[eid,bex])))

b3p <- intersect(bmdm.conds,names(CSSs.tc.3prime))
x11()
gridPlotCSS(eid, CSSs.tc.3prime[b3p], data.matrix = dm.3prime, nx=5, ny=7, ylim=c(0,max(dm.3prime[eid,bex])))

dex  <- intersect(dc.conds,names(CSSs.tc.exon))
x11()
gridPlotCSS(eid, CSSs.tc.exon[dex], data.matrix = dm.exon, nx=5, ny=5, ylim =c(0,max(dm.exon[eid,bex])))


### Which DC stims were also done for BMDM


for ( cond in dc.conds ){
  query <-  CSSs.tc[[cond]][c("Strain","Stimulus 1","Stimulus 2")]
  v <- inListSoft(query, CSSs.tc.exon[bex])
  if ( length(v) > 0 ) {
    cat("***",cond, "***\n")
    cat(names(CSSs.tc.exon[bex[v]]),"\n")
  }
}

for ( cond in dc.conds ){
  query <-  CSSs.tc[[cond]][c("Strain","Stimulus 1","Stimulus 2")]
  v <- inListSoft(query, CSSs.tc.3prime[b3p])
  if ( length(v) > 0 ) {
    cat("***",cond, "***\n")
    cat(names(CSSs.tc.3prime[b3p[v]]),"\n")
  }
}

