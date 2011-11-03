
dev.set(2)  
gridPlotCSS(gene.eid["Il17rb"],CSSs.tc.3prime,data.matrix=dm.3prime,nx=7,ny=7)
##x11()
dev.set(3)
gridPlotCSS(gene.eid["Il23r"],CSSs.tc.exon,data.matrix=dm.exon,nx=10,ny=10)

###########
########## Test plots for the above
##########

  
gridPlotCSS(gene.eid["Il12b"],CSSs.tc.3prime[inc],data.matrix=dm.3prime,nx=3,ny=4)
x11()
gridPlotCSS(gene.eid["Il12b"],CSSs.tc.exon[inc],data.matrix=dm.exon,nx=3,ny=4)


for ( co in inc ){
  cat("==========", co,"============\n")
  cat("3 prime:",CSSs.tc.3prime[[co]][["Time 1"]],"\n")
  cat("Exon:",CSSs.tc.exon[[co]][["Time 1"]],"\n")  
}

  
eid <- gene.eid["Arg1"]

v <- binarizeCSSTs(eid,CSSs.tc.exon,data.matrix=dm.exon) 
gridPlotCSS(eid,CSSs.tc.exon,data.matrix=dm.exon,nx=10,ny=8,labvec=v)

eid <- gene.eid["Chi3l3"]

v <- binarizeCSSTs(eid,CSSs.tc.3prime,data.matrix=dm.3prime) 
 
gridPlotCSS(eid,CSSs.tc.3prime[so],data.matrix=dm.3prime,nx=6,ny=5,labvec)

####

both <- intersect(names(CSSs.tc.3prime),names(CSSs.tc.exon))
gridPlotCSS(eid,CSSs.tc.exon[both],data.matrix=dm.exon,nx=3,ny=2)
x11()
gridPlotCSS(eid,CSSs.tc.3prime[both],data.matrix=dm.3prime,nx=3,ny=2)



### Make choices for which arrays to get repsentative values
CSSs.tc["BMDM_Bl6_AC-LDL__Female"] <- CSSs.tc.3prime["BMDM_Bl6_AC-LDL__Female"] ## Both at 24 hrs, but most LDLs were done for 3' (?)
CSSs.tc["BMDM_Bl6_Ifn beta__Female"] <- CSSs.tc.3prime["BMDM_Bl6_Ifn beta__Female"] ## more time points. Tmax is 24hrs for exon, though.
CSSs.tc["BMDM_Bl6_Ifn gamma__Female"] <- CSSs.tc.exon["BMDM_Bl6_Ifn gamma__Female"] ## Tough call, see below
## CSSs.tc.3prime[["BMDM_Bl6_Ifn gamma__Female"]][["Time 1"]] = 0   60  120  240  480 1440
## CSSs.tc.exon[["BMDM_Bl6_Ifn gamma__Female"]][["Time 1"]] == 0  240  720 1020 1440
CSSs.tc["BMDM_Bl6_LPS__Female"] <- CSSs.tc.3prime["BMDM_Bl6_LPS__Female"] ## Longer coverage
CSSs.tc["BMDM_Bl6_PAM3_Poly IC_Female"] <- CSSs.tc.3prime["BMDM_Bl6_PAM3_Poly IC_Female"] ## more coverage
CSSs.tc["BMDM_Bl6_PAM3__Female"] <- CSSs.tc.3prime["BMDM_Bl6_PAM3__Female"] ## Longer
CSSs.tc["BMDM_Bl6_Poly IC__Female"] <- CSSs.tc.3prime["BMDM_Bl6_Poly IC__Female"] ## Longer
CSSs.tc["BMDM_Bl6_salmonella.wt__Female"] <- CSSs.tc.exon["BMDM_Bl6_salmonella.wt__Female"] ## longer TC
CSSs.tc["BMDM_Cebp/d null_LPS__Female"] <- CSSs.tc.exon["BMDM_Cebp/d null_LPS__Female"] ## Longer tc
CSSs.tc["BMDM_Cebp/d null_PAM3__Female"] <- CSSs.tc.exon["BMDM_Cebp/d null_PAM3__Female"] ## Longer tc
CSSs.tc["BMDM_Trif null_PAM3_Poly IC_Female"] <- CSSs.tc.exon["BMDM_Trif null_PAM3_Poly IC_Female"] ## Not sure here (both are at 12hrs only)
