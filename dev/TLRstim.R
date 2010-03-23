
## Run the beginning of CombineArrayPlatforms
## but do not combine

dm.lps.exon <- dm.exon[,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]][["DM Column"]]]
dm.lps.3prime <- dm.3prime[,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]][["DM Column"]]]

save(file="lps.RData",list=c("dm.lps.exon","dm.lps.3prime"))
