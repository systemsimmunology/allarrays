
unique(unlist(lapply(CSSs.tc.exon,"[[","Strain")))


mutants.3prime <- setdiff(unique(unlist(lapply(CSSs.tc.3prime,"[[","Strain"))),"Bl6")



for ( strain in mutants.3prime ){
  eid <- gene.eid[mt$GeneSymbols[which(strain==mt$Strain)]] ## See below to get this to work
  conds <- names(which(unlist(lapply(CSSs.tc.3prime,"[[","Strain"))==strain))
  cat(strain,gene.symbol[eid],"No. TCs",length(conds),"\n")
  cat(conds,"\n")
}

gridPlotCSS(eid,CSSs.tc.3prime[conds],data.matrix=dm.3prime,nx=3,ny=3)
x11()
gridPlotCSS(eid,CSSs.tc.3prime,data.matrix=dm.3prime,nx=4,ny=5)


mutants.exon <- setdiff(unique(unlist(lapply(CSSs.tc.exon,"[[","Strain"))),"Bl6")

for ( strain in mutants.exon ){
  eid <- gene.eid[mt$GeneSymbols[which(strain==mt$Strain)]] ## See below to get this to work
  conds <- names(which(unlist(lapply(CSSs.tc.exon,"[[","Strain"))==strain))
  cat(strain,gene.symbol[eid],"No. TCs",length(conds),"\n")
  ##cat(conds,"\n")
}




ind <- 2
strain <- mutants.exon[ind]
eid <- gene.eid[mt$GeneSymbols[which(strain==mt$Strain)]] ## See below to get this to work
conds <- names(which(unlist(lapply(CSSs.tc.exon,"[[","Strain"))==strain))
cat(strain,"No. conds",length(conds),"\n")


gridPlotCSS(eid,CSSs.tc.exon[conds],data.matrix=dm.exon,nx=3,ny=3)
x11()
gridPlotCSS(eid,CSSs.tc.exon,data.matrix=dm.exon,nx=8,ny=10)


###  Table of mutants
mt <- read.table("MutantEffects.tsv",as.is=TRUE,sep='\t',header=TRUE)
## Strains
mt$Strain
## works for all but Trif-Myd88
gene.eid[mt$GeneSymbols]
## Can be used to map directly to gene symbols or eids in the above
