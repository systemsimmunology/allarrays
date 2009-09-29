

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




