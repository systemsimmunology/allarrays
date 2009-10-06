
dMat <- 1-cor(t(lps.ratios.w.zero[targSet,]))
dr <- dist(dMat)

hc <- hclust( dr , "complete" )
x11()
plot(hc, main="complete")

dMat.exon <- 1-cor(dm.exon[eids,])
dr.exon <- dist(dMat.exon)
hc.exon <- hclust( dr.exon , "complete" )
x11()
plot(hc.exon, main="complete")

dMat.3prime <- 1-cor(dm.3prime[eids,])
dr.3prime <- dist(dMat.3prime)
hc.3prime <- hclust( dr.3prime , "complete" )
x11()
plot(hc.3prime, main="complete")
