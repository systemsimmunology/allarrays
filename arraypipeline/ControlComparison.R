
##go to tmp directory. Turn off rm -f * in scripts.
##grep 0000 $AA/auxfiles/ThreePrimeMasterFile.tsv | awk '{OFS="\t"; print $3,$2}' > tempfile
##Edit to yield only triplicated 
##mv tempfile unstim-groupings.tsv
##../VERASAMpipeline.sh $AA/data/20100426.curated.3prime/emat.tsv unstim-groupings.tsv BMDM_Bl6_LPS_0000___Female

ls <- read.matrix("~/allarrays/arraypipeline/tmp/expression.lambdas")

collect <- numeric()
for ( n in 1:ncol(ls)) {
  collect <- c(collect, quantile(ls[,n],0.95))
}

mean(collect) ## 26.61275

