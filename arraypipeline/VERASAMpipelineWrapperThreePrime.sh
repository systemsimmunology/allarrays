#!/bin/bash
##Important!! run in a tmp directory
conds=("LPS" "PAM2" "PAM3" "R848" "Poly-IC" "CpG" "MyD88")
zeros=("BMDM_Bl6_LPS_0000___Female" "BMDM_Bl6_PAM2_0000___Male" "BMDM_Bl6_PAM3_0000___Male" "BMDM_Bl6_R848_0000___Female" "BMDM_Bl6_Poly-IC_0000___Female" "BMDM_Bl6_CpG_0000___Female" "BMDM_Myd88-null_LPS_0000___Female")

cnt=${#conds[@]}
for (( i = 0 ; i < cnt ; i++ ))
do
    cond=${conds[$i]}
    zero=${zeros[$i]}
    echo "Processing [$i]: $cond. Zero: $zero"

    grep ^$cond $AA/auxfiles/ThreePrimeMasterFile.tsv  | awk '{OFS="\t"; print $3,$2}' > groupings.tsv
    $AA/arraypipeline/VERASAMpipeline.sh $AA/data/20100426.curated.3prime/emat.tsv groupings.tsv $zero
    cp expression.mus /Users/thorsson/allarrays/data/20100426.curated.3prime/$cond.mus
    cp expression.log10_ratios /Users/thorsson/allarrays/data/20100426.curated.3prime/$cond.log10_ratios
    cp expression.lambdas /Users/thorsson/allarrays/data/20100426.curated.3prime/$cond.lambdas      
    cp expression.muStandardErrors /Users/thorsson/allarrays/data/20100426.curated.3prime/$cond.muStandardErrors
    rm -f *

done
