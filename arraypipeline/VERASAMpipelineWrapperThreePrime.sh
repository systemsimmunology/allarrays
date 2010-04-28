#!/bin/bash
##Important!! run in a tmp directory
stims=("PAM2" "PAM3" "R848" "Poly-IC" "CpG")
zeros=("BMDM_Bl6_PAM2_0000___Male" "BMDM_Bl6_PAM3_0000___Male" "BMDM_Bl6_R848_0000___Female" "BMDM_Bl6_Poly-IC_0000___Female" "BMDM_Bl6_CpG_0000___Female")

cnt=${#stims[@]}
echo "Number of elements: $cnt"
for (( i = 0 ; i < cnt ; i++ ))
do
    stim=${stims[$i]}
    zero=${zeros[$i]}
    echo "Processing [$i]: $stim. Zero: $zero"

    grep ^$stim $AA/auxfiles/try.tsv | awk '{OFS="\t"; print $3,$2}' > groupings.tsv
    $AA/arraypipeline/VERASAMpipeline.sh $AA/data/20100426.curated.3prime/emat.tsv groupings.tsv $zero
    cp expression.mus /Users/thorsson/allarrays/data/20100426.curated.3prime/$stim.mus
    cp expression.log10_ratios /Users/thorsson/allarrays/data/20100426.curated.3prime/$stim.log10_ratios
    cp expression.lambdas /Users/thorsson/allarrays/data/20100426.curated.3prime/$stim.lambdas      
    cp expression.muStandardErrors /Users/thorsson/allarrays/data/20100426.curated.3prime/$stim.muStandardErrors
    rm -f *

done
