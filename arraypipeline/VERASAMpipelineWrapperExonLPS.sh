#!/bin/bash
##Important!! run in a tmp directory
cond="LPS"
zero="BMDM_Bl6_LPS_0000___Female"
$AA/arraypipeline/VERASAMpipeline.sh $AA/data/20121001.curated.exon/emat.tsv $AA/auxfiles/LPSexonCelFiles-groupings.tsv $zero
cp expression.mus $AA/data/20121001.curated.exon/$cond.mus
cp expression.log10_ratios $AA/data/20121001.curated.exon/$cond.log10_ratios
cp expression.lambdas $AA/data/20121001.curated.exon/$cond.lambdas      
cp expression.muStandardErrors $AA/data/20121001.curated.exon/$cond.muStandardErrors

