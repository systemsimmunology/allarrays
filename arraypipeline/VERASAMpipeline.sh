#!/bin/bash

#
#  Run VERA and SAM 
#

#WRONGARGS=1
#if [ $# != 3 ]
#then
#    echo "Usage: `basename $0` <Expression matrix with replicates> <Group File> <Output Folder>" >&2
#    exit $WRONGARGS
#fi

BINDIR=$HOME/bin
SCRIPTDIR=$AA/arraypipeline
VERAOCDIR=$HOME/VERAandSAM/oneChannel
SAMSTDERRDIR=$HOME/VERAandSAM/SAM_standard_error
WORKDIR=./tmp

#groupfile=$HOME/allarrays/auxfiles/LPSthreePrimeCelFiles-groupings.tsv
groupfile=$HOME/allarrays/auxfiles/grouptemp.tsv
#metafile=$HOME/allarrays/auxfiles/LPSthreePrimeCelFiles-GroupMetaData.tsv
metafile=$HOME/allarrays/auxfiles/metatemp.tsv
expfile=$HOME/allarrays/data/20100407.curated.3prime/emat.tsv
## create all auxiliary files
$SCRIPTDIR/createHeaders.py $groupfile

groups=`awk '{print $1}' $metafile | tail -n +3`
ref="BMDM_Bl6_LPS_0000___Female"

echo '****************Obtaining Error Model for ' $ref ' *****************'
$BINDIR/keepColumns.py $expfile $ref.headers > $ref.m1
$SCRIPTDIR/twoPower.py $ref.m1 $ref.dataColumnHeaders > $ref.m2
$VERAOCDIR/VERA_OC $ref.m2 $ref.mod
for gr in $groups; do
    echo '****************Obtaining Error Model for ' $gr ' *****************'
    $BINDIR/keepColumns.py $expfile $gr.headers > $gr.m1
    $SCRIPTDIR/twoPower.py $gr.m1 $gr.dataColumnHeaders > $gr.m2
    $VERAOCDIR/VERA_OC $gr.m2 $gr.mod
done


## Create input files for SAM
echo '************ Creating input files for SAM *******************'
## Create paired data file 
$BINDIR/keepColumns.py $ref.m2 $ref.XtoYheaders > $ref.Y
paste $ref.m2 $ref.Y > $ref.m3
for gr in $groups; do
    $BINDIR/keepColumns.py $gr.m2 $gr.XtoYheaders > $gr.Y
    paste $gr.m2 $gr.Y > $gr.m3
done

## Create model files for SAM (mock two-channel model)
for gr in $groups; do
    ##ofile=`$gr + 'vs' + $ref + '.mod'`
    $SCRIPTDIR/createSAMmodfile.py $gr.mod $ref.mod > "$gr"vs"$ref".mod
done

echo '****************Significance Testing*****************'
for gr in $groups; do
    $SAMSTDERRDIR/SAM $gr.m3 "$gr"vs"$ref".mod $gr.sig
done

awk '{print $1}' $expfile > ps 

echo '****************Constructing Intensity Matrix*****************'
cp ps expression.mus
for gr in $groups; do
    $BINDIR/keepOneColumn.py $gr.sig mu_X $gr > t1
    paste expression.mus t1 > t2
    mv -f t2 expression.mus
done

echo '****************Constructing Ratio Matrix*****************'
cp ps expression.log10_ratios
for gr in $groups; do
    $BINDIR/keepOneColumn.py $gr.sig muRATIO $gr > t1
    paste expression.log10_ratios t1 > t2
    mv -f t2 expression.log10_ratios
done

echo '****************Constructing Lambda Matrix*****************'
cp ps expression.lambdas
for gr in $groups; do
    $BINDIR/keepOneColumn.py $gr.sig lambda $gr > t1
    paste expression.lambdas t1 > t2
    mv -f t2 expression.lambdas
done

echo '****************Constructing Standard Error Matrix*****************'
cp ps expression.muStandardErrors
for gr in $groups; do
    $BINDIR/keepOneColumn.py $gr.sig StandErr_X $gr > t1
    paste expression.muStandardErrors t1 > t2
    mv -f t2 expression.muStandardErrors
done

## Clean up
rm -f *.headers *.dataColumnHeaders *.XtoYheaders *.Y
rm -f *.m1 *.m2 *.m3
rm -f ps t1 t2 

