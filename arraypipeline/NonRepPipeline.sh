#!/bin/bash
#
#  Run simple single column extraction for condition sets with no replicates
#
#WRONGARGS=1
if [ $# != 1 ]
then
    echo "Usage: `basename $0` <Expression matrix with replicates> " >&2
    exit $WRONGARGS
fi

BINDIR=$HOME/bin
SCRIPTDIR=$AA/arraypipeline

expfile=$1
##
##
##

conds=("MyD88-PAM3" "TRIF-LPS" "TRIF-PAM2")

cnt=${#conds[@]}
for (( i = 0 ; i < cnt ; i++ ))
do
    cond=${conds[$i]}
    grep ^$cond $AA/auxfiles/ThreePrimeMasterFile.tsv > condfile
    echo -e 'ProbesetID\tProbesetID' > headerfile
    awk  '{OFS="\t";print $3,$2}' condfile >> headerfile ## columns to extract from expression mat
    awk  '{OFS="\t";print $2}' headerfile | tail -n +2 > tpfile ## columns for twoPower
    $BINDIR/keepColumns.py $expfile headerfile > t1
    $SCRIPTDIR/twoPower.py t1 tpfile > $cond.mus
done
rm -f headerfile tpfile t1 condfile
