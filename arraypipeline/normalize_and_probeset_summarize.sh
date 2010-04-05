#!/bin/bash

WRONGARGS=1
if [ $# != 3 ]
then
  echo "Usage: `basename $0` <CEL File list with full paths> <Destination directory> <Temporary Working directory (Created and Deleted)>" >&2
  exit $WRONGARGS
fi

cellist=$1
OUTDIR=$2
TEMPDIR=$3

ORIGDIR=`echo $PWD`
#TEMPDIR="./tmp"

mkdir $TEMPDIR
cd $TEMPDIR

cfiles=( $(cat $cellist ) )
nc=${#cfiles[@]}

##
## Create symbolic links and file list without path
##
for (( i=0 ; i<$nc ; i++ )); do
    cfile_i=${cfiles[$i]}
    cflocal=`basename $cfile_i`
    echo $cflocal >> filelist
    ln -s $cfile_i .
done

R --vanilla --slave --args filelist . < ../normalize_and_probeset_summarize.R

## Move files out of Temp dir and clean up
mv emat.RData $OUTDIR
mv emat.tsv $OUTDIR
mv filelist $OUTDIR/celfiles_no_path

cd $ORIGDIR
rm -fr $TEMPDIR
