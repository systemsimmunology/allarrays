#!/bin/bash

WRONGARGS=1
if [ $# != 4 ]
then
  echo "Usage: `basename $0` <CEL File list with full paths> <cdfName> <Destination directory> <Temporary Working directory (Created and Deleted)>" >&2
  echo "cdfName e.g. Mouse4302_Mm_ENTREZG, MoEx10stv1_Mm_ENTREZG" >&2

  exit $WRONGARGS
fi


cellist=$1
cdfName=$2
OUTDIR=$3
TEMPDIR=$4

ORIGDIR=`echo $PWD`

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

R --vanilla --slave --args filelist . $cdfName  < $AA/arraypipeline/normalize_and_probeset_summarize.R

## Move files out of Temp dir and clean up
mv emat.RData $OUTDIR
mv emat.tsv $OUTDIR
mv filelist $OUTDIR/celfiles_no_path

cd $ORIGDIR
rm -fr $TEMPDIR
