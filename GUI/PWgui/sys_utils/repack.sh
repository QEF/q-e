#!/bin/sh -f

HERE=$(pwd)

if [ $# -ne 3 ]
    then
    echo "Usage: $0 src_filehead dst_filehead dst_dirname"
    exit 1
fi

# arguments ...
SRC_FILEHEAD=$1
DST_FILEHEAD=$2
DST_DIR=$3

# files ...
SRC_TAR=$SRC_FILEHEAD.tar
DST_TGZ=$DST_FILEHEAD.tgz

# directory ...
TMPDIR=/tmp/$DST_DIR

if [ -d $TMPDIR ]
    then 
    rm -rf $TMPDIR
fi
mkdir $TMPDIR

if [ -f $SRC_TAR ]; 
    then
    cp $SRC_TAR $TMPDIR/
else
    echo "File: $SRC_TAR does not exists !!!"
    exit 1
fi


cd $TMPDIR 
tar xvf $SRC_TAR; rm $SRC_TAR

cd $TMPDIR/..
tar zcvf $DST_TGZ $DST_DIR/
mv $DST_TGZ $HERE/
rm -rf $TMPDIR
