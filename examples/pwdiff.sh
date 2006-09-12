#!/bin/sh

# pwdiff.sh -- script for checking outputs of PWscf examples
# checking is done in three steps: preprocess, diff against reference
# data, postprocess
# this way "false positives" are eliminated, that is, differences that
# don't mean that something went wrong

# pre/postprocessing scripts must be in the same directory as this one
dir=`echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
prediff=$dir/prediff.awk
postdiff=$dir/postdiff.awk
if test ! -f "$prediff"
then
    echo error: file prediff.awk not found
    exit -1
fi
if test ! -f "$postdiff"
then
    echo error: file postdiff.awk not found
    exit -1
fi

# check that files exist
if test ! -f "$1"
then
    echo $0: file $1 does not exist
    exit -1
fi
if test ! -f "$2"
then
    echo $0: file $2 does not exist
    exit -1
fi

# preprocess
awk -f $prediff $1 > pwdiff.sh.tmp1
awk -f $prediff $2 > pwdiff.sh.tmp2

# uncomment to debug
# cp pwdiff.sh.tmp1 pwdiff.sh.tmp6
# cp pwdiff.sh.tmp2 pwdiff.sh.tmp7

rm -f pwdiff.sh.tmp5
offset1=0
offset2=0
check1=1
check2=1

# check preprocessed files
while test $check1 -gt 0 -a $check2 -gt 0
do
    # look for next checkpoints
    check1=`grep -n "@CHECKPOINT@" pwdiff.sh.tmp1 | head -1 | sed 's/:.*//'`
    check2=`grep -n "@CHECKPOINT@" pwdiff.sh.tmp2 | head -1 | sed 's/:.*//'`
    eof1=0
    eof2=0
    if test "$check1" = ""
    then
	check1=`wc pwdiff.sh.tmp1 | awk '{print $1}'`
	eof1=1
    fi
    if test "$check2" = ""
    then
	check2=`wc pwdiff.sh.tmp2 | awk '{print $1}'`
	eof2=1
    fi
    # some OS do not like "head -0"
    if test $check1 -gt 0
    then
        head -$check1 pwdiff.sh.tmp1 > pwdiff.sh.tmp3
    else
        touch pwdiff.sh.tmp3
    fi
    if test $check2 -gt 0
    then
        head -$check2 pwdiff.sh.tmp2 > pwdiff.sh.tmp4
    else
        touch pwdiff.sh.tmp4
    fi
    # diff up to next checkpoints, then postprocess
    diff -wib pwdiff.sh.tmp3 pwdiff.sh.tmp4 \
	   | awk -f $postdiff o1=$offset1 o2=$offset2 >> pwdiff.sh.tmp5

    # discard processed part
    sed "1,${check1}d" pwdiff.sh.tmp1 > pwdiff.sh.tmp3
    sed "1,${check2}d" pwdiff.sh.tmp2 > pwdiff.sh.tmp4
    mv -f pwdiff.sh.tmp3 pwdiff.sh.tmp1
    mv -f pwdiff.sh.tmp4 pwdiff.sh.tmp2
    offset1=`expr $offset1 + $check1 + $eof1 - 1`
    offset2=`expr $offset2 + $check2 + $eof2 - 1`

done

# return success if there's no output, failure otherwise (as diff does)
rvalue=`wc pwdiff.sh.tmp5 | awk '{print ($1 == 0) ? 0 : 1}'`
cat pwdiff.sh.tmp5
rm -f pwdiff.sh.tmp[1-5]
exit $rvalue;
