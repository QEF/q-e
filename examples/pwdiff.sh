#!/bin/sh

# check that files exist
if ! test -f "$1"
then
    echo $0: file $1 does not exist
    exit 1
fi
if ! test -f "$2"
then
    echo $0: file $2 does not exist
    exit 1
fi

# preprocess
awk -f prediff.awk $1 > pwdiff.sh.tmp1
awk -f prediff.awk $2 > pwdiff.sh.tmp2

# uncomment to debug
# cp pwdiff.sh.tmp1 pwdiff.sh.tmp6
# cp pwdiff.sh.tmp2 pwdiff.sh.tmp7

rm -f pwdiff.sh.tmp5
offset1=0
offset2=0
check1=1
check2=1

while test $check1 -gt 0 || test $check2 -gt 0
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

    head -$check1 pwdiff.sh.tmp1 > pwdiff.sh.tmp3
    head -$check2 pwdiff.sh.tmp2 > pwdiff.sh.tmp4

    # diff up to next checkpoint, then postprocess
    diff pwdiff.sh.tmp3 pwdiff.sh.tmp4 \
	| awk -f postdiff.awk o1=$offset1 o2=$offset2 >> pwdiff.sh.tmp5

    # discard processed part
    sed "1,${check1}d" pwdiff.sh.tmp1 > pwdiff.sh.tmp3
    sed "1,${check2}d" pwdiff.sh.tmp2 > pwdiff.sh.tmp4
    mv -f pwdiff.sh.tmp3 pwdiff.sh.tmp1
    mv -f pwdiff.sh.tmp4 pwdiff.sh.tmp2
    offset1=`expr $offset1 + $check1 + $eof1 - 1`
    offset2=`expr $offset2 + $check2 + $eof2 - 1`
done

# return success if there's no output, failure otherwise (just like diff)
rvalue=`wc pwdiff.sh.tmp5 | awk '{print ($1 == 0) ? 0 : 1}'`
cat pwdiff.sh.tmp5
rm -f pwdiff.sh.tmp[1-5]
exit $rvalue;
