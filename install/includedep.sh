#!/bin/sh 
# includedep.sh -- script that computes dependencies on preprocessor includes

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# files whose dependencies must be computed
sources=`echo *.c *.f90 |
sed 's/\*\.c//g
     s/\*\.f90//g'`       # remove the "*.c" and "*.f90" that remain
#                         # when there are no such files
if test "$sources" = " " ; then exit ; fi

# files that may be included
# extra directories may be specified on the command line
includes=`echo *.h`
for dir in $*
do
    includes="$includes `echo $dir/*.h`"
done
includes=`echo $includes |
sed 's/[^ ]*\*\.h//g'`     # remove the "dir/*.h" that remain
#                          # when there are no such files

# create list of include dependencies
# each line is of the form:
# file_name.o : @include_file.h@
egrep -H '^ *# *include *"' $sources |  # look for #include "..." statements
#                                    #   ignore #include <...> ones
sed 's/f90:/o /
     s/c:/o /
     s/# *include *//
     s/\"/ /g' |                     # replace extension, insert space
#                                    #   remove '# include' statements
#                                    #   remove quotes 
awk '{print $1 " : @" $2 "@"}' |     # create dependency entry
sort | uniq > includedep.tmp1        # remove duplicates

# create list of available include files
# for each file, create a line of the form:
# s/@file_name@/pathname/g
echo $includes | tr " " "\n" |
sed 's/\//\\\//g
     s/.*\/\([^/]*\)/\1 &/' |         # escape slashes
awk '{print "s/@" $1 "@/" $2 "/" }' > includedep.tmp2

# replace file names with pathnames
# by applying the file of substitution patterns just created
sed -f includedep.tmp2 includedep.tmp1

rm -f includedep.tmp1 includedep.tmp2 # remove temporary files
