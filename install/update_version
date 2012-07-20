#!/bin/sh

# check if svn info available (do not write anything)
svn info 2> /dev/null > /dev/null

if [ $? = 0 ] ; then 
# svn info available: get svn revision
   svn_rev=$(svnversion -n)
else
# svn info available: revert to no info
   svn_rev=unknown
fi

# write svn into file version_tmp.f90
cat version.f90.in | sed 's/unknown/'$svn_rev'/' > version.f90.tmp

# check if a previous version.f90 file exists
if test -f version.f90 ; then

# version.f90 existing: check if new and previous files differ
   diff -wib version.f90.tmp version.f90  2> /dev/null > /dev/null

   if [ $? = 1 ] ; then 
# they differ: update file version.f90
      mv version.f90.tmp version.f90
   fi
# do not update if files are the same (prevents useless recompilation)

else

# file version.f90 not existing: create one
   mv version.f90.tmp version.f90

fi

