# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_MAKEDEPS], [

if test ! $topdir -ef $topbuilddir ; then
   $topdir/install/makedeps.sh $topbuilddir
else
   $topdir/install/makedeps.sh
fi

]
)
