#!/bin/bash
#
#IMPORTANT NOTE: run the script from dev-tools/

export ESPRESSO_ROOT=${PWD}/..

cat > ignorelist << END
*_tmp.f90
*.mod
*.o
*.a
*.x
END

cat > dirlist << END
upftools
upftools/HGH2QE
XSpectra/src
LR_Modules
COUPLE/src
COUPLE/examples
LAXlib
PP/src
PP/simple_transport/src
TDDFPT/src
TDDFPT/tools
FFTXlib
atomic/src
QHA/SRC
QHA/Debye
EPW/src
NEB/src
CPV/src
GWW/bse
GWW/gww
GWW/pw4gww
GWW/head
PHonon/FD
PHonon/Gamma
PHonon/PH
PW/src
PW/tools
PlotPhon/SRC
PWCOND/src
GUI/PWgui/external/src/
clib
END

for x in `cat dirlist`
do
svn propset svn:ignore -F ./ignorelist ${ESPRESSO_ROOT}/${x}
done

# Special case for "Modules"
cat > ignorelist << END
version.f90
version.f90.tmp
*_tmp.f90
*.mod
*.o
*.a
*.x
END

svn propset svn:ignore -F ./ignorelist ${ESPRESSO_ROOT}/Modules

# Cleanup
rm ignorelist dirlist
