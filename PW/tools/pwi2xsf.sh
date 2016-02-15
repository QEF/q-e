#!/bin/sh
#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#                                                                           #
# Copyright (c) 2004 by Anton Kokalj                                        #
#############################################################################

#------------------------------------------------------------------------
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#------------------------------------------------------------------------
# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

#
# pwi2xsf.sh: PW-input to XSF converison
#
# Usage: pwi2xsf [-r] pw-input-file
#
# Last major rewrite by Tone Kokalj on Mon Feb  9 12:48:10 CET 2004
#                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

if [ "$#" -lt 1 ]; then
    echo "
Usage: pwi2xsf.sh [-r] pw-input

Option for PWscf version < 1.2:

-r ... one must specify i.e. ityp->nat conversion, and the corresponding
       data are writen to file nuclei.charges. The -r flag deletes this
       file.
"
    exit 1
fi

r=0
if [ "$1" = "-r" ]; then
    r=1
    shift
fi

#
# check if we have OLD or NEW PW.X input format
#
new_format1=`grep 'ATOMIC_POSITIONS' $1`
new_format2=`grep -i '&system' $1`

if [ "$new_format1" != ""  -a  "$new_format2" != "" ]; then
    # we have NEW PW.X input format
    #
cat $1 | awk 'BEGIN {RS=",";} {print $0}' | awk '
BEGIN {
  calculation="";
  num_of_images="";
  nml_end=0;
  nml_end_string="";
}

toupper($0) ~ /&SYSTEM/          { print; }

/=/ { 
  if ( toupper($1) ~ /^IBRAV($|=)|^CELLDM\([1-6]\)($|=)|^NAT($|=)|^A($|=)|^B($|=)|^C($|=)|^COSAB($|=)|^COSAC($|=)|^COSBC($|=)/ ) { print; } 
  
  if ( toupper($1) ~ /^CALCULATION($|=)/ ) { calculation=toupper($0); }
  
  if ( toupper($1) ~ /^STRING_METHOD($|=)/ ) { calculation="calculation=" toupper($3); }

  if ( toupper($1) ~ /^NUM_OF_IMAGES($|=)/ ) { num_of_images=toupper($0); }
}

/ATOMIC_POSITIONS|CELL_PARAMETERS|FIRST_IMAGE/ {
  if ( !nml_end) {
     # first finish the namelist
     nml_end=1;
     if (calculation != "")   print calculation;
     if (num_of_images != "") print num_of_images;
     print nml_end_string;
  }
  # now print the current record
  print_line=1;
  print toupper($0);
  next;
}


toupper($0) ~ /&END|^\/|^ \// { 
  nml_end_string=$0;
}

/a*/ {
  if ( print_line == 1 ) {
    print toupper($0);
  }
}'> pw.$$
else 
    # we have OLD PW.X input format	
    echo "
------------------------------------------------------------------------
   ERROR: This is NOT a PW-input or an input for an older PW version
------------------------------------------------------------------------
"
    exit 1
fi

#
# execute $PWI2XSF fortran program and print the XSF file
#
if test -f $XCRYSDEN_TOPDIR/bin/$PWI2XSF ; then
    $XCRYSDEN_TOPDIR/bin/pwi2xsf < pw.$$ | tee pwi2xsf.xsf_out
else
    pwi2xsf.x < pw.$$ | tee pwi2xsf.xsf_out
fi
rm -f pw.$$

if [ "$r" -eq 1 ]; then
    if [ -f nuclei.charges ]; then rm nuclei.charges; fi
fi

exit 0
