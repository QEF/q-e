#!/bin/sh
#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwi2xsf.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#------------------------------------------------------------------------
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# pwi2xsf.sh: PW-input to XSF conversion
#
# Usage: pwi2xsf [-r] pw-input-file
#
#------------------------------------------------------------------------



pwError() {
    echo "
 ========================================================================
    $1
 ========================================================================
" 1>&2
    if [ "$2" -ge 0 ]; then
	exit $2
    fi
}

# --------------------------------------------------------------------------
# FUNCTION:   pwNucleiCharges --
#
# Purpose:    ityp->nat conversion data 
#
# Usage:      pwNucleiCharges pw_input|pw_output outfile
#
# Side efect: creates nuclei.charges file
# --------------------------------------------------------------------------

pwNucleiCharges() {
    #
    # if file nuclei.charges does not exists prompt for ityp->nat conversion !!
    #

    if [ \( "$1" = "" \) -o \( "$2" = "" \) ]; then
	pwError "Usage:  pwNucleiCharges  pw_input|pw_output  outfile" 1
    fi
    
    # do we have PW-INPUT or PW-OUTPUT file ???
    
    if [ "`cat $1 | egrep '&input|&system'`" != "" ]; then
	# it is PW-INPUT
	ntyp=`cat "$1" | awk '{gsub(",","\n"); print}' | grep ntyp \
		| awk '{split($0,a,"=|,"); print a[2];}'`
    else
	# PW-OUTPUT
	ntyp=`cat "$1" | grep 'number of atomic types' | \
	    head -1 | awk '{print $NF}'`
	#echo 'NTYP=$ntyp'
	if [ "$ntyp" = "" ]; then
	    # some older PWSCF versions didn't have "number of atomic
	    # types" printout -> user will have to make nuclei.charges
	    # file by himself/herself !!!
	    pwError "This is either non PW-output file or is a PW-output file 
    produced with some old PWSCF version" -1
	    echo -n "How many ityp->nat replacements ? " 1>&2
	    read ntyp
	fi
    fi		     

    if [ ! -f nuclei.charges ]; then
	echo -n "Please enter $ntyp ityp->nat replacements !!! " 1>&2
	echo $ntyp > nuclei.charges
	i=0
	while [ $i -lt "$ntyp" ]
	do
	    i=`echo "$i + 1"|bc`
	    echo "" 1>&2
	    echo "Replacement #${i}: ityp->nat" 1>&2
	    echo -n "ityp[$i]=$i; nat[$i]=" 1>&2; read nat
	    echo "$i $nat" >> nuclei.charges
	done
    fi
	
    cat nuclei.charges > "$2"
}


# ------------------------------------------------------------------------
#  MAIN
# ------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
    echo "
Usage: pwi2xsf.sh [-r] pw-input

-r ... one must specify i.e. ityp->nat conversion, and the corresponding
       data are written to file nuclei.charges. The -r flag deletes this
       file.
" 1>&2
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
  end=0;
}

toupper($0) ~ /&SYSTEM/          { sysnml=1; print; }
toupper($0) ~ /IBRAV|CELLDM|NAT/ { print; }

/ATOMIC_SPECIES|ATOMIC_POSITIONS|K_POINTS|CELL_PARAMETERS/ {
  print_line=1;
  print toupper($0);
  next;
}

toupper($0) ~ /&END|^\/|^ \// { 
  if ( sysnml == 1 ) { 
    print; 
    sysnml=0; 
  } 
}

/a*/ {
  if ( print_line == 1 ) {
    print toupper($0);
  }
}'> pw.$$
    PWI2XSF=pwi2xsf_new
else 
    # we have OLD PW.X input format	

    pwNucleiCharges $1 /dev/null

    cat $1 | awk 'BEGIN {RS=",";} {print}' | awk '
BEGIN {
    end=0;
}
toupper($0) ~ /&INPUT|CELLDM|NAT|LTAUCRY/ { print; }
toupper($0) ~ /IBRAV/ { 
    print;
    split($0,a,"=");
    split(a[1],b,",");
    ibrav = b[1];
}
toupper($0) ~ /&END|^\/|^ \// { end=1; }
/a*/ {
    if ( end == 1 ) print;
}' > pw.$$
    PWI2XSF=pwi2xsf
fi
#
# execute $PWI2XSF fortran program and print the XSF file
#
if [ \( "$XCRYSDEN_TOPDIR" != "" \) -a \( -x $XCRYSDEN_TOPDIR/bin/$PWI2XSF \) ]; then
    $XCRYSDEN_TOPDIR/bin/$PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
else
    $PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
fi
rm -f pw.$$

if [ "$r" -eq 1 ]; then
    if [ -f nuclei.charges ]; then rm nuclei.charges; fi
fi

exit 0
