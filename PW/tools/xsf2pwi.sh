#!/bin/sh

# Usage: xsf2pwi.sh [-c] XSF2-file
# Purpose: convert XSF2 file to PW.X input syntax

coor_only=0
if test x$1 = x"-c"; then
    # coor only option specified
    coor_only=1
    shift
fi

head=prefix
if test $# -lt 1; then
    input=-
else
    input=$1
    head=${input%.xsf*}
fi


cat $input | awk -v coor_only=$coor_only '
BEGIN {
  f=1.0;
  bohr=0.529177;
}
/PRIMVEC/ { 
  if ( $2 != "bohr" ) {
    f = 1.0 / bohr;  
  }

  if (!coor_only) {
     print "CELL_PARAMETERS cubic";
     getline; printf "   %12.6f %12.6f %12.6f\n", $1*f, $2*f, $3*f; 
     getline; printf "   %12.6f %12.6f %12.6f\n", $1*f, $2*f, $3*f; 
     getline; printf "   %12.6f %12.6f %12.6f\n", $1*f, $2*f, $3*f; 
     print "";
  }
}
/PRIMCOORD/ {
  if ( NF < 2 ) {
    unit="angstrom";
  } else {
    unit=$2;
  }
  print "ATOMIC_POSITIONS ", unit;
  getline;
  nat=$1;
  for (i=0; i<nat; i++) 
   {
     getline; print;
   }
   print "";
}
'
