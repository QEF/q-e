BEGIN {nr=0; nat=0; nline=0; print}
{ if ($3=="atoms/cell" && nr==0) {nat=$5 };
  if ($1=="lattice" && $2=="parameter" && nr==0 ) {alat= $5*0.529177}
  if ($1=="Search") {nr=NR; nline=nat+2; print "frame ",$7,$9,$10; print " "}; 
  if (NR-nr>=2 && NR-nr<nline)  print $2*alat,$3*alat,$4*alat 
  if (NR-nr==nline) {print;nline=0}
}
