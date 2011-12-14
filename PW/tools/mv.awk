BEGIN {nr=0; nat=0; nline=0; nframe=0; label=""; print}
{ if ($3=="atoms/cell" && nr==0) {nat=$5 };
  if ($1=="lattice" && $2=="parameter" && nr==0 ) {alat= $5*0.529177}
  if ($1=="Search") {label="BFGS Search" };
  if ($1=="Final" && $2=="estimate") {label=$0};
  if ($1=="ATOMIC_POSITIONS") {nframe=nframe+1; nr=NR; nline=nat; print "frame ",nframe,label ; print " "}; 
  if (NR-nr>=1 && NR-nr<=nline)  print $2*alat,$3*alat,$4*alat 
  if (NR-nr==nline) {print " ";nline=0;label=""}
}
