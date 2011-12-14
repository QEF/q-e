BEGIN {nr=0; nrs=0; nat=0; nstep=0; 
       print "* XBS file created by pawk.bs "; print "";

print ""
print "* the following are AUXILIARY lines defining the bonds as "
print "* bonds spec1 spec2 dmin dmax bondthickness grayscale (white==1.0)"
print "* "
print "* bonds S H 0.1 0.6  0.0500     1.0"
print ""

       }
{ if ($3=="atoms/cell" && nr==0) {nat=$5};
  if ($1=="lattice" && $2=="parameter" && nr==0 ) {alat= $5*0.529177}
  if ($1=="a(1)" && nr==0) \
                 {print"* it might be useful to duplicate as follows" ;
                  print "* dup ",$4*alat,$5*alat,$6*alat}
  if ($1=="a(2)" && nr==0) {print "* dup ",$4*alat,$5*alat,$6*alat}
  if ($1=="a(3)" && nr==0) {print "* dup ",$4*alat,$5*alat,$6*alat;
                            print " "}
  if ($1=="atomic" && $2=="species" && nrs==0 ) \
     {nrs=NR+nat+1
      print "* the following are MANDATORY lines defining the atomic species as"
      print "* spec  name   radius grayscale (white==1.0) "
      print "* "}
  if (NR<=nrs) {if (NF==0) {print ""; nrs=-1}
                if (NF>0 && $1!="atomic") \
                   printf ( "spec %2s %6.2f %4.2f \n", \
                            $1, 0.4, 1.0/$2 ) }
  if ($1=="site" && nr==0 ) \
     {nr=NR
      print "* the following are MANDATORY lines defining the atomic positions"
      print "* atom  name  x y z  dummyinteger"
      print "* "}
  if (NR-nr>0 && NR-nr<=nat && nr>0) \
     printf ( "atom %2s %10.7f %10.7f %10.7f %3d \n", \
                    $2,  $(NF-3)*alat, $(NF-2)*alat, $(NF-1)*alat, $1) 
 }
END{
print ""
#print "* the following are AUXILIARY lines defining the bonds as "
#print "* bonds spec1 spec2 dmin dmax bondthickness grayscale (white==1.0)"
#print "* "
#print "* bonds S H 0.1 0.6  0.0500     1.0"
print ""
print "tmat  0.000  1.000  0.000 0.000  0.000  1.000 1.000  0.000  0.000"
print ""
print "dist    50.000"
print "inc      1.000"
print "scale   50.000"
print "rfac 1.00"
print "bfac 1.00"
print "pos    -10.000    -100.00"
print "switches 1 0 1 0 0 1 1 0 0"
}
