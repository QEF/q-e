# Description

This example show how to use pw.x and kcw.x to compute the 
electronic structure of FCC Silicon with and without SOC. 

Different simulations are compared to test the non-collinear implementation 
 * collinear     w/o spin         (nspin = 1) 
 * collinear     w   spin         (nspin = 2) 
 * non collinear w/o spin w/o SOC (nspin = 4, domag=F, lspinorb=F) 
 * non collinear w   spin w/o SOC (nspin = 4, domag=T, lspinorb=F)
 * non collinear w/o spin w   SOC (nspin = 4, domag=F, lspinorb=T) 
 * non collinear w   spin w   SOC (nspin = 4, domag=T, lspinorb=T)

Without SOC, the non-collinear calculations w/o and w spin (domag=F, domag=T) 
should match with the reference collinear calculations w/o and w spin (nspin=1, nspin=2)
This is checked by the script compare.gnu that compares the KI band structure for the 
collinear and non-collinear-without-SOC (with and without spin degrees of freedom) cases.

