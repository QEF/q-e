!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvloc_of_g (lloc, lmax, numeric, mesh, msh, rab, r, &
 vnl, cc, alpc, nlc, nnl, zp, aps, alps, tpiba2, ngl, gl, omega, &
 dvloc)
!----------------------------------------------------------------------
!
#include "machine.h"
use parameters
use allocate 
implicit none  
!
!    first the dummy variables
!
integer :: nlc, nnl, ngl, lloc, lmax, mesh, msh  
                             ! input: analytic, number of erf functions
                             ! input: analytic, number of gaussian funct
                             ! input: the number of shell of G vectors
                             ! input: the non local part which becomes l
                             ! input: the maximum non local angular mome
                             ! input: numeric, the dimensions of the mes
                             ! number of mesh points for radial integrat

real(kind=DP) :: cc (2), alpc (2), alps (3, 0:3), aps (6, 0:3), &
 zp, rab (mesh), r (mesh), vnl (mesh), tpiba2, omega, gl (ngl), &
 dvloc (ngl)
                             ! input: analytic, c of the erf functions
                             ! input: analytic, alpha of the erf
                             ! input: analytic, alpha of the gaussians
                             ! input: analytic, a and b of the gaussians
                             ! input: valence pseudocharge
                             ! input: numeric, the derivative of mesh po
                             ! input: numeric, the mesh points
                             ! input: numeric, the pseudo on the radial
                             ! input: 2 pi / alat
                             ! input: the volume of the unit cell
                             ! input: the moduli of g vectors for each s
                             ! output: the fourier transform dVloc/dG
logical :: numeric  
                             ! input: if true the pseudo is numeric
!
!   here the parameters
!

real(kind=DP) :: pi, fpi, e2, eps  
                             ! pi
                             ! four times pi
                             ! the square of the charge
                             ! a small number
parameter (pi = 3.14159265358979d0, fpi = 4.d0 * pi, e2 = 2.d0, &
 eps = 1.d-8)
!
!    and the local variables
!

real(kind=DP) :: vlcp, g2a
real(kind=DP), pointer ::  aux (:), aux1 (:)
real(kind=DP) ::  erf, gx  
                             !
                             ! auxiliary variables for faster computatio
                             ! auxiliary variable for numerical integrat
                             ! auxiliary variable for numerical integrat
                             ! the erf function
                             ! the modulus of g vectors

integer :: i, igl, igl0, l  
                             ! counter on erf functions or gaussians
                             ! counter on g shells vectors
                             ! first shells with g != 0
                             ! the angular momentum

! the  G=0 component is not computed
if (gl (1) .lt.eps) then  
   dvloc (1) = 0.0d0  
   igl0 = 2  
else  
   igl0 = 1  

endif  

if (.not.numeric) then  
   do igl = igl0, ngl  
   dvloc (igl) = 0.d0  

   enddo  
   do i = 1, nlc  
   do igl = igl0, ngl  
   g2a = gl (igl) * tpiba2 / 4.d0 / alpc (i)  
   vlcp = fpi / omega * zp * e2 * cc (i) * exp ( - g2a) * (g2a + &
    1.d0) / (gl (igl) * tpiba2) **2
   dvloc (igl) = dvloc (igl) + vlcp  
   enddo  


   enddo  
! Add the local part l=lloc term (only if l < lmax)
   l = lloc  
   if (l.le.lmax) then  
      do i = 1, nnl  
      do igl = igl0, ngl  
      g2a = gl (igl) * tpiba2 / 4.d0 / alps (i, l)  
      vlcp = - e2 * (pi / alps (i, l) ) **1.5 * exp ( - g2a) &
       * g2a * (aps (i, l) + aps (i + 3, l) * (2.5d0 - g2a) &
       / alps (i, l) ) / gl (igl) / tpiba2 / omega
      dvloc (igl) = dvloc (igl) + vlcp  
      enddo  
      enddo  

   endif  


else  
! Pseudopotentials in numerical form (Vnl(lloc) contain the local part)
! In order to perform the Fourier transform, a term erf(r)/r is
! subtracted in real space and added again in G space
   call mallocate(aux, mesh)  

   call mallocate(aux1, mesh)  
!
!   This is the part of the integrand function
!   indipendent of |G| in real space
!
   do i = 1, msh  
   aux1 (i) = r (i) * vnl (i) + zp * e2 * erf (r (i) )  

   enddo  
   do igl = igl0, ngl  


   gx = sqrt (gl (igl) * tpiba2)  
!
!    and here we perform the integral, after multiplying for the |G|
!    dependent  part
!
! DV(g)/Dg = Integral of r (Dj_0(gr)/Dg) V(r) dr
   do i = 1, msh  
   aux (i) = aux1 (i) * (r (i) * cos (gx * r (i) ) / gx - sin (gx &
    * r (i) ) / gx**2)
   enddo  


   call simpson (msh, aux, rab, vlcp)  
! DV(g^2)/Dg^2 = (DV(g)/Dg)/2g


   vlcp = fpi / omega / 2.0 / gx * vlcp  
! subtract the long-range term
   g2a = gl (igl) * tpiba2 / 4.d0  
   vlcp = vlcp + fpi / omega * zp * e2 * exp ( - g2a) * (g2a + &
    1.d0) / (gl (igl) * tpiba2) **2
   dvloc (igl) = vlcp  

   enddo  
   call mfree (aux1)  
   call mfree (aux)  

endif  
return  
end subroutine dvloc_of_g

