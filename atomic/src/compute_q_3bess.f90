!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine compute_q_3bess(ldip,lam,ik,chir,phi_out,ecutrho)
  !--------------------------------------------------------------------------
  !
  !   This routine computes the phi_out function by pseudizing the
  !   chir function with a linear combination of three Bessel functions
  !   multiplied by r**2. In input it receives the point
  !   ik where the cut is done, the angular momentum lam of the 
  !   bessel functions and the function chir.
  !   Phi_out has the same ldip dipole moment of chir.
  !
  !      
  use kinds, only : DP
  use radial_grids, only: ndmx
  use ld1inc, only: grid
  implicit none

  integer ::    &
       ldip,    & ! input: the order of the dipole
       lam,     & ! input: the angular momentum
       ik         ! input: the point corresponding to rc

  real(DP) :: &
       xc(8)      ! output: the coefficients of the Bessel functions

  real(DP) ::         &
       chir(ndmx),    &   ! input: the all-electron function
       phi_out(ndmx)      ! output: the phi function
  !
  real(DP) ::  &
       ecutrho,& ! the expected cut-off on the charge density for this q
       fae,    & ! the value of the all-electron function
       f1ae,   & ! its first derivative
       f2ae,   & ! the second derivative
       dip       ! the norm of the function

  integer ::    &
       n, nst, nc

  real(DP) :: &
       gi(ndmx), j1(ndmx,3), jnor, cm(3), bm(3), delta, gam

  real(DP), external :: deriv_7pts, deriv2_7pts, int_0_inf_dr

  integer ::  &
       iok,   &  ! flag
       nbes      ! number of Bessel functions to be used

  nbes = 3
  !
  nst=lam+2+ldip
  !
  !   compute the first and second derivative of input function at r(ik)
  !
  fae=chir(ik)
  f1ae=deriv_7pts(chir,ik,grid%r(ik),grid%dx)
  f2ae=deriv2_7pts(chir,ik,grid%r(ik),grid%dx)
  !
  !   compute the ldip dipole moment of the input function
  !
  do n=1,ik+1
     gi(n)=chir(n) * grid%r(n)**ldip  
  enddo
  dip=int_0_inf_dr(gi,grid,ik,nst)
  !
  !   RRKJ: the pseudo-wavefunction is written as an expansion into 3  
  !         spherical Bessel functions for r < r(ik)
  !   find q_i with the correct log derivatives
  !   
  call find_qi(f1ae/fae,xc(nbes+1),ik,ldip,nbes,2,iok)
 
  if (iok.ne.0) &
     call errore('compute_q_3bess', 'problem with the q_i coefficients', 1)
  !
  !   compute the Bessel functions and multiply by r**2
  !
  do nc=1,nbes
     call sph_bes(ik+5,grid%r,xc(nbes+nc),ldip,j1(1,nc))
     jnor=j1(ik,nc)*grid%r2(ik)
     do n=1,ik+5
        j1(n,nc)=j1(n,nc)*grid%r2(n)*chir(ik)/jnor
     enddo
  enddo
  !
  !    compute the bm functions (second derivative of the j1)
  !    and the ldip dipole moment of the Bessel function (cm)
  !
  do nc=1, nbes
     bm(nc)=deriv2_7pts(j1(1,nc),ik,grid%r(ik),grid%dx)
     do n=1,ik
        gi(n)=j1(n,nc)*grid%r(n)**ldip
     enddo
     cm(nc)=int_0_inf_dr(gi,grid,ik,nst)
  enddo
  !
  !    solve the linear system to find the coefficients
  !
  gam=(bm(3)-bm(1))/(bm(2)-bm(1))
  delta=(f2ae-bm(1))/(bm(2)-bm(1))
  !   
  xc(3)= (dip-cm(1)+delta*(cm(1)-cm(2)))/(gam*(cm(1)-cm(2))+cm(3)-cm(1))
  xc(2)=-xc(3)*gam+delta
  xc(1)=1.0_dp-xc(2)-xc(3)
  !
  !  Set the function for r<=r(ik)
  !
  do n=1,ik
     phi_out(n)= xc(1)*j1(n,1) + xc(2)*j1(n,2) + xc(3)*j1(n,3)
  enddo
  !
  !  for r > r(ik) the function does not change
  !
  do n=ik+1,grid%mesh
     phi_out(n)= chir(n)
  enddo
  ecutrho=2.0_dp*xc(6)**2

  return
end subroutine compute_q_3bess
