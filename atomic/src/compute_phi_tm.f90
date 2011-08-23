!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine compute_phi_tm(lam,ik,chir,phi_out,iflag,xc,e,els_in)
  !--------------------------------------------------------------------------
  !
  !     This routine computes the phi functions by pseudizing the
  !     all_electron chi functions. In input it receives, the point
  !     ik where the cut is done, the angular momentum lam,
  !     and the all electron wavefunction
  !
  !
  !      
  use io_global, only : stdout
  use kinds, only : DP
  use constants, only: pi
  use radial_grids, only: ndmx
  use ld1inc, only: grid, vpot
  implicit none

  integer ::    &
       lam,  &   ! input: the angular momentum
       ik,   &   ! input: the point corresponding to rc
       iflag     ! input: if 1 print the message

  real(DP) :: &
       e         ! input: the energy of the level

  real(DP) ::       &
       xc(8),       &    ! output: parameters of the fit
       chir(ndmx),   &    ! input: the all-electron function
       phi_out(ndmx)      ! output: the phi function

  character(len=2) :: els_in  ! input: the label of the state
  !
  real(DP) :: &
       faenor    ! the norm of the function

  integer ::    &
       n, nst

  real(DP) :: &
       gi(ndmx), &
       cn(6), c2

  real(DP), external :: deriv_7pts, deriv2_7pts, int_0_inf_dr, pr
  !
  !
  !   compute the norm of the all-electron function
  !
  nst=(lam+1)*2
  do n=1,ik+1
     gi(n)=chir(n)**2  
  enddo
  faenor=int_0_inf_dr(gi,grid,ik,nst)
  !
  !
  ! TM: the pseudo-wavefunction is written as polynomial times exponential
  !
  call find_coefficients &
          (ik, chir, e, grid%r, grid%dx, faenor, vpot, cn, c2, lam, grid%mesh)
  do n=1,ik
     phi_out(n)=sign(grid%r(n)**(lam+1)*exp(pr(cn,c2,grid%r(n))),chir(ik))
  end do
  xc(1:6) = cn(:)
  xc(7) = c2
  !
  !      for r > r(ik) the pseudo and all-electron psi(r) coincide
  !
  do n=ik+1,grid%mesh
     phi_out(n)= chir(n)
  enddo
  !
  !      check for the norm of the pseudo wavefunction
  !
!  do n=1,ik
!     gi(n)=phi_out(n)**2
!  enddo
!  psnor=int_0_inf_dr(gi,grid,ik,nst)
! write(stdout,'(5x," AE norm = ",f6.3,"  PS norm = ",f6.3)') faenor, psnor

  if (iflag == 1) then
     write(stdout,120) els_in, grid%r(ik)
120  format (/ /5x, ' Wfc  ',a3,'  rcut=',f6.3, &
       '  Using Troullier-Martins method ')
  end if

  return
end subroutine compute_phi_tm
