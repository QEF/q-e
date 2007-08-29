!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine all_electron(ild)
  !---------------------------------------------------------------
  !
  !  this routine is a driver to an all-electron calculation
  !  with the parameters given in input
  !
  !
  use kinds, only : DP
  use radial_grids, only: ndmx
  use ld1inc, only: isic, grid, zeta, rho, enne, vpot, vxt, enl, &
                     deld, encl, etot, ecxc, evxt, ehrt, epseu, ekin, &
                     vnl, vh, lsd, nspin, nlcc, vdw, nn, ll, oc, nwf, &
                     zed, zval, vxc, exc, excgga, v0, verbosity
  implicit none

  logical :: ild    ! if true compute log der
  !
  !    compute an initial estimate of the potential
  !
  call starting_potential(ndmx,grid%mesh,zval,zed,nwf,oc,nn,ll,grid%r,&
       enl,v0,vxt,vpot,enne,nspin)
  !
  !     solve the eigenvalue self-consistent equation
  !
  call scf
  !
  !  compute total energy
  !
  call elsd (zed,grid,rho,vxt,vh,vxc,exc,excgga,nwf,nspin,enl,oc,    &
             etot,ekin,encl,ehrt,ecxc,evxt)
  !
  if (verbosity=='high') call elsd_highv( )
  !
  !   add sic correction if needed
  !
  if(isic /= 0) call esic  
  !
  !   print results
  !
  call write_results 
  !
  !  compute logarithmic derivative
  !
  if (deld > 0.0_DP .and. ild) call lderiv
  !
  ! compute C6 coefficient if required
  !
  if (vdw) then
     call c6_tfvw ( grid%mesh, zed, grid, rho(1,1) )
     call c6_dft  ( grid%mesh, zed, grid )
  end if
  !
  return
  !
end subroutine all_electron
