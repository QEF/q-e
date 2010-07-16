!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE all_electron(ild,ic)
  !---------------------------------------------------------------
  !
  !  this routine is a driver to an all-electron calculation
  !  with the parameters given in input
  !
  !
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: ndmx
  USE ld1inc, ONLY: isic, grid, zeta, rho, enne, vpot, vxt, enl, &
                     deld, encl, etot, ecxc, evxt, ehrt, epseu, ekin, &
                     vnl, vh, lsd, nspin, nlcc, vdw, nn, ll, oc, nwf, &
                     zed, zval, vxc, exc, excgga, v0, verbosity, &
                     relpert, evel, edar, eso, vsic, vsicnew, vhn1, egc
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ic   ! counter on configurations
  LOGICAL :: ild    ! if true compute log der
  !
  !    compute an initial estimate of the potential
  !
  CALL starting_potential (ndmx, grid%mesh, zval, zed, nwf, oc, nn, ll,&
                           grid%r,enl, v0, vxt, vpot, enne, nspin )
  !
  ! allocate variables for SIC, if needed
  !
  IF (isic /= 0) THEN
     ALLOCATE(vsic(ndmx,nwf), vsicnew(ndmx), vhn1(ndmx), egc(ndmx))
     vsic=0.0_dp
  ENDIF
  !
  !     solve the eigenvalue self-consistent equation
  !
  CALL scf(ic)
  !
  !   compute relativistic corrections to the eigenvalues
  !
  IF ( relpert ) CALL compute_relpert(evel,edar,eso)
  !
  !  compute total energy
  !
  CALL elsd (zed,grid,rho,vxt,vh,vxc,exc,excgga,nwf,nspin,enl,oc,    &
             etot,ekin,encl,ehrt,ecxc,evxt)
  !
  IF (verbosity=='high') CALL elsd_highv(ic)
  !
  !   add sic correction if needed
  !
  IF(isic /= 0) CALL esic
  !
  !   print results
  !
  CALL write_results
  !
  !  compute logarithmic derivative
  !
  IF (deld > 0.0_DP .and. ild) CALL lderiv
  !
  ! compute C6 coefficient if required
  !
  IF (vdw) THEN
     CALL c6_tfvw ( grid%mesh, zed, grid, rho(1,1) )
     CALL c6_dft  ( grid%mesh, zed, grid )
  ENDIF
  !
  IF (isic /= 0) THEN
     DEALLOCATE(egc, vhn1, vsicnew, vsic)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE all_electron
