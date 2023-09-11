
!
! Copyright (C) 2005-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE tpssmeta(nnr, nspin,grho,rho,kedtau,etxc)
  !     ===================
  !--------------------------------------------------------------------
  !! Calculates XC energy and potential for metaGGA case.
  !
  use kinds,   only: dp
  use xc_lib,  only: xclib_set_threshold, xc_metagcx
  !
  IMPLICIT NONE
  !
  ! input
  integer :: nspin, nnr
  real(dp) :: grho(3,nnr,nspin), rho(nnr,nspin), kedtau(nnr,nspin)
  ! output: excrho: exc * rho ;  E_xc = \int excrho(r) d_r
  ! output: rhor:   contains the exchange-correlation potential
  real(dp) :: etxc
  REAL(dp) :: zeta, rh, grh2
  INTEGER :: k, ipol, is
  REAL(dp), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10
  !
  etxc = 0.d0
  ! calculate the gradient of rho+rho_core in real space
  !
  call exch_corr_meta() !HK/MCA
  !
  RETURN
  !
  !
 contains
  !
  !
  subroutine  exch_corr_meta()
    !
    implicit none
    !
    INTEGER :: np
    REAL(DP), ALLOCATABLE :: sx(:), v1x(:,:), v2x(:,:), v3x(:,:), &
                             sc(:), v1c(:,:), v2c(:,:,:), v3c(:,:)
    !
    np=1
    if (nspin==2) np=3
    !
    allocate( sx(nnr), v1x(nnr,nspin), v2x(nnr,nspin),    v3x(nnr,nspin) )
    allocate( sc(nnr), v1c(nnr,nspin), v2c(np,nnr,nspin), v3c(nnr,nspin) )
    !
    if (nspin==1) then
      !
      call xclib_set_threshold( 'mgga', epsr, epsg, epsr )
      !
      call xc_metagcx( nnr, 1, np, rho, grho, kedtau, sx, sc, &
                       v1x, v2x, v3x, v1c, v2c, v3c )
      !
      rho(:,1) = v1x(:,1) + v1c(:,1)
      kedtau(:,1) = (v3x(:,1) + v3c(:,1)) * 0.5d0
      ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
      do ipol = 1, 3  
         grho(ipol,:,1) = (v2x(:,1) + v2c(1,:,1))*grho(ipol,:,1) 
      enddo
      !
    else
      !
      call xclib_set_threshold( 'mgga', epsr )
      !
      call xc_metagcx( nnr, 2, np, rho, grho, kedtau, sx, sc, &
                       v1x, v2x, v3x, v1c, v2c, v3c )
      !
      rho(:,1) = v1x(:,1) + v1c(:,1)
      rho(:,2) = v1x(:,2) + v1c(:,2)
      !
      ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
      !
      do ipol = 1, 3  
        grho(ipol,:,1) = (v2x(:,1)*grho(ipol,:,1) + v2c(ipol,:,1))
        grho(ipol,:,2) = (v2x(:,2)*grho(ipol,:,2) + v2c(ipol,:,2)) 
      enddo
      !
      kedtau(:,1) = (v3x(:,1) + v3c(:,1)) * 0.5d0
      kedtau(:,2) = (v3x(:,2) + v3c(:,2)) * 0.5d0
      !
    endif
    etxc = etxc + SUM( sx(:) + sc(:) )
    !
    deallocate( sx, v1x, v2x, v3x )
    deallocate( sc, v1c, v2c, v3c )
    !
    !
    return
    !
  end subroutine exch_corr_meta
  !
END SUBROUTINE tpssmeta

!-----------------------------------------------------------------------
