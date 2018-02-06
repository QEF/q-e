
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
  use kinds, only: dp
  use funct, only: tau_xc, tau_xc_spin, tau_xc_array, tau_xc_array_spin, get_meta
  IMPLICIT NONE
  !
  ! input
  integer nspin , nnr
  real(dp)  grho(3,nnr,nspin), rho(nnr,nspin),kedtau(nnr,nspin)
  ! output: excrho: exc * rho ;  E_xc = \int excrho(r) d_r
  ! output: rhor:   contains the exchange-correlation potential
  real(dp)  etxc
  REAL(dp) :: zeta, rh, grh2
  INTEGER :: k, ipol, is
  REAL(dp), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10
  INTEGER :: imeta
  etxc = 0.d0
  ! calculate the gradient of rho+rho_core in real space
  imeta = get_meta()
  if (imeta.eq.5.or.imeta.eq.6.or.imeta.eq.7) then
    call exch_corr_meta_array_mode() !HK/MCA: currently only implmented for SCAN
  else
    call exch_corr_meta_scalar_mode() !HK/MCA: compatibility for the original implementation
  end if
  RETURN
contains

  subroutine  exch_corr_meta_array_mode()
    implicit none
    real(dp)  :: grho_(3,nnr,nspin) !MCA/HK : store grho only in nspin=2
    REAL(dp)  :: arho(nnr), segno(nnr), grho2 (nnr),                 &
      & sx(nnr), sc(nnr),                                   &
      & v1x(nnr,nspin), v2x(nnr,nspin*2-1), v3x(nnr,nspin), & !MCA/HK
      & v1c(nnr,nspin), v2c(nnr,nspin*2-1), v3c(nnr,nspin)    !MCA/HK
    IF (nspin == 1) THEN
      !
      !$omp parallel do
      do k = 1, nnr
        !
        grho2(k) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
        arho(k)  = ABS (rho (k,1) )
        segno(k) = SIGN (1.d0, rho (k,1) )
        !
      end do !k
      !$omp end parallel do
      !
      CALL tau_xc_array (nnr,arho,grho2,kedtau,sx,sc,v1x,v2x,v3x,v1c,v2c,v3c)
      !
      ! store potentials
      !
      rho (:, 1)  = ( v1x(:,1) + v1c(:,1) )
      kedtau(:,1) = ( v3x(:,1) + v3c(:,1) ) *0.5_dp
      !
      ! v2 contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
      !
      DO ipol = 1, 3  
        grho(ipol,:,1) =  ( v2x(:,1) + v2c(:,1) )*grho (ipol,:,1) 
      ENDDO
      !
    ELSE
      !
      !MCA/HK: only SCAN is available 
      CALL tau_xc_array_spin (nnr, rho, grho, kedtau, sx, sc, v1x, v2x, &
        & v3x, v1c, v2c, v3c)
      !
      ! MCA/HK : store grho to compute v2x cross terms
      !
      grho_ = grho
      !
      DO is = 1,nspin
        !
        rho(:, is) =  v1x(:,is) + v1c(:,is)
        !
        DO ipol = 1, 3  !MCA/HK: second line is the cross term
          grho(ipol,:,is) = ( v2x(:,2*is-1) + v2c(:,2*is-1) ) * grho(ipol,:,is) &
            & + 0.5_dp * ( v2x(:,2) + v2c(:,2) ) * grho_(ipol,:,MOD(is,2)+1) 
        ENDDO
        !
        kedtau(:,is)=  ( v3x(:,is) + v3c(:,is) ) *0.5d0
        !
        segno = 1.0 !MCA: not the most efficient way
        !
      ENDDO
      !
    ENDIF
    !
    ! compute exc energy contribution from the current process
    !
    etxc = 0.0_dp
    !$omp parallel do reduction(+:etxc)
    do k = 1, nnr
      etxc = etxc +  (sx(k) + sc(k)) * segno(k)
    end do !k
    !$omp end parallel do
    return
  end subroutine exch_corr_meta_array_mode

  subroutine  exch_corr_meta_scalar_mode()
    implicit none
    REAL(dp) :: grho2 (2), sx, sc, v1x, v2x, v3x,v1c, v2c, v3c, &
      v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw ,v2cup(3),v2cdw(3), &
      v3xup, v3xdw,grhoup(3),grhodw(3),v3cup, v3cdw, segno, arho, atau
    DO k = 1, nnr
      DO is = 1, nspin
        grho2 (is) = grho(1,k, is)**2 + grho(2,k,is)**2 + grho(3,k,is)**2
      ENDDO
      IF (nspin == 1) THEN
        !
        !    This is the spin-unpolarised case
        !
        arho = ABS (rho (k, 1) )
        segno = SIGN (1.d0, rho (k, 1) )
        atau = kedtau(k,1)
        IF (arho.GT.epsr.AND.grho2 (1) .GT.epsg.AND.ABS(atau).GT.epsr) THEN
          CALL tau_xc (arho, grho2(1), atau, sx, sc, &
            v1x, v2x, v3x, v1c, v2c, v3c)
          rho (k, 1) =  (v1x + v1c )
          kedtau(k,1)=  (v3x + v3c) *0.5d0
          ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
          DO ipol = 1, 3  
            grho(ipol,k,1) =  (v2x + v2c)*grho (ipol,k,1) 
          ENDDO
          etxc = etxc +  (sx + sc) * segno 
        ELSE  
          DO ipol = 1, 3  
            grho (ipol,k,1) = 0.d0  
          ENDDO
          kedtau(k,1)=0.d0
        ENDIF
      ELSE
        !
        !    spin-polarised case
        !
        !CALL tpsscx_spin(rho (k, 1), rho (k, 2), grho2 (1), grho2 (2), &
        !     kedtau(k,1),kedtau(k,2),sx, &
        !     v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
        rh = rho (k, 1) + rho (k, 2)
        IF (rh.GT.epsr) THEN
          !zeta = (rho (k, 1) - rho (k, 2) ) / rh
          DO ipol=1,3
            grhoup(ipol)=grho(ipol,k,1)
            grhodw(ipol)=grho(ipol,k,2)
          END DO
          ! atau=kedtau(k,1)+kedtau(k,2)
          call tau_xc_spin (rho(k,1), rho(k,2), grhoup, grhodw, &
            kedtau(k,1), kedtau(k,2), sx, sc, v1xup, v1xdw, v2xup, v2xdw, &
            v3xup, v3xdw, v1cup, v1cdw, v2cup, v2cdw,&
            v3cup, v3cdw)
          !CALL tpsscc_spin(rh,zeta,grhoup,grhodw, &
          !     atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3c)
        ELSE
          sx = 0.d0  
          sc = 0.d0  
          v1xup = 0.d0  
          v1xdw = 0.d0  
          v2xup=0.d0
          v2xdw=0.d0
          v3xup=0.d0
          v3xdw=0.d0
          v1cup = 0.d0  
          v1cdw = 0.d0  
          v2cup=0.d0
          v2cdw=0.d0
          v3cup=0.d0
          v3cdw=0.d0
          !
        ENDIF
        !
        ! first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        rho(k, 1) =  (v1xup + v1cup)
        rho(k, 2) =  (v1xdw + v1cdw) 
        !
        ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        DO ipol = 1, 3  
          grho(ipol,k,1) = (v2xup*grho(ipol,k,1) + v2cup(ipol))
          grho(ipol,k,2) = (v2xdw*grho(ipol,k,2) + v2cdw(ipol)) 
        ENDDO
        kedtau(k,1)=  (v3xup + v3cup) *0.5d0
        kedtau(k,2)=  (v3xdw + v3cdw) *0.5d0
        etxc = etxc +  (sx + sc)
      ENDIF
    ENDDO
    return
  end subroutine exch_corr_meta_scalar_mode
END SUBROUTINE tpssmeta

!-----------------------------------------------------------------------
