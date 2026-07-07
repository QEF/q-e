!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE v_of_rho( rho, rho_core, rhog_core, tau_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, v )
  !----------------------------------------------------------------------------
  !! This routine computes the Hartree and Exchange and Correlation
  !! potential and energies which corresponds to a given charge density
  !! The XC potential is computed in real space, while the
  !! Hartree potential is computed in reciprocal space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE ions_base,        ONLY : nat, tau
  USE ldaU,             ONLY : lda_plus_u
  USE xc_lib,           ONLY : xclib_dft_is
  USE scf,              ONLY : scf_type
  USE cell_base,        ONLY : alat
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : ts_vdw, mbd_vdw, sic
  USE tsvdw_module,     ONLY : tsvdw_calculate, UtsvdW
  USE libmbd_interface, ONLY : mbd_interface
  USE sic_mod,          ONLY : add_vsic
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  TYPE(scf_type), INTENT(INOUT) :: v
  !! the scf (Hxc) potential
  !=================> NB: NOTE that in F90 derived data type must be INOUT and 
  !=================> not just OUT because otherwise their allocatable or pointer
  !=================> components are NOT defined 
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(IN) :: tau_core(dfftp%nnr)
  !! the kinetic energy density of the core charge
  REAL(DP), INTENT(OUT) :: vtxc
  !! the integral V_xc * rho
  REAL(DP), INTENT(OUT) :: etxc
  !! the E_xc energy
  REAL(DP), INTENT(OUT) :: ehart
  !! the hartree energy
  REAL(DP), INTENT(OUT) :: eth
  !! the hubbard energy
  REAL(DP), INTENT(OUT) :: charge
  !! the integral of the charge
  REAL(DP), INTENT(INOUT) :: etotefield
  !! electric field energy - inout due to the screwed logic of add_efield
  !
  INTEGER :: is, ir
  !
  CALL start_clock( 'v_of_rho' )
  !
  ! ... calculate exchange-correlation potential
  !
  !$acc data copyin(rho,v)
  !$acc data copyin(rho%of_r,rho%of_g,rho_core,rhog_core) copyout(v%of_r,v%kin_r)
  !
  IF (xclib_dft_is('meta')) THEN
     CALL v_xc_meta( rho, rho_core, rhog_core, tau_core, etxc, vtxc, v%of_r, v%kin_r )
  ELSE
     CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  ENDIF
  !
  ! ... add a magnetic field  (if any)
  !
  CALL add_bfield( v%of_r, rho%of_r )
  !
  ! ... calculate hartree potential
  !
  CALL v_h( rho%of_g(:,1), ehart, charge, v%of_r )
  !
  !$acc end data
  !$acc end data
  !
  ! ... DFT+U(+V): build up (extended) Hubbard potential 
  !
  IF (lda_plus_u) CALL v_hubbard ( noncolin, rho, v, eth )
  !
  ! ... add an electric field
  ! 
  DO is = 1, nspin_lsda
     CALL add_efield(v%of_r(1,is), etotefield, rho%of_r(:,1), .false. )
  END DO
  !
  ! ... add Tkatchenko-Scheffler potential (factor 2: Ha -> Ry)
  !
  IF (ts_vdw .or. mbd_vdw) THEN
     CALL tsvdw_calculate(tau*alat,rho%of_r(:,1))
     DO is = 1, nspin_lsda
        DO ir=1,dfftp%nnr
           v%of_r(ir,is)=v%of_r(ir,is)+2.0d0*UtsvdW(ir)
        END DO
     END DO
  END IF
  !
  IF (mbd_vdw) THEN
    call mbd_interface() ! self-consistent but only up to TS level
  END IF
  !
  IF (sic) CALL add_vsic(rho, rho_core, rhog_core, tau_core, v)
  !
  CALL stop_clock( 'v_of_rho' )
  !
  RETURN
  !
END SUBROUTINE v_of_rho
!
!
!----------------------------------------------------------------------------
SUBROUTINE v_xc_meta( rho, rho_core, rhog_core, tau_core, etxc, vtxc, v, kedtaur )
  !----------------------------------------------------------------------------
  !! Exchange-Correlation potential (meta) Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g, ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE funct,            ONLY : dft_is_nonlocc, nlc
  USE scf,              ONLY : scf_type
  USE xc_lib,           ONLY : xc_metagcx, xclib_get_ID
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge in real space
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(IN) :: tau_core(dfftp%nnr)
  !! the kinetic energy density of the core charge
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
  !! V_xc potential
  REAL(DP), INTENT(INOUT) :: kedtaur(dfftp%nnr,nspin)
  !! local K energy density                     
  REAL(DP), INTENT(INOUT) :: vtxc
  !! integral V_xc * rho
  REAL(DP), INTENT(INOUT) :: etxc
  !! E_xc energy
  !
  ! ... local variables
  !
  REAL(DP) :: zeta, rh, sgn_is
  REAL(DP) :: etxc0, vtxc0
  INTEGER  :: k, ipol, is, np, dfftp_nnr
  LOGICAL  :: lda_gga_terms
  !
  REAL(DP), ALLOCATABLE :: ex(:), ec(:), v0(:,:)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v3x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:,:), v3c(:,:)
  !
  REAL(DP) :: fac, rhoneg1, rhoneg2
  REAL(DP), DIMENSION(2) :: grho2
  REAL(DP), DIMENSION(3) :: grhoup, grhodw
  !
  REAL(DP), ALLOCATABLE :: h(:,:,:), dh(:)
  REAL(DP), ALLOCATABLE :: rho_updw(:,:), rhotot(:,:), grho(:,:,:), tau(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:)
  REAL(DP), PARAMETER :: eps12 = 1.0d-12, zero=0._dp
  !
  CALL start_clock( 'v_xc_meta' )
  !
  etxc = zero
  vtxc = zero
  rhoneg1 = zero ; rhoneg2 = zero
  fac = 1.D0 / DBLE( nspin )
  np = 1
  IF (nspin==2) np=3
  dfftp_nnr = dfftp%nnr !to avoid unnecessary copies in acc loop
  !
  !$acc data copyin( rho ) copyout( kedtaur, v )
  !
  ALLOCATE( grho(3,dfftp%nnr,nspin) )
  ALLOCATE( h(3,dfftp%nnr,nspin) )
  ALLOCATE( rhogsum(ngm), tau(dfftp%nnr,nspin) )
  !$acc data create( tau, grho, h )
  !
  ALLOCATE( ex(dfftp%nnr), ec(dfftp%nnr) )
  ALLOCATE( v1x(dfftp%nnr,nspin), v2x(dfftp%nnr,nspin)   , v3x(dfftp%nnr,nspin) )
  ALLOCATE( v1c(dfftp%nnr,nspin), v2c(np,dfftp%nnr,nspin), v3c(dfftp%nnr,nspin) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  ! ... in LSDA case rhogsum is in (up,down) format
  !
  !$acc data create( rhogsum ) copyin( rhog_core, rho%of_g, rho%kin_r )
  DO is = 1, nspin
     !
     sgn_is = (-1.d0)**(is+1)
     !
     !$acc parallel loop present(rho)
     DO k = 1, ngm
       rhogsum(k) = fac*rhog_core(k) + ( rho%of_g(k,1) + sgn_is*rho%of_g(k,nspin) )*0.5D0
     ENDDO
     !
     CALL fft_gradient_g2r( dfftp, rhogsum, g, grho(:,:,is) ) 
     !
  ENDDO
  !
  !$acc parallel loop collapse(2) present(rho)
  DO is = 1, nspin
    DO k = 1, dfftp_nnr
      tau(k,is) = rho%kin_r(k,is)/e2 + fac*tau_core(k)
    ENDDO
  ENDDO
  ! Not sure if the above expression for the kinetic energy density is correct, needs checking...
  ! rho_kin%r(k,is) is in up/down format while rho%of_g(k,is) is in total/magnetization format.

  !
  !$acc end data
  DEALLOCATE( rhogsum )
  !
  !$acc data copyin( rho%of_r )
  !$acc data create( ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  IF (nspin == 1) THEN
    !
    ALLOCATE( rhotot(dfftp%nnr, 1) )
    !$acc data create( rhotot )
    !$acc parallel loop present( rho_core, rho, rho%of_r )
    DO k = 1, dfftp_nnr
      rhotot(k,1) = rho%of_r(k,1) + rho_core(k)
    ENDDO
    !
    CALL xc_metagcx( dfftp_nnr, 1, np, rhotot, grho, tau, ex, ec, &
                     v1x, v2x, v3x, v1c, v2c, v3c, gpu_args_=.TRUE. )
    !
    !$acc parallel loop reduction(+:etxc,vtxc,rhoneg1,rhoneg2) present(rho)
    DO k = 1, dfftp_nnr
       !
       v(k,1) = (v1x(k,1)+v1c(k,1)) * e2
       !
       ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
       DO ipol = 1, 3
         h(ipol,k,1) = (v2x(k,1)+v2c(1,k,1)) * grho(ipol,k,1) * e2
       ENDDO
       !
       kedtaur(k,1) = (v3x(k,1)+v3c(k,1)) * 0.5d0 * e2
       !
       etxc = etxc + (ex(k)+ec(k)) * e2
       vtxc = vtxc + (v1x(k,1)+v1c(k,1)) * e2 * ABS(rho%of_r(k,1))
       !
       IF (rho%of_r(k,1) < zero) rhoneg1 = rhoneg1-rho%of_r(k,1)
       !
    ENDDO
    !
    !$acc end data
    DEALLOCATE( rhotot )
    !
  ELSE
    !
    ALLOCATE( rho_updw(dfftp%nnr,2) )
    !$acc data create( rho_updw )
    !
    !$acc parallel loop present(rho, rho_core)
    DO k = 1, dfftp_nnr
        rho_updw(k,1) = ( rho%of_r(k,1) + rho%of_r(k,2) ) * 0.5d0 + fac*rho_core(k)
        rho_updw(k,2) = ( rho%of_r(k,1) - rho%of_r(k,2) ) * 0.5d0 + fac*rho_core(k)
    ENDDO
    !
    CALL xc_metagcx( dfftp_nnr, 2, np, rho_updw, grho, tau, ex, ec, &
                     v1x, v2x, v3x, v1c, v2c, v3c, gpu_args_=.TRUE. )
    !
    ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
    !
    !$acc parallel loop reduction(+:etxc,vtxc,rhoneg1,rhoneg2)
    DO k = 1, dfftp_nnr
       !
       v(k,1) = (v1x(k,1) + v1c(k,1)) * e2
       v(k,2) = (v1x(k,2) + v1c(k,2)) * e2
       !
       ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
       !
       DO ipol = 1, 3
         h(ipol,k,1) = (v2x(k,1) * grho(ipol,k,1) + v2c(ipol,k,1)) * e2
         h(ipol,k,2) = (v2x(k,2) * grho(ipol,k,2) + v2c(ipol,k,2)) * e2
       ENDDO
       !
       kedtaur(k,1) = (v3x(k,1) + v3c(k,1)) * 0.5d0 * e2
       kedtaur(k,2) = (v3x(k,2) + v3c(k,2)) * 0.5d0 * e2
       !
       etxc = etxc + (ex(k)+ec(k)) * e2
       vtxc = vtxc + (v1x(k,1)+v1c(k,1)) * ABS(rho_updw(k,1) - fac*rho_core(k)) * e2 + &
                     (v1x(k,2)+v1c(k,2)) * ABS(rho_updw(k,2) - fac*rho_core(k)) * e2
       !
       IF ( rho_updw(k,1) < 0.d0 ) rhoneg1 = rhoneg1 - rho_updw(k,1)
       IF ( rho_updw(k,2) < 0.d0 ) rhoneg2 = rhoneg2 - rho_updw(k,2)
       !
    ENDDO
    !
    !$acc end data
    DEALLOCATE( rho_updw )
    !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( ex, ec )
  DEALLOCATE( v1x, v2x, v3x )
  DEALLOCATE( v1c, v2c, v3c )
  !
  ALLOCATE( dh( dfftp%nnr ) )
  !$acc data create( dh )
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin
     CALL fft_graddot( dfftp, h(1,1,is), g, dh )
     !
     sgn_is = (-1.d0)**(is+1)
     !
     !$acc parallel loop reduction(+:vtxc) present(rho)
     DO k = 1, dfftp_nnr
       v(k,is) = v(k,is) - dh(k)
       vtxc = vtxc - dh(k) * ( rho%of_r(k,1) + sgn_is*rho%of_r(k,nspin) )*0.5D0
     ENDDO
  ENDDO
  !
  !$acc end data
  DEALLOCATE( dh )
  !
  !$acc end data
  !$acc end data
  !$acc end data
  !
  CALL mp_sum( rhoneg1, intra_bgrp_comm )
  CALL mp_sum( rhoneg2, intra_bgrp_comm )
  !
  rhoneg1 = rhoneg1 * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  rhoneg2 = rhoneg2 * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ((rhoneg1 > eps8) .OR. (rhoneg2 > eps8)) THEN
    write (stdout, '(/,5x, "negative rho (up,down): ", 2es10.3)') rhoneg1, rhoneg2
  ENDIF
  !
  vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) 
  etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
  !
  CALL mp_sum( vtxc, intra_bgrp_comm )
  CALL mp_sum( etxc, intra_bgrp_comm )
  !
  !
  ! ... calculate and add LDA+GGA terms separately, if needed (not standard)
  !
  lda_gga_terms = (xclib_get_ID('LDA','EXCH') + xclib_get_ID('LDA','CORR') + &
                   xclib_get_ID('GGA','EXCH') + xclib_get_ID('GGA','CORR')) /= 0
  !
  IF ( lda_gga_terms ) THEN
    ALLOCATE(v0(dfftp%nnr,nspin))
    !
    CALL v_xc( rho, rho_core, rhog_core, etxc0, vtxc0, v0 )
    !
    etxc = etxc + etxc0
    vtxc = vtxc + vtxc0
    v = v + v0
    !
    DEALLOCATE(v0)
  ENDIF
  !
  DEALLOCATE( tau, grho )
  DEALLOCATE( h )
  !
  CALL stop_clock( 'v_xc_meta' )
  !
  RETURN
  !
END SUBROUTINE v_xc_meta
!
!------------------------------------------------------------------------------
SUBROUTINE v_xc( rho, rho_core, rhog_core, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !! Exchange-Correlation potential from charge density.
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE noncollin_module, ONLY : domag
  USE funct,            ONLY : nlc, dft_is_nonlocc
  USE scf,              ONLY : scf_type
  USE xc_lib,           ONLY : xc
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin)
  !! \(V_{xc}\) potential
  REAL(DP), INTENT(OUT) :: vtxc
  !! integral \(V_{xc}\cdot\text{rho}\)
  REAL(DP), INTENT(OUT) :: etxc
  !! \(E_{xc}\) energy
  !
  ! ... local variables
  !
  REAL(DP) :: rhoneg(2), vs
  REAL(DP) :: rhoup2, rhodw2, rhoneg1, rhoneg2
  REAL(DP) :: arho, amag, vtxc24
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: vx(:,:), vc(:,:)
  ! In order:
    ! the absolute value of the total charge
    ! the absolute value of the magnetization
    ! zeta = amag / arhox
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir, ipol, dfftp_nnr
    ! counter on mesh points
    ! counter on polarization components
    ! number of mesh points (=dfftp%nnr)
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag    = 1.D-20
  !
  CALL start_clock( 'v_xc' )
  !
  dfftp_nnr = dfftp%nnr !to avoid unnecessary copies in acc loop
  !
  etxc = 0.D0 ;  rhoneg1 = 0.D0
  vtxc = 0.D0 ;  rhoneg2 = 0.D0
  !
  !$acc data copyin( rho_core, rhog_core, rho ) copyout( v )
  !$acc data copyin( rho%of_r, rho%of_g )
  !
  ALLOCATE( ex(dfftp%nnr), vx(dfftp%nnr,nspin) )
  ALLOCATE( ec(dfftp%nnr), vc(dfftp%nnr,nspin) )
  !$acc data create( ex, ec, vx, vc )
  !
  !$acc parallel loop
  DO ir = 1, dfftp_nnr
    rho%of_r(ir,1) = rho%of_r(ir,1) + rho_core(ir)
  ENDDO
  !
  IF ( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
     ! ... spin-unpolarized case
     !
     CALL xc( dfftp_nnr, 1, 1, rho%of_r, ex, ec, vx, vc, gpu_args_=.TRUE. )
     !
     !$acc parallel loop reduction(+:etxc,vtxc,rhoneg1) present(rho)
     DO ir = 1, dfftp_nnr
        v(ir,1) = e2*( vx(ir,1) + vc(ir,1) )
        etxc = etxc + e2*( ex(ir) + ec(ir) )*rho%of_r(ir,1)
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        vtxc = vtxc + v(ir,1)*rho%of_r(ir,1)
        IF (rho%of_r(ir,1) < 0.D0) rhoneg1 = rhoneg1-rho%of_r(ir,1)
     ENDDO
     !
     !
  ELSEIF ( nspin == 2 ) THEN
     ! ... spin-polarized case
     !
     CALL xc( dfftp_nnr, 2, 2, rho%of_r, ex, ec, vx, vc, gpu_args_=.TRUE. )
     !
     !$acc parallel loop reduction(+:etxc,vtxc,rhoneg1,rhoneg2) &
     !$acc&              present(rho)
     DO ir = 1, dfftp_nnr
        v(ir,1) = e2*( vx(ir,1) + vc(ir,1) )
        v(ir,2) = e2*( vx(ir,2) + vc(ir,2) )
        etxc = etxc + e2*( (ex(ir) + ec(ir))*rho%of_r(ir,1) )
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        vtxc = vtxc + ( ( v(ir,1) + v(ir,2) )*rho%of_r(ir,1) + &
                        ( v(ir,1) - v(ir,2) )*rho%of_r(ir,2) )*0.5d0
        !
        rhoup2 = rho%of_r(ir,1)+rho%of_r(ir,2)
        rhodw2 = rho%of_r(ir,1)-rho%of_r(ir,2)
        IF (rhoup2 < 0.d0) rhoneg1 = rhoneg1 - rhoup2*0.5d0
        IF (rhodw2 < 0.d0) rhoneg2 = rhoneg2 - rhodw2*0.5d0
     ENDDO
     !
   ELSEIF ( nspin == 4 ) THEN
      ! ... noncollinear case
      !
      CALL xc( dfftp_nnr, 4, 2, rho%of_r, ex, ec, vx, vc, gpu_args_=.TRUE. )
      !
      !$acc parallel loop reduction(+:etxc,vtxc,rhoneg1,rhoneg2) present(rho)
      DO ir = 1, dfftp_nnr
         arho = ABS( rho%of_r(ir,1) )
         IF ( arho < vanishing_charge ) THEN
           v(ir,1) = 0.d0 ;  v(ir,2) = 0.d0
           v(ir,3) = 0.d0 ;  v(ir,4) = 0.d0
           CYCLE
         ENDIF
         vs = 0.5D0*( vx(ir,1) + vc(ir,1) - vx(ir,2) - vc(ir,2) )
         v(ir,1) = e2*( 0.5D0*( vx(ir,1) + vc(ir,1) + vx(ir,2) + vc(ir,2) ) )
         !
         amag = SQRT( rho%of_r(ir,2)**2 + rho%of_r(ir,3)**2 + rho%of_r(ir,4)**2 )
         IF ( amag > vanishing_mag ) THEN
            v(ir,2) = e2 * vs * rho%of_r(ir,2) / amag
            v(ir,3) = e2 * vs * rho%of_r(ir,3) / amag
            v(ir,4) = e2 * vs * rho%of_r(ir,4) / amag
            vtxc24 = v(ir,2) * rho%of_r(ir,2) + v(ir,3) * rho%of_r(ir,3) + &
                     v(ir,4) * rho%of_r(ir,4)
         ELSE
            v(ir,2) = 0.d0 ;  v(ir,3) = 0.d0 ;  v(ir,4) = 0.d0
            vtxc24 = 0.d0
         ENDIF
         etxc = etxc + e2*( ex(ir) + ec(ir) ) * arho
         !
         rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
         IF ( rho%of_r(ir,1) < 0.D0 )  rhoneg1 = rhoneg1 - rho%of_r(ir,1)
         IF (   amag / arho  > 1.D0 )  rhoneg2 = rhoneg2 + 1.D0/omega
         vtxc = vtxc + vtxc24 + v(ir,1) * rho%of_r(ir,1)
      ENDDO
      !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( ex, vx )
  DEALLOCATE( ec, vc )
  !
  CALL mp_sum(  rhoneg1 , intra_bgrp_comm )
  CALL mp_sum(  rhoneg2 , intra_bgrp_comm )
  !
  rhoneg1 = rhoneg1 * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  rhoneg2 = rhoneg2 * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( rhoneg1 > eps8 .OR. rhoneg2 > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2ES10.3)') rhoneg1, rhoneg2
  !
  ! ... energy terms, local-density contribution
  !
  vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  ! ... add gradient corrections (if any)
  !
  CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, v )
  !
  ! ... to avoid NaN in some rare cases (see summations in subroutine delta_e)
  IF ( nspin==4 .AND. .NOT.domag ) THEN
     !$acc kernels
     v(:,2:nspin) = 0.D0
     !$acc end kernels
  ENDIF
  !
  ! ... add non local corrections (if any)
  !
  IF ( dft_is_nonlocc() ) THEN
    !$acc update self(v)
    CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
    !$acc update device(v)
  ENDIF
  !
  !$acc end data
  !$acc end data
  !
  CALL mp_sum(  vtxc , intra_bgrp_comm )
  CALL mp_sum(  etxc , intra_bgrp_comm )
  !
  CALL stop_clock( 'v_xc' )
  !
  RETURN
  !
END SUBROUTINE v_xc
!
!----------------------------------------------------------------------------
SUBROUTINE v_h( rhog, ehart, charge, v )
  !----------------------------------------------------------------------------
  !! Hartree potential VH(r) from n(G)
  !
  USE constants,         ONLY : fpi, e2
  USE kinds,             ONLY : DP
  USE fft_base,          ONLY : dfftp
  USE fft_rho,           ONLY : rho_g2r
  USE gvect,             ONLY : ngm, gg, gstart
  USE lsda_mod,          ONLY : nspin
  USE cell_base,         ONLY : omega, tpiba2
  USE control_flags,     ONLY : gamma_only
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  USE esm,               ONLY : do_comp_esm, esm_hartree, esm_bc
  USE Coul_cut_2D,       ONLY : do_cutoff_2D, cutoff_hartree
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: rhog(ngm)
  !! the charge density in reciprocal space
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
  !! Hartree potential
  REAL(DP), INTENT(OUT) :: ehart
  !! Hartree energy
  REAL(DP), INTENT(OUT) :: charge
  !
  !  ... local variables
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:), vh(:)
  REAL(DP)              :: rgtot_re, rgtot_im, eh_corr
  INTEGER               :: is, ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), rgtot(:), vaux(:)
  INTEGER               :: nt
  !
  CALL start_clock( 'v_h' )
  !
  !$acc data copyin(rhog) copy(v)
  !
  ALLOCATE( aux(dfftp%nnr), aux1(2,ngm), vh(dfftp%nnr) )
  !$acc data create(aux,aux1,vh)
  !
  charge = 0.D0
  !
  IF ( gstart == 2 ) THEN
     !
     charge = omega*REAL( rhog(1) )
     !
  ENDIF
  !
  CALL mp_sum( charge, intra_bgrp_comm )
  !
  ! ... calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  IF ( do_comp_esm .AND. ( esm_bc .NE. 'pbc' ) ) THEN
     !
     ! ... calculate modified Hartree potential for ESM
     !
     CALL esm_hartree( rhog, ehart, aux )
     !$acc update device(aux)
     !
  ELSE
     !
     ehart = 0.D0
     !$acc kernels
     aux1(:,:) = 0.D0
     !$acc end kernels
     !
     IF (do_cutoff_2D) THEN  !TS
        CALL cutoff_hartree(rhog(:), aux1, ehart)
     ELSE
#if defined(_OPENACC)
        !$acc parallel loop reduction(+:ehart)
#elif defined(__OPENMP)
        !$omp parallel do private( fac, rgtot_re, rgtot_im ), reduction(+:ehart)
#endif
        DO ig = gstart, ngm
           !
           fac = 1.D0 / gg(ig) 
           !
           rgtot_re = REAL(  rhog(ig) )
           rgtot_im = AIMAG( rhog(ig) )
           !
           ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
           !
           aux1(1,ig) = rgtot_re * fac
           aux1(2,ig) = rgtot_im * fac
           !
        ENDDO
#if defined(__OPENMP)
        !$omp end parallel do
#endif
     ENDIF
     !
     fac = e2 * fpi / tpiba2
     !
     ehart = ehart * fac
     !
     !$acc kernels
     aux1(:,:) = aux1(:,:) * fac
     !$acc end kernels
     !
     IF ( gamma_only ) THEN
        !
        ehart = ehart * omega
        !
     ELSE
        !
        ehart = ehart * 0.5D0 * omega
        !
     ENDIF
     !
     IF (do_comp_mt) THEN
        ALLOCATE( vaux(ngm), rgtot(ngm) )
        !$acc data create(vaux,rgtot)
        !$acc kernels
        rgtot(:) = rhog(:)
        !$acc end kernels
        CALL wg_corr_h( omega, ngm, rgtot, vaux, eh_corr )
        !$acc kernels
        aux1(1,1:ngm) = aux1(1,1:ngm) + REAL( vaux(1:ngm))
        aux1(2,1:ngm) = aux1(2,1:ngm) + AIMAG(vaux(1:ngm))
        !$acc end kernels
        ehart = ehart + eh_corr
        !$acc end data
        DEALLOCATE( rgtot, vaux )
     ENDIF
     !
     CALL mp_sum( ehart, intra_bgrp_comm )
     !
     !$acc kernels
     aux(1:ngm) = CMPLX( aux1(1,1:ngm), aux1(2,1:ngm), KIND=DP )
     !$acc end kernels
     !
  ENDIF
  !
  ! ... transform Hartree potential to real space
  !
  CALL rho_g2r( dfftp, aux, vh )
  !
  ! ... add Hartree potential to the xc potential
  !
  IF ( nspin == 4 ) THEN
     !
     !$acc kernels
     v(:,1) = v(:,1) + vh(:)
     !$acc end kernels
     !
  ELSE
     !
     DO is = 1, nspin
        !
        !$acc kernels
        v(:,is) = v(:,is) + vh(:)
        !$acc end kernels
        !
     ENDDO
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( aux, aux1, vh )
  !
  !$acc end data
  !
  CALL stop_clock( 'v_h' )
  !
  RETURN
  !
END SUBROUTINE v_h
!
!----------------------------------------------------------------------------
SUBROUTINE v_h_without_esm( rhog, ehart, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential VH(r) from n(G), with do_comp_esm = .FALSE.
  !
  USE kinds,    ONLY : DP
  USE fft_base, ONLY : dfftp
  USE gvect,    ONLY : ngm
  USE lsda_mod, ONLY : nspin
  USE esm,      ONLY : do_comp_esm
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)    :: rhog(ngm)
  REAL(DP),    INTENT(INOUT) :: v(dfftp%nnr,nspin)
  REAL(DP),    INTENT(OUT)   :: ehart, charge
  !
  LOGICAL :: do_comp_esm_org
  !
  do_comp_esm_org = do_comp_esm
  do_comp_esm = .FALSE.
  !
  CALL v_h( rhog, ehart, charge, v )
  !
  do_comp_esm = do_comp_esm_org
  !
END SUBROUTINE v_h_without_esm
!----------------------------------------------------------------------------
SUBROUTINE v_h_of_rho_r( rhor, ehart, charge, v )
  !----------------------------------------------------------------------------
  !! Hartree potential VH(r) from a density in R space n(r) 
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft
  USE lsda_mod,        ONLY : nspin
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rhor( dfftp%nnr )
  REAL( DP ), INTENT(INOUT)  :: v( dfftp%nnr )
  REAL( DP ), INTENT(OUT)    :: ehart, charge
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhog( : )
  COMPLEX( DP ), ALLOCATABLE :: aux( : )
  REAL( DP ), ALLOCATABLE :: vaux(:,:)
  INTEGER :: is
  !
  ! ... bring the (unsymmetrized) rho(r) to G-space (use aux as work array)
  !
  ALLOCATE( rhog( dfftp%ngm ) )
  ALLOCATE( aux( dfftp%nnr ) )
  aux = CMPLX(rhor,0.D0,kind=dp)
  CALL fwfft ('Rho', aux, dfftp)
  rhog(:) = aux(dfftp%nl(:))
  DEALLOCATE( aux )
  !
  ! ... compute VH(r) from n(G)
  !
  ALLOCATE( vaux( dfftp%nnr, nspin ) )
  vaux = 0.D0
  CALL v_h( rhog, ehart, charge, vaux )
  v(:) = v(:) + vaux(:,1)
  !
  DEALLOCATE( rhog )
  DEALLOCATE( vaux )
  !
  RETURN
  !
END SUBROUTINE v_h_of_rho_r
!----------------------------------------------------------------------------
SUBROUTINE gradv_h_of_rho_r( rho, gradv )
  !----------------------------------------------------------------------------
  !! Gradient of Hartree potential in R space from a total 
  !! (spinless) density in R space n(r)
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE constants,       ONLY : fpi, e2
  USE control_flags,   ONLY : gamma_only
  USE cell_base,       ONLY : tpiba, omega
  USE gvect,           ONLY : ngm, gg, gstart, g
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rho(dfftp%nnr)
  REAL( DP ), INTENT(OUT)    :: gradv(3, dfftp%nnr)
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhoaux(:)
  COMPLEX( DP ), ALLOCATABLE :: gaux(:)
  COMPLEX( DP ), ALLOCATABLE :: rgtot(:), vaux(:)
  REAL( DP )                 :: fac, eh_corr
  INTEGER                    :: ig, ipol
  !
  ! ... Bring rho to G space
  !
  ALLOCATE( rhoaux( dfftp%nnr ) )
  rhoaux( : ) = CMPLX( rho( : ), 0.D0, KIND=dp ) 
  !
  CALL fwfft( 'Rho', rhoaux, dfftp )
  !
  ! ... Compute total potential in G space
  !
  ALLOCATE( gaux( dfftp%nnr ) )
  !
  DO ipol = 1, 3
    !
    gaux(:) = (0.0_dp,0.0_dp)
    !
    DO ig = gstart, ngm
      !
      fac = g(ipol,ig) / gg(ig)
      gaux(dfftp%nl(ig)) = CMPLX(-AIMAG(rhoaux(dfftp%nl(ig))),REAL(rhoaux(dfftp%nl(ig))),KIND=dp)*fac 
      !
    ENDDO
    !
    ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of 
    !  V = e2 * fpi divided by the 2\pi/a factor missing in G  
    !
    fac = e2 * fpi / tpiba
    gaux = gaux * fac 
    !
    ! ...add martyna-tuckerman correction, if needed
    ! 
    IF (do_comp_mt) THEN
       ALLOCATE( vaux( ngm ), rgtot(ngm) )
       rgtot(1:ngm) = rhoaux(dfftp%nl(1:ngm))
       CALL wg_corr_h( omega, ngm, rgtot, vaux, eh_corr )
       DO ig = gstart, ngm
         fac = g(ipol,ig) * tpiba
         gaux(dfftp%nl(ig)) = gaux(dfftp%nl(ig)) + CMPLX(-AIMAG(vaux(ig)),REAL(vaux(ig)),kind=dp)*fac 
       END DO
       DEALLOCATE( rgtot, vaux )
    ENDIF
    !
    IF ( gamma_only ) THEN
      !
      gaux(dfftp%nlm(:)) = &
        CMPLX( REAL( gaux(dfftp%nl(:)) ), -AIMAG( gaux(dfftp%nl(:)) ) ,kind=DP)
       !
    END IF
    !
    ! ... bring back to R-space, (\grad_ipol a)(r) ...
    !
    CALL invfft( 'Rho', gaux, dfftp )
    !
    gradv(ipol,:) = REAL( gaux(:) )
    !
  ENDDO
  !
  DEALLOCATE(gaux)
  !
  DEALLOCATE(rhoaux)
  !
  RETURN
  !
END SUBROUTINE gradv_h_of_rho_r
