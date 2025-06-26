!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE v_of_rho( rho, rho_core, rhog_core, &
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
  USE ldaU,             ONLY : lda_plus_u, lda_plus_u_kind, ldmx_b, &
                               nsg, v_nsg, Hubbard_l, Hubbard_lmax, apply_U, orbital_resolved 
  USE xc_lib,           ONLY : xclib_dft_is
  USE scf,              ONLY : scf_type
  USE cell_base,        ONLY : alat
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : ts_vdw, mbd_vdw, sic
  USE tsvdw_module,     ONLY : tsvdw_calculate, UtsvdW
  USE libmbd_interface, ONLY : mbd_interface
  USE sic_mod,          ONLY : add_vsic
#if defined (__OSCDFT)
  USE plugin_flags,     ONLY : use_oscdft
  USE oscdft_base,      ONLY : oscdft_ctx
  USE oscdft_functions, ONLY : oscdft_v_constraint
#endif  
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
  REAL(DP) :: eth1
  !! the hubbard energy coming from the background states
  REAL(DP), INTENT(INOUT) :: etotefield
  !! electric field energy - inout due to the screwed logic of add_efield
  !
  INTEGER :: is, ir
  !
  CALL start_clock( 'v_of_rho' )
  !
  ! ... calculate exchange-correlation potential
  !
  IF (xclib_dft_is('meta')) THEN
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
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
  ! ... DFT+U(+V): build up (extended) Hubbard potential 
  !
  IF (lda_plus_u) THEN
     !
     IF (lda_plus_u_kind == 0) THEN
        !
        ! DFT+U (simplified)
        !
        IF (noncolin) THEN
           IF (orbital_resolved) THEN        
              CALL v_hubbard_resolved_nc (rho%ns_nc, v%ns_nc, eth)
           ELSE
              CALL v_hubbard_nc (rho%ns_nc, v%ns_nc, eth)
           ENDIF
        ELSE
           IF (orbital_resolved) THEN
              CALL v_hubbard_resolved(rho%ns, v%ns, eth)
           ELSE
              CALL v_hubbard (rho%ns, v%ns, eth)
           ENDIF
        ENDIF   
        !
        ! Background
        IF (ldmx_b.GT.0) THEN
           CALL v_hubbard_b (rho%nsb, v%nsb, eth1)
           eth = eth + eth1
        ENDIF
        !
     ELSEIF (lda_plus_u_kind == 1) THEN
        !
        ! DFT+U (full)
        !
        IF (noncolin) THEN
           CALL v_hubbard_full_nc (rho%ns_nc, v%ns_nc, eth)
        ELSE
           CALL v_hubbard_full (rho%ns, v%ns, eth)
        ENDIF
        !
     ELSEIF (lda_plus_u_kind == 2) THEN
        !
        ! DFT+U+V (simplified)
        !
        IF (noncolin) THEN
           CALL v_hubbard_extended_nc (nsg, v_nsg, eth)
        ELSE
           CALL v_hubbard_extended (nsg, v_nsg, eth)
        ENDIF
     ELSE
        !
        CALL errore('v_of_rho', 'Not allowed value of lda_plus_u_kind',1)
        !
     ENDIF
     !
  ENDIF
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (lda_plus_u) THEN
        IF (lda_plus_u_kind == 0) THEN
           CALL oscdft_v_constraint (oscdft_ctx, Hubbard_lmax, Hubbard_l, rho%ns, v%ns, eth)
        ELSEIF (lda_plus_u_kind == 2) THEN
           CALL oscdft_v_constraint_extended (nsg, v_nsg, eth)
        ENDIF
     ENDIF
  ENDIF
#endif
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
  IF (sic) CALL add_vsic(rho, rho_core, rhog_core, v)
  !
  CALL stop_clock( 'v_of_rho' )
  !
  RETURN
  !
END SUBROUTINE v_of_rho
!
!
!----------------------------------------------------------------------------
SUBROUTINE v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
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
  REAL(DP), ALLOCATABLE :: rho_updw(:,:), grho(:,:,:), tau(:,:)
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
      tau(k,is) = rho%kin_r(k,is)/e2
    ENDDO
  ENDDO
  !
  !$acc end data
  DEALLOCATE( rhogsum )
  !
  !$acc data copyin( rho%of_r )
  !$acc data create( ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  IF (nspin == 1) THEN
    !
    CALL xc_metagcx( dfftp_nnr, 1, np, rho%of_r, grho, tau, ex, ec, &
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
  ELSE
    !
    ALLOCATE( rho_updw(dfftp%nnr,2) )
    !$acc data create( rho_updw )
    !
    !$acc parallel loop present(rho)
    DO k = 1, dfftp_nnr  
        rho_updw(k,1) = ( rho%of_r(k,1) + rho%of_r(k,2) ) * 0.5d0
        rho_updw(k,2) = ( rho%of_r(k,1) - rho%of_r(k,2) ) * 0.5d0
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
       vtxc = vtxc + (v1x(k,1)+v1c(k,1)) * ABS(rho_updw(k,1)) * e2 + &
                     (v1x(k,2)+v1c(k,2)) * ABS(rho_updw(k,2)) * e2
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
  !$acc end data
  !$acc end data
  !
  ! ... to avoid NaN in some rare cases (see summations in subroutine delta_e)
  IF ( nspin==4 .AND. .NOT.domag ) v(:,2:nspin) = 0.D0
  !
  ! ... add non local corrections (if any)
  !
  IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
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
  USE Coul_cut_2D,       ONLY : do_cutoff_2D, cutoff_2D, cutoff_hartree  
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
  ALLOCATE( aux(dfftp%nnr), aux1(2,ngm), vh(dfftp%nnr) )
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
     !
  ELSE
     !
     ehart     = 0.D0
     aux1(:,:) = 0.D0
     !
     IF (do_cutoff_2D) THEN  !TS
        CALL cutoff_hartree(rhog(:), aux1, ehart)
     ELSE
!$omp parallel do private( fac, rgtot_re, rgtot_im ), reduction(+:ehart)
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
!$omp end parallel do
     ENDIF
     !
     fac = e2 * fpi / tpiba2
     !
     ehart = ehart * fac
     !
     aux1 = aux1 * fac
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
        rgtot(:) = rhog(:)
        CALL wg_corr_h( omega, ngm, rgtot, vaux, eh_corr )
        aux1(1,1:ngm) = aux1(1,1:ngm) + REAL( vaux(1:ngm))
        aux1(2,1:ngm) = aux1(2,1:ngm) + AIMAG(vaux(1:ngm))
        ehart = ehart + eh_corr
        DEALLOCATE( rgtot, vaux )
     ENDIF
     !
     CALL mp_sum( ehart, intra_bgrp_comm )
     !
     aux(1:ngm) = CMPLX( aux1(1,1:ngm), aux1(2,1:ngm), KIND=DP )
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
     v(:,1) = v(:,1) + vh(:)
     !
  ELSE
     !
     DO is = 1, nspin
        !
        v(:,is) = v(:,is) + vh(:)
        !
     ENDDO
     !
  ENDIF
  !
  DEALLOCATE( aux, aux1, vh )
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
!
!-----------------------------------------------------------------------
SUBROUTINE v_hubbard( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !! DFT+U+J0: B. Himmetoglu et al., Phys. Rev. B 84, 115108 (2011).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_alpha, Hubbard_J0, Hubbard_beta
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity, dfpt_hub
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  COMPLEX(DP) :: check
  !
  !  ... local variables
  !
  REAL(DP) :: effU, sgn(2) 
  INTEGER  :: is, isop, na, nt, m1, m2
  !
  eth    = 0.d0
  sgn(1) =  1.d0  
  sgn(2) = -1.d0
  !
  v_hub(:,:,:,:) = 0.d0
  check = (0.0,0.0)
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
        !
        IF (Hubbard_J0(nt) /= 0.d0) THEN
           effU = Hubbard_U(nt) - Hubbard_J0(nt)
        ELSE
           effU = Hubbard_U(nt)
        ENDIF 
        ! 
        DO is = 1, nspin
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              eth = eth + ( Hubbard_alpha(nt) + 0.5D0*effU )*ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                                   Hubbard_alpha(nt)  + 0.5D0*effU
               check = check + Hubbard_alpha(nt)  + 0.5D0*effU
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 eth = eth - 0.5D0 * effU * ns(m2,m1,is,na)* ns(m1,m2,is,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                        effU * ns(m2,m1,is,na)
                 check = check -effU * ns(m2,m1,is,na)
              ENDDO
           ENDDO
        ENDDO
        !
     ENDIF
     !
     IF (Hubbard_J0(nt) /= 0.d0 .OR. Hubbard_beta(nt) /= 0.d0) THEN
        !
        DO is = 1, nspin
           isop = 1
           IF ( nspin == 2 .AND. is == 1) isop = 2
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              eth = eth + sgn(is)*Hubbard_beta(nt) * ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + sgn(is)*Hubbard_beta(nt)
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 eth = eth + 0.5D0*Hubbard_J0(nt)*ns(m2,m1,is,na)*ns(m1,m2,isop,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + Hubbard_J0(nt) * &
                                                            ns(m2,m1,isop,na)
              ENDDO
           ENDDO
        ENDDO
        !
     END IF
     !   
  ENDDO
  !
  IF (nspin==1) eth = 2.d0 * eth
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
     !write(stdout,*) "check coll U", check
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
SUBROUTINE v_hubbard_nc( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy in the \textbf{noncollinear formulation}.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !! DFT+U+J0: B. Himmetoglu et al., Phys. Rev. B 84, 115108 (2011).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_alpha, Hubbard_J0, Hubbard_beta
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity, dfpt_hub
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  COMPLEX(DP) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  !
  !  ... local variables
  !
  INTEGER  :: is, is1, na, nt, m1, m2
  !
  eth    = 0.d0
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt) /= 0.d0) THEN
        DO is = 1, nspin
           !
           IF (is == 2) THEN
            is1 = 3
           ELSEIF (is == 3) THEN
            is1 = 2
           ELSE
            is1 = is
           ENDIF
           !
           ! Non spin-flip contribution
           ! (diagonal [spin indexes] occupancy matrices)
           IF (is1 == is) THEN
              !     
              ! diagonal part [spin indexes]     
              DO m1 = 1, 2*Hubbard_l(nt) + 1
                 ! Hubbard energy
                 eth = eth + ( Hubbard_alpha(nt) + 0.5D0*Hubbard_U(nt) )&
                             * ns(m1,m1,is,na)
                 ! Hubbard potential
                 v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + & 
                                      Hubbard_alpha(nt) + 0.5D0*Hubbard_U(nt)
                 ! 
                 ! NON-diagonal part [spin indexes]
                 DO m2 = 1, 2 * Hubbard_l(nt) + 1
                    ! Hubbard energy
                    eth = eth - 0.5D0 * Hubbard_U(nt) * ns(m1,m2,is,na)*ns(m2,m1,is,na)
                    ! Hubbard potential
                    v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                         Hubbard_U(nt) * ns(m2,m1,is,na)
                 ENDDO
              ENDDO  
           !
           ! Spin-flip contribution
           ! (NON-diagonal [spin indexes] occupancy matrices)   
           ELSE
              DO m1 = 1, 2*Hubbard_l(nt) + 1
                 DO m2 = 1, 2 * Hubbard_l(nt) + 1 
                    ! Hubbard energy
                    eth = eth - 0.5D0 * Hubbard_U(nt) &
                          * ns(m1,m2,is,na)*ns(m2,m1,is1,na)
                    ! Hubbard potential
                    v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                         Hubbard_U(nt) * ns(m2,m1,is1,na)
                 ENDDO
              ENDDO         
           ENDIF
        ENDDO  
     ENDIF
     !
  ENDDO
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_nc
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
SUBROUTINE v_hubbard_b (ns, v_hub, eth)
  !-------------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy for background states.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_J0, Hubbard_beta, Hubbard_U2,  &
                                   ldim_back, ldmx_b, Hubbard_alpha_back, &
                                   is_hubbard_back
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity, dfpt_hub
  USE io_global,            ONLY : stdout

  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(ldmx_b,ldmx_b,nspin,nat)
  REAL(DP), INTENT(OUT) :: v_hub(ldmx_b,ldmx_b,nspin,nat)
  REAL(DP), INTENT(OUT) :: eth
  REAL(DP) :: effU
  INTEGER :: is, is1, na, nt, m1, m2, m3, m4
  !
  eth = 0.d0
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_J0(nt).NE.0.d0) &
          CALL errore('v_hubbard_b', 'J0 is not supported in DFT+U with multiple channels per atomic type',1)
     !
     IF (Hubbard_beta(nt).NE.0.d0) &
     CALL errore('v_hubbard_b', 'Hubbard_beta is not supported in DFT+U with multiple channels per atomic type',1) 
     !
     IF (is_hubbard_back(nt)) THEN
        !
        effU = Hubbard_U2(nt)
        !
        DO is = 1, nspin
           !
           DO m1 = 1, ldim_back(nt) 
              !
              eth = eth + ( Hubbard_alpha_back(nt) + 0.5D0 * effU ) * &
                              ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                            ( Hubbard_alpha_back(nt) + 0.5D0 * effU )
              !
              DO m2 = 1, ldim_back(nt)
                 !
                 eth = eth - 0.5D0 * effU * &
                            ns(m2,m1,is,na)* ns(m1,m2,is,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                       effU * ns(m2,m1,is,na)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF (nspin.EQ.1) eth = 2.d0 * eth
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,'(/5x,"HUBBARD BACKGROUND ENERGY = ",f9.4,1x," (Ry)")') eth
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_b
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE v_hubbard_full( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy.
  !! DFT+U(+J) : Formulation by Liechtenstein et al., Phys. Rev. B 52, R5467 (1995).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_J, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  !
  !  ... local variables
  !
  REAL(DP) :: n_tot, n_spin, eth_dc, eth_u, mag2
  INTEGER  :: is, isop, is1, na, nt, m1, m2, m3, m4
  REAL(DP), ALLOCATABLE :: u_matrix(:,:,:,:)
  !
  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )
  !
  eth    = 0.d0
  eth_dc = 0.d0
  eth_u  = 0.d0
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt)/=0.d0) THEN
        !
        ! Initialize U(m1,m2,m3,m4) matrix 
        !
        CALL hubbard_matrix( Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                               Hubbard_J(1,nt), u_matrix )
        !
        ! Total N and M^2 for DC (double counting) term
        !
        n_tot = 0.d0
        !
        DO is = 1, nspin
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              n_tot = n_tot + ns(m1,m1,is,na)
           ENDDO
        ENDDO
        !
        IF (nspin==1) n_tot = 2.d0 * n_tot
        !
        mag2  = 0.d0
        ! 
        IF (nspin==2) THEN
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              mag2 = mag2 + ns(m1,m1,1,na) - ns(m1,m1,2,na)
           ENDDO
        ENDIF
        mag2  = mag2**2
        !
        ! Hubbard energy: DC term
        !
        eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*n_tot*(n_tot-1.d0) - &
                                  Hubbard_J(1,nt)*n_tot*(0.5d0*n_tot-1.d0) - &
                                  0.5d0*Hubbard_J(1,nt)*mag2 )
        !
        DO is = 1, nspin
           !
           ! n_spin = up/down N
           !
           n_spin = 0.d0
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              n_spin = n_spin + ns(m1,m1,is,na)
           ENDDO
           !
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              !
              ! Hubbard potential: DC contribution  
              !
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_spin + &
                         0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot
              !
              ! +U contributions 
              !
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 DO m3 = 1, 2 * Hubbard_l(nt) + 1
                    DO m4 = 1, 2 * Hubbard_l(nt) + 1
                       !
                       DO is1 = 1, nspin
                          v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + (MOD(nspin,2)+1) * &
                                               u_matrix(m1,m3,m2,m4) * ns(m3,m4,is1,na)
                       ENDDO
                       !
                       v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                            u_matrix(m1,m3,m4,m2) * ns(m3,m4,is,na)
                       !
                       eth_u = eth_u + 0.5d0*( ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) ) * &
                                       ns(m1,m3,is,na)*ns(m2,m4,is,na)+u_matrix(m1,m2,m3,m4)   * &
                                       ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )
                    ENDDO ! m4
                 ENDDO ! m3
              ENDDO ! m2
              !
           ENDDO ! m1
           !
        ENDDO ! is
        !
     ENDIF
     !
  ENDDO ! na
  ! 
  IF (nspin==1) eth_u = 2.d0 * eth_u
  eth = eth_u - eth_dc
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGIES (dc, U, total) ",3f9.4,1x," (Ry)")') eth_dc, eth_u, eth
  ENDIF
  !
  DEALLOCATE (u_matrix)
  !
  RETURN
  !
END SUBROUTINE v_hubbard_full
!---------------------------------------------------------------

!---------------------------------------------------------------
SUBROUTINE v_hubbard_full_nc( ns, v_hub, eth )
  !-------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy (noncollinear case).
  !! DFT+U(+J) : Formulation by Liechtenstein et al., Phys. Rev. B 52, R5467 (1995).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE noncollin_module,     ONLY : noncolin
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_J, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  COMPLEX(DP) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP) :: eth
  !! Hubbard matrix
  !
  !  ... local variables
  !
  REAL(DP) :: eth_dc, eth_noflip, eth_flip, mx, my, mz, mag2
  INTEGER :: is, is1, js, i, j, na, nt, m1, m2, m3, m4
  COMPLEX(DP) :: n_tot, n_aux
  REAL(DP), ALLOCATABLE :: u_matrix(:,:,:,:)
  !
  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )
  !
  eth        = 0.d0
  eth_dc     = 0.d0  
  eth_noflip = 0.d0  
  eth_flip   = 0.d0  
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat  
     !
     nt = ityp (na)  
     !
     IF (Hubbard_U(nt) /= 0.d0) THEN  
        !
        ! Initialize U(m1,m2,m3,m4) matrix 
        !
        CALL hubbard_matrix( Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                             Hubbard_J(1,nt), u_matrix )
        !
        ! Total N and M^2 for DC (double counting) term
        !
        n_tot = 0.d0
        mx    = 0.d0
        my    = 0.d0
        mz    = 0.d0
        DO m1 = 1, 2*Hubbard_l(nt)+1
          n_tot = n_tot + ns(m1,m1,1,na) + ns(m1,m1,4,na)
          mx = mx + DBLE( ns(m1, m1, 2, na) + ns(m1, m1, 3, na) )
          my = my + 2.d0 * AIMAG( ns(m1, m1, 2, na) )
          mz = mz + DBLE( ns(m1, m1, 1, na) - ns(m1, m1, 4, na) )
        ENDDO  
        mag2 = mx**2 + my**2 + mz**2  
        !
        ! Hubbard energy: DC term
        !
        mx = REAL(n_tot)
        eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*mx*(mx-1.d0) - &
                                  Hubbard_J(1,nt)*mx*(0.5d0*mx-1.d0) - &
                                  0.5d0*Hubbard_J(1,nt)*mag2 )   
        !
        DO is = 1, nspin  
           !
           IF (is == 2) THEN
            is1 = 3
           ELSEIF (is == 3) THEN
            is1 = 2
           ELSE
            is1 = is
           ENDIF
           !
           ! Hubbard energy:
           !
           IF (is1 == is) THEN
             !
             ! Non spin-flip contribution
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                DO m3 = 1, 2*Hubbard_l(nt)+1
                 DO m4 = 1, 2*Hubbard_l(nt)+1
                   eth_noflip = eth_noflip + 0.5d0*(                            &
                              ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) )*  & 
                              ns(m1,m3,is,na)*ns(m2,m4,is,na) +                 &
                      u_matrix(m1,m2,m3,m4)*ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )
                 ENDDO
                ENDDO
              ENDDO
             ENDDO
             !
           ELSE
             ! 
             ! Spin-flip contribution
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
               DO m3 = 1, 2*Hubbard_l(nt)+1
                DO m4 = 1, 2*Hubbard_l(nt)+1
                   eth_flip = eth_flip - 0.5d0*u_matrix(m1,m2,m4,m3)* &
                                     ns(m1,m3,is,na)*ns(m2,m4,is1,na) 
                ENDDO
               ENDDO
              ENDDO
             ENDDO
             !
           ENDIF
           !
           ! Hubbard potential: non spin-flip contribution 
           !
           IF (is1 == is) THEN
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
               DO m3 = 1, 2*Hubbard_l(nt)+1
                DO m4 = 1, 2*Hubbard_l(nt)+1
                  v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                    u_matrix(m1,m3,m2,m4)*( ns(m3,m4,1,na)+ns(m3,m4,4,na) ) 
                ENDDO 
               ENDDO
              ENDDO
             ENDDO
             !
           ENDIF
           !
           ! n_aux = /sum_{i} n_{i,i}^{sigma2, sigma1} for DC term
           !
           n_aux = 0.d0
           DO m1 = 1, 2*Hubbard_l(nt)+1
              n_aux = n_aux + ns(m1,m1,is1,na)  
           ENDDO
           !
           DO m1 = 1, 2*Hubbard_l(nt)+1
             ! 
             ! Hubbard potential: DC contribution  
             !
             v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_aux
             !
             IF (is1 == is) THEN
                v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                     0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot  
             ENDIF
             !
             ! Hubbard potential: spin-flip contribution
             !
             DO m2 = 1, 2*Hubbard_l(nt)+1  
              DO m3 = 1, 2*Hubbard_l(nt)+1
               DO m4 = 1, 2*Hubbard_l(nt)+1
                  v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                             u_matrix(m1,m3,m4,m2) * ns(m3,m4,is1,na) 
               ENDDO 
              ENDDO
             ENDDO
             !
           ENDDO
           !
        ENDDO ! is
        !
     ENDIF
     !
  ENDDO ! na
  !
  eth = eth_noflip + eth_flip - eth_dc
  !
  ! Hubbard energies
  !
  IF ( iverbosity > 0 ) THEN
    WRITE(stdout,*) '--- in v_hubbard ---'
    WRITE(stdout,'("Hub. E (dc, noflip, flip, total) ",4f9.4)') &
                                 eth_dc, eth_noflip, eth_flip, eth 
    WRITE(stdout,*) '-------'
  ENDIF
  !
  DEALLOCATE (u_matrix)
  !
  RETURN
  !
END SUBROUTINE v_hubbard_full_nc
!----------------------------------------------------------------------------

!------------------------------------------------------------------------------------
SUBROUTINE v_hubbard_extended (nsg, v_hub, eth)
  !-----------------------------------------------------------------------------------
  !
  !! Computes extended Hubbard potential and Hubbard energy.
  !! DFT+U+V: Simplified rotationally-invariant formulation by
  !! V.L. Campo Jr and M. Cococcioni, J. Phys.: Condens. Matter 22, 055602 (2010).
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp
  USE ldaU,              ONLY : Hubbard_l, Hubbard_alpha, Hubbard_J0, Hubbard_beta,   &
                                ldim_u, ldmx_tot, max_num_neighbors, at_sc, neighood, &
                                Hubbard_V, Hubbard_alpha_back, is_hubbard, is_hubbard_back
  USE lsda_mod,          ONLY : nspin
  USE control_flags,     ONLY : iverbosity, dfpt_hub
  USE io_global,         ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
  COMPLEX(DP), INTENT(OUT) :: v_hub(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
  REAL(DP),    INTENT(OUT) :: eth
  ! 
  ! Local variables 
  !
  INTEGER :: is, isop, na, na1, na2, nt, nt1, nt2, m1, m2, viz, equiv_na2, i_type
  COMPLEX(DP) :: check, check_en
  INTEGER, EXTERNAL :: type_interaction, find_viz
  !
  eth  = 0.d0
  v_hub(:,:,:,:,:) = (0.d0, 0.d0)
  check = (0.d0, 0.d0)
  check_en = (0.d0, 0.d0)
  !
  DO na1 = 1, nat
     !
     nt1 = ityp(na1)
     !
     IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
        !
        DO is = 1, nspin
           !
           DO viz = 1, neighood(na1)%num_neigh
              !
              na2 = neighood(na1)%neigh(viz)
              equiv_na2 = at_sc(na2)%at
              nt2 = ityp(equiv_na2)
              !
              IF ((is_hubbard(nt2).OR.is_hubbard_back(nt2)) .AND. &
                  (Hubbard_V(na1,na2,1).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,2).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,3).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,4).NE.0.d0) ) THEN
                  !
                  ! For both standard and background states of a center atom
                  DO m1 = 1, ldim_u(nt1)
                     ! For both standard and background states of the neighbor atom
                     DO m2 = 1, ldim_u(nt2)
                        !
                        i_type = type_interaction(na1,m1,equiv_na2,m2)
                        !
                        v_hub(m2,m1,viz,na1,is) = &
                           - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,i_type)
                        check = check  - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,i_type)

                        !
                        eth = eth - nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                    * Hubbard_V(na1,na2,i_type) * 0.5d0
                        check_en = check_en -nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                    * Hubbard_V(na1,na2,i_type) * 0.5d0
                        !
                     ENDDO
                  ENDDO
                  !
                  IF ( na1.EQ.na2 ) THEN
                     !
                     na = find_viz(na1,na1)
                     !
                     ! This is the diagonal term (like in the DFT+U only case)
                     ! 
                     DO m1 = 1, ldim_u(nt1)
                        !
                        i_type = type_interaction(na1,m1,equiv_na2,m1)
                        !
                        v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) &
                                       + Hubbard_V(na1,na1,i_type) * 0.5d0
                        check = check  + Hubbard_V(na1,na1,i_type) * 0.5d0
                        ! 
                        eth = eth + nsg(m1,m1,na,na1,is) &
                                       * Hubbard_V(na1,na1,i_type) * 0.5d0
                        check_en = check_en +nsg(m1,m1,na,na1,is) &
                                       * Hubbard_V(na1,na1,i_type) * 0.5d0
                        !
                     ENDDO
                     !
                     ! Hubbard_J0 (only on-site)
                     !
                     IF ( nspin.EQ.2 .AND. &
                          (Hubbard_J0(nt1).NE.0.d0 .OR. Hubbard_beta(nt1).NE.0.d0) ) THEN
                          !
                          IF (is.EQ.1) THEN
                             isop = 2
                          ELSE
                             isop = 1
                          ENDIF
                          !
                          DO m1 = 1, 2*Hubbard_l(nt1)+1
                             !
                             v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) - &
                                                      Hubbard_J0(nt1) * 0.5d0
                             eth = eth - 0.5d0 * Hubbard_J0(nt1) * nsg(m1,m1,na,na1,is)
                             !
                             IF (is.EQ.1) THEN
                                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + &
                                                         Hubbard_beta(nt1)
                                eth = eth + Hubbard_beta(nt1) * nsg(m1,m1,na,na1,is)
                             ELSE
                                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) - &
                                                         Hubbard_beta(nt1)
                                eth = eth - Hubbard_beta(nt1) * nsg(m1,m1,na,na1,is)
                             ENDIF
                             !
                             DO m2 = 1, 2*Hubbard_l(nt1)+1
                                v_hub(m2,m1,na,na1,is) = v_hub(m2,m1,na,na1,is) + &
                                   CONJG(nsg(m2,m1,na,na1,is) + nsg(m2,m1,na,na1,isop)) * &
                                   Hubbard_J0(nt1)
                                eth = eth + nsg(m2,m1,na,na1,is) *    &
                                   CONJG(nsg(m2,m1,na,na1,is) + nsg(m2,m1,na,na1,isop)) * &
                                   Hubbard_J0(nt1) * 0.5d0
                             ENDDO
                             !
                          ENDDO ! m1
                          !
                     ENDIF
                     !
                  ENDIF
                  !
              ENDIF
              !
           ENDDO ! viz
           !
        ENDDO ! is
        !  
     ENDIF
     !
     ! Hubbard_alpha or Hubbard_alpha_back
     !
     IF ( ldim_u(nt1).GT.0 .AND. &
          (Hubbard_alpha(nt1).NE.0.0d0 .OR. Hubbard_alpha_back(nt1).NE.0.0d0) ) THEN
          !
          na = find_viz(na1,na1)
          !
          DO is = 1, nspin
             !
             DO m1 = 1, 2*Hubbard_l(nt1)+1
                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + Hubbard_alpha(nt1)
                eth = eth + nsg(m1,m1,na,na1,is)* Hubbard_alpha(nt1)
             ENDDO
             !
             IF ( ldim_u(nt1).GT.2*Hubbard_l(nt1)+1 .AND. &
                  Hubbard_alpha_back(nt1).NE.0.0d0 ) THEN
                ! Background states
                DO m1 = 2*Hubbard_l(nt1)+2, ldim_u(nt1)
                   v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + Hubbard_alpha_back(nt1)
                   eth = eth + nsg(m1,m1,na,na1,is)* Hubbard_alpha_back(nt1)
                ENDDO
             ENDIF
             !
          ENDDO
          !
     ENDIF
     !
  ENDDO ! na1
  !
  IF (nspin.EQ.1) eth = eth * 2.d0
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
     !write(stdout,*) "check col UV",  check
     !write(stdout,*) "check_en col UV",  check_en
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_extended
!---------------------------------------------------------------------

SUBROUTINE v_hubbard_extended_nc (nsg, v_hub, eth)
   !-----------------------------------------------------------------------------------
   !
   !! Computes extended Hubbard potential and Hubbard energy.
   !! DFT+U+V: Simplified rotationally-invariant formulation by
   !! V.L. Campo Jr and M. Cococcioni, J. Phys.: Condens. Matter 22, 055602 (2010).
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE ldaU,              ONLY : Hubbard_l, Hubbard_alpha, Hubbard_J0, Hubbard_beta,   &
                                 ldim_u, ldmx_tot, max_num_neighbors, at_sc, neighood, &
                                 Hubbard_V, Hubbard_alpha_back, is_hubbard, is_hubbard_back
   USE lsda_mod,          ONLY : nspin
   USE control_flags,     ONLY : iverbosity, dfpt_hub
   USE io_global,         ONLY : stdout
   USE noncollin_module,  ONLY : npol
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   COMPLEX(DP), INTENT(OUT) :: v_hub(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   REAL(DP),    INTENT(OUT) :: eth
   COMPLEX(DP) :: check, check_en
   ! 
   ! Local variables 
   !
   INTEGER :: is, is1, isop, na, na1, na2, nt, nt1, nt2, m1, m2, viz, equiv_na2
   INTEGER, EXTERNAL :: find_viz
   !
   eth  = 0.d0
   v_hub(:,:,:,:,:) = (0.d0, 0.d0)
   check = (0.0,0.0)
   check_en = (0.0,0.0)
   !
   !write(stdout,*) nsg
   DO na1 = 1, nat
      !
      nt1 = ityp(na1)
      !
      IF ( is_hubbard(nt1) ) THEN
         !
         DO is = 1, nspin
            !
            IF (is == 2) THEN
               is1 = 3
            ELSEIF (is == 3) THEN
               is1 = 2
            ELSE
               is1 = is
            ENDIF
            DO viz = 1, neighood(na1)%num_neigh
               !
               na2 = neighood(na1)%neigh(viz)
               equiv_na2 = at_sc(na2)%at
               nt2 = ityp(equiv_na2)
               !
               IF (is_hubbard(nt2) .AND. &
                   (Hubbard_V(na1,na2,1).NE.0.d0) ) THEN
                   !
                   ! Here no need to use is1: complex conjugation is enough
                   ! For both standard and background states of a center atom
                   DO m1 = 1, ldim_u(nt1)
                      ! For both standard and background states of the neighbor atom
                      DO m2 = 1, ldim_u(nt2)
                         !
                         v_hub(m2,m1,viz,na1,is) = - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,1)
                        check = check - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,1)
                         !
                         eth = eth - nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                     * Hubbard_V(na1,na2,1) * 0.5d0
                        check_en = check_en - nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                     * Hubbard_V(na1,na2,1) * 0.5d0
                         !
                      ENDDO
                   ENDDO
                   !
                   IF ( na1.EQ.na2 .AND. is1.EQ.is) THEN
                      !
                      na = find_viz(na1,na1)
                      !
                      ! This is the diagonal term (like in the DFT+U only case)
                      ! 
                      DO m1 = 1, ldim_u(nt1)
                         !
                         v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) &
                                        + Hubbard_V(na1,na1,1) * 0.5d0
                        check = check + Hubbard_V(na1,na1,1) * 0.5d0
                         ! 
                         eth = eth + nsg(m1,m1,na,na1,is) &
                                        * Hubbard_V(na1,na1,1) * 0.5d0
                        check_en = check_en + nsg(m1,m1,na,na1,is) &
                                        * Hubbard_V(na1,na1,1) * 0.5d0
                         !
                      ENDDO
                      !
                   ENDIF
                   !
               ENDIF
               !
            ENDDO ! viz
            !
         ENDDO ! is
         !  
      ENDIF
      !
      ! Hubbard_alpha
      !!
      IF ( ldim_u(nt1).GT.0 .AND. (Hubbard_alpha(nt1).NE.0.0d0 ) ) THEN
         !
         na = find_viz(na1,na1)
         !
         DO is = 1, npol
            !
            DO m1 = 1, 2*Hubbard_l(nt1)+1
               v_hub(m1,m1,na,na1,is**2) = v_hub(m1,m1,na,na1,is**2) + Hubbard_alpha(nt1)
               eth = eth + nsg(m1,m1,na,na1,is**2)* Hubbard_alpha(nt1)
            ENDDO
            !
         ENDDO
         !
      ENDIF
      !
   ENDDO ! na1
   !
   !
   IF (nspin.EQ.1) eth = eth * 2.d0
   !
   ! Hubbard energy
   !
   IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
      WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
   ENDIF
   !
   RETURN
   !
 END SUBROUTINE v_hubbard_extended_nc
!
!-----------------------------------------------------------------------
SUBROUTINE v_hubbard_resolved( ns, v_hub, eth )
!---------------------------------------------------------------------
!
!! Computes Hubbard potential and Hubbard energy
!! for a manifold of selected spin- or magnetic quantum orbitals.
!! The Hubbard potential is first calculated in the diagonal representation
!! based on the eigenvalues of the occupation matrix and then 
!! re-rotated using the eigenvectors in order to retain compatiblity with
!! other parts of the code.
!! See Macke et al., arXiv:2312.13580 (2023).
!! Uses the simplified rotationally-invariant formulation by
!! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_Um, &
                                 Hubbard_alpha_m, lambda_ns, &
                                 eigenvecs_ref, order_um, apply_U, hub_pot_fix
USE lsda_mod,             ONLY : nspin
USE constants,            ONLY : eps16, RYTOEV
USE control_flags,        ONLY : iverbosity, dfpt_hub
USE io_global,            ONLY : stdout
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! occupation matrix
REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! Hubbard potential
REAL(DP), INTENT(OUT) :: eth
!! Hubbard energy
!
!  ... local variables
!
! Hubbard potential in the diagonal representation
REAL(DP)                 :: v_hub_diag(2*Hubbard_lmax+1,nspin,nat), temp
REAL(DP)                 :: effU, effalpha
! eigenvectors of the ns occupation matrix in the current iteration
COMPLEX(DP)              :: eigenvecs_current(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin)
!
INTEGER                  :: is, na, nt, m1, m2, m3, ldim, m_order
! the ordering vector for the orbital-tracking routine
INTEGER                  :: order(2*Hubbard_lmax+1)
LOGICAL                  :: is_first
IF (.NOT. ALLOCATED(order_um)) THEN 
  IF (nspin == 2 ) THEN 
     ALLOCATE(order_um(2*Hubbard_lmax+1,nspin, nat))  
  ELSE 
     ALLOCATE(order_um(2*Hubbard_lmax+1,1,nat)) 
  END IF 
  order_um = 0
END IF 

!
!
eth    = 0.d0
lambda_ns(:,:,:) = 0.d0
! orbital occupations (=eigenvalues of rho%ns)
!
v_hub(:,:,:,:) = 0.d0
v_hub_diag(:,:,:) = 0.d0 
!
IF ( .NOT. apply_U ) RETURN
! Do not apply corrections before the eigenstates have stabilized
! The eigenstates are considered stable a) when starting from a
! converged charge density or b) when the treshold of iterative 
! diagonalization gets small enough (typically within 3-6 iterations).
! This is controlled by electrons.f90.
! 
DO na = 1, nat
   !
   nt = ityp (na)
   !
   IF ( ANY(ABS(Hubbard_Um(:,:,nt)) .GT. eps16) .OR. &
        ANY(ABS(Hubbard_alpha_m(:,:,nt)) .GT. eps16) ) THEN       
      !
      ldim = 2 * Hubbard_l(nt) + 1
      eigenvecs_current(:,:,:) = CMPLX(0.d0,0.d0, kind=dp)
      is_first = ALL(eigenvecs_ref(:,:,:,na) == eigenvecs_current)  

      !
      effU = 0.0
      effalpha = 0.0
      !
      CALL diag_ns( ldim, ns(1:ldim,1:ldim,:,na), lambda_ns(1:ldim,:,na), &
                    eigenvecs_current(1:ldim,1:ldim,:) )
      !
      DO is = 1, nspin
         !
         ! sort eigenvectors with respect to the (reference) order established in eigvecs_first
         IF (is_first)  THEN  
             order(1:ldim) = order_um(1:ldim,is,na)
             IF (ALL(order(1:ldim)==0)) order(1:ldim) = [(m1,m1=1,ldim)] 
             DO m1 =1, ldim 
               eigenvecs_ref(1:ldim, order(m1), is, na) = eigenvecs_current(1:ldim, m1, is) 
             END DO 
         END IF 
         order(:) = 0
         CALL order_eigenvecs( order(1:ldim), eigenvecs_current(1:ldim,1:ldim,is), &
                                 eigenvecs_ref(1:ldim,1:ldim,is,na), ldim )
         !
         
         order_um(1:ldim,is,na) = order(1:ldim) 
         DO m1 = 1, ldim
            !
            ! calculate Hubbard potential and -energy
            ! using the eigenvalues and their ordering
            m_order = order(m1)
            effU = Hubbard_Um(m_order,is,nt)
            effalpha = Hubbard_alpha_m(m_order,is,nt)
            !
            ! linear Hubbard U terms:
            eth = eth + ( effalpha + 0.5D0*effU ) * lambda_ns(m1,is,na)
            v_hub_diag(m1,is,na) = v_hub_diag(m1,is,na) + &
                                       effalpha + 0.5D0*effU
            ! quadratic Hubbard U terms:
            eth = eth - 0.5D0 * effU * lambda_ns(m1,is,na) * lambda_ns(m1,is,na)
            v_hub_diag(m1,is,na) = v_hub_diag(m1,is,na) - &
                                       effU * lambda_ns(m1,is,na)
            !
            ! The following can be eliminated once code is approved
#if defined(__DEBUG)
            WRITE( stdout,'(/5x,"v_of_rho:")')
            WRITE( stdout,'(/5x,"AT:",i2," ,m_order:",i2,", is: ",i2,", LAMBDA:",f7.5," U:",f7.5)') &
                  na,m_order,is,lambda_ns(m1,is,na),effU*RYTOEV
            IF ( effalpha /= 0.d0) &
               WRITE( stdout,'(/5x,"m_order:",i2,", is: ",i2,", LAMBDA:",f8.5," alpha:",f8.5)') &
                  m_order,is,lambda_ns(m1,is,na),effalpha*RYTOEV
            WRITE( stdout,'(/5x,"Hubbard potential of this state (v_hub_diag) = ",f8.5)') &
                  v_hub_diag(m1,is,na)
            WRITE( stdout,'(/5x,"Cumulative running total of Hubbard energies (eth) = ",f8.5)') eth
#endif
            !
         ENDDO
         !
         ! backrotation of v_hub_diag to the non-diagonal v_hub
         DO m1 = 1, ldim
            DO m2 = 1, ldim
               temp = CMPLX(0.d0,0.d0, kind=dp)
               DO m3 = 1, ldim
                  temp = temp + CONJG(eigenvecs_current(m1,m3,is))* &
                     v_hub_diag(m3,is,na)*eigenvecs_current(m2,m3,is)
               ENDDO
               v_hub(m1,m2,is,na) = DBLE(temp)
            ENDDO
         ENDDO
      ENDDO ! is
      !
      IF ( ALL(eigenvecs_ref(:,:,:,na) .EQ. 0.d0) ) THEN
         !
         ! if this routine is executed for the first time,
         ! save the current eigenvecs as reference eigenvecs
         eigenvecs_ref(1:ldim,1:ldim,:,na) = eigenvecs_current(1:ldim,1:ldim,:)
      ENDIF
      !
   ENDIF
   !
ENDDO ! nt
!
IF (nspin==1) eth = 2.d0 * eth
!
! print Hubbard energy
!
IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
   WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.5,1x," (Ry)")') eth
ENDIF
!
RETURN
!
END SUBROUTINE v_hubbard_resolved
!----------------------------------------------------------------------------
SUBROUTINE v_hubbard_resolved_nc( ns, v_hub, eth )
!----------------------------------------------------------------------------
!
!! Computes Hubbard potential and Hubbard energy
!! for a manifold of selected orbitals (noncollinear formulation).
!! The Hubbard potential is first calculated in the diagonal representation
!! based on the eigenvalues of the occupation matrix and then 
!! re-rotated using the eigenvectors in order to retain compatiblity with
!! other parts of the code.
!! Uses the simplified rotationally-invariant formulation by
!! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_Um_nc, &
                                 Hubbard_alpha_m_nc, lambda_ns, order_um,&
                                 eigenvecs_ref, apply_U, hub_pot_fix
USE lsda_mod,             ONLY : nspin
USE constants,            ONLY : eps16, RYTOEV
USE control_flags,        ONLY : iverbosity, dfpt_hub
USE io_global,            ONLY : stdout
!
IMPLICIT NONE
!
COMPLEX(DP), INTENT(IN) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! occupation matrix
COMPLEX(DP), INTENT(OUT):: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! Hubbard potential
REAL(DP), INTENT(OUT) :: eth
!! Hubbard energy
!
!  ... local variables
!
REAL(DP)                 :: effU, effalpha
COMPLEX(DP)              :: v_hub_diag(4*Hubbard_lmax+2,nat), temp
COMPLEX(DP)              :: v_hub_temp(4*Hubbard_lmax+2,4*Hubbard_lmax+2)
COMPLEX(DP)              :: eigenvecs_current(4*Hubbard_lmax+2,4*Hubbard_lmax+2)
!
INTEGER                  :: is, na, nt, m1, m2, m3, ldim, m_order
!
INTEGER                  :: order(4*Hubbard_lmax+2) 
LOGICAL                  :: is_first
!! ordering vector
!
IF (.NOT. ALLOCATED(order_um)) THEN 
  ALLOCATE(order_um(4*Hubbard_lmax+2,1,nat)) 
  order_um = 0
END IF 

eth    = 0.d0
lambda_ns(:,:,:) = 0.d0
! orbital occupations (=eigenvalues of rho%ns)
!
v_hub(:,:,:,:) = 0.d0
v_hub_diag(:,:) = 0.d0 
! Hubbard potential in the diagonal representation
!
IF ( .NOT. apply_U ) RETURN
! Do not apply corrections before the eigenstates have stabilized
! The eigenstates are considered stable a) when starting from a
! converged charge density or b) when the treshold of iterative 
! diagonalization gets small enough (typically within 3-6 iterations)
! this is checked in electrons.f90
! 
DO na = 1, nat
   !
   nt = ityp (na)
   !
   IF ( ANY(ABS(Hubbard_Um_nc(:,nt)) .GT. eps16) .OR. &
        ANY(ABS(Hubbard_alpha_m_nc(:,nt)) .GT. eps16) ) THEN       
      !
      ldim = 2 * Hubbard_l(nt) + 1
      eigenvecs_current(:,:) = CMPLX(0.d0,0.d0, kind=dp)
      is_first = ALL(eigenvecs_ref(:,1,:,na) == eigenvecs_current)  
      !
      effU = 0.0
      effalpha = 0.0
      !
      CALL diag_ns_nc( ldim, ns(1:ldim,1:ldim,:,na), lambda_ns(1:2*ldim,1,na), &
                     eigenvecs_current(1:2*ldim,1:2*ldim) )
                              !
      ! sort eigenvectors with respect to the (reference) order established in eigvecs_first
      IF (is_first)  THEN  
        order(1:2*ldim) = order_um(1:2*ldim,1,na)
        IF (ALL(order(1:2*ldim)==0)) order(1:2*ldim) = [(m1,m1=1,2*ldim)] 
        DO m1 =1, 2*ldim 
          eigenvecs_ref(1:2*ldim, order(m1), 1, na) = eigenvecs_current(1:2*ldim, m1) 
        END DO 
      END IF 
      order(:) = 0
      !
      CALL order_eigenvecs( order(1:2*ldim), eigenvecs_current(1:2*ldim,1:2*ldim), &
                                 eigenvecs_ref(1:2*ldim,1:2*ldim,1,na), 2*ldim )
      !
      ! No need to iterate over is (all done in diag_ns_nc)
      order_um(1:2*ldim,1,na) = order(1:2*ldim) 
      DO m1 = 1, 2*ldim
         !
         ! calculate Hubbard potential and -energy
         ! using the eigenvalues and their ordering
         m_order = order(m1)
         effU = Hubbard_Um_nc(m_order,nt)
         effalpha = Hubbard_alpha_m_nc(m_order,nt)
         !
         ! linear Hubbard U terms:
         eth = eth + ( effalpha + 0.5D0*effU ) * lambda_ns(m1,1,na)
         v_hub_diag(m1,na) = v_hub_diag(m1,na) + &
                                    effalpha + 0.5D0*effU
         ! quadratic Hubbard U terms:
         eth = eth - 0.5D0 * effU * lambda_ns(m1,1,na) * lambda_ns(m1,1,na)
         v_hub_diag(m1,na) = v_hub_diag(m1,na) - &
                                    effU * lambda_ns(m1,1,na)
         !
         ! The following block can be eliminated once code is approved
#if defined(__DEBUG)
         WRITE( stdout,'(/5x,"v_of_rho (NC):")')
         WRITE( stdout,'(/5x,"AT:",i2," ,m_order:",i2,", LAMBDA:",f7.5," U:",f7.5)') &
               na,m_order,lambda_ns(m1,1,na),effU*RYTOEV
         IF ( effalpha /= 0.d0) &
            WRITE( stdout,'(/5x,"m_order:",i2,", LAMBDA:",f8.5," alpha:",f8.5)') &
               m_order,lambda_ns(m1,1,na),effalpha*RYTOEV
         WRITE( stdout,'(/5x,"Hubbard potential of this state (v_hub_diag) = ",f8.5)') &
               v_hub_diag(m1,na)
         WRITE( stdout,'(/5x,"Cumulative running total of Hubbard energies (eth) = ",f8.5)') eth
#endif
         !
      ENDDO
         !
      ! backrotation of v_hub_diag to the non-diagonal v_hub
      ! first, go back to the 2ldim*2ldim matrix
      v_hub_temp(:,:) = 0.d0
      !
      DO m1 = 1, 2*ldim
         DO m2 = 1, 2*ldim
            temp = CMPLX(0.d0,0.d0, kind=dp)
            DO m3 = 1, 2*ldim
               temp = temp + CONJG(eigenvecs_current(m1,m3))* &
                  v_hub_diag(m3,na)*eigenvecs_current(m2,m3)
            ENDDO
            v_hub_temp(m1,m2) = temp
         ENDDO
      ENDDO
      ! now, sort the different quadrants of v_hub_temp into the actual
      ! v_hub. upper left quadrant is spin 1, upper right spin 2 etc.
      !
      DO m1 = 1, ldim
         DO m2 = 1,ldim
            v_hub(m1,m2,1,na) = v_hub_temp(m1,m2)
            v_hub(m1,m2,2,na) = v_hub_temp(m1,ldim+m2)
            v_hub(m1,m2,3,na) = v_hub_temp(ldim+m1,m2)
            v_hub(m1,m2,4,na) = v_hub_temp(ldim+m1,ldim+m2)
         ENDDO
      ENDDO
      !
      IF ( ALL(eigenvecs_ref(:,:,1,na) .EQ. 0.d0) ) THEN
         ! if this routine is executed for the first time,
         ! save the current eigenvecs as reference eigenvecs
         eigenvecs_ref(1:2*ldim,1:2*ldim,1,na) = &
                     & eigenvecs_current(1:2*ldim,1:2*ldim)
      ENDIF
      !
   ENDIF
   !
ENDDO ! nt
!
! print Hubbard energy
!
IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
   WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.5,1x," (Ry)")') eth
ENDIF
!
RETURN
!
END SUBROUTINE v_hubbard_resolved_nc
!----------------------------------------------------------------------------
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
