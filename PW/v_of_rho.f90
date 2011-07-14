!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE v_of_rho( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the Hartree and Exchange and Correlation
  ! ... potential and energies which corresponds to a given charge density
  ! ... The XC potential is computed in real space, while the
  ! ... Hartree potential is computed in reciprocal space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : noncolin
  USE ions_base,        ONLY : nat
  USE ldaU,             ONLY : lda_plus_U, Hubbard_lmax, Hubbard_l, &
                               Hubbard_U, Hubbard_alpha
  USE funct,            ONLY : dft_is_meta
  USE scf,              ONLY : scf_type
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(IN) :: rho  ! the valence charge
  TYPE(scf_type), INTENT(INOUT) :: v ! the scf (Hxc) potential 
  !!!!!!!!!!!!!!!!! NB: NOTE that in F90 derived data type must be INOUT and 
  !!!!!!!!!!!!!!!!! not just OUT because otherwise their allocatable or pointer
  !!!!!!!!!!!!!!!!! components are NOT defined !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: vtxc, etxc, ehart, eth, charge
    ! the integral V_xc * rho
    ! the E_xc energy
    ! the hartree energy
    ! the hubbard energy
    ! the integral of the charge
  REAL(DP), INTENT(INOUT) :: etotefield
    ! electric field energy - inout due to the screwed logic of add_efield
  ! ! 
  INTEGER :: is
  !
  CALL start_clock( 'v_of_rho' )
  !
  ! ... calculate exchange-correlation potential
  !
  if (dft_is_meta()) then
     call v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
  else
     CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  endif
  !
  ! ... add a magnetic field  (if any)
  !
  CALL add_bfield( v%of_r, rho%of_r )
  !
  ! ... calculate hartree potential
  !
  CALL v_h( rho%of_g, ehart, charge, v%of_r )
  !
  if (lda_plus_u) call v_Hubbard(rho%ns,v%ns,eth)
  !
  ! ... add an electric field
  ! 
  DO is = 1, nspin
     CALL add_efield(v%of_r(1,is), etotefield, rho%of_r, .false. )
  END DO
  !
  CALL stop_clock( 'v_of_rho' )
  !
  RETURN
  !
END SUBROUTINE v_of_rho
!
SUBROUTINE v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
  !----------------------------------------------------------------------------
  !
  ! ... Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : xc, xc_spin, get_igcx, get_igcc
  USE scf,              ONLY : scf_type
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin), kedtaur(dfftp%nnr,nspin), vtxc, etxc
    ! V_xc potential
    ! local K energy density 
    ! integral V_xc * rho
    ! E_xc energy
  !
  ! ... local variables
  !
  CALL start_clock( 'v_xc_meta' )
  !
  etxc   = 0.D0
  vtxc   = 0.D0
  v(:,:) = 0.D0
  !
!  IF (get_igcx()==7.AND.get_igcc()==6) THEN
     call v_xc_tpss( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
!  ELSE
!     CALL errore('v_xc_meta','wrong igcx and/or igcc',1)
!  ENDIF
  CALL stop_clock( 'v_xc_meta' )
  RETURN
END SUBROUTINE v_xc_meta
!
SUBROUTINE v_xc_tpss( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
  !     ===================
  !--------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE gvect,            ONLY : g,nl,ngm
  USE fft_base,         ONLY : dfftp
  USE scf,              ONLY : scf_type
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega, alat
  USE constants,        ONLY : e2
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum

  IMPLICIT NONE
  !
  ! input
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP),INTENT(IN) :: rho_core(dfftp%nnr)
  COMPLEX(DP),INTENT(IN) :: rhog_core(ngm)
  REAL(DP),INTENT(OUT) :: etxc, vtxc, v(dfftp%nnr,nspin), kedtaur(dfftp%nnr,nspin)
!  integer nspin , nnr
!  real(8)  grho(nnr,3,nspin), rho(nnr,nspin),kedtau(nnr,nspin)
  ! output: excrho: exc * rho ;  E_xc = \int excrho(r) d_r
  ! output: rhor:   contains the exchange-correlation potential
  REAL(DP) :: zeta, rh
  INTEGER :: k, ipol, is
  REAL(DP) :: grho2 (2), sx, sc, v1x, v2x, v3x,v1c, v2c, v3c, &
       v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw ,v2cup(3),v2cdw(3), &
       v3xup, v3xdw,grhoup(3),grhodw(3),&
       segno, arho, atau, fac
  !
  REAL(DP),    ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP),    ALLOCATABLE :: rhoout(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:,:)
  REAL(DP), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10
  !
  ALLOCATE (grho(3,dfftp%nnr,nspin))
  ALLOCATE (h(3,dfftp%nnr,nspin))
  ALLOCATE (rhoout(dfftp%nnr,nspin))
  ALLOCATE (rhogsum(ngm,nspin))
  !
  vtxc = 0.d0
  etxc = 0.d0
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  rhoout(:,1:nspin)=rho%of_r(:,1:nspin)
  rhogsum(:,1:nspin)=rho%of_g(:,1:nspin)
  fac = 1.D0 / DBLE( nspin )
  !
  DO is = 1, nspin
     !
     rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)
     rhogsum(:,is) = fac * rhog_core(:) + rhogsum(:,is)
     !
     CALL gradrho( dfftp%nnr, rhogsum(1,is), ngm, g, nl, grho(1,1,is) )
     !
  END DO
  !
  DO k = 1, dfftp%nnr
     DO is = 1, nspin
        grho2 (is) = grho(1,k, is)**2 + grho(2,k,is)**2 + grho(3,k, is)**2
     ENDDO
     IF (nspin == 1) THEN
        !
        !    This is the spin-unpolarised case
        !
        arho = ABS (rho%of_r (k, 1) )
        segno = SIGN (1.d0, rho%of_r (k, 1) )
        atau = rho%kin_r(k,1) / e2  ! kinetic energy density in Hartree
        IF (arho.GT.epsr.AND.grho2 (1) .GT.epsg.AND.ABS(atau).GT.epsr) THEN
           CALL tpsscxc (arho, grho2(1),atau,sx, sc, v1x, v2x,v3x,v1c, v2c,v3c)
!           if (mod(k,100).eq.0) then
!             write(6,*) 'PON k=',k
!             write(6,*) ' arho,atau=',arho,atau
!             write(6,*) ' sx,sc=',sx,sc
!             write(6,*) ' v1x,v2x,v3c=',v1x,v2x,v3x
!             write(6,*) ' v1c,v2c,v3c=',v1c,v2c,v3c
!           endif
           v(k, 1) =  (v1x + v1c )*e2
           kedtaur(k,1)=  (v3x + v3c) * 0.5d0 * e2 
           ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
           h(:,k,1) =  (v2x + v2c)*grho (:,k,1) *e2
           etxc = etxc +  (sx + sc) * segno *e2
           vtxc = vtxc + (v1x+v1c)*e2*arho + 1.5d0*kedtaur(k,1)*rho%kin_r(k,1)
        ELSE  
           h (:, k, 1) = 0.d0  
           kedtaur(k,1)=0.d0
        ENDIF
     ELSE
        !
        !    spin-polarised case
        !
        CALL tpsscx_spin(rho%of_r(k, 1), rho%of_r(k, 2), grho2 (1), grho2 (2), &
             rho%kin_r(k,1)/e2,rho%kin_r(k,2)/e2,sx, &
             v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
        rh = rho%of_r(k, 1) + rho%of_r(k, 2)
        IF (rh.GT.epsr) THEN
           zeta = (rho%of_r(k, 1) - rho%of_r(k, 2) ) / rh
           DO ipol=1,3
              grhoup(ipol)=grho(ipol,k,1)
              grhodw(ipol)=grho(ipol,k,2)
           END DO
           atau = (rho%kin_r(k,1)+rho%kin_r(k,2)) / e2 ! KE-density in Hartree
           CALL tpsscc_spin(rh,zeta,grhoup,grhodw, &
                atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3c)
        ELSE
           sc = 0.d0  
           v1cup = 0.d0  
           v1cdw = 0.d0  
           v2cup=0.d0
           v2cdw=0.d0
           v3c=0.d0
           !
        ENDIF
        !
        ! first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        v(k, 1) =  (v1xup + v1cup)*e2
        v(k, 2) =  (v1xdw + v1cdw)*e2
        !
        ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        h(:,k,1) = (v2xup*grho(:,k,1) + v2cup(:))*e2
        h(:,k,2) = (v2xdw*grho(:,k,2) + v2cdw(:)) *e2
        kedtaur(k,1)=  (v3xup + v3c) * 0.5d0 * e2
        kedtaur(k,2)=  (v3xdw + v3c) * 0.5d0 * e2
        etxc = etxc +  (sx + sc)*e2
        vtxc = vtxc + (v1xup+v1cup+v1xdw+v1cdw)*e2*rh
     ENDIF
  ENDDO
  !
  ALLOCATE( dh( dfftp%nnr ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin
     !
     CALL grad_dot( dfftp%nnr, h(1,1,is), ngm, g, nl, alat, dh )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     rhoout(:,is)=rhoout(:,is)-fac*rho_core(:)
     vtxc = vtxc - SUM( dh(:) * rhoout(:,is) )
     !
  END DO
  DEALLOCATE(dh)
  !
  vtxc = omega * (vtxc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 ))
  etxc = omega * etxc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  !
  CALL mp_sum(  vtxc , intra_pool_comm )
  CALL mp_sum(  etxc , intra_pool_comm )
  DEALLOCATE(grho)
  DEALLOCATE(h)
  DEALLOCATE(rhoout)
  DEALLOCATE(rhogsum)
  !
  RETURN
  !
END SUBROUTINE v_xc_tpss
!----------------------------------------------------------------------------
SUBROUTINE v_xc( rho, rho_core, rhog_core, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !
  ! ... Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : xc, xc_spin
  USE scf,              ONLY : scf_type
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum

  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin), vtxc, etxc
    ! V_xc potential
    ! integral V_xc * rho
    ! E_xc energy
  !
  ! ... local variables
  !
  REAL(DP) :: rhox, arhox, zeta, amag, vs, ex, ec, vx(2), vc(2), rhoneg(2)
    ! the total charge in each point
    ! the absolute value of the charge
    ! the absolute value of the charge
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir, ipol
    ! counter on mesh points
    ! counter on nspin
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag    = 1.D-20
  !
  !
  CALL start_clock( 'v_xc' )
  !
  etxc   = 0.D0
  vtxc   = 0.D0
  v(:,:) = 0.D0
  rhoneg = 0.D0
  !
  IF ( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
     !
     ! ... spin-unpolarized case
     !
!$omp parallel do private( rhox, arhox, ex, ec, vx, vc ), &
!$omp             reduction(+:etxc,vtxc), reduction(-:rhoneg)
     DO ir = 1, dfftp%nnr
        !
        rhox = rho%of_r(ir,1) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           CALL xc( arhox, ex, ec, vx(1), vc(1) )
           !
           v(ir,1) = e2*( vx(1) + vc(1) )
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           !
           vtxc = vtxc + v(ir,1) * rho%of_r(ir,1)
           !
        ENDIF
        !
        IF ( rho%of_r(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        !
     END DO
!$omp end parallel do
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ! ... spin-polarized case
     !
!$omp parallel do private( rhox, arhox, zeta, ex, ec, vx, vc ), &
!$omp             reduction(+:etxc,vtxc), reduction(-:rhoneg)
     DO ir = 1, dfftp%nnr
        !
        rhox = rho%of_r(ir,1) + rho%of_r(ir,2) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = ( rho%of_r(ir,1) - rho%of_r(ir,2) ) / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) zeta = SIGN( 1.D0, zeta )
           !
           IF ( rho%of_r(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
           IF ( rho%of_r(ir,2) < 0.D0 ) rhoneg(2) = rhoneg(2) - rho%of_r(ir,2)
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           v(ir,:) = e2*( vx(:) + vc(:) )
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           !
           vtxc = vtxc + v(ir,1) * rho%of_r(ir,1) + v(ir,2) * rho%of_r(ir,2)
           !
        END IF
        !
     END DO
!$omp end parallel do
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     ! ... noncolinear case
     !
     DO ir = 1,dfftp%nnr
        !
        amag = SQRT( rho%of_r(ir,2)**2 + rho%of_r(ir,3)**2 + rho%of_r(ir,4)**2 )
        !
        rhox = rho%of_r(ir,1) + rho_core(ir)
        !
        IF ( rho%of_r(ir,1) < 0.D0 )  rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = amag / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) THEN
              !
              rhoneg(2) = rhoneg(2) + 1.D0 / omega
              !
              zeta = SIGN( 1.D0, zeta )
              !
           END IF
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           vs = 0.5D0*( vx(1) + vc(1) - vx(2) - vc(2) )
           !
           v(ir,1) = e2*( 0.5D0*( vx(1) + vc(1) + vx(2) + vc(2 ) ) )
           !
           IF ( amag > vanishing_mag ) THEN
              !
              DO ipol = 2, 4
                 !
                 v(ir,ipol) = e2 * vs * rho%of_r(ir,ipol) / amag
                 !
                 vtxc = vtxc + v(ir,ipol) * rho%of_r(ir,ipol)
                 !
              END DO
              !
           END IF
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           vtxc = vtxc + v(ir,1) * rho%of_r(ir,1)
           !
        END IF
        !
     END DO
     !
  END IF
  !
  CALL mp_sum(  rhoneg , intra_pool_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2E10.3)') rhoneg
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
  ! ... add non local corrections (if any)
  !
  CALL nonloccorr(rho%of_r, rho_core, etxc, vtxc, v)
  !
  CALL mp_sum(  vtxc , intra_pool_comm )
  CALL mp_sum(  etxc , intra_pool_comm )
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
  !
  ! ... Hartree potential VH(r) from n(G)
  !
  USE constants, ONLY : fpi, e2
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : invfft
  USE gvect,     ONLY : nl, nlm, ngm, gg, gstart
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : omega, tpiba2
  USE control_flags, ONLY : gamma_only
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  USE esm,       ONLY: do_comp_esm, esm_hartree, esm_bc
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: rhog(ngm,nspin)
  REAL(DP),  INTENT(INOUT) :: v(dfftp%nnr,nspin)
  REAL(DP),    INTENT(OUT) :: ehart, charge
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:)
  REAL(DP)              :: rgtot_re, rgtot_im, eh_corr
  INTEGER               :: is, ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), rgtot(:), vaux(:)
  INTEGER               :: nt
  !
  CALL start_clock( 'v_h' )
  !
  ALLOCATE( aux( dfftp%nnr ), aux1( 2, ngm ) )
  charge = 0.D0
  !
  IF ( gstart == 2 ) THEN
     !
     charge = omega*REAL( rhog(1,1) )
     !
     IF ( nspin == 2 ) charge = charge + omega*REAL( rhog(1,2) )
     !
  END IF
  !
  CALL mp_sum(  charge , intra_pool_comm )
  !
  ! ... calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
     !
     ! ... calculate modified Hartree potential for ESM
     !
     CALL esm_hartree (rhog, ehart, aux)
     !
  ELSE
     !
     ehart     = 0.D0
     aux1(:,:) = 0.D0
     !
!$omp parallel do private( fac, rgtot_re, rgtot_im ), reduction(+:ehart)
     DO ig = gstart, ngm
        !
        fac = 1.D0 / gg(ig)
        !
        rgtot_re = REAL(  rhog(ig,1) )
        rgtot_im = AIMAG( rhog(ig,1) )
        !
        IF ( nspin == 2 ) THEN
           !
           rgtot_re = rgtot_re + REAL(  rhog(ig,2) )
           rgtot_im = rgtot_im + AIMAG( rhog(ig,2) )
           !
        END IF
        !
        ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
        !
        aux1(1,ig) = rgtot_re * fac
        aux1(2,ig) = rgtot_im * fac
        !
     ENDDO
!$omp end parallel do
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
     END IF
     ! 
     if (do_comp_mt) then
        ALLOCATE( vaux( ngm ), rgtot(ngm) )
        rgtot(:) = rhog(:,1)
        if (nspin==2) rgtot(:) = rgtot(:) + rhog(:,2)
        CALL wg_corr_h (omega, ngm, rgtot, vaux, eh_corr)
        aux1(1,1:ngm) = aux1(1,1:ngm) + REAL( vaux(1:ngm))
        aux1(2,1:ngm) = aux1(2,1:ngm) + AIMAG(vaux(1:ngm))
        ehart = ehart + eh_corr
        DEALLOCATE( rgtot, vaux )
     end if
     !
     CALL mp_sum(  ehart , intra_pool_comm )
     ! 
     aux(:) = 0.D0
     !
     aux(nl(1:ngm)) = CMPLX ( aux1(1,1:ngm), aux1(2,1:ngm), KIND=dp )
     !
     IF ( gamma_only ) THEN
        !
        aux(nlm(1:ngm)) = CMPLX ( aux1(1,1:ngm), -aux1(2,1:ngm), KIND=dp )
        !
     END IF
  END IF
  !
  ! ... transform hartree potential to real space
  !
  CALL invfft ('Dense', aux, dfftp)
  !
  ! ... add hartree potential to the xc potential
  !
  IF ( nspin == 4 ) THEN
     !
     v(:,1) = v(:,1) + DBLE (aux(:))
     !
  ELSE
     !
     DO is = 1, nspin
        !
        v(:,is) = v(:,is) + DBLE (aux(:))
        !
     END DO
     !
  END IF
  !
  DEALLOCATE( aux, aux1 )
  !
  CALL stop_clock( 'v_h' )
  !
  RETURN
  !
END SUBROUTINE v_h
!
!-----------------------------------------------------------------------
SUBROUTINE v_hubbard(ns, v_hub, eth)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin

  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) 
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) 
  REAL(DP), INTENT(OUT) :: eth
  INTEGER :: is, na, nt, m1, m2
  !
  ! Now the contribution to the total energy is computed. 
  !
  eth = 0.d0  
  v_hub(:,:,:,:) = 0.d0
  DO na = 1, nat  
     nt = ityp (na)  
     IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN  
        DO is = 1, nspin  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1  
              eth = eth + ( Hubbard_alpha(nt) + 0.5D0 * Hubbard_U(nt) ) * &
                            ns(m1,m1,is,na) 
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                          ( Hubbard_alpha(nt) + 0.5D0 * Hubbard_U(nt) ) 
              DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                 eth = eth - 0.5D0 * Hubbard_U(nt) * &
                                     ns(m2,m1,is,na)* ns(m1,m2,is,na) 
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                      Hubbard_U(nt) * ns(m2,m1,is,na) 
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  IF (nspin.EQ.1) eth = 2.d0 * eth

  RETURN

END SUBROUTINE v_hubbard
#ifdef __SOLVENT
  !----------------------------------------------------------------------------
  SUBROUTINE v_h_of_rho_r( rho, v )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential in R space from a total (spinless) density 
  !     in R space. 
  !     THIS SUBROUTINE IS ONLY USED BY THE SOLVENT MODULE
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE control_flags,   ONLY : gamma_only
  USE constants,       ONLY : fpi, e2
  USE cell_base,       ONLY : tpiba2
  USE gvect,           ONLY : nl, ngm, nlm, gg, gstart
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rho( dfftp%nnr )
  REAL( DP ), INTENT(OUT)    :: v( dfftp%nnr )
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhoaux( : )
  COMPLEX( DP ), ALLOCATABLE :: vaux( : )
  REAL( DP )                 :: fac
  INTEGER                    :: ig
  !
  ! ... Bring rho to G space
  !
  ALLOCATE( rhoaux( dfftp%nnr ) )
  rhoaux( : ) = CMPLX(rho( : ),0.D0,kind=dp) 
  !
  CALL fwfft('Dense', rhoaux, dfftp)
  !
  ! ... Compute total potential in G space
  !
  ALLOCATE( vaux( dfftp%nnr ) )
  vaux( : ) = CMPLX(0.D0,0.D0,kind=dp)
  !
  DO ig = gstart, ngm
    !
    fac = 1.D0 / gg(ig)
    vaux(nl(ig)) = CMPLX(REAL(rhoaux(nl(ig))),AIMAG(rhoaux(nl(ig))),kind=dp) * fac
    !
  ENDDO
  !
  DEALLOCATE(rhoaux)
  !
  IF ( gamma_only ) THEN
     !
     vaux(nlm(1:ngm)) = CMPLX(REAL( vaux(nl(:))),-AIMAG(vaux(nl(:))),kind=DP)
     !
  END IF
  !
  fac = e2 * fpi / tpiba2
  vaux = vaux * fac 
  !
  ! ... Bring V to R space
  !
  CALL invfft('Dense', vaux, dfftp)
  !
  v(:) = REAL(vaux(:))
  !
  DEALLOCATE( vaux )
  !
  RETURN
  !
  END SUBROUTINE v_h_of_rho_r
#endif
#ifdef __SOLVENT
  !----------------------------------------------------------------------------
  SUBROUTINE gradv_h_of_rho_r( rho, gradv )
  !----------------------------------------------------------------------------
  !
  ! ... Gradient of Hartree potential in R space from a total 
  !     (spinless) density in R space
  !     THIS SUBROUTINE IS ONLY USED BY THE SOLVENT MODULE
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE constants,       ONLY : fpi, e2
  USE control_flags,   ONLY : gamma_only
  USE cell_base,       ONLY : tpiba
  USE gvect,           ONLY : nl, ngm, nlm, gg, gstart, g
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rho( dfftp%nnr )
  REAL( DP ), INTENT(OUT)    :: gradv( 3, dfftp%nnr )
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhoaux( : )
  COMPLEX( DP ), ALLOCATABLE :: gaux( : )
  REAL( DP )                 :: fac
  INTEGER                    :: ig, ipol
  !
  ! ... Bring rho to G space
  !
  ALLOCATE( rhoaux( dfftp%nnr ) )
  rhoaux( : ) = CMPLX( rho( : ), 0.D0 ) 
  !
  CALL fwfft('Dense', rhoaux, dfftp)
  !
  ! ... Compute total potential in G space
  !
  ALLOCATE( gaux( dfftp%nnr ) )
  !
  DO ipol = 1, 3
    !
    gaux(:) = CMPLX(0.d0,0.d0,kind=dp)
    !
    DO ig = gstart, ngm
      !
      fac = g(ipol,ig) / gg(ig)
      gaux(nl(ig)) = CMPLX(-AIMAG(rhoaux(nl(ig))),REAL(rhoaux(nl(ig))),kind=dp) * fac 
      !
    END DO
    !
    IF ( gamma_only ) THEN
      !
      gaux(nlm(:)) = &
        CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
       !
    END IF
    !
    ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of 
    !  V = e2 * fpi divided by the 2\pi/a factor missing in G  
    !
    fac = e2 * fpi / tpiba
    gaux = gaux * fac 
    !
    ! ... bring back to R-space, (\grad_ipol a)(r) ...
    !
    CALL invfft ('Dense', gaux, dfftp)
    !
    gradv(ipol,:) = REAL( gaux(:) )
    !
  ENDDO
  !
  DEALLOCATE(gaux)
  !
  RETURN
  !
  END SUBROUTINE gradv_h_of_rho_r
#endif
