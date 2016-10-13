!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
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
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE ions_base,        ONLY : nat, tau
  USE ldaU,             ONLY : lda_plus_U 
  USE funct,            ONLY : dft_is_meta
  USE scf,              ONLY : scf_type
  USE cell_base,        ONLY : alat
  USE control_flags,    ONLY : ts_vdw
  USE tsvdw_module,     ONLY : tsvdw_calculate, UtsvdW
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
  INTEGER :: is, ir
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
  ! ... LDA+U: build up Hubbard potential 
  !
  if (lda_plus_u) then
     if(noncolin) then
        call v_hubbard_nc(rho%ns_nc,v%ns_nc,eth)
     else
        call v_hubbard(rho%ns,v%ns,eth)
     endif
  endif
  !
  ! ... add an electric field
  ! 
  DO is = 1, nspin_lsda
     CALL add_efield(v%of_r(1,is), etotefield, rho%of_r, .false. )
  END DO
  !
  ! ... add Tkatchenko-Scheffler potential (factor 2: Ha -> Ry)
  ! 
  IF (ts_vdw) THEN
     CALL tsvdw_calculate(tau*alat,rho%of_r)
     DO is = 1, nspin_lsda
        DO ir=1,dfftp%nnr
           v%of_r(ir,is)=v%of_r(ir,is)+2.0d0*UtsvdW(ir)
        END DO
     END DO
  END IF
  !
  CALL stop_clock( 'v_of_rho' )
  !
  RETURN
  !
END SUBROUTINE v_of_rho
!----------------------------------------------------------------------------
SUBROUTINE v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
  !----------------------------------------------------------------------------
  !
  ! ... Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g, nl,ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega, alat
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : xc, xc_spin, tau_xc, tau_xc_spin, get_meta
  USE scf,              ONLY : scf_type
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    ! the core charge in real space
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin), kedtaur(dfftp%nnr,nspin), &
                           vtxc, etxc
    ! v:      V_xc potential
    ! kedtau: local K energy density 
    ! vtxc:   integral V_xc * rho
    ! etxc:   E_xc energy
    !
    ! ... local variables
    !
  REAL(DP) :: zeta, rh
  INTEGER  :: k, ipol, is
  REAL(DP) :: ex, ec, v1x, v2x, v3x,v1c, v2c, v3c,                     &
  &           v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw,                &
  &           v3xup, v3xdw,v3cup, v3cdw,                               &
  &           arho, atau, fac, rhoup, rhodw, ggrho2, tauup,taudw          
       
  REAL(DP), DIMENSION(2)   ::    grho2, rhoneg
  REAL(DP), DIMENSION(3)   ::    grhoup, grhodw, v2cup, v2cdw
  
  !
  REAL(DP),    ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP),    ALLOCATABLE :: rhoout(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:,:)
  REAL(DP), PARAMETER      :: eps12 = 1.0d-12, zero=0._dp
  !
  !----------------------------------------------------------------------------
  !
  !
  CALL start_clock( 'v_xc_meta' )
  !
  !
  etxc      = zero
  vtxc      = zero
  v(:,:)    = zero
  rhoneg(:) = zero
  !
  !
  ALLOCATE (grho(3,dfftp%nnr,nspin))
  ALLOCATE (h(3,dfftp%nnr,nspin))
  ALLOCATE (rhoout(dfftp%nnr,nspin))
  ALLOCATE (rhogsum(ngm,nspin))
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
  do k = 1, dfftp%nnr
  
     do is = 1, nspin
        grho2 (is) = grho(1,k, is)**2 + grho(2,k,is)**2 + grho(3,k, is)**2
     end do
     
     if (nspin == 1) then
        !
        !    This is the spin-unpolarised case
        !
        arho = ABS (rho%of_r (k, 1) )

        atau = rho%kin_r(k,1) / e2  ! kinetic energy density in Hartree
        
        if ( (arho > eps8) .and. (grho2 (1) > eps12) .and. &
                                 (abs(atau) > eps8)) then
           
           call tau_xc (arho, grho2(1),atau, ex, ec, v1x, v2x,v3x,v1c, v2c,v3c)
           
           v(k, 1) =  (v1x + v1c )*e2 
           
           ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
           h(:,k,1) =  (v2x + v2c)*grho (:,k,1) *e2
           
           kedtaur(k,1)=  (v3x + v3c) * 0.5d0 * e2
           
           etxc = etxc +  (ex + ec) *e2 !* segno
           vtxc = vtxc + (v1x+v1c)*e2*arho
           
        else  
           h (:, k, 1) = zero  
           kedtaur(k,1)= zero
        end if
        
        if (rho%of_r (k, 1) < zero ) rhoneg(1) = rhoneg(1) - rho%of_r (k, 1)
        
     else
        !
        !    spin-polarised case
        !
        rhoup=rho%of_r(k, 1)
        rhodw=rho%of_r(k, 2)
        
        rh   = rhoup + rhodw 
        
        do ipol=1,3
            grhoup(ipol)=grho(ipol,k,1)
            grhodw(ipol)=grho(ipol,k,2)
        end do
        
        ggrho2  = ( grho2 (1) + grho2 (2) ) * 4._dp
        
        tauup = rho%kin_r(k,1) / e2
        taudw = rho%kin_r(k,2) / e2
        atau  = tauup + taudw

        if ((rh > eps8) .and. (ggrho2 > eps12) .and. (abs(atau) > eps8) ) then
                
        call tau_xc_spin (rhoup, rhodw, grhoup, grhodw, tauup, taudw, ex, ec, &
                      v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw, v1cup, v1cdw, &
                      v2cup, v2cdw, v3cup, v3cdw ) 
          !
          ! first term of the gradient correction : D(rho*Exc)/D(rho)
          !
          v(k, 1) =  (v1xup + v1cup)*e2
          v(k, 2) =  (v1xdw + v1cdw)*e2
          !
          ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
          !
          if (get_meta()==1) then  ! tpss functional
            !
            h(:,k,1) = (v2xup * grhoup(:) + v2cup(:)) * e2
            h(:,k,2) = (v2xdw * grhodw(:) + v2cdw(:)) * e2
            !
          else
            !
            h(:,k,1) = (v2xup + v2cup(1)) * grhoup(:) * e2
            h(:,k,2) = (v2xdw + v2cdw(1)) * grhodw(:) * e2
            !
          end if
          !
          kedtaur(k,1)=  (v3xup + v3cup) * 0.5d0 * e2
          kedtaur(k,2)=  (v3xdw + v3cdw) * 0.5d0 * e2
          !
          etxc = etxc + (ex + ec) * e2
          vtxc = vtxc + (v1xup+v1cup+v1xdw+v1cdw) * e2 * rh
          !
        else
          h(:,k,1) = zero
          h(:,k,2) = zero
          !
          kedtaur(k,1)=  zero
          kedtaur(k,2)=  zero

        end if
        
        if (rho%of_r (k, 1) < zero ) rhoneg(1) = rhoneg(1) - rho%of_r (k, 1)
        if (rho%of_r (k, 2) < zero ) rhoneg(2) = rhoneg(2) - rho%of_r (k, 2)
        
     end if
  end do
  !
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
  call mp_sum ( rhoneg, intra_bgrp_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  if ((rhoneg(1) > eps8) .or. (rhoneg(2) > eps8)) then
    write (stdout, '(/,5x, "negative rho (up,down): ", 2es10.3)') rhoneg(:)
  end if
  !
  vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) 
  etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  CALL mp_sum(  vtxc , intra_bgrp_comm )
  CALL mp_sum(  etxc , intra_bgrp_comm )
  !
  DEALLOCATE(grho)
  DEALLOCATE(h)
  DEALLOCATE(rhoout)
  DEALLOCATE(rhogsum)
  !
  RETURN
  !
END SUBROUTINE v_xc_meta
!
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
  USE funct,            ONLY : xc, xc_spin, nlc, dft_is_nonlocc
  USE scf,              ONLY : scf_type
  USE mp_bands,         ONLY : intra_bgrp_comm
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
           vtxc = vtxc + ( v(ir,1)*rho%of_r(ir,1) + v(ir,2)*rho%of_r(ir,2) )
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
  CALL mp_sum(  rhoneg , intra_bgrp_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2ES10.3)') rhoneg
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
  USE mp_bands,  ONLY: intra_bgrp_comm
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
  CALL mp_sum(  charge , intra_bgrp_comm )
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
     CALL mp_sum(  ehart , intra_bgrp_comm )
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
  !
  ! Computes Hubbard potential and Hubbard energy
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
                                   Hubbard_J0, Hubbard_beta
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout

  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) 
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) 
  REAL(DP), INTENT(OUT) :: eth
  REAL(DP) :: n_tot, n_spin, eth_dc, eth_u, mag2, effU
  INTEGER :: is, isop, is1, na, nt, m1, m2, m3, m4
  REAL(DP),    ALLOCATABLE :: u_matrix(:,:,:,:)

  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )

  eth    = 0.d0
  eth_dc = 0.d0
  eth_u  = 0.d0

  v_hub(:,:,:,:) = 0.d0

  if (lda_plus_u_kind.eq.0) then

    DO na = 1, nat
       nt = ityp (na)
       IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN
          IF (Hubbard_J0(nt).NE.0.d0) THEN
             effU = Hubbard_U(nt) - Hubbard_J0(nt)
          ELSE
             effU = Hubbard_U(nt)
          END IF  
          DO is = 1, nspin
             DO m1 = 1, 2 * Hubbard_l(nt) + 1
                eth = eth + ( Hubbard_alpha(nt) + 0.5D0 * effU ) * &
                              ns(m1,m1,is,na)
                v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                            ( Hubbard_alpha(nt) + 0.5D0 * effU )
                DO m2 = 1, 2 * Hubbard_l(nt) + 1
                   eth = eth - 0.5D0 * effU * &
                                       ns(m2,m1,is,na)* ns(m1,m2,is,na)
                   v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                        effU * ns(m2,m1,is,na)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       IF (Hubbard_J0(nt).NE.0.d0 .OR. Hubbard_beta(nt).NE.0.d0) THEN
          DO is=1, nspin
             IF (is .eq. 2) THEN
                isop = 1
             ELSE
                isop = 2
             END IF
             DO m1 = 1, 2 * Hubbard_l(nt) + 1
                IF ( is .eq. 1) THEN
                   eth = eth + Hubbard_beta(nt) * ns(m1,m1,is,na)
                   v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_beta(nt)
                   DO m2 = 1, 2 * Hubbard_l(nt) + 1
                      eth = eth + 0.5D0 * Hubbard_J0(nt) * &
                            ns(m2,m1,is,na)* ns(m1,m2,isop,na)
                      v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                                           Hubbard_J0(nt) * ns(m2,m1,isop,na)
                   END DO
                ELSE IF (is .eq. 2) THEN
                   eth = eth - Hubbard_beta(nt) * ns(m1,m1,is,na)
                   v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) - Hubbard_beta(nt)
                   DO m2 = 1, 2 * Hubbard_l(nt) + 1
                      eth = eth + 0.5D0 * Hubbard_J0(nt) * &
                            ns(m2,m1,is,na) * ns(m1,m2,isop,na)
                      v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                                           Hubbard_J0(nt) * ns(m2,m1,isop,na)
                   END DO
                END IF
             END DO
          END DO
       END IF
        
    END DO

    IF (nspin.EQ.1) eth = 2.d0 * eth

!-- output of hubbard energies:
    IF ( iverbosity > 0 ) THEN
      write(stdout,*) '--- in v_hubbard ---'
      write(stdout,'("Hubbard energy ",f9.4)') eth
      write(stdout,*) '-------'
    ENDIF
!--

  else

    DO na = 1, nat
       nt = ityp (na)
       IF (Hubbard_U(nt).NE.0.d0) THEN

!       initialize U(m1,m2,m3,m4) matrix 
          call hubbard_matrix (Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                               Hubbard_J(1,nt), u_matrix)

!---      total N and M^2 for DC (double counting) term
          n_tot = 0.d0
          do is = 1, nspin
            do m1 = 1, 2 * Hubbard_l(nt) + 1
              n_tot = n_tot + ns(m1,m1,is,na)
            enddo
          enddo
          if (nspin.eq.1) n_tot = 2.d0 * n_tot

          mag2  = 0.d0
          if (nspin.eq.2) then
            do m1 = 1, 2 * Hubbard_l(nt) + 1
              mag2 = mag2 + ns(m1,m1,1,na) - ns(m1,m1,2,na)
            enddo
          endif
          mag2  = mag2**2
!---

!---      hubbard energy: DC term

          eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*n_tot*(n_tot-1.d0) -       &
                                    Hubbard_J(1,nt)*n_tot*(0.5d0*n_tot-1.d0) - &
                                    0.5d0*Hubbard_J(1,nt)*mag2 )
!--
          DO is = 1, nspin

!---        n_spin = up/down N

            n_spin = 0.d0
            do m1 = 1, 2 * Hubbard_l(nt) + 1
              n_spin = n_spin + ns(m1,m1,is,na)
            enddo
!---

            DO m1 = 1, 2 * Hubbard_l(nt) + 1

!             hubbard potential: DC contribution  

              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_spin + &
                      0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot

!             +U contributions 

              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                do m3 = 1, 2 * Hubbard_l(nt) + 1
                  do m4 = 1, 2 * Hubbard_l(nt) + 1

                    if (nspin.eq.1) then
                      v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                                           2.d0*u_matrix(m1,m3,m2,m4)*ns(m3,m4,is,na)
                    else
                      do is1 = 1, nspin
                         v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                                          u_matrix(m1,m3,m2,m4)*ns(m3,m4,is1,na)
                      enddo
                    endif

                    v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                      u_matrix(m1,m3,m4,m2) * ns(m3,m4,is,na)

                    eth_u = eth_u + 0.5d0*(                            &
                              ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) )*  &
                              ns(m1,m3,is,na)*ns(m2,m4,is,na)                +  &
                       u_matrix(m1,m2,m3,m4)*ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )

                  enddo
                enddo
              ENDDO
            ENDDO

          ENDDO

       endif
    enddo

    if (nspin.eq.1) eth_u = 2.d0 * eth_u
    eth = eth_u - eth_dc

!-- output of hubbard energies:
    IF ( iverbosity > 0 ) THEN
      write(stdout,*) '--- in v_hubbard ---'
      write(stdout,'("Hubbard energies (dc, U, total) ",3f9.4)') eth_dc, eth_u, eth
      write(stdout,*) '-------'
    ENDIF
!--

  endif

  DEALLOCATE (u_matrix)
  RETURN

END SUBROUTINE v_hubbard
!-------------------------------------

!-------------------------------------
SUBROUTINE v_hubbard_nc(ns, v_hub, eth)
  !
  ! Noncollinear version of v_hubbard.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_J, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout

  IMPLICIT NONE
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) 
  COMPLEX(DP) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) 
  REAL(DP) :: eth, eth_dc, eth_noflip, eth_flip, psum, mx, my, mz, mag2
  
  INTEGER :: is, is1, js, i, j, na, nt, m1, m2, m3, m4
  COMPLEX(DP) :: n_tot, n_aux
  REAL(DP),    ALLOCATABLE :: u_matrix(:,:,:,:)

  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )

  eth        = 0.d0
  eth_dc     = 0.d0  
  eth_noflip = 0.d0  
  eth_flip   = 0.d0  

  v_hub(:,:,:,:) = 0.d0

  DO na = 1, nat  
     nt = ityp (na)  
     IF (Hubbard_U(nt).NE.0.d0) THEN  

!       initialize U(m1,m2,m3,m4) matrix 
        call hubbard_matrix (Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                             Hubbard_J(1,nt), u_matrix)

!---    total N and M^2 for DC (double counting) term
        n_tot = 0.d0
        mx    = 0.d0
        my    = 0.d0
        mz    = 0.d0
        do m1 = 1, 2 * Hubbard_l(nt) + 1
          n_tot = n_tot + ns(m1,m1,1,na) + ns(m1,m1,4,na)
          mx = mx + DBLE( ns(m1, m1, 2, na) + ns(m1, m1, 3, na) )
          my = my + 2.d0 * AIMAG( ns(m1, m1, 2, na) )
          mz = mz + DBLE( ns(m1, m1, 1, na) - ns(m1, m1, 4, na) )
        enddo  
        mag2 = mx**2 + my**2 + mz**2  
!---

!---    hubbard energy: DC term
        mx = REAL(n_tot)
        eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*mx*(mx-1.d0) -       &
                                  Hubbard_J(1,nt)*mx*(0.5d0*mx-1.d0) - &
                                  0.5d0*Hubbard_J(1,nt)*mag2 )   
!--
        DO is = 1, nspin  
           if (is.eq.2) then
            is1 = 3
           elseif (is.eq.3) then
            is1 = 2
           else
            is1 = is
           endif   

!---       hubbard energy:
           if (is1.eq.is) then

!           non spin-flip contribution
            DO m1 = 1, 2 * Hubbard_l(nt) + 1
             DO m2 = 1, 2 * Hubbard_l(nt) + 1
               do m3 = 1, 2 * Hubbard_l(nt) + 1
                do m4 = 1, 2 * Hubbard_l(nt) + 1

                  eth_noflip = eth_noflip + 0.5d0*(                            &
                             ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) )*  & 
                             ns(m1,m3,is,na)*ns(m2,m4,is,na)                +  &
                     u_matrix(m1,m2,m3,m4)*ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )

                enddo
               enddo
              ENDDO
             ENDDO

           else
!           spin-flip contribution
            DO m1 = 1, 2 * Hubbard_l(nt) + 1
             DO m2 = 1, 2 * Hubbard_l(nt) + 1
               do m3 = 1, 2 * Hubbard_l(nt) + 1
                do m4 = 1, 2 * Hubbard_l(nt) + 1

                  eth_flip = eth_flip - 0.5d0*u_matrix(m1,m2,m4,m3)*       &
                                    ns(m1,m3,is,na)*ns(m2,m4,is1,na) 

                enddo
               enddo
              ENDDO
             ENDDO

           endif
!---

!---       hubbard potential: non spin-flip contribution 
           if (is1.eq.is) then
            DO m1 = 1, 2 * Hubbard_l(nt) + 1
             DO m2 = 1, 2 * Hubbard_l(nt) + 1  

               do m3 = 1, 2 * Hubbard_l(nt) + 1
                do m4 = 1, 2 * Hubbard_l(nt) + 1
                  v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                    u_matrix(m1,m3,m2,m4)*( ns(m3,m4,1,na)+ns(m3,m4,4,na) ) 
                enddo 
               enddo

              ENDDO
             ENDDO
           endif
!---
            
!---       n_aux = /sum_{i} n_{i,i}^{sigma2, sigma1} for DC term
           n_aux = 0.d0
           do m1 = 1, 2 * Hubbard_l(nt) + 1
             n_aux = n_aux + ns(m1,m1,is1,na)  
           enddo
!---

           DO m1 = 1, 2 * Hubbard_l(nt) + 1  

!---          hubbard potential: DC contribution  
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_aux
              if (is1.eq.is) then
                v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                 0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot  
              endif
!---
              
!---          hubbard potential: spin-flip contribution
              DO m2 = 1, 2 * Hubbard_l(nt) + 1  

               do m3 = 1, 2 * Hubbard_l(nt) + 1
                do m4 = 1, 2 * Hubbard_l(nt) + 1
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                      u_matrix(m1,m3,m4,m2) * ns(m3,m4,is1,na) 
                enddo 
               enddo

              ENDDO
!---
           ENDDO
        ENDDO

     ENDIF
  ENDDO

  eth = eth_noflip + eth_flip - eth_dc

!-- output of hubbard energies:
  IF ( iverbosity > 0 ) THEN
    write(stdout,*) '--- in v_hubbard ---'
    write(stdout,'("Hub. E (dc, noflip, flip, total) ",4f9.4)') &
                                 eth_dc, eth_noflip, eth_flip, eth 
    write(stdout,*) '-------'
  ENDIF
!--

  DEALLOCATE (u_matrix)
  RETURN
END SUBROUTINE v_hubbard_nc
!-------------------------------------------

!----------------------------------------------------------------------------
SUBROUTINE v_h_of_rho_r( rhor, ehart, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential VH(r) from a density in R space n(r) 
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft
  USE gvect,           ONLY : nl, ngm
  USE lsda_mod,        ONLY : nspin
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rhor( dfftp%nnr, nspin )
  REAL( DP ), INTENT(INOUT)  :: v( dfftp%nnr, nspin )
  REAL( DP ), INTENT(OUT)    :: ehart, charge
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhog( : , : )
  COMPLEX( DP ), ALLOCATABLE :: aux( : )
  INTEGER :: is
  !
  ! ... bring the (unsymmetrized) rho(r) to G-space (use aux as work array)
  !
  ALLOCATE( rhog( ngm, nspin ) )
  ALLOCATE( aux( dfftp%nnr ) )
  DO is = 1, nspin
     aux(:) = CMPLX(rhor( : , is ),0.D0,kind=dp) 
     CALL fwfft ('Dense', aux, dfftp)
     rhog(:,is) = aux(nl(:))
  END DO
  DEALLOCATE( aux )
  !
  ! ... compute VH(r) from n(G) 
  !
  CALL v_h( rhog, ehart, charge, v )
  DEALLOCATE( rhog )
  !
  RETURN
  !
END SUBROUTINE v_h_of_rho_r
!----------------------------------------------------------------------------
SUBROUTINE gradv_h_of_rho_r( rho, gradv )
  !----------------------------------------------------------------------------
  !
  ! ... Gradient of Hartree potential in R space from a total 
  !     (spinless) density in R space n(r)
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE constants,       ONLY : fpi, e2
  USE control_flags,   ONLY : gamma_only
  USE cell_base,       ONLY : tpiba, omega
  USE gvect,           ONLY : nl, ngm, nlm, gg, gstart, g
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
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
  COMPLEX( DP ), ALLOCATABLE :: rgtot(:), vaux(:)
  REAL( DP )                 :: fac, eh_corr
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
    ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of 
    !  V = e2 * fpi divided by the 2\pi/a factor missing in G  
    !
    fac = e2 * fpi / tpiba
    gaux = gaux * fac 
    !
    ! ...add martyna-tuckerman correction, if needed
    ! 
    if (do_comp_mt) then
       ALLOCATE( vaux( ngm ), rgtot(ngm) )
       rgtot(1:ngm) = rhoaux(nl(1:ngm))
       CALL wg_corr_h (omega, ngm, rgtot, vaux, eh_corr)
       DO ig = gstart, ngm
         fac = g(ipol,ig) * tpiba
         gaux(nl(ig)) = gaux(nl(ig)) + CMPLX(-AIMAG(vaux(ig)),REAL(vaux(ig)),kind=dp)*fac 
       END DO
       DEALLOCATE( rgtot, vaux )
    end if
    !
    IF ( gamma_only ) THEN
      !
      gaux(nlm(:)) = &
        CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
       !
    END IF
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
  DEALLOCATE(rhoaux)
  !
  RETURN
  !
END SUBROUTINE gradv_h_of_rho_r
