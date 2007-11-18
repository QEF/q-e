!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
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
  USE gvect,            ONLY : nrxx, ngm
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
  REAL(DP), INTENT(IN) :: rho_core(nrxx)
    ! input: the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: vtxc, etxc, ehart, eth, charge, etotefield
    ! output: the integral V_xc * rho
    ! output: the E_xc energy
    ! output: the hartree energy
    ! output: the integral of the charge
  !
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
     !
     CALL add_efield( rho%of_r, v%of_r(1,is), etotefield, 0 )
     !
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
  USE gvect,            ONLY : nr1, nr2, nr3, nrxx, ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : xc, xc_spin, get_igcx, get_igcc
  USE scf,              ONLY : scf_type
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(nrxx)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(nrxx,nspin), kedtaur(nrxx,nspin), vtxc, etxc
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
  CALL start_clock( 'v_xc_meta' )
  !
  etxc   = 0.D0
  vtxc   = 0.D0
  v(:,:) = 0.D0
  rhoneg = 0.D0
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
!  use gvecp, only: ng => ngm
  USE kinds,            ONLY : DP
  USE gvect,            ONLY : nrxx, nrx1,nrx2,nrx3,nr1,nr2,nr3, &
                               g,nl,ngm
  USE gsmooth,          ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE scf,              ONLY : scf_type
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega, alat
  USE constants,        ONLY : e2
  IMPLICIT NONE
  !
  ! input
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP),INTENT(IN) :: rho_core(nrxx)
  COMPLEX(DP),INTENT(IN) :: rhog_core(ngm)
  REAL(DP),INTENT(OUT) :: etxc, vtxc, v(nrxx,nspin), kedtaur(nrxx,nspin)
!  integer nspin , nnr
!  real(8)  grho(nnr,3,nspin), rho(nnr,nspin),kedtau(nnr,nspin)
  ! output: excrho: exc * rho ;  E_xc = \int excrho(r) d_r
  ! output: rhor:   contains the exchange-correlation potential
  REAL(DP) :: zeta, rh, grh2
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
  ALLOCATE (grho(3,nrxx,nspin))
  ALLOCATE (h(3,nrxx,nspin))
  ALLOCATE (rhoout(nrxx,nspin))
  ALLOCATE (rhogsum(ngm,nspin))
  !
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
     CALL gradrho( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                   rhogsum(1,is), ngm, g, nl, grho(1,1,is) )
     !
  END DO
  !
  DO k = 1, nrxx
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
  ALLOCATE( dh( nrxx ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin
     !
     CALL grad_dot( nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                    nrxx, h(1,1,is), ngm, g, nl, alat, dh )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     rhoout(:,is)=rhoout(:,is)-fac*rho_core(:)
     vtxc = vtxc - SUM( dh(:) * rhoout(:,is) )
     !
  END DO
  DEALLOCATE(dh)
  !
  vtxc = omega * (vtxc / ( nr1 * nr2 * nr3 ))
  etxc = omega * etxc / ( nr1 * nr2 * nr3 )
  !
  CALL reduce( 1, vtxc )
  CALL reduce( 1, etxc )
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
  USE gvect,            ONLY : nr1, nr2, nr3, nrxx, ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : xc, xc_spin
  USE scf,              ONLY : scf_type
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(nrxx)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(nrxx,nspin), vtxc, etxc
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
     DO ir = 1, nrxx
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
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ! ... spin-polarized case
     !
     DO ir = 1, nrxx
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
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     ! ... noncolinear case
     !
     DO ir = 1,nrxx
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
  CALL reduce( 2, rhoneg )
  !
  rhoneg(:) = rhoneg(:) * omega / ( nr1*nr2*nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2E10.3)') rhoneg
  !
  ! ... energy terms, local-density contribution
  !
  vtxc = omega * vtxc / ( nr1*nr2*nr3 )
  etxc = omega * etxc / ( nr1*nr2*nr3 )
  !
  ! ... add gradient corrections (if any)
  !
  CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, v )
  !
  CALL reduce( 1, vtxc )
  CALL reduce( 1, etxc )
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
  USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                        nl, nlm, ngm, gg, gstart
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : omega
  USE wvfct,     ONLY : gamma_only
  USE cell_base, ONLY : tpiba2
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: rhog(ngm,nspin)
  REAL(DP),  INTENT(INOUT) :: v(nrxx,nspin)
  REAL(DP),    INTENT(OUT) :: ehart, charge
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux(:,:), aux1(:,:)
  REAL(DP)              :: rgtot_re, rgtot_im
  INTEGER               :: is, ig
  !
  CALL start_clock( 'v_h' )
  !
  ALLOCATE( aux( 2, nrxx ), aux1( 2, ngm ) )
  !
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
  CALL reduce( 1, charge )
  !
  ! ... calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  ehart     = 0.D0
  aux1(:,:) = 0.D0
  !
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
  CALL reduce( 1, ehart )
  ! 
  aux(:,:) = 0.D0
  !
  aux(:,nl(1:ngm)) = aux1(:,1:ngm)
  !
  IF ( gamma_only ) THEN
     !
     aux(1,nlm(1:ngm)) =   aux1(1,1:ngm)
     aux(2,nlm(1:ngm)) = - aux1(2,1:ngm)
     !
  END IF
  !
  ! ... transform hartree potential to real space
  !
  CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  ! ... add hartree potential to the xc potential
  !
  IF ( nspin == 4 ) THEN
     !
     v(:,1) = v(:,1) + aux(1,:)
     !
  ELSE
     !
     DO is = 1, nspin
        !
        v(:,is) = v(:,is) + aux(1,:)
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
