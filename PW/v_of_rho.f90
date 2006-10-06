!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE v_of_rho( rho, rhog, rho_core, rhog_core, &
                     ehart, etxc, vtxc, etotefield, charge, v )
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
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho(nrxx,nspin), rho_core(nrxx)
    ! input: the valence charge
    ! input: the core charge
  COMPLEX(DP), INTENT(IN) :: rhog(ngm,nspin), rhog_core(ngm)
    ! input: the valence charge in reciprocal space
    ! input: the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: vtxc, etxc, ehart, charge, etotefield, v(nrxx,nspin) 
    ! output: the integral V_xc * rho
    ! output: the E_xc energy
    ! output: the hartree energy
    ! output: the integral of the charge
    ! output: the H+xc_up  potential
  !
  INTEGER :: is
  !
  CALL start_clock( 'v_of_rho' )
  !
  ! ... calculate exchange-correlation potential
  !
  CALL v_xc( rho, rhog, rho_core, rhog_core, etxc, vtxc, v )
  !
  ! ... add a magnetic field 
  !
  IF ( noncolin ) CALL add_bfield( v, rho )
  !
  ! ... calculate hartree potential
  !
  CALL v_h( rhog, ehart, charge, v )
  !
  ! ... add an electric field
  ! 
  DO is = 1, nspin
     !
     CALL add_efield( rho, v(1,is), etotefield, 0 )
     !
  END DO
  !
  CALL stop_clock( 'v_of_rho' )
  !
  RETURN
  !
END SUBROUTINE v_of_rho
!
!----------------------------------------------------------------------------
SUBROUTINE v_xc( rho, rhog, rho_core, rhog_core, etxc, vtxc, v )
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
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho(nrxx,nspin), rho_core(nrxx)
    ! the valence charge
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog(ngm,nspin), rhog_core(ngm)
    ! input: the valence charge in reciprocal space
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
        rhox = rho(ir,1) + rho_core(ir)
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
           vtxc = vtxc + v(ir,1) * rho(ir,1)
           !
        ENDIF
        !
        IF ( rho(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho(ir,1)
        !
     END DO
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ! ... spin-polarized case
     !
     DO ir = 1, nrxx
        !
        rhox = rho(ir,1) + rho(ir,2) + rho_core(ir)
        !
        arhox = ABS( rhox )
        !
        IF ( arhox > vanishing_charge ) THEN
           !
           zeta = ( rho(ir,1) - rho(ir,2) ) / arhox
           !
           IF ( ABS( zeta ) > 1.D0 ) zeta = SIGN( 1.D0, zeta )
           !
           IF ( rho(ir,1) < 0.D0 ) rhoneg(1) = rhoneg(1) - rho(ir,1)
           IF ( rho(ir,2) < 0.D0 ) rhoneg(2) = rhoneg(2) - rho(ir,2)
           !
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           !
           v(ir,:) = e2*( vx(:) + vc(:) )
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           !
           vtxc = vtxc + v(ir,1) * rho(ir,1) + v(ir,2) * rho(ir,2)
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
        amag = SQRT( rho(ir,2)**2 + rho(ir,3)**2 + rho(ir,4)**2 )
        !
        rhox = rho(ir,1) + rho_core(ir)
        !
        IF ( rho(ir,1) < 0.D0 )  rhoneg(1) = rhoneg(1) - rho(ir,1)
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
                 v(ir,ipol) = e2 * vs * rho(ir,ipol) / amag
                 !
                 vtxc = vtxc + v(ir,ipol) * rho(ir,ipol)
                 !
              END DO
              !
           END IF
           !
           etxc = etxc + e2*( ex + ec ) * rhox
           vtxc = vtxc + v(ir,1) * rho(ir,1)
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
  CALL gradcorr( rho, rhog, rho_core, rhog_core, etxc, vtxc, v )
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
  REAL(DP),    INTENT(OUT) :: v(nrxx,nspin), ehart, charge
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux(:,:), aux1(:,:)
  REAL(DP)              :: rgtot_re, rgtot_im
  INTEGER               :: is, ig
  !
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
