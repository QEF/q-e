!
! Copyright (C) 2007-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... original code written by Giovanni Cantele and Paolo Cazzato
! ... adapted to work in the parallel case by Carlo Sbraccia
! ... code for the calculation of the vacuum level written by Carlo Sbraccia
!
!#define _PRINT_ON_FILE
!
!---------------------------------------------------------------------------
SUBROUTINE makov_payne( etot )
  !---------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE ions_base, ONLY : nat, tau, ityp, zv
  USE cell_base, ONLY : at, bg
#ifdef __SOLVENT
  USE solvent_base, ONLY : do_solvent, pol_dipole, pol_quadrupole
#endif
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: etot
  !
  INTEGER  :: ia
  REAL(DP) :: x0(3), zvtot, qq
  REAL(DP) :: e_dipole(0:3), e_quadrupole
  !
  ! ... x0 is the center of charge of the system
  !
  zvtot = 0.D0
  x0(:) = 0.D0
  !
  DO ia = 1, nat
     !
     zvtot = zvtot + zv(ityp(ia))
     !
     x0(:) = x0(:) + tau(:,ia)*zv(ityp(ia))
     !
  END DO
  !
  x0(:) = x0(:) / zvtot
  !
  CALL compute_e_dipole( x0, e_dipole, e_quadrupole )
#ifdef __SOLVENT
  !
  IF ( do_solvent ) CALL compute_pol_dipole( x0, pol_dipole, pol_quadrupole )
#endif  
  !
  CALL write_dipole( etot, x0, e_dipole, e_quadrupole, qq )
  !
  CALL vacuum_level( x0, zvtot )
  !
  RETURN
  !
END SUBROUTINE makov_payne
!
!---------------------------------------------------------------------------
SUBROUTINE compute_e_dipole( x0, e_dipole, e_quadrupole )
  !---------------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, tau
  USE cell_base,  ONLY : at, bg, omega, alat
  USE scf,        ONLY : rho
  USE lsda_mod,   ONLY : nspin
  USE fft_base,   ONLY : dfftp
  USE mp_global,  ONLY : me_pool, intra_pool_comm, &
                         me_bgrp, intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: x0(3)
  REAL(DP), INTENT(OUT) :: e_dipole(0:3), e_quadrupole
  !
  REAL(DP) :: r(3), rhoir
  INTEGER  :: i, j, k, ip, ir, index, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  !
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
  e_dipole(:)  = 0.D0
  e_quadrupole = 0.D0
  !
  index0 = 0
  !
#if defined (__PARA)
  !
#ifdef __BANDS
  DO i = 1, me_bgrp
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
#else
  DO i = 1, me_pool
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
#endif
  !
#endif
  !
  DO ir = 1, dfftp%nnr
     !
     ! ... three dimensional indexes
     !
     index = index0 + ir - 1
     k     = index / (dfftp%nr1x*dfftp%nr2x)
     index = index - (dfftp%nr1x*dfftp%nr2x)*k
     j     = index / dfftp%nr1x
     index = index - dfftp%nr1x*j
     i     = index
     !
     DO ip = 1, 3
        r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
     END DO
     !
     r(:) = r(:) - x0(:)
     !
     ! ... minimum image convenction
     !
     CALL cryst_to_cart( 1, r, bg, -1 )
     !
     r(:) = r(:) - ANINT( r(:) )
     !
     CALL cryst_to_cart( 1, r, at, 1 )
     !
     rhoir = rho%of_r(ir,1)
     !
     IF ( nspin == 2 ) rhoir = rhoir + rho%of_r(ir,2)
     !
     ! ... dipole(0) = charge density
     !
     e_dipole(0) = e_dipole(0) + rhoir
     !
     DO ip = 1, 3
        !
        e_dipole(ip) = e_dipole(ip) + rhoir*r(ip)
        e_quadrupole = e_quadrupole + rhoir*r(ip)**2
        !
     END DO
     !
  END DO
  !
#ifdef __BANDS
  CALL mp_sum(  e_dipole(0:3) , intra_bgrp_comm )
  CALL mp_sum(  e_quadrupole  , intra_bgrp_comm )
#else
  CALL mp_sum(  e_dipole(0:3) , intra_pool_comm )
  CALL mp_sum(  e_quadrupole  , intra_pool_comm )
#endif
  !
  e_dipole(0) = e_dipole(0)*omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  DO ip = 1, 3
     e_dipole(ip) = e_dipole(ip)*omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) * alat
  END DO
  !
  e_quadrupole = e_quadrupole*omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) * alat**2
  !
  RETURN
  !
END SUBROUTINE compute_e_dipole
!
!---------------------------------------------------------------------------
SUBROUTINE write_dipole( etot, x0, dipole_el, quadrupole_el, qq )
  !---------------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : e2, pi, rytoev, au_debye
  USE ions_base,  ONLY : nat, ityp, tau, zv
  USE cell_base,  ONLY : at, bg, omega, alat, ibrav
  USE io_global,  ONLY : ionode
#ifdef __SOLVENT
  USE solvent_base, ONLY : do_solvent, pol_dipole, pol_quadrupole
#endif
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: etot
  REAL(DP), INTENT(IN)  :: x0(3)
  REAL(DP), INTENT(IN)  :: dipole_el(0:3), quadrupole_el
  REAL(DP), INTENT(OUT) :: qq
  !
  REAL(DP) :: dipole_ion(3), quadrupole_ion, dipole(3), quadrupole
  REAL(DP) :: zvia, zvtot
  REAL(DP) :: corr1, corr2, aa, bb
  INTEGER  :: ia, ip
  !
  ! ... Note that the definition of the Madelung constant used here
  ! ... differs from the "traditional" one found in the literature. See
  ! ... Lento, Mozos, Nieminen, J. Phys.: Condens. Matter 14 (2002), 2637-2645
  !
  REAL(DP), PARAMETER :: madelung(3) = (/ 2.8373D0, 2.8883D0, 2.885D0 /)
  !
  !
  IF ( .NOT. ionode ) RETURN
  !
  ! ... compute ion dipole moments
  !
  zvtot          = 0.D0
  dipole_ion     = 0.D0
  quadrupole_ion = 0.D0
  !
  DO ia = 1, nat
     !
     zvia = zv(ityp(ia))
     !
     zvtot = zvtot + zvia
     !
     DO ip = 1, 3
        !
        dipole_ion(ip) = dipole_ion(ip) + &
                         zvia*( tau(ip,ia) - x0(ip) )*alat
        quadrupole_ion = quadrupole_ion + &
                         zvia*( ( tau(ip,ia) - x0(ip) )*alat )**2
        !
     END DO
  END DO
  !
  ! ... compute ionic+electronic total charge, dipole and quadrupole moments
  !
  qq = -dipole_el(0) + zvtot
  !
  dipole(:)  = -dipole_el(1:3) + dipole_ion(:)
  quadrupole = -quadrupole_el  + quadrupole_ion
#ifdef __SOLVENT
  IF ( do_solvent ) THEN 
    qq = qq - pol_dipole(0)
    dipole(:)  = dipole(:) - pol_dipole(1:3) 
    quadrupole = quadrupole - pol_quadrupole  
  ENDIF
#endif
  !
  WRITE( stdout, '(/5X,"charge density inside the ", &
       &               "Wigner-Seitz cell:",3F14.8," el.")' ) dipole_el(0)
  !
  WRITE( stdout, &
         '(/5X,"reference position (x0):",5X,3F14.8," bohr")' ) x0(:)*alat
  !
  ! ... A positive dipole goes from the - charge to the + charge.
  !
  WRITE( stdout, '(/5X,"Dipole moments (with respect to x0):")' )
  WRITE( stdout, '( 5X,"Elect",3F9.4," au (Ha),",3F9.4," Debye")' ) &
      (-dipole_el(ip), ip = 1, 3), (-dipole_el(ip)*au_debye, ip = 1, 3 )
  WRITE( stdout, '( 5X,"Ionic",3F9.4," au (Ha),", 3F9.4," Debye")' ) &
      ( dipole_ion(ip),ip = 1, 3), ( dipole_ion(ip)*au_debye,ip = 1, 3 )
#ifdef __SOLVENT
  IF ( do_solvent ) &
    WRITE( stdout, '( 5X,"Diele",3F9.4," au (Ha),", 3F9.4," Debye")' ) &
      (-pol_dipole(ip),ip = 1, 3), (-pol_dipole(ip)*au_debye,ip = 1, 3 )
#endif
  WRITE( stdout, '( 5X,"Total",3F9.4," au (Ha),", 3F9.4," Debye")' ) &
      ( dipole(ip),    ip = 1, 3), ( dipole(ip)*au_debye,    ip = 1, 3 )
  !
  ! ... print the electronic, ionic and total quadrupole moments
  !
  WRITE( stdout, '(/5X,"Electrons quadrupole moment",F20.8," a.u. (Ha)")' )  &
      -quadrupole_el
  WRITE( stdout, '( 5X,"     Ions quadrupole moment",F20.8," a.u. (Ha)")' ) &
      quadrupole_ion
#ifdef __SOLVENT
  IF ( do_solvent ) &
    WRITE( stdout, '( 5X,"Dielectr. quadrupole moment",F20.8," a.u. (Ha)")' )  &
      -pol_quadrupole
#endif
  WRITE( stdout, '( 5X,"    Total quadrupole moment",F20.8," a.u. (Ha)")' ) &
      quadrupole
  !
  IF ( ibrav < 1 .OR. ibrav > 3 ) THEN
     call errore(' write_dipole', &
                  'Makov-Payne correction defined only for cubic lattices', 1)
     !
  END IF
  !
  ! ... Makov-Payne correction, PRB 51, 4014 (1995)
  ! ... Note that Eq. 15 has the wrong sign for the quadrupole term
  !
  corr1 = - madelung(ibrav) / alat * qq**2 / 2.0D0 * e2
  !
  aa = quadrupole
  bb = dipole(1)**2 + dipole(2)**2 + dipole(3)**2
  !
  corr2 = ( 2.D0 / 3.D0 * pi )*( qq*aa - bb ) / alat**3 * e2
  !
  ! ... print the Makov-Payne correction
  !
  WRITE( stdout, '(/,5X,"*********    MAKOV-PAYNE CORRECTION    *********")' )
  WRITE( stdout, &
         '(/5X,"Makov-Payne correction with Madelung constant = ",F8.4)' ) &
      madelung(ibrav)
  !
  WRITE( stdout,'(/5X,"Makov-Payne correction ",F14.8," Ry = ",F6.3, &
       &              " eV (1st order, 1/a0)")'   ) -corr1, -corr1*rytoev
  WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
       &              " eV (2nd order, 1/a0^3)")' ) -corr2, -corr2*rytoev
  WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
       &              " eV (total)")' ) -corr1-corr2, (-corr1-corr2)*rytoev
  !
  WRITE( stdout,'(/"!    Total+Makov-Payne energy  = ",F16.8," Ry")' ) &
      etot - corr1 - corr2
  !
  RETURN
  !
END SUBROUTINE write_dipole
!
!---------------------------------------------------------------------------
SUBROUTINE vacuum_level( x0, zion )
  !---------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE io_files,  ONLY : prefix
  USE constants, ONLY : e2, pi, tpi, fpi, rytoev, eps32
  USE gvect,     ONLY : g, gg, ngm, gstart, igtongl
  USE scf,       ONLY : rho
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : at, alat, tpiba, tpiba2
  USE ions_base, ONLY : nsp
  USE vlocal,    ONLY : strf, vloc
  USE mp_global, ONLY : intra_pool_comm, intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE control_flags, ONLY : gamma_only
  USE basic_algebra_routines, ONLY : norm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: x0(3)
  REAL(DP), INTENT(IN) :: zion
  !
  INTEGER                  :: i, ir, ig, first_point
  REAL(DP)                 :: r, dr, rmax, rg, phase, sinxx
  COMPLEX(DP), ALLOCATABLE :: vg(:)
  COMPLEX(DP)              :: vgig, qgig
  REAL(DP)                 :: vsph, qsph, qqr
  REAL(DP)                 :: absg, qq, vol, fac, rgtot_re, rgtot_im
  INTEGER,  PARAMETER      :: npts = 100
  REAL(DP), PARAMETER      :: x(3) = (/ 0.5D0, 0.0D0, 0.0D0 /), &
                              y(3) = (/ 0.0D0, 0.5D0, 0.0D0 /), &
                              z(3) = (/ 0.0D0, 0.0D0, 0.5D0 /)
  !
  !
  IF ( .NOT.gamma_only ) RETURN
  !
  rmax = norm( MATMUL( at(:,:), x(:) ) )
  !
  rmax = MIN( rmax, norm( MATMUL( at(:,:), y(:) ) ) )
  rmax = MIN( rmax, norm( MATMUL( at(:,:), z(:) ) ) )
  !
  rmax = rmax*alat
  !
  dr = rmax / DBLE( npts )
  !
  ALLOCATE( vg( ngm ) )
  !
  ! ... the local ionic potential
  !
  vg(:) = ( 0.D0, 0.D0 )
  !
  DO i = 1, nsp
     DO ig = 1, ngm
        vg(ig) = vg(ig) + vloc(igtongl(ig),i)*strf(ig,i)
     END DO
  END DO
  !
  ! ... add the hartree potential in G-space (NB: V(G=0)=0 )
  !
  DO ig = gstart, ngm
     !
     fac = e2*fpi / ( tpiba2*gg(ig) )
     !
     rgtot_re = REAL(  rho%of_g(ig,1) )
     rgtot_im = AIMAG( rho%of_g(ig,1) )
     !
     IF ( nspin == 2 ) THEN
        !
        rgtot_re = rgtot_re + REAL(  rho%of_g(ig,2) )
        rgtot_im = rgtot_im + AIMAG( rho%of_g(ig,2) )
        !
     END IF
     !
     vg(ig) = vg(ig) + CMPLX( rgtot_re, rgtot_im ,kind=DP)*fac
     !
  END DO
  !
  first_point = npts
  !
#if defined _PRINT_ON_FILE
  !
  first_point = 1
  !
  IF ( ionode ) THEN
     !
     OPEN( UNIT = 123, FILE = TRIM( prefix ) // ".E_vac.dat" )
     !
     WRITE( 123, '("# estimate of the vacuum level as a function of r")' )
     WRITE( 123, '("#",/,"#",8X,"r (bohr)", &
                  &8X,"E_vac (eV)",6X,"integrated charge")' )
     !
  END IF
  !
#endif
  !
  DO ir = first_point, npts
     !
     ! ... r is in atomic units
     !
     r = dr*ir
     !
     vol = ( 4.D0 / 3.D0 )*pi*(r*r*r)
     !
     vsph = ( 0.D0, 0.D0 )
     qsph = ( 0.D0, 0.D0 )
     qqr  = ( 0.D0, 0.D0 )
     !
     DO ig = 1, ngm
        !
        ! ... g vectors are in units of 2pi / alat :
        ! ... to go to atomic units g must be multiplied by 2pi / alat
        !
        absg = tpiba*SQRT( gg(ig) )
        !
        rg = r*absg
        !
        IF ( r == 0.D0 .AND. absg /= 0 ) THEN
           !
           sinxx = 1.D0
           !
        ELSE IF ( absg == 0 ) THEN
           !
           sinxx = 0.5D0
           !
        ELSE
           !
           sinxx = SIN( rg ) / rg
           !
        END IF
        !
        vgig = vg(ig)
        qgig = rho%of_g(ig,1)
        !
        IF ( nspin == 2 ) qgig = qgig + rho%of_g(ig,2)
        !
        ! ... add the phase factor corresponding to the translation of the
        ! ... origin by x0 (notice that x0 is in alat units)
        !
        phase = tpi*( g(1,ig)*x0(1) + g(2,ig)*x0(2) + g(3,ig)*x0(3) )
        !
        vgig = vgig*CMPLX( COS( phase ), SIN( phase ) ,kind=DP)
        qgig = qgig*CMPLX( COS( phase ), SIN( phase ) ,kind=DP)
        !
        vsph = vsph + 2.D0*REAL( vgig )*sinxx
        qsph = qsph + 2.D0*REAL( qgig )*sinxx
        !
        IF ( absg /= 0.D0 ) THEN
           !
           qqr = qqr + 2.D0*REAL( qgig )* &
                       ( fpi / absg**3 )*( SIN( rg ) - rg*COS( rg ) )
           !
        ELSE
           !
           qqr = qqr + REAL( qgig )*vol
           !
        END IF
        !
     END DO
     !
#ifdef __BANDS
     CALL mp_sum( vsph, intra_bgrp_comm )
     CALL mp_sum( qsph, intra_bgrp_comm )
     CALL mp_sum( qqr,  intra_bgrp_comm )
#else
     CALL mp_sum( vsph, intra_pool_comm )
     CALL mp_sum( qsph, intra_pool_comm )
     CALL mp_sum( qqr,  intra_pool_comm )
#endif
     !
     qq = ( zion - qqr )
     !
#if defined _PRINT_ON_FILE
     IF ( ionode ) &
        WRITE( 123, '(3(2X,F16.8))' ) r, ( vsph + e2*qq / r )*rytoev, qqr
#endif
     !
  END DO
  !
#if defined _PRINT_ON_FILE
  IF ( ionode ) CLOSE( UNIT = 123 )
#endif
  !
  WRITE( stdout, '(5X,"Corrected vacuum level    = ",F16.8," eV")' ) &
      ( vsph + e2*qq / rmax )*rytoev
  !
  DEALLOCATE( vg )
  !
  RETURN
  !
END SUBROUTINE vacuum_level

#ifdef __SOLVENT
!---------------------------------------------------------------------------
SUBROUTINE compute_pol_dipole( x0, pol_dipole, pol_quadrupole )
  !---------------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, tau
  USE cell_base,  ONLY : at, bg, omega, alat
  USE solvent_base, ONLY : epsinfty, rhopol
  USE fft_base,   ONLY : dfftp
  USE mp_global,  ONLY : me_pool, intra_pool_comm, &
                         me_bgrp, intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: x0(3)
  REAL(DP), INTENT(OUT) :: pol_dipole(0:3), pol_quadrupole
  !
  REAL(DP) :: r(3), rhoir
  INTEGER  :: i, j, k, ip, ir, index, index0
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  !
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
  pol_dipole(:)  = 0.D0
  pol_quadrupole = 0.D0
  !
  IF ( epsinfty .LE. 1.D0 ) RETURN
  !
  index0 = 0
  !
#if defined (__PARA)
  !
#ifdef __BANDS
  DO i = 1, me_bgrp
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
#else
  DO i = 1, me_pool
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
#endif
  !
#endif
  !
  DO ir = 1, dfftp%nnr
     !
     ! ... three dimensional indexes
     !
     index = index0 + ir - 1
     k     = index / (dfftp%nr1x*dfftp%nr2x)
     index = index - (dfftp%nr1x*dfftp%nr2x)*k
     j     = index / dfftp%nr1x
     index = index - dfftp%nr1x*j
     i     = index
     !
     DO ip = 1, 3
        r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
     END DO
     !
     r(:) = r(:) - x0(:)
     !
     ! ... minimum image convenction
     !
     CALL cryst_to_cart( 1, r, bg, -1 )
     !
     r(:) = r(:) - ANINT( r(:) )
     !
     CALL cryst_to_cart( 1, r, at, 1 )
     !
     rhoir = rhopol(ir)
     !
     ! ... dipole(0) = charge density
     !
     pol_dipole(0) = pol_dipole(0) + rhoir
     !
     DO ip = 1, 3
        !
        pol_dipole(ip) = pol_dipole(ip) + rhoir*r(ip)
        pol_quadrupole = pol_quadrupole + rhoir*r(ip)**2
        !
     END DO
     !
  END DO
  !
#ifdef __BANDS
  CALL mp_sum(  pol_dipole(0:3) , intra_bgrp_comm )
  CALL mp_sum(  pol_quadrupole  , intra_bgrp_comm )
#else
  CALL mp_sum(  pol_dipole(0:3) , intra_pool_comm )
  CALL mp_sum(  pol_quadrupole  , intra_pool_comm )
#endif
  !
  pol_dipole(0) = pol_dipole(0)*omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  DO ip = 1, 3
     pol_dipole(ip) = pol_dipole(ip)*omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) * alat
  END DO
  !
  pol_quadrupole = pol_quadrupole*omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) * alat**2
  !
  RETURN
  !
END SUBROUTINE compute_pol_dipole
!
#endif
