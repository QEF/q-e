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
MODULE makovpayne
IMPLICIT NONE
CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE makov_payne( etot, tau, rhor, rhog, strf, vloc, gamma_only, etot_in_hartree, &
                        output_in_hartree, vacuum_level )
  !---------------------------------------------------------------------------
  !
  !USE vlocal,    ONLY : strf, vloc
  !USE control_flags, ONLY : gamma_only
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE ions_base, ONLY : nat, ityp, zv
  USE fft_base,  ONLY : dfftp
  USE plugin_flags
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: etot
  REAL(DP), INTENT(IN) :: tau(:,:)
  !! atomic positions in cell coordinates
  REAL(DP), INTENT(IN) :: rhor(:)
  !! total density r space
  COMPLEX(DP), INTENT(IN) :: rhog(:)
  !! total density g space
  COMPLEX(DP), INTENT(IN) :: strf(:,:) 
  !! the structure factor
  REAL(DP), INTENT(IN) :: vloc(:,:)
  !! the local potential for each atom type
  LOGICAL, INTENT(IN) :: gamma_only, etot_in_hartree, output_in_hartree, &
                         vacuum_level
  !
  INTEGER  :: ia
  REAL(DP) :: x0(3), zvtot, qq
  REAL(DP) :: e_dipole(0:3), e_quadrupole(3)
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
  CALL compute_dipole( dfftp%nnr, rhor, x0, e_dipole, e_quadrupole)
  !
  IF (etot_in_hartree) THEN
     CALL write_dipole( etot*2.0_dp, tau, x0, e_dipole, e_quadrupole, qq, &
                        output_in_hartree )
  ELSE
     CALL write_dipole( etot, tau, x0, e_dipole, e_quadrupole, qq,&
                        output_in_hartree )
  ENDIF
  !
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
  !
  IF (vacuum_level) &
     CALL compute_vacuum_level( x0, zvtot, rhog, strf, vloc, gamma_only )
  !
  RETURN
  !
END SUBROUTINE makov_payne
!
!---------------------------------------------------------------------------
SUBROUTINE write_dipole( etot, tau, x0, dipole_el, quadrupole_el, qq, &
                output_in_hartree )
  !---------------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : e2, pi, rytoev, au_debye
  USE ions_base,  ONLY : nat, ityp, zv
  USE cell_base,  ONLY : at, omega, alat, ibrav
  USE io_global,  ONLY : ionode
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: etot
  REAL(DP), INTENT(IN)  :: tau(:,:)
  !! atomic positions in cell coordinates
  REAL(DP), INTENT(IN)  :: x0(3)
  REAL(DP), INTENT(IN)  :: dipole_el(0:3), quadrupole_el(3)
  REAL(DP), INTENT(OUT) :: qq
  LOGICAL, INTENT(IN)   :: output_in_hartree
  !
  REAL(DP) :: dipole_ion(3), quadrupole_ion(3), dipole(3), quadrupole(3)
  REAL(DP) :: zvia, zvtot
  REAL(DP) :: corr1, corr2, aa, bb
  INTEGER  :: ia, ip, ibrav_mp
  INTEGER, EXTERNAL :: at2ibrav
  !
  ! Note that the definition of the Madelung constant used here:
  !   Lento, Mozos, Nieminen, J. Phys.: Condens. Matter 14 (2002), 2637-2645
  ! differs from the "traditional" one found in the literature, e.g.:
  !   Leslie and Gillam https://doi.org/10.1088/0022-3719/18/5/005,
  !   Dabo et al. at https://doi.org/10.1103/PhysRevB.77.115139:
  ! because different definitions of the length parameter L are adopted
  !
  REAL(DP), PARAMETER :: madelung(3) = (/ 2.8373D0, 2.8883D0, 2.8885D0 /)
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
        quadrupole_ion(ip) = quadrupole_ion(ip) + &
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
  WRITE( stdout, '( 5X,"Total",3F9.4," au (Ha),", 3F9.4," Debye")' ) &
      ( dipole(ip),    ip = 1, 3), ( dipole(ip)*au_debye,    ip = 1, 3 )
  !
  ! ... print the electronic, ionic and total quadrupole moments
  !
  WRITE( stdout, '(/5X,"Electrons quadrupole moment",F20.8," a.u. (Ha)")' )  &
      -SUM(quadrupole_el(:))
  WRITE( stdout, '( 5X,"     Ions quadrupole moment",F20.8," a.u. (Ha)")' ) &
      SUM(quadrupole_ion(:))
  WRITE( stdout, '( 5X,"    Total quadrupole moment",F20.8," a.u. (Ha)")' ) &
      SUM(quadrupole(:))
  !
  ibrav_mp = ibrav
  IF ( ibrav .EQ. 0 ) ibrav_mp = at2ibrav(at(:, 1), at(:, 2), at(:, 3))
  IF ( ibrav_mp < 1 .OR. ibrav_mp > 3 ) THEN
     call errore(' write_dipole', &
                  'Makov-Payne correction defined only for cubic lattices', 1)
     !
  END IF
  !
  ! ... Makov-Payne correction, PRB 51, 4014 (1995)
  ! ... Note that Eq. 15 has the wrong sign for the quadrupole term
  !
  corr1 = - madelung(ibrav_mp) / alat * qq**2 / 2.0D0 * e2
  !
  aa = SUM(quadrupole(:))
  bb = dipole(1)**2 + dipole(2)**2 + dipole(3)**2
  !
  corr2 = ( 2.D0 / 3.D0 * pi )*( qq*aa - bb ) / alat**3 * e2
  !
  ! ... print the Makov-Payne correction
  !
  WRITE( stdout, '(/,5X,"*********    MAKOV-PAYNE CORRECTION    *********")' )
  WRITE( stdout, &
         '(/5X,"Makov-Payne correction with Madelung constant = ",F8.4)' ) &
      madelung(ibrav_mp)
  !
  IF (output_in_hartree) THEN
      WRITE( stdout,'(/5X,"Makov-Payne correction ",F14.8," Ha = ",F6.3, &
           &              " eV (1st order, 1/a0)")'   ) -corr1/2.0_dp, -corr1*rytoev
      WRITE( stdout,'( 5X,"                       ",F14.8," Ha = ",F6.3, &
           &              " eV (2nd order, 1/a0^3)")' ) -corr2/2.0_dp, -corr2*rytoev
      WRITE( stdout,'( 5X,"                       ",F14.8," Ha = ",F6.3, &
           &              " eV (total)")' ) (-corr1-corr2)/2.0_dp, (-corr1-corr2)*rytoev
      !
      WRITE( stdout,'(/"!    Total+Makov-Payne energy  = ",F16.8," Ha")' ) &
          (etot - corr1 - corr2 )/2.0_dp
  ELSE
      WRITE( stdout,'(/5X,"Makov-Payne correction ",F14.8," Ry = ",F6.3, &
           &              " eV (1st order, 1/a0)")'   ) -corr1, -corr1*rytoev
      WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
           &              " eV (2nd order, 1/a0^3)")' ) -corr2, -corr2*rytoev
      WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
           &              " eV (total)")' ) -corr1-corr2, (-corr1-corr2)*rytoev
      !
      WRITE( stdout,'(/"!    Total+Makov-Payne energy  = ",F16.8," Ry")' ) &
          etot - corr1 - corr2 
  ENDIF
  !
  RETURN
  !
END SUBROUTINE write_dipole
!
!---------------------------------------------------------------------------
SUBROUTINE compute_vacuum_level( x0, zion, rhog, strf, vloc, gamma_only )
  !---------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE io_files,  ONLY : prefix
  USE constants, ONLY : e2, pi, tpi, fpi, rytoev, eps32
  USE gvect,     ONLY : g, gg, ngm, gstart, igtongl
  USE cell_base, ONLY : at, alat, tpiba, tpiba2
  USE ions_base, ONLY : nsp
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE basic_algebra_routines, ONLY : norm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: x0(3)
  REAL(DP), INTENT(IN) :: zion
  COMPLEX(DP), INTENT(IN) :: rhog(:)
  !! total density
  COMPLEX(DP), INTENT(IN) :: strf(:,:) 
  !! the structure factor
  REAL(DP), INTENT(IN) :: vloc(:,:)
  !! the local potential for each atom type
  LOGICAL, INTENT(IN) :: gamma_only
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
     rgtot_re = REAL(  rhog(ig) )
     rgtot_im = AIMAG( rhog(ig) )
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
        qgig = rhog(ig)
        !
        ! ... add the phase factor corresponding to the translation of the
        ! ... origin by x0 (notice that x0 is in alat units)
        !
        phase = tpi*( g(1,ig)*x0(1) + g(2,ig)*x0(2) + g(3,ig)*x0(3) )
        !
        vgig = vgig*CMPLX( COS( phase ), SIN( phase ) ,kind=DP)
        qgig = qgig*CMPLX( COS( phase ), SIN( phase ) ,kind=DP)
        !
        ! ... vsph is the spherical average of the periodic electrostatic  
        ! ... potential on a sphere of radius r centered in x0 
        ! ... so this should be the monopole term in the potential  
        !
        vsph = vsph + 2.D0*REAL( vgig )*sinxx
        qsph = qsph + 2.D0*REAL( qgig )*sinxx
        !
        IF ( absg /= 0.D0 ) THEN
           !
           ! ... qqr is the integral of the electronic charge in the sphere
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
     CALL mp_sum( vsph, intra_bgrp_comm )
     CALL mp_sum( qsph, intra_bgrp_comm )
     CALL mp_sum( qqr,  intra_bgrp_comm )
     !
     ! ... qq is therefore the total (positive) charge  of the system
     !
     qq = ( zion - qqr )
     !
     ! ... that by Gauss theorem gives a monopole average potential on the 
     ! ... sphere seen by electrons of: - qq e2 / r
     ! ... so  (vsph + qq e2/r) should be the shift of the isolated molecule 
     ! ... monopole wrt the periodic potential
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
  ! ... one should see (if a range of r's are computed) that this corrections 
  ! ... should become a constant when the charge density of the molecule is decayed
  !
  WRITE( stdout, '(5X,"Corrected vacuum level    = ",F16.8," eV")' ) &
      ( vsph + e2*qq / rmax )*rytoev
  !
  DEALLOCATE( vg )
  !
  RETURN
  !
END SUBROUTINE compute_vacuum_level
END MODULE
