!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE report_mag
      !------------------------------------------------------------------------
      !! This subroutine prints out information about the local magnetization
      !! and/or charge, integrated around the atomic positions at points which
      !! are calculated in make_pointlists.
      !
      USE kinds,            ONLY: DP
      USE ions_base,        ONLY: nat, tau, ityp
      USE io_global,        ONLY: stdout
      use constants,        ONLY: pi
      USE scf,              ONLY: rho
      USE noncollin_module, ONLY: noncolin, mcons, i_cons
      USE lsda_mod,         ONLY: nspin
      !
      IMPLICIT NONE
      !
      REAL(DP) :: theta, phi, norm, norm1
      INTEGER :: ipol, iat
      REAL(DP) :: r1_loc(nat), m1_loc(nspin-1,nat)
      !
      ! get_local integrates on the previously determined points
      !
      CALL get_locals( r1_loc, m1_loc, rho%of_r )
      !
      IF (nspin == 2) THEN
         WRITE( stdout, * )
         WRITE( stdout, '(5X,"Magnetic moment per site:")' )
         DO iat = 1, nat
            WRITE(stdout,1020) iat, r1_loc(iat), m1_loc(1,iat), mcons(1,ityp(iat))
         END DO
         !
      ELSE IF (noncolin) THEN
         DO iat = 1, nat
            !
            ! norm is the length of the magnetic moment vector
            !
            norm = DSQRT( m1_loc(1,iat)**2+m1_loc(2,iat)**2+m1_loc(3,iat)**2 )
            !
            ! norm1 is the length of the projection of the mm vector into
            ! the xy plane
            !
            norm1 = DSQRT( m1_loc(1,iat)**2+m1_loc(2,iat)**2 )
            !
            ! calculate the polar angles of the magnetic moment
            !
            IF (norm > 1.d-10) THEN
               theta = ACOS(m1_loc(3,iat)/norm)
               IF (norm1 > 1.d-10) THEN
                  phi = ACOS(m1_loc(1,iat)/norm1)
                  IF (m1_loc(2,iat) < 0.d0) phi = - phi
               ELSE
                  phi = 2.d0*pi
               ENDIF
            ELSE
               theta = 2.d0*pi
               phi = 2.d0*pi
            ENDIF
            !
            ! go to degrees
            theta = theta*180.d0/pi
            phi = phi*180.d0/pi
            !
            WRITE( stdout,1010)
            WRITE( stdout,1011) iat,(tau(ipol,iat),ipol=1,3)
            WRITE( stdout,1014) r1_loc (iat)
            WRITE( stdout,1012) (m1_loc(ipol,iat),ipol=1,3)
            WRITE( stdout,1018) (m1_loc(ipol,iat)/r1_loc(iat),ipol=1,3)
            WRITE( stdout,1013) norm,theta,phi
            IF (i_cons==1) THEN
               WRITE( stdout,1015) (mcons(ipol,ityp(iat)),ipol=1,3)
            ELSEIF (i_cons==2) THEN
               WRITE( stdout,1017) 180.d0 * ACOS(mcons(3,ityp(iat)))/pi
            ENDIF
            WRITE( stdout,1010)
         ENDDO
         !
      END IF
      !
      !
 1010 FORMAT (/,1x,78('='))
 1011 FORMAT (5x,'atom number ',i4,' relative position : ',3f9.4)
 1012 FORMAT (5x,'magnetization :      ',3f12.6)
 1013 FORMAT (5x,'polar coord.: r, theta, phi [deg] : ',3f12.6)
 1014 FORMAT (5x,'charge : ',f12.6)
 1018 FORMAT (5x,'magnetization/charge:',3f12.6)
 1015 FORMAT (5x,'constrained moment : ',3f12.6) 
 1017 FORMAT (5x,'constrained theta [deg] : ',f12.6) 
 1020 FORMAT (5x,'atom: ',i4,4X,'charge: ',F9.4,4X,'magn: ',F9.4,4X,'constr: ',f9.4)
      !
END SUBROUTINE report_mag
