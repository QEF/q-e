!
! Copyright (C) 2005-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE report_mag(save_locals)
      !------------------------------------------------------------------------
      !! This subroutine prints out information about the local magnetization
      !! and/or charge, integrated around the atomic positions at points which
      !! are calculated in make_pointlists. Uses an optional variable in input
      !! either declare this interface or import if from pwcom module. 
      !
      USE kinds,            ONLY: DP
      USE ions_base,        ONLY: nat, tau, ityp
      USE io_global,        ONLY: stdout
      use constants,        ONLY: pi
      USE scf,              ONLY: rho
      USE noncollin_module, ONLY: noncolin, mcons, i_cons, r_m
      USE lsda_mod,         ONLY: nspin, local_charges, local_mag
      !
      IMPLICIT NONE
      LOGICAL,OPTIONAL,INTENT(IN) :: save_locals 
      !! if .TRUE. locals are saved  in two local_charges and local_mag of lsda_mod 
      !
      REAL(DP) :: theta, phi, norm, norm1
      INTEGER :: ipol, na, nt
      REAL(DP) :: r1_loc(nat), m1_loc(nspin-1,nat)
      !
      ! get_local integrates on the previously determined points
      !
      CALL get_locals( r1_loc, m1_loc, rho%of_r ) 
      IF (PRESENT(save_locals)) THEN
        IF (save_locals) THEN  
          IF (ALLOCATED (local_charges)) DEALLOCATE (local_charges)
          IF (ALLOCATED (local_mag))     DEALLOCATE (local_mag)
          ALLOCATE (local_charges, SOURCE = r1_loc)
          ALLOCATE (local_mag, SOURCE =  m1_loc)
        END IF
      END IF  
      !
      IF (nspin == 2) THEN
         WRITE( stdout, * )
         WRITE( stdout, '(5X,"Magnetic moment per site ", &
              &" (integrated on atomic sphere of radius R)")' )
         DO na = 1, nat
            nt = ityp(na)
            IF (i_cons > 0) THEN
               WRITE(stdout,1020) na, r_m(nt), r1_loc(na), m1_loc(1,na), mcons(1,nt)
            ELSE
               WRITE(stdout,1021) na, r_m(nt), r1_loc(na), m1_loc(1,na)
            END IF
         END DO
         !
      ELSE IF (noncolin) THEN
         DO na = 1, nat
            !
            ! norm is the length of the magnetic moment vector
            !
            norm = SQRT( m1_loc(1,na)**2+m1_loc(2,na)**2+m1_loc(3,na)**2 )
            !
            ! norm1 is the length of the projection of the mm vector into
            ! the xy plane
            !
            norm1 = SQRT( m1_loc(1,na)**2+m1_loc(2,na)**2 )
            !
            ! calculate the polar angles of the magnetic moment
            !
            IF (norm > 1.d-10) THEN
               theta = ACOS(m1_loc(3,na)/norm)
               IF (norm1 > 1.d-10) THEN
                  phi = ACOS(m1_loc(1,na)/norm1)
                  IF (m1_loc(2,na) < 0.d0) phi = - phi
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
            WRITE( stdout,1011) na,(tau(ipol,na),ipol=1,3)
            WRITE( stdout,1014) r1_loc (na), r_m(ityp(na))
            WRITE( stdout,1012) (m1_loc(ipol,na),ipol=1,3)
            WRITE( stdout,1018) (m1_loc(ipol,na)/r1_loc(na),ipol=1,3)
            WRITE( stdout,1013) norm,theta,phi
            IF (i_cons==1) THEN
               WRITE( stdout,1015) (mcons(ipol,ityp(na)),ipol=1,3)
            ELSEIF (i_cons==2) THEN
               WRITE( stdout,1017) 180.d0 * ACOS(mcons(3,ityp(na)))/pi
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
 1014 FORMAT (5x,'charge : ',f12.6,'  (integrated on a sphere of radius ',F5.3,')')
 1018 FORMAT (5x,'magnetization/charge:',3f12.6)
 1015 FORMAT (5x,'constrained moment : ',3f12.6) 
 1017 FORMAT (5x,'constrained theta [deg] : ',f12.6) 
1020  FORMAT (5x,'atom',i4,' (R=',F5.3,')  charge=',F8.4,'  magn=',F8.4,&
           & '   constr=',F8.4)
1021  FORMAT (5x,'atom',i4,' (R=',F5.3,')  charge=',F8.4,'  magn=',F8.4)
      !
END SUBROUTINE report_mag
