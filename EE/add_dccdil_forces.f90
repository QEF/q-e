!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!--------------------------------------------------------------------
      SUBROUTINE add_dccdil_forces( pot, force_vcorr, &
                                    nr1, nr2, nr3, nrx1,nrx2, nrx3, nrxx )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dpot = nabla pot
      ! ... to be subtracted from the Poisson-Boltzmann equation
      !
      USE kinds,      ONLY: DP
      USE io_global,  ONLY: stdout     
      USE cell_base,  ONLY: alat, omega, at
      USE ions_base,  ONLY : nat, ityp, ntyp => nsp, zv, tau
      USE constants,  ONLY: e2, fpi
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nr1, nr2, nr3, nrx1,nrx2,nrx3,nrxx
      REAL( DP )                  :: pot( nrxx )
      REAL( DP )                  :: force_vcorr(3, nat )
      !
      REAL( DP )                  :: force_gouy(3, nat )
      REAL( DP )                  :: dvdtao(3, nat )
      REAL( DP )                  :: delta1, delta2, delta3, l1, l2, l3 
      REAL( DP )                  :: sumfor
      REAL( DP )                  :: delta1n
      REAL( DP )                  :: delta2n
      REAL( DP )                  :: delta3n
      !
      !
      REAL( DP ), ALLOCATABLE     :: gradpot( :, : )
      REAL( DP ), ALLOCATABLE     :: gradpot5p( :, : )
      !
      INTEGER :: ipol, ir, m1, m2, m3, ion_pos, na
      INTEGER, EXTERNAL      :: compindex
      !
      ALLOCATE( gradpot( 3, nrxx ) )
      ALLOCATE( gradpot5p( 3, nrxx ) )
      !
      ! ... Initializes variables
      !
      l1 = alat * at( 1, 1 )
      l2 = alat * at( 2, 2 )
      l3 = alat * at( 3, 3 )

      delta1 = l1 / DBLE( nr1 )
      delta2 = l2 / DBLE( nr2 )
      delta3 = l3 / DBLE( nr3 )
!
      delta1n = alat * at( 1, 1 ) / DBLE( nr1 )
      delta2n = alat * at( 2, 2 ) / DBLE( nr2 )
      delta3n = alat * at( 3, 3 ) / DBLE( nr3 )

      gradpot(:,:) = 0.D0
      gradpot5p(:,:) = 0.D0

      !
      force_gouy(:,:) = 0.D0
      force_vcorr(:,:) = 0.D0
      dvdtao(:,:) = 0.D0

      call dvdr_tao(dvdtao, nr1, nr2, nr3, nrx1, nrx2, nrx3 )

      DO na = 1,nat

         m1 = INT(tau(1, na)*alat/delta1) + 1
         m2 = INT(tau(2, na)*alat/delta2) + 1
         m3 = INT(tau(3, na)*alat/delta3) + 1

        DO ipol = 1,3      
!             
           ion_pos = compindex (m1, m2, m3, nr1, nr2, nr3)

           force_vcorr(ipol, na)  = force_vcorr(ipol, na)   &
                          + zv(ityp( na ))* dvdtao(ipol,na)

        END DO
      END DO

      sumfor = 0.D0

      DO na = 1, nat
         
        sumfor = sumfor + force_vcorr(1,na)**2 + &
                          force_vcorr(2,na)**2 + &
                          force_vcorr(3,na)**2

      END DO

      sumfor = SQRT( sumfor )


      DEALLOCATE( gradpot )
      DEALLOCATE( gradpot5p )
      !
      RETURN
      !

9035  FORMAT(5X,'ATOM ',I3,' type ',I2,'   force = ',3F14.8)

!--------------------------------------------------------------------
      END SUBROUTINE add_dccdil_forces
!--------------------------------------------------------------------
