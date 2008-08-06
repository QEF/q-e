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
      SUBROUTINE add_dccdil_forces( vcomp, force_vcorr, &
                                    nr1, nr2, nr3, nrx1,nrx2, nrx3, nrxx )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dpot = nabla pot
      ! ... to be subtracted from the Poisson-Boltzmann equation
      !
      USE kinds,      ONLY: DP
      USE io_global,  ONLY: stdout     
      USE cell_base,  ONLY: alat, omega, at
      USE ions_base,  ONLY: nat, ityp, ntyp => nsp, zv, tau
      USE constants,  ONLY: e2, fpi
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nr1, nr2, nr3, nrx1,nrx2,nrx3,nrxx
      REAL( DP )                  :: vcomp( nrx1*nrx2*nrx3 )
      REAL( DP )                  :: force_vcorr(3, nat )
      !
      REAL( DP )                  :: dvdtao(3, nat )
      REAL( DP )                  :: sumfor
      !
      INTEGER :: na 
      !
      ! ... Initializes variables
      !
      force_vcorr(:,:) = 0.D0
      dvdtao(:,:) = 0.D0

      call dvdr_tao(dvdtao, vcomp, nr1, nr2, nr3, nrx1, nrx2, nrx3 )

      DO na = 1,nat
         force_vcorr(1:3, na)  = force_vcorr(1:3, na)   &
                               + zv(ityp( na )) * dvdtao(1:3,na)
      END DO

      sumfor = 0.D0

      DO na = 1, nat
        sumfor = sumfor + force_vcorr(1,na)**2 + &
                          force_vcorr(2,na)**2 + &
                          force_vcorr(3,na)**2
      END DO

      sumfor = SQRT( sumfor )

      RETURN
      !
9035  FORMAT(5X,'ATOM ',I3,' type ',I2,'   force = ',3F14.8)

!--------------------------------------------------------------------
      END SUBROUTINE add_dccdil_forces
!--------------------------------------------------------------------
