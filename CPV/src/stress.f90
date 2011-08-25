!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
   SUBROUTINE pstress_conv( de3x3, de6, ainv )
!------------------------------------------------------------------------------!

      USE kinds,         ONLY: DP
      USE mp_global,     ONLY: intra_bgrp_comm
      USE mp,            ONLY: mp_sum
      USE stress_param,  ONLY: alpha, beta

      IMPLICIT NONE

      REAL(DP) :: de3x3(3,3)
      REAL(DP), INTENT(IN) :: de6(6)
      REAL(DP), INTENT(IN) :: ainv(3,3)
      REAL(DP) :: tmp(3,3)

      INTEGER  :: k

      DO k = 1, 6
        tmp( alpha(k), beta(k)  ) = de6(k)
        tmp( beta(k),  alpha(k) ) = tmp(alpha(k),beta(k))
      END DO

      de3x3 = MATMUL( tmp(:,:), TRANSPOSE( ainv(:,:) ) )

      CALL mp_sum( de3x3, intra_bgrp_comm )


      RETURN
   END SUBROUTINE 



!------------------------------------------------------------------------------!
   SUBROUTINE pseudo_stress_x( deps, epseu, gagb, sfac, dvps, rhoeg, omega )
!------------------------------------------------------------------------------!
      !
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: nsp
      USE gvecs,              ONLY: ngms
      USE electrons_base,     ONLY: nspin
      USE stress_param,       ONLY: dalbe
      USE cp_interfaces,      ONLY: stress_local

      IMPLICIT NONE

      REAL(DP),     INTENT(IN)  :: omega
      REAL(DP),     INTENT(OUT) :: deps(:)
      REAL(DP),     INTENT(IN)  :: gagb(:,:)
      COMPLEX(DP),  INTENT(IN)  :: rhoeg(:,:)
      COMPLEX(DP),  INTENT(IN)  :: sfac(:,:)
      REAL(DP),     INTENT(IN)  :: dvps(:,:)
      REAL(DP),     INTENT(IN)  :: epseu

      INTEGER     :: k
      COMPLEX(DP) :: rhets, depst(6)
      COMPLEX(DP), ALLOCATABLE :: rhoe( : )
      COMPLEX(DP), ALLOCATABLE :: drhoe( :, : )
      !
      ALLOCATE( drhoe( ngms, 6 ), rhoe( ngms ) )

      rhoe( 1:ngms ) = rhoeg( 1:ngms, 1 )
      IF( nspin > 1 ) rhoe( 1:ngms ) = rhoe( 1:ngms ) + rhoeg( 1:ngms, 2 )

      DO k = 1, 6
         drhoe( 1:ngms, k ) = - rhoe( 1:ngms ) * dalbe( k )
      END DO

      CALL stress_local( deps, epseu, gagb, sfac, rhoe, drhoe, omega )

      DEALLOCATE( drhoe, rhoe )

      RETURN
   END SUBROUTINE pseudo_stress_x



!------------------------------------------------------------------------------!
   SUBROUTINE stress_local_x( deps, epseu, gagb, sfac, rhoe, drhoe, omega )
!------------------------------------------------------------------------------!
      !
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: nsp
      USE gvect, ONLY: gstart
      USE gvecs,              ONLY: ngms
      USE electrons_base,     ONLY: nspin
      USE local_pseudo,       ONLY: vps, dvps

      IMPLICIT NONE

      REAL(DP),     INTENT(IN)  :: omega
      REAL(DP),     INTENT(OUT) :: deps(:)
      REAL(DP),     INTENT(IN)  :: gagb(:,:)
      COMPLEX(DP),  INTENT(IN)  :: rhoe(:)
      COMPLEX(DP),  INTENT(IN)  :: drhoe(:,:)
      COMPLEX(DP),  INTENT(IN)  :: sfac(:,:)
      REAL(DP),     INTENT(IN)  :: epseu

      INTEGER     :: ig,k,is, ispin
      COMPLEX(DP) :: dsvp, svp, depst(6)
      REAL(DP)    :: wz
      !
      depst = (0.d0,0.d0)

      wz = 2.0d0

      DO ig = gstart, ngms
         svp = 0.0d0
         DO is = 1, nsp
            svp = svp + sfac( ig, is ) * vps( ig, is )
         END DO
         depst = depst + wz * CONJG( drhoe( ig, : ) ) * svp
      END DO
      IF( gstart == 2 ) THEN
         svp = 0.0d0
         DO is = 1, nsp
            svp = svp + sfac( 1, is ) * vps( 1, is )
         END DO
         depst = depst + CONJG( drhoe( 1, : ) ) * svp
      END IF

      DO ig = gstart, ngms
         dsvp = 0.0d0
         DO is = 1, nsp
            dsvp = dsvp + sfac( ig, is ) * dvps( ig, is )
         END DO
         DO k = 1, 6
            depst( k ) = depst( k ) - wz * 2.0d0 * CONJG( rhoe( ig ) ) * dsvp * gagb( k, ig )
         END DO
      END DO

      deps = omega * DBLE( depst )

      RETURN
   END SUBROUTINE stress_local_x




!------------------------------------------------------------------------------!
   SUBROUTINE stress_kin_x( dekin, c0_bgrp, occ_bgrp ) 
!------------------------------------------------------------------------------!

!  this routine computes the kinetic energy contribution to the stress 
!  tensor
!
!  dekin(:) = - 2 (sum over i) f(i) * 
!    ( (sum over g) gagb(:,g) CONJG( c0(g,i) ) c0(g,i) )
!                       

      USE kinds,              ONLY: DP
      USE gvecw,              ONLY: q2sigma, ecfixed, qcutz, ngw
      USE constants,          ONLY: pi
      USE gvect, ONLY: gstart, gg, g
      USE cell_base,          ONLY: tpiba2
      USE electrons_base,     ONLY: nspin, iupdwn_bgrp, nupdwn_bgrp
      USE stress_param,       ONLY: alpha, beta

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP),    INTENT(OUT) :: dekin(:)
      COMPLEX(DP), INTENT(IN)  :: c0_bgrp(:,:)
      REAL(DP),    INTENT(IN)  :: occ_bgrp(:)

! ... declare other variables
      REAL(DP)  :: sk(6), scg, efac
      REAL(DP), ALLOCATABLE :: arg(:)
      INTEGER    :: ib, ig, ispin, iwfc

! ... end of declarations
!  ----------------------------------------------

      dekin = 0.0_DP
      ALLOCATE( arg( ngw ) ) 

      efac = 2.0d0 * qcutz / q2sigma / SQRT(pi)
      IF( efac > 0.0d0 ) THEN
        DO ig = gstart, ngw
          arg(ig) = 1.0d0 + efac * exp( -( ( tpiba2 *gg(ig) - ecfixed ) / q2sigma )**2 )
        END DO
      ELSE
        arg = 1.0d0
      END IF

      ! ... compute kinetic energy contribution

      DO ispin = 1, nspin
        DO ib = 1, nupdwn_bgrp( ispin )
          sk = 0.0_DP
          iwfc = ib + iupdwn_bgrp( ispin ) - 1
          DO ig = gstart, ngw
            scg = arg(ig) * CONJG( c0_bgrp( ig, iwfc ) ) * c0_bgrp( ig, iwfc )
            sk(1)  = sk(1) + scg * g( alpha( 1 ), ig ) * g( beta( 1 ), ig )
            sk(2)  = sk(2) + scg * g( alpha( 2 ), ig ) * g( beta( 2 ), ig )
            sk(3)  = sk(3) + scg * g( alpha( 3 ), ig ) * g( beta( 3 ), ig )
            sk(4)  = sk(4) + scg * g( alpha( 4 ), ig ) * g( beta( 4 ), ig )
            sk(5)  = sk(5) + scg * g( alpha( 5 ), ig ) * g( beta( 5 ), ig )
            sk(6)  = sk(6) + scg * g( alpha( 6 ), ig ) * g( beta( 6 ), ig )
          END DO
          dekin = dekin  + occ_bgrp( iwfc ) * sk * tpiba2
        END DO
      END DO
      dekin = - 2.0_DP * dekin
      DEALLOCATE(arg) 
      RETURN
   END SUBROUTINE stress_kin_x




!------------------------------------------------------------------------------!
   SUBROUTINE add_drhoph_x( drhot, sfac, gagb )
!------------------------------------------------------------------------------!
      !
      USE kinds,        ONLY: DP
      USE gvecs,        ONLY: ngms
      USE ions_base,    ONLY: nsp, rcmax
      USE local_pseudo, ONLY: rhops
      USE stress_param, ONLY: dalbe
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(INOUT) :: drhot( :, : )
      COMPLEX(DP), INTENT(IN) :: sfac( :, : )
      REAL(DP),    INTENT(IN) :: gagb( :, : )
      !
      INTEGER     :: ij, is, ig
      COMPLEX(DP) :: drhop
      !
      DO ij = 1, 6
         IF( dalbe( ij ) > 0.0d0 ) THEN
            DO is = 1, nsp
               DO ig = 1, ngms
                  drhot(ig,ij) = drhot(ig,ij) - sfac(ig,is)*rhops(ig,is)
               ENDDO
            END DO
         END IF
      END DO
      DO ig = 1, ngms
         drhop = 0.0d0
         DO is = 1, nsp
           drhop = drhop - sfac( ig, is ) * rhops(ig,is) * rcmax(is)**2 * 0.5D0
         END DO
         DO ij = 1, 6
             drhot(ig,ij) = drhot(ig,ij) - drhop * gagb( ij, ig )
         END DO
      END DO
      RETURN
   END SUBROUTINE add_drhoph_x




!------------------------------------------------------------------------------!
   SUBROUTINE stress_har_x(deht, ehr, sfac, rhoeg, gagb, omega ) 
!------------------------------------------------------------------------------!

      use kinds,              only: DP
      use ions_base,          only: nsp, rcmax
      use mp_global,          ONLY: me_bgrp, root_bgrp
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2
      USE gvect, ONLY: gstart
      USE gvecs,              ONLY: ngms
      USE gvect,              ONLY: ngm
      USE local_pseudo,       ONLY: rhops
      USE electrons_base,     ONLY: nspin
      USE stress_param,       ONLY: dalbe
      USE cp_interfaces,      ONLY: add_drhoph, stress_hartree

      IMPLICIT NONE

      REAL(DP),    INTENT(OUT) :: DEHT(:)
      REAL(DP),    INTENT(IN)  :: omega, EHR, gagb(:,:)
      COMPLEX(DP), INTENT(IN)  :: RHOEG(:,:)
      COMPLEX(DP), INTENT(IN)  :: sfac(:,:)

      COMPLEX(DP)    DEHC(6)
      COMPLEX(DP)    RHOP,DRHOP
      COMPLEX(DP)    RHET,RHOG,RHETS,RHOGS
      COMPLEX(DP)    CFACT
      COMPLEX(DP), ALLOCATABLE :: rhot(:), drhot(:,:)
      REAL(DP)       hgm1

      INTEGER       ig, is, k, ispin


      ALLOCATE( rhot( ngm ) )
      ALLOCATE( drhot( ngm, 6 ) )

      ! sum up spin components
      !
      DO ig = gstart, ngm
         rhot( ig ) = rhoeg( ig, 1 )
         IF( nspin > 1 ) rhot( ig ) = rhot( ig ) + rhoeg( ig, 2 )
      END DO
      !
      ! add Ionic pseudo charges  rho_I
      !
      DO is = 1, nsp
         DO ig = gstart, ngms
            rhot( ig ) = rhot( ig ) + sfac( ig, is ) * rhops( ig, is )
         END DO
      END DO

      ! add drho_e / dh
      !
      DO k = 1, 6
         IF( dalbe( k ) > 0.0d0 ) THEN
            drhot( :, k ) = - rhoeg( :, 1 )
            IF( nspin > 1 ) drhot( :, k ) = drhot( :, k ) + rhoeg( :, 2 )
         ELSE
            drhot( :, k ) = 0.0d0
         END IF
      END DO

      ! add drho_I / dh
      !
      CALL add_drhoph( drhot, sfac, gagb )

      CALL stress_hartree(deht, ehr, sfac, rhot, drhot, gagb, omega ) 

      DEALLOCATE( rhot, drhot )

      RETURN
   END SUBROUTINE stress_har_x




!------------------------------------------------------------------------------!
   SUBROUTINE stress_hartree_x(deht, ehr, sfac, rhot, drhot, gagb, omega ) 
!------------------------------------------------------------------------------!

      ! This subroutine computes: d E_hartree / dh  =
      !   E_hartree * h^t + 
      !   4pi omega rho_t * CONJG( rho_t ) / G^2 / G^2 * G_alpha * G_beta +
      !   4pi omega Re{ CONJG( rho_t ) * drho_t / G^2 }
      ! where:
      !   rho_t  = rho_e + rho_I
      !   drho_t = d rho_t / dh = -rho_e + d rho_hard / dh  + d rho_I / dh

      use kinds,              only: DP
      use ions_base,          only: nsp, rcmax
      use mp_global,          ONLY: me_bgrp, root_bgrp
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2
      USE gvect, ONLY: gstart, gg
      USE gvecs,              ONLY: ngms
      USE gvect,              ONLY: ngm
      USE local_pseudo,       ONLY: rhops
      USE electrons_base,     ONLY: nspin
      USE stress_param,       ONLY: dalbe

      IMPLICIT NONE

      REAL(DP),    INTENT(OUT) :: DEHT(:)
      REAL(DP),    INTENT(IN)  :: omega, EHR, gagb(:,:)
      COMPLEX(DP) :: rhot(:)  ! total charge: Sum_spin ( rho_e + rho_I )
      COMPLEX(DP) :: drhot(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)

      COMPLEX(DP)    DEHC(6)
      COMPLEX(DP)    CFACT
      REAL(DP), ALLOCATABLE :: hgm1( : )
      REAL(DP)    :: wz

      INTEGER       ig, is, k, iss

      DEHC  = (0.D0,0.D0)
      DEHT  = 0.D0

      wz = 2.0d0

      ALLOCATE( hgm1( ngm ) )

      hgm1( 1 ) = 0.0d0
      DO ig = gstart, ngm
         hgm1( ig ) = 1.D0 / gg(ig) / tpiba2
      END DO

      ! Add term  rho_t * CONJG( rho_t ) / G^2 * G_alpha * G_beta / G^2

      DO ig = gstart, ngm
         cfact = rhot( ig ) * CONJG( rhot( ig ) ) * hgm1( ig ) ** 2 
         dehc = dehc + cfact * gagb(:,ig)
      END DO

      ! Add term  2 * Re{ CONJG( rho_t ) * drho_t / G^2 }

      DO ig = gstart, ngm
         DO k = 1, 6
            dehc( k ) = dehc( k ) +  rhot( ig ) * CONJG( drhot( ig, k ) ) * hgm1( ig )
         END DO
      END DO

      ! term:  E_h * h^t

      if ( me_bgrp == root_bgrp ) then
        deht = wz * fpi * omega * DBLE(dehc) + ehr * dalbe
      else
        deht = wz * fpi * omega * DBLE(dehc)
      end if

      DEALLOCATE( hgm1 )

      RETURN
   END SUBROUTINE stress_hartree_x



!------------------------------------------------------------------------------!
      SUBROUTINE stress_debug(dekin, deht, dexc, desr, deps, denl, htm1)
!------------------------------------------------------------------------------!

        USE kinds,        ONLY: DP
        USE io_global,    ONLY: stdout
        USE mp_global,    ONLY: intra_bgrp_comm
        USE mp,           ONLY: mp_sum
        USE stress_param, ONLY: alpha, beta

        IMPLICIT NONE

        REAL(DP) :: dekin(6), deht(6), dexc(6), desr(6), deps(6), denl(6)
        REAL(DP) :: detot(6), htm1(3,3)
        REAL(DP) :: detmp(3,3)

        INTEGER :: k, i, j

        detot = dekin + deht + dexc + desr + deps + denl

        DO k=1,6
          detmp(alpha(k),beta(k)) = detot(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(tot)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = dekin(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(kin)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = deht(k) + desr(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(electrostatic)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = deht(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(h)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = desr(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(sr)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = deps(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(ps)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = denl(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(nl)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = dexc(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_bgrp_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(xc)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

      RETURN
      END SUBROUTINE stress_debug




!------------------------------------------------------------------------------!
      SUBROUTINE compute_gagb_x( gagb, g, ngm, tpiba2 )
!------------------------------------------------------------------------------!

         ! ... compute G_alpha * G_beta  

         USE kinds,        ONLY: DP
         USE stress_param, ONLY: alpha, beta

         IMPLICIT NONE

         INTEGER,  INTENT(IN)  :: ngm
         REAL(DP), INTENT(IN)  :: g(:,:)
         REAL(DP), INTENT(OUT) :: gagb(:,:)
         REAL(DP), INTENT(IN)  :: tpiba2

         INTEGER :: k, ig
      
!$omp parallel do default(shared), private(k)
         DO ig = 1, ngm          
            DO k = 1, 6
               gagb( k, ig ) = g( alpha( k ), ig ) * g( beta( k ), ig ) * tpiba2
            END DO
         END DO

      END SUBROUTINE compute_gagb_x
