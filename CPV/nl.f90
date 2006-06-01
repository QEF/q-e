!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!----------------------------------------------
 MODULE nl
!----------------------------------------------

        USE kinds
        USE spherical_harmonics
        USE nl_base

        IMPLICIT NONE
        SAVE

        PRIVATE


        PUBLIC :: nlrh_m

!----------------------------------------------
 CONTAINS
!----------------------------------------------



   REAL(DP) FUNCTION nlrh_m( c0, cdesc, tforce, atoms, bec, becdr, eigr )

      !  this routine computes:
      !  Kleinman-Bylander pseudopotential terms (see nlsm1)
      !  enl: nonlocal potential contribution to total energy (see ene_nl)
      !  nonlocal potential contribution to forces on ions, see nlsm2
      !
      ! ... include modules

      USE wave_types,        ONLY: wave_descriptor
      USE pseudo_projector,  ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags,     ONLY: force_pairing
      USE pseudopotential,   ONLY: nspnl
      USE electrons_base,    ONLY: iupdwn, nupdwn
      USE uspp,              ONLY: becsum, nkb

      IMPLICIT NONE

      ! ... declare subroutine arguments

      COMPLEX(DP)                           :: eigr(:,:)     ! exp(i G dot r)
      COMPLEX(DP),           INTENT(INOUT)  :: c0(:,:,:)     ! wave functions
      TYPE (wave_descriptor), INTENT(IN)    :: cdesc         ! wave functions descriptor
      LOGICAL,               INTENT(IN)     :: tforce        ! if .TRUE. compute forces on ions
      TYPE(atoms_type),      INTENT(INOUT)  :: atoms         ! ions structure
      REAL(DP)                              :: bec(:,:)
      REAL(DP)                              :: becdr(:,:,:)

      REAL(DP)    :: ennl
      EXTERNAL     :: ennl

      ! ... declare other variables
      !
      INTEGER      :: iss, iss_wfc, i, j
      REAL(DP)    :: etmp
      REAL(DP), ALLOCATABLE :: btmp( :, :, : )
      REAL(DP), ALLOCATABLE :: fion( :, : )

!  end of declarations
!  ----------------------------------------------

      DO iss = 1, cdesc%nspin
         !
         iss_wfc = iss
         IF( force_pairing ) iss_wfc = 1
         !
         CALL nlsm1 ( cdesc%nbl( iss ), 1, nspnl, eigr(1,1),    &
                      c0( 1, 1, iss_wfc ), bec(1, iupdwn( iss ) ) )
         !
         IF( tforce ) THEN
            !
            ALLOCATE( btmp( nkb, nupdwn( iss ), 3 ) ) 
            !
            CALL nlsm2( cdesc%ngwl, nkb, nupdwn( iss ), eigr(1,1), &
                        c0( 1, 1, iss_wfc ), btmp( 1, 1, 1 ), .false. )
            !
            DO i = 1, 3
               DO j = iupdwn( iss ), iupdwn( iss ) + nupdwn( iss ) - 1
                  becdr( :, j , i ) = btmp( :, j - iupdwn( iss ) + 1, i ) 
               END DO
            END DO
            !
            DEALLOCATE( btmp )
            !
         END IF
         !
      END DO
      
      nlrh_m = ennl( becsum, bec )

      IF( tforce ) THEN
         !
         CALL force_nl( atoms%for, bec, becdr )
         !
      END IF
      !
      RETURN
   END FUNCTION nlrh_m

!  ----------------------------------------------

!  ----------------------------------------------


   SUBROUTINE nlsm2_s( ispin, wnl, atoms, eigr, c, cdesc, g2, gx, dfnl, kk)

      !  this routine computes the derivatives of the Kleinman-Bylander
      !  factors fnl, to be used for Hellmann-Feynman forces evaluation
      !

      ! ... declare modules
      USE pseudopotential,   ONLY: nspnl
      USE wave_types,        ONLY: wave_descriptor
      USE pseudo_projector,  ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE cell_base,         ONLY: tpiba
      USE uspp_param,        only: nh, lmaxkb
      USE uspp,              only: nhtol, nhtolm, indv

      IMPLICIT   NONE

      ! ... declare subroutine arguments
      COMPLEX(DP), INTENT(IN) :: c(:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (projector), INTENT(OUT) :: dfnl
      TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
      REAL(DP), INTENT(IN) :: wnl(:,:,:), g2(:), gx(:,:)
      COMPLEX(DP), INTENT(IN) :: eigr(:,:)
      INTEGER, INTENT(IN) :: kk
      INTEGER, INTENT(IN) :: ispin

      ! ... declare other variables
      REAL(DP), ALLOCATABLE :: gwork(:,:)
      INTEGER :: is, ia, igh, isa, ig, iss, ll, l, m, ngw, nb, lda, ldw, ldf
      INTEGER :: iy, ih, iv
      COMPLEX(DP), ALLOCATABLE :: auxc(:,:), gxtmp(:)
      COMPLEX(DP), PARAMETER :: ONE  = (1.0d0,0.0d0), ZERO = (0.0d0,0.0d0)
      ! ... (-i) * i^l
      COMPLEX(DP), PARAMETER :: csign(0:3) = (/ (0.0d0,-1.0d0), &
        (1.0d0,0.0d0), (0.0d0,1.0d0), (-1.0d0,0.0d0) /)


      ngw = cdesc%ngwl
      nb  = cdesc%nbl( ispin )
      IF( cdesc%gamma ) THEN
        lda = 2*SIZE(c, 1)
        ldw = 2*SIZE(c, 1)
        ldf = SIZE(dfnl%r, 1) * SIZE(dfnl%r, 2)
        dfnl%r = 0.0d0
      ELSE
        lda = SIZE(c, 1)
        ldw = SIZE(c, 1)
        ldf = SIZE(dfnl%c, 1) * SIZE(dfnl%c, 2)
        dfnl%c = 0.0d0
      END IF

      ALLOCATE(gwork(ngw, (lmaxkb+1)**2 ), gxtmp(ngw))

      CALL ylmr2( (lmaxkb+1)**2, ngw, gx, g2, gwork )
      !
      DO iy = 1, (lmaxkb+1)**2 
        gwork(1:ngw,iy) = tpiba * gx(kk,1:ngw) * gwork(1:ngw,iy)
      END DO

      iss = 1
      SPECS: DO is = 1, nspnl
        ALLOCATE(auxc(ngw,atoms%na(is)))
        LM: DO ih = 1, nh( is )
          iv  = indv  ( ih, is )
          iy  = nhtolm( ih, is )
          ll  = nhtol ( ih, is ) + 1
          l   = ll - 1
          igh = ih
          ! write( 6, * ) 'DEBUG = ', SUM( wnl( :, iv, is ) ), SUM( gwork( :, iy ) )
          gxtmp(1:ngw) = csign(l) * wnl(1:ngw,iv,is) * gwork(1:ngw,iy)
          DO ia = 1, atoms%na(is)
            auxc(1:ngw,ia) = gxtmp(1:ngw) * eigr(1:ngw,iss + ia - 1)
          END DO
          IF( cdesc%gamma ) THEN
             CALL DGEMM('T', 'N', atoms%na(is), nb, 2*ngw, 1.0d0, auxc(1,1), 2*ngw, &
                c(1,1), ldw, 0.0d0, dfnl%r(iss,igh,1), ldf)
           ELSE
             CALL ZGEMM('C', 'N', atoms%na(is), nb, ngw, one, auxc(1,1), ngw, &
                c(1,1), ldw, zero, dfnl%c(iss,igh,1), ldf)
           END IF
        END DO LM
        DEALLOCATE(auxc)
        iss = iss + atoms%na(is)
      END DO SPECS
      IF( cdesc%gamma ) CALL DSCAL(size(dfnl%r),2.0d0,dfnl%r(1,1,1),1)
      !write( 6, * ) 'DEBUG ==== ', SUM( dfnl%r )
      DEALLOCATE(gwork, gxtmp)
      RETURN
   END SUBROUTINE nlsm2_s


!----------------------------------------------
 END MODULE nl
!----------------------------------------------
