!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

MODULE guess

! ...   declare modules
        USE kinds
        USE parallel_toolkit, ONLY: rep_matmul_drv
        USE dspev_module,     ONLY: diagonalize

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(DP), ALLOCATABLE :: rho_save( :, : )
   
! ...   declare module-scope variables
        LOGICAL :: tguess

        PUBLIC :: guess_setup, guessc0, guess_closeup

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE guess_setup(tguess_inp)

          LOGICAL, INTENT(IN) :: tguess_inp
          tguess = tguess_inp

          RETURN
        END SUBROUTINE guess_setup


!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE guessc0( tk, bec, c0, cm, cdesc )

!  this subroutine updates the wavefunctions, leaving the new wave
!  functions in the KS base
!  ----------------------------------------------

! ...     declare modules
          USE mp_global,        ONLY : nproc_image, me_image, intra_image_comm
          USE wave_types,       ONLY : wave_descriptor
          USE control_flags,    ONLY : force_pairing
          USE uspp,             ONLY : vkb, nkb

          IMPLICIT NONE

! ...     declare subroutine arguments
          COMPLEX(DP), INTENT(INOUT) ::  c0(:,:)
          COMPLEX(DP), INTENT(INOUT) ::  cm(:,:)
          REAL(DP), INTENT(INOUT) ::  bec(:,:)
          TYPE (wave_descriptor), INTENT(IN) ::  cdesc
          LOGICAL, INTENT(IN) :: tk

! ...     declare other variables
          COMPLEX(DP) :: ctemp( cdesc%ngwl )
          REAL(DP)    :: alpha

          REAL(DP), ALLOCATABLE :: uu(:,:) 
          REAL(DP), ALLOCATABLE :: a(:,:)  
          REAL(DP), ALLOCATABLE :: ap(:,:) 
          COMPLEX(DP), ALLOCATABLE :: cuu(:,:) 
          COMPLEX(DP), ALLOCATABLE :: ca(:,:)  
          COMPLEX(DP), ALLOCATABLE :: cap(:,:) 

          COMPLEX(DP), ALLOCATABLE :: crot(:,:) 
          REAL(DP),   ALLOCATABLE :: aloc(:,:)
          COMPLEX(DP), ALLOCATABLE :: caloc(:,:)
          REAL(DP),   ALLOCATABLE :: evloc(:,:)
          COMPLEX(DP), ALLOCATABLE :: cevloc(:,:)
          REAL(DP),   ALLOCATABLE :: e(:)

          REAL(DP)   costh2 ( cdesc%ngwl )
          REAL(DP)   costemp( cdesc%ngwl )

          INTEGER jl, i,j,k,ig,h,n,ngw,nrl,ik,nk

! ...     end of declarations
!  ----------------------------------------------


          IF( force_pairing ) &
            CALL errore( ' guess ', ' force_pairing not yet implemented ', 1 )

          IF( cdesc%nspin > 1 ) &
            CALL errore( ' guess ', ' guess with spin not yet implemented ', 1 )

          n   = cdesc%nbl( 1 )
          ngw = cdesc%ngwl
          nk  = 1

          IF( tguess ) THEN

! ...       uu(i,j)=<cm_i|c0_j>

              ALLOCATE(uu(n,n))
              ALLOCATE(a(n,n))
              ALLOCATE(ap(n,n))
              ALLOCATE(crot(ngw,n))

              CALL ucalc(cm(:,:),c0(:,:),ngw,cdesc%gzero,n,uu)
              CALL rep_matmul_drv('T','N',n,n,n,1.0d0,uu,n,uu,n,0.0d0,a,n,intra_image_comm)
              CALL diagonalize(1,a,SIZE(a,1),costemp,ap,SIZE(ap,1),n,nproc_image,me_image)
              DO j=1,n
                DO i=1,n
                  a(i,j)=ap(i,n-j+1)
                END DO
              END DO
              DO i=1,n
                costh2(i)=1.0d0/sqrt(costemp(n-i+1))
              END DO
              CALL rep_matmul_drv('N','N',n,n,n,1.0d0,uu,n,a,n,0.0d0,ap,n,intra_image_comm)
              DO j=1,n
                DO i=1,n
                  ap(i,j)=ap(i,j) * costh2(i)
                END DO
              END DO
              crot = (0.d0,0.d0)
              DO i=1,n
                DO j=1,n
                  CALL daxpy(2*ngw,a(j,i),c0(1,j),1,crot(1,i),1)
                END DO
              END DO
              c0(:,:) = crot
              crot = (0.d0,0.d0)
              DO i=1,n
                DO j=1,n
                  CALL daxpy(2*ngw,ap(i,j),cm(1,j),1,crot(1,i),1)
                END DO
              END DO
              cm(:,:) = crot

              DEALLOCATE(crot)
              DEALLOCATE(ap)
              DEALLOCATE(a)
              DEALLOCATE(uu)


            DO ik = 1, nk
              DO i=1,n
                ctemp(:) = 2.d0*c0(:,i) - cm(:,i)
                cm(:,i) = c0(:,i)
                c0(:,i) = ctemp(:)
              END DO
            END DO

          ELSE

            cm = c0

          END IF

          CALL gram( vkb, bec, nkb, c0(1,1), SIZE(c0,1), cdesc%nbt( 1 ) )

          RETURN
          END SUBROUTINE guessc0
      
!  ----------------------------------------------
!  ----------------------------------------------
          SUBROUTINE guessrho(rho, cm, c0, cdesc, fi, ht ) 

!  (describe briefly what this routine does...)
!  ----------------------------------------------

             USE cell_base, only: boxdimensions
             use cp_interfaces, only: rhoofr
             USE wave_types

! ...        declare subroutine argument
             REAL(DP), INTENT(OUT) :: rho(:,:)
             COMPLEX(DP), INTENT(IN) :: c0(:,:), cm(:,:)
             TYPE (wave_descriptor), INTENT(IN) :: cdesc
             TYPE (boxdimensions), INTENT(IN) :: ht
             REAL(DP), INTENT(IN) :: fi(:)

! ...        declare other variables
             REAL(DP), ALLOCATABLE :: rho0( :, : )
             REAL(DP) :: edum, dedum(6)

             LOGICAL, SAVE :: tfirst = .TRUE.
             INTEGER :: ispin, nspin

! ...        end of declarations
!  ----------------------------------------------

             nspin = SIZE( rho, 2 )

             IF( tfirst ) THEN
               ALLOCATE( rho_save( SIZE( rho, 1 ), nspin ) )
               CALL rhoofr( 1, .false., cm, fi, rho_save, ht%deth, edum, dedum )
               tfirst = .FALSE.
             END IF

             ALLOCATE( rho0( SIZE( rho, 1 ), nspin ) )
             CALL rhoofr( 1, .false., c0, fi, rho0, ht%deth, edum, dedum )

             rho = 2.0d0 * rho0 - rho_save

             rho_save = rho0

             deallocate(rho0)

             RETURN
          END SUBROUTINE guessrho

!  ----------------------------------------------


!  ----------------------------------------------
      SUBROUTINE ucalc(a,b,ngw,gzero,n,lambda)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ... declare modules
      USE mp,        ONLY: mp_sum
      USE mp_global, ONLY: intra_image_comm

      IMPLICIT NONE

! ... declare subroutine arguments
      INTEGER n,ngw
      LOGICAL gzero 
      COMPLEX(DP) a(ngw,*),b(ngw,*)
      REAL(DP) lambda(n,n)

! ... declare other variables
      REAL(DP), ALLOCATABLE :: tmp(:,:)

      INTEGER i,j,jp1,jp2

! ... end of declarations
!  ----------------------------------------------

      ALLOCATE(tmp(n,2*ngw))

      DO i=1,n
        DO j=1,ngw
          jp1 = j + j - 1
          jp2 = j + j
          tmp(i,jp1) =  DBLE(a(j,i))
          tmp(i,jp2) = AIMAG(a(j,i))
        END DO
      END DO

      CALL DGEMM('N','N',n,n,2*ngw,2.0d0,tmp,n,b,2*ngw,0.0d0,lambda,n)
      IF(gzero) THEN
        CALL DGER(n,n,-1.0d0,a,2*ngw,b,2*ngw,lambda,n)
      END IF

      CALL mp_sum(lambda,intra_image_comm)

      DEALLOCATE(tmp)

      RETURN
      END SUBROUTINE ucalc

!  ----------------------------------------------

      SUBROUTINE guess_closeup()
        IF( ALLOCATED( rho_save ) ) DEALLOCATE( rho_save )
      RETURN
      END SUBROUTINE guess_closeup

!  ----------------------------------------------
!  ----------------------------------------------

END MODULE guess

