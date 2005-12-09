!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
#if defined __T3E
#  define daxpy saxpy
#  define zaxpy caxpy
#endif

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sun Nov 21 11:19:43 MET 1999
!  ----------------------------------------------
!  BEGIN manual

      MODULE guess

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE guess_setup(tguess_inp)
!  SUBROUTINE guessc0(tk,c0,cm)
!  SUBROUTINE guessrho(rho,cm,c0,occ,ht) 
!  SUBROUTINE ucalc_kp(a,b,ngw,gzero,n,lambda)
!  SUBROUTINE ucalc(a,b,ngw,gzero,n,lambda)
!  ----------------------------------------------
!  END manual

! ...   declare modules
        USE kinds
        USE parallel_toolkit, ONLY: matmulp, cmatmulp, &
          diagonalize, cdiagonalize

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(DP), ALLOCATABLE :: rho_save( :, :, :, : )
   
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
          USE mp_global,        ONLY : nproc, mpime, group
          USE wave_types,       ONLY : wave_descriptor
          USE control_flags,    ONLY : force_pairing
          USE uspp,             ONLY : vkb, nkb

          IMPLICIT NONE

! ...     declare subroutine arguments
          COMPLEX(DP), INTENT(INOUT) ::  c0(:,:,:,:)
          COMPLEX(DP), INTENT(INOUT) ::  cm(:,:,:,:)
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
          nk  = cdesc%nkl

          IF( tguess ) THEN

! ...       uu(i,j)=<cm_i|c0_j>
            IF(tk) THEN

              ALLOCATE(cuu(n,n))
              ALLOCATE(ca(n,n))
              ALLOCATE(cap(n,n))
              ALLOCATE(crot(ngw,n))

              DO ik = 1, nk

                CALL ucalc_kp(cm(:,:,ik,1),c0(:,:,ik,1),ngw,cdesc%gzero,n,cuu)
                CALL cmatmulp('N','C',cuu,cuu,ca,n)
                CALL cdiagonalize(1,ca,costemp,cap,n,nproc,mpime)
                DO j=1,n
                  DO i=1,n
                    ca(i,j)=cap(i,n-j+1)
                  END DO
                END DO
                DO i=1,n
                  costh2(i)=1.0d0/sqrt(costemp(n-i+1))
                END DO
                CALL cmatmulp('N','N',cuu,ca,cap,n)
                DO j=1,n
                  DO i=1,n
                    cap(i,j)=cap(i,j) * costh2(i)
                  END DO
                END DO
                crot = (0.d0,0.d0)
                DO i=1,n
                  DO j=1,n
                    CALL zaxpy(ngw,ca(j,i),c0(1,j,ik,1),1,crot(1,i),1)
                  END DO
                END DO
                c0(:,:,ik,1) = crot
                crot = (0.d0,0.d0)
                DO i=1,n
                  DO j=1,n
                    CALL zaxpy(ngw,cap(i,j),cm(1,j,ik,1),1,crot(1,i),1)
                  END DO
                END DO
                cm(:,:,ik,1) = crot

              END DO

              DEALLOCATE(crot)
              DEALLOCATE(cap)
              DEALLOCATE(ca)
              DEALLOCATE(cuu)

            ELSE

              ALLOCATE(uu(n,n))
              ALLOCATE(a(n,n))
              ALLOCATE(ap(n,n))
              ALLOCATE(crot(ngw,n))

              CALL ucalc(cm(:,:,1,1),c0(:,:,1,1),ngw,cdesc%gzero,n,uu)
              CALL matmulp('T','N',uu,uu,a,n)
              CALL diagonalize(1,a,costemp,ap,n,nproc,mpime)
              DO j=1,n
                DO i=1,n
                  a(i,j)=ap(i,n-j+1)
                END DO
              END DO
              DO i=1,n
                costh2(i)=1.0d0/sqrt(costemp(n-i+1))
              END DO
              CALL matmulp('N','N',uu,a,ap,n)
              DO j=1,n
                DO i=1,n
                  ap(i,j)=ap(i,j) * costh2(i)
                END DO
              END DO
              crot = (0.d0,0.d0)
              DO i=1,n
                DO j=1,n
                  CALL daxpy(2*ngw,a(j,i),c0(1,j,1,1),1,crot(1,i),1)
                END DO
              END DO
              c0(:,:,1,1) = crot
              crot = (0.d0,0.d0)
              DO i=1,n
                DO j=1,n
                  CALL daxpy(2*ngw,ap(i,j),cm(1,j,1,1),1,crot(1,i),1)
                END DO
              END DO
              cm(:,:,1,1) = crot

              DEALLOCATE(crot)
              DEALLOCATE(ap)
              DEALLOCATE(a)
              DEALLOCATE(uu)

            END IF

            DO ik = 1, nk
              DO i=1,n
                ctemp(:) = 2.d0*c0(:,i,ik,1) - cm(:,i,ik,1)
                cm(:,i,ik,1) = c0(:,i,ik,1)
                c0(:,i,ik,1) = ctemp(:)
              END DO
            END DO

          ELSE

            cm = c0

          END IF

          CALL gram( vkb, bec, nkb, c0(1,1,1,1), SIZE(c0,1), cdesc%nbt( 1 ) )

          RETURN
          END SUBROUTINE guessc0
      
!  ----------------------------------------------
!  ----------------------------------------------
          SUBROUTINE guessrho(rho, desc, cm, c0, cdesc, occ, ht ) 

!  (describe briefly what this routine does...)
!  ----------------------------------------------

             USE cell_module, only: boxdimensions
             use charge_density, only: rhoofr
             use brillouin, only: kpoints, kp
             USE wave_types
             USE parameters, ONLY: nspinx
             USE charge_types, ONLY: charge_descriptor

! ...        declare subroutine argument
             REAL(DP), INTENT(OUT) :: rho(:,:,:,:)
             TYPE (charge_descriptor), INTENT(IN) :: desc
             COMPLEX(DP), INTENT(IN) :: c0(:,:,:,:), cm(:,:,:,:)
             TYPE (wave_descriptor), INTENT(IN) :: cdesc
             TYPE (boxdimensions), INTENT(IN) :: ht
             REAL(DP), INTENT(IN) :: occ(:,:,:)

! ...        declare other variables
             REAL(DP), ALLOCATABLE :: rho0( :, :, :, : )

             LOGICAL, SAVE :: tfirst = .TRUE.
             INTEGER :: ispin, nspin, nx, ny, nz

! ...        end of declarations
!  ----------------------------------------------

             nx = SIZE( rho, 1 )
             ny = SIZE( rho, 2 )
             nz = SIZE( rho, 3 )
             nspin = SIZE( rho, 4 )

             IF( tfirst ) THEN
               ALLOCATE( rho_save( nx, ny, nz, nspin ) )
               CALL rhoofr( 1, cm, cdesc, occ, rho_save, desc, ht)
               tfirst = .FALSE.
             END IF

             ALLOCATE( rho0( nx, ny, nz, nspin ) )
             CALL rhoofr( 1, c0, cdesc, occ, rho0, desc, ht)

             rho = 2.0d0 * rho0 - rho_save

             rho_save = rho0

             deallocate(rho0)

             RETURN
          END SUBROUTINE guessrho

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE ucalc_kp(a,b,ngw,gzero,n,lambda)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ...   declare modules
        USE mp

        IMPLICIT NONE

! ...   declare subroutine arguments
        INTEGER, INTENT(IN) :: n,ngw
        LOGICAL, INTENT(IN) :: gzero 
        COMPLEX(DP), INTENT(IN) :: a(ngw,*),b(ngw,*)
        COMPLEX(DP), INTENT(OUT) :: lambda(n,n)

! ...   declare other variables
        INTEGER   nproc,mpime,gid

! ...   end of declarations
!  ----------------------------------------------

        CALL mp_env(nproc,mpime,gid)
        CALL DGEMM('C','N',n,n,ngw,1.0d0,a,n,b,ngw,0.0d0,lambda,n)
        CALL mp_sum(lambda,gid)

        RETURN
      END SUBROUTINE ucalc_kp

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE ucalc(a,b,ngw,gzero,n,lambda)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ... declare modules
      USE mp

      IMPLICIT NONE

! ... declare subroutine arguments
      INTEGER n,ngw
      LOGICAL gzero 
      COMPLEX(DP) a(ngw,*),b(ngw,*)
      REAL(DP) lambda(n,n)

! ... declare other variables
      REAL(DP), ALLOCATABLE :: tmp(:,:)

      INTEGER i,j,jp1,jp2
      INTEGER nproc,mpime,gid

! ... end of declarations
!  ----------------------------------------------

      CALL mp_env(nproc,mpime,gid)
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

      CALL mp_sum(lambda,gid)

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

