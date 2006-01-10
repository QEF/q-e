!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


!=----------------------------------------------------------------------------=!
   MODULE wave_functions
!=----------------------------------------------------------------------------=!

! ...   include modules
        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

          PUBLIC :: crot, proj, fixwave
          INTERFACE crot
            MODULE PROCEDURE crot_kp, crot_gamma
          END INTERFACE
          INTERFACE proj
            MODULE PROCEDURE proj_kp, proj_gamma, proj2
          END INTERFACE
          INTERFACE fixwave
            MODULE PROCEDURE fixwave_s, fixwave_v, fixwave_m
          END INTERFACE

          PUBLIC :: cp_kinetic_energy
          PUBLIC :: update_wave_functions, wave_rand_init

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!


   SUBROUTINE fixwave_s ( ispin, c, cdesc, kmask )

      USE wave_types, ONLY: wave_descriptor

      IMPLICIT NONE

      COMPLEX(DP), INTENT(INOUT) :: c(:,:)
      INTEGER, INTENT(IN) :: ispin
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      REAL(DP), INTENT(IN) :: kmask(:)
      INTEGER :: i

        IF( .NOT. cdesc%gamma ) THEN

          IF( SIZE( c, 1 ) /= SIZE( kmask ) ) &
            CALL errore( ' fixwave_s ', ' wrong dimensions ', 3 )

          DO i = 1, cdesc%nbl( ispin )
            c(:,i) = c(:,i) * kmask(:)
          END DO

        ELSE 

          IF( cdesc%gzero ) THEN
            DO i = 1, cdesc%nbl( ispin )
              c( 1, i ) = DBLE( c( 1, i ) )
            END DO
          END IF

        END IF

      RETURN
   END SUBROUTINE fixwave_s

!=----------------------------------------------------------------------------=!

   SUBROUTINE fixwave_v ( ispin, c, cdesc, kmask )
      USE wave_types, ONLY: wave_descriptor
      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: c(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      REAL(DP), INTENT(IN) :: kmask(:,:)
      INTEGER, INTENT(IN) :: ispin
      INTEGER :: i
        DO i = 1, cdesc%nkl
          CALL fixwave_s ( ispin, c(:,:,i), cdesc, kmask(:,i) )
        END DO
      RETURN
   END SUBROUTINE fixwave_v

!=----------------------------------------------------------------------------=!

   SUBROUTINE fixwave_m ( c, cdesc, kmask )
      USE wave_types, ONLY: wave_descriptor
      USE control_flags, ONLY: force_pairing
      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: c(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      REAL(DP), INTENT(IN) :: kmask(:,:)
      INTEGER :: i, j, nspin
      !
      nspin = cdesc%nspin
      IF( force_pairing ) nspin = 1
      !
      DO j = 1, nspin
        DO i = 1, cdesc%nkl
          CALL fixwave_s ( j, c(:,:,i,j), cdesc, kmask(:,i) )
        END DO
      END DO
      RETURN
   END SUBROUTINE fixwave_m

!=----------------------------------------------------------------------------=!

   REAL(DP) FUNCTION cp_kinetic_energy( ispin, cp, cm, cdesc, pmss, delt)

!  (describe briefly what this routine does...)
!  if ekinc_fp will hold the full electron kinetic energy (paired and unpaired) and
!  the function returns the paired electrons' kinetic energy only
!  ----------------------------------------------

! ... declare modules
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY:  group
      USE brillouin, ONLY: kpoints, kp
      USE wave_types, ONLY: wave_descriptor
      USE wave_base, ONLY: wave_speed2

      IMPLICIT NONE

! ... declare subroutine arguments
      COMPLEX(DP), INTENT(IN) :: cp(:,:,:), cm(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      INTEGER, INTENT( IN ) :: ispin
      REAL(DP), INTENT(IN) :: delt
      REAL(DP) :: pmss(:)

! ... declare other variables
      COMPLEX(DP) speed
      REAL(DP)  ekinc, ekinct, dt2, fact
      INTEGER    ib, j, ik

! ... end of declarations
!  ----------------------------------------------

      ekinct  = 0.d0
      dt2     = delt * delt 

      DO ik = 1, cdesc%nkl 

        ekinc  = 0.d0
        fact   = 1.0d0
        IF( cdesc%gamma .AND. cdesc%gzero ) fact =  0.5d0

        DO ib = 1, cdesc%nbl( ispin )
          ekinc = ekinc + wave_speed2( cp(:,ib,ik),  cm(:,ib,ik), pmss, fact )
        END DO

        IF( cdesc%gamma ) ekinc = ekinc * 2.0d0

        ekinct = ekinct + kp%weight(ik) * ekinc

      END DO

      CALL mp_sum( ekinct, group )

      cp_kinetic_energy = ekinct / (4.0d0 * dt2)

      RETURN
   END FUNCTION cp_kinetic_energy

!=----------------------------------------------------------------------------=!

   SUBROUTINE update_wave_functions(cm, c0, cp, cdesc)

      USE energies, ONLY: dft_energy_type
      USE wave_types, ONLY: wave_descriptor
      USE control_flags, ONLY: force_pairing

      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN) :: cp(:,:,:,:)
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:,:,:)
      COMPLEX(DP), INTENT(OUT) :: cm(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc

      INTEGER :: ispin, ik, nspin

      nspin = cdesc%nspin
      IF( force_pairing ) nspin = 1
      
      DO ispin = 1, nspin
        DO ik = 1, cdesc%nkl
          cm(:,:,ik,ispin) = c0(:,:,ik,ispin)
          c0(:,:,ik,ispin) = cp(:,:,ik,ispin)
        END DO
      END DO

      RETURN
   END SUBROUTINE update_wave_functions

!=----------------------------------------------------------------------------=!

   SUBROUTINE crot_gamma ( ispin, c0, cdesc, lambda, eig )

!  this routine rotates the wave functions to the Kohn-Sham base
!  it works with a block-like distributed matrix
!  of the Lagrange multipliers ( lambda ).
!  no replicated data are used, allowing scalability for large problems.
!  the layout of lambda is as follows :
!
!  (PE 0)                 (PE 1)               ..  (PE NPE-1)
!  lambda(1      ,1:nx)   lambda(2      ,1:nx) ..  lambda(NPE      ,1:nx)
!  lambda(1+  NPE,1:nx)   lambda(2+  NPE,1:nx) ..  lambda(NPE+  NPE,1:nx)
!  lambda(1+2*NPE,1:nx)   lambda(2+2*NPE,1:nx) ..  lambda(NPE+2*NPE,1:nx)
!
!  distributes lambda's rows across processors with a blocking factor
!  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row NPROC+1 to PE 1 and
!  so on).
!  nrl = local number of rows
!  ----------------------------------------------

! ... declare modules
      USE mp, ONLY: mp_bcast
      USE mp_global, ONLY: nproc, mpime, group
      USE wave_types, ONLY: wave_descriptor
      USE parallel_toolkit, ONLY: pdspev_drv, dspev_drv

      IMPLICIT NONE

! ... declare subroutine arguments
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      INTEGER, INTENT(IN) :: ispin
      REAL(DP) :: lambda(:,:)
      REAL(DP) :: eig(:)

! ... declare other variables
      INTEGER   ::  nx, ngw, nrl
      COMPLEX(DP), ALLOCATABLE :: c0rot(:,:)
      REAL(DP), ALLOCATABLE :: uu(:,:), vv(:,:)
      INTEGER   :: i, j, k, ip
      INTEGER   :: jl, nrl_ip

! ... end of declarations
!  ----------------------------------------------

      nx  = cdesc%nbl( ispin )
  
      IF( nx < 1 ) THEN
        RETURN
      END IF

      ngw = cdesc%ngwl
      nrl = SIZE(lambda, 1)
      ALLOCATE(uu(nrl,nx))
      ALLOCATE(vv(nrl,nx))
      ALLOCATE(c0rot(ngw,nx))

      c0rot = 0.0d0
      uu    = lambda

      CALL pdspev_drv( 'V', uu, nrl, eig, vv, nrl, nrl, nx, nproc, mpime)

      DEALLOCATE(uu)

      DO ip = 1, nproc

        nrl_ip = nx/nproc
        IF((ip-1).LT.mod(nx,nproc)) THEN
          nrl_ip = nrl_ip + 1
        END IF

        ALLOCATE(uu(nrl_ip,nx))
        IF(mpime.EQ.(ip-1)) THEN
          uu = vv
        END IF
        CALL mp_bcast(uu, (ip-1), group)

        j      = ip
        DO jl = 1, nrl_ip
          DO i = 1, nx
            CALL DAXPY(2*ngw,uu(jl,i),c0(1,j),1,c0rot(1,i),1)
          END DO
          j = j + nproc
        END DO
        DEALLOCATE(uu)

      END DO

      c0(:,:) = c0rot(:,:)

      DEALLOCATE(vv)
      DEALLOCATE(c0rot)

      RETURN
   END SUBROUTINE crot_gamma

!=----------------------------------------------------------------------------=!

   SUBROUTINE crot_kp ( ispin, ik, c0, cdesc, lambda, eig )

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ... declare modules
      USE mp, ONLY: mp_bcast
      USE mp_global, ONLY: nproc, mpime, group
      USE wave_types, ONLY: wave_descriptor
      USE parallel_toolkit, ONLY: pzhpev_drv, zhpev_drv

      IMPLICIT   NONE

! ... declare subroutine arguments
      INTEGER, INTENT(IN) :: ik
      INTEGER, INTENT(IN) :: ispin
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      COMPLEX(DP)  :: lambda(:,:)
      REAL(DP)      :: eig(:)

! ... declare other variables
      INTEGER   ngw, nx
      COMPLEX(DP), ALLOCATABLE :: c0rot(:,:)
      COMPLEX(DP), ALLOCATABLE :: vv(:,:)
      COMPLEX(DP), ALLOCATABLE :: uu(:,:)
      INTEGER   i,j,jl,nrl,ip,nrl_ip

! ... end of declarations
!  ----------------------------------------------

        nx  = cdesc%nbl( ispin )

        IF( nx < 1 ) THEN
          RETURN
        END IF

        ngw = cdesc%ngwl 
        nrl = SIZE(lambda,1)

        ALLOCATE( vv(nrl, nx), c0rot(ngw, nx) )
        c0rot = (0.d0,0.d0)

        ALLOCATE(uu(nrl,nx))
        uu    = lambda
        CALL pzhpev_drv( 'V', uu, nrl, eig, vv, nrl, nrl, nx, nproc, mpime)
        DEALLOCATE(uu)

        DO ip = 1,nproc
          j = ip
          nrl_ip = nx/nproc
          IF((ip-1).LT.mod(nx,nproc)) THEN
            nrl_ip = nrl_ip + 1
          END IF
          ALLOCATE(uu(nrl_ip,nx))
          IF(mpime.EQ.(ip-1)) THEN
            uu = vv
          END IF
          CALL mp_bcast(uu, (ip-1), group)
          DO jl=1,nrl_ip
            DO i=1,nx
              CALL ZAXPY(ngw,uu(jl,i),c0(1,j,ik),1,c0rot(1,i),1)
            END DO
            j = j + nproc
          END DO
          DEALLOCATE(uu)
        END DO

        c0(:,:,ik) = c0rot(:,:)

        DEALLOCATE(c0rot, vv)

        RETURN
   END SUBROUTINE crot_kp

!=----------------------------------------------------------------------------=!

   SUBROUTINE proj_gamma( ispin, a, adesc, b, bdesc, lambda)

!  projection A=A-SUM{B}<B|A>B
!  no replicated data are used, allowing scalability for large problems.
!  The layout of lambda is as follows :
!
!  (PE 0)                 (PE 1)               ..  (PE NPE-1)
!  lambda(1      ,1:nx)   lambda(2      ,1:nx) ..  lambda(NPE      ,1:nx)
!  lambda(1+  NPE,1:nx)   lambda(2+  NPE,1:nx) ..  lambda(NPE+  NPE,1:nx)
!  lambda(1+2*NPE,1:nx)   lambda(2+2*NPE,1:nx) ..  lambda(NPE+2*NPE,1:nx)
!
!  distribute lambda's rows across processors with a blocking factor
!  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row NPROC+1 to PE 1 and so on).
!  ----------------------------------------------

! ...   declare modules
        USE mp_global, ONLY: nproc,mpime,group
        USE wave_types, ONLY: wave_descriptor
        USE wave_base, ONLY: dotp

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(INOUT) :: a(:,:), b(:,:)
        TYPE (wave_descriptor), INTENT(IN) :: adesc, bdesc
        REAL(DP), OPTIONAL :: lambda(:,:)
        INTEGER, INTENT( IN ) :: ispin

! ...   declare other variables
        REAL(DP), ALLOCATABLE :: ee(:)
        INTEGER :: i, j, ngwc, jl
        INTEGER :: nstate_a, nstate_b
        COMPLEX(DP) :: alp

! ... end of declarations
!  ----------------------------------------------

        ngwc     = adesc%ngwl
        nstate_a = adesc%nbl( ispin )
        nstate_b = bdesc%nbl( ispin )

        IF( nstate_b < 1 ) THEN
          RETURN
        END IF

        ALLOCATE( ee( nstate_b ) )
        DO i = 1, nstate_a
          DO j = 1, nstate_b
            ee(j) = -dotp(adesc%gzero, ngwc, b(:,j), a(:,i))
          END DO
          IF( PRESENT(lambda) ) THEN
            IF( MOD( (i-1), nproc ) == mpime ) THEN
              DO j = 1, MIN( SIZE( lambda, 2 ), SIZE( ee ) )
                lambda( (i-1) / nproc + 1, j ) = ee(j)
              END DO
            END IF
          END IF
          DO j = 1, nstate_b
            alp = CMPLX(ee(j),0.0d0)
            CALL ZAXPY(ngwc,alp,b(1,j),1,a(1,i),1)
          END DO
        END DO
        DEALLOCATE(ee)

        RETURN
   END SUBROUTINE proj_gamma

!=----------------------------------------------------------------------------=!

   SUBROUTINE proj2( ispin, a, adesc, b, bdesc, c, cdesc)

!  projection A=A-SUM{B}<B|A>B-SUM{C}<C|A>
!
!  ----------------------------------------------

! ...   declare modules
        USE mp_global, ONLY: nproc,mpime,group
        USE wave_types, ONLY: wave_descriptor
        USE wave_base, ONLY: dotp

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(INOUT) :: a(:,:), b(:,:), c(:,:)
        TYPE (wave_descriptor), INTENT(IN) :: adesc, bdesc, cdesc
        INTEGER, INTENT( IN ) :: ispin

! ...   declare other variables
        COMPLEX(DP), ALLOCATABLE :: ee(:)
        INTEGER :: i, j, ngwc, jl
        INTEGER :: nstate_a, nstate_b, nstate_c

! ... end of declarations
!  ----------------------------------------------

        ngwc     = adesc%ngwl
        nstate_a = adesc%nbl( ispin )
        nstate_b = bdesc%nbl( ispin )
        nstate_c = cdesc%nbl( ispin )

        ALLOCATE( ee( MAX( nstate_b, nstate_c, 1 ) ) )

        DO i = 1, nstate_a
          DO j = 1, nstate_b
            IF( adesc%gamma ) THEN
              ee(j) = -dotp(adesc%gzero, ngwc, b(:,j), a(:,i))
            ELSE
              ee(j) = -dotp(ngwc, b(:,j), a(:,i))
            END IF
          END DO
! ...     a(:,i) = a(:,i) - (sum over j) e(i,j) b(:,j)
          DO j = 1, nstate_b
            CALL ZAXPY(ngwc, ee(j), b(1,j), 1, a(1,i), 1)
          END DO
        END DO

        DO i = 1, nstate_a
          DO j = 1, nstate_c
            IF( adesc%gamma ) THEN
              ee(j) = -dotp(adesc%gzero, ngwc, c(:,j), a(:,i))
            ELSE
              ee(j) = -dotp(ngwc, c(:,j), a(:,i))
            END IF
          END DO
          DO j = 1, nstate_c
            CALL ZAXPY(ngwc, ee(j), c(1,j), 1, a(1,i), 1)
          END DO
        END DO
        DEALLOCATE(ee)
        RETURN
   END SUBROUTINE proj2

!=----------------------------------------------------------------------------=!

   SUBROUTINE proj_kp( ispin, ik, a, adesc, b, bdesc, lambda)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ...   declare modules
        USE mp_global, ONLY: mpime, nproc
        USE wave_types, ONLY: wave_descriptor
        USE wave_base, ONLY: dotp

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(INOUT) :: a(:,:,:), b(:,:,:)
        TYPE (wave_descriptor), INTENT(IN) :: adesc, bdesc
        COMPLEX(DP), OPTIONAL  :: lambda(:,:)
        INTEGER, INTENT(IN) :: ik
        INTEGER, INTENT(IN) :: ispin

! ...   declare other variables
        COMPLEX(DP), ALLOCATABLE :: ee(:)
        INTEGER      :: ngwc, i, j
        INTEGER      :: nstate_a, nstate_b

! ... end of declarations
!  ----------------------------------------------

        ngwc     = adesc%ngwl
        nstate_a = adesc%nbl( ispin )
        nstate_b = bdesc%nbl( ispin )

        IF( nstate_b < 1 ) THEN
          RETURN
        END IF

! ...   lambda(i,j) = b(:,i,ik) dot a(:,j,ik)
        ALLOCATE(ee(nstate_b))
        DO i = 1, nstate_a
          DO j = 1, nstate_b
            ee(j) = -dotp(ngwc, b(:,j,ik), a(:,i,ik))
          END DO
          IF( PRESENT( lambda ) ) THEN
            IF(mod((i-1),nproc).EQ.mpime) THEN
              DO j = 1, MIN( SIZE( lambda, 2 ), SIZE( ee ) )
                lambda((i-1)/nproc+1,j) = ee(j)
              END DO
            END IF
          END IF
! ...     a(:,i,ik) = a(:,i,ik) - (sum over j) lambda(i,j) b(:,j,ik)
          DO j = 1, nstate_b
            CALL ZAXPY(ngwc, ee(j), b(1,j,ik), 1, a(1,i,ik), 1)
          END DO
        END DO
        DEALLOCATE(ee)

        RETURN
   END SUBROUTINE proj_kp



   SUBROUTINE wave_rand_init( cm )

!  this routine sets the initial wavefunctions at random
!  ----------------------------------------------

! ... declare modules
      USE mp, ONLY: mp_sum
      USE mp_wave, ONLY: splitwf
      USE mp_global, ONLY: mpime, nproc, root
      USE reciprocal_vectors, ONLY: ig_l2g, ngw, ngwt, gzero
      USE io_base, ONLY: stdout
      
      IMPLICIT NONE

! ... declare module-scope variables

! ... declare subroutine arguments 
      COMPLEX(DP), INTENT(OUT) :: cm(:,:)
      
      REAL(DP) :: rranf
      EXTERNAL rranf

! ... declare other variables
      INTEGER :: ntest, ig, ib
      REAL(DP) ::  rranf1, rranf2, ampre
      COMPLEX(DP), ALLOCATABLE :: pwt( : )

! ... end of declarations
!  ----------------------------------------------

! 
! ... Check array dimensions
      IF( SIZE( cm, 1 ) < ngw ) THEN 
        CALL errore(' wave_rand_init ', ' wrong dimensions ', 3)
      END IF

! ... Reset them to zero
!
      cm = 0.0d0

! ... initialize the wave functions in such a way that the values
! ... of the components are independent on the number of processors
!

      ampre = 0.01d0
      ALLOCATE( pwt( ngwt ) )

      ntest = ngwt / 4
      IF( ntest < SIZE( cm, 2 ) ) THEN
         ntest = ngwt
      END IF
      !
      ! ... assign random values to wave functions
      !
      DO ib = 1, SIZE( cm, 2 )
        pwt( : ) = 0.0d0
        DO ig = 3, ntest
          rranf1 = 0.5d0 - rranf()
          rranf2 = rranf()
          pwt( ig ) = ampre * CMPLX(rranf1, rranf2)
        END DO
        CALL splitwf ( cm( :, ib ), pwt, ngw, ig_l2g, mpime, nproc, 0 )
      END DO
      IF ( gzero ) THEN
        cm( 1, : ) = (0.0d0, 0.0d0)
      END IF

      DEALLOCATE( pwt )

      RETURN
    END SUBROUTINE wave_rand_init


!=----------------------------------------------------------------------------=!
   END MODULE wave_functions
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   MODULE wave_constrains
!=----------------------------------------------------------------------------=!

     ! ...   include modules
     USE kinds

     IMPLICIT NONE
     SAVE

     PRIVATE

     PUBLIC :: interpolate_lambda, update_lambda

     INTERFACE update_lambda
       MODULE PROCEDURE update_rlambda, update_clambda
     END INTERFACE

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

     SUBROUTINE interpolate_lambda( lambdap, lambda, lambdam )
       IMPLICIT NONE
       REAL(DP) :: lambdap(:,:), lambda(:,:), lambdam(:,:) 
       !
       ! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
       !
       lambdap(:,:) = 2.d0*lambda(:,:)-lambdam(:,:)
       lambdam(:,:)=lambda (:,:)
       lambda (:,:)=lambdap(:,:)
       RETURN
     END SUBROUTINE interpolate_lambda


     SUBROUTINE update_rlambda( i, lambda, c0, cdesc, c2 )
       USE electrons_module, ONLY: ib_owner, ib_local
       USE mp_global, ONLY: mpime
       USE mp, ONLY: mp_sum
       USE wave_base, ONLY: hpsi
       USE wave_types, ONLY: wave_descriptor
       IMPLICIT NONE
       REAL(DP) :: lambda(:,:)
       COMPLEX(DP) :: c0(:,:), c2(:)
       TYPE (wave_descriptor), INTENT(IN) :: cdesc
       INTEGER :: i
       !
       REAL(DP), ALLOCATABLE :: prod(:)
       INTEGER :: ibl
       !
       ALLOCATE( prod( SIZE( c0, 2 ) ) )
       prod = hpsi( cdesc%gzero, c0(:,:), c2 )
       CALL mp_sum( prod )
       IF( mpime == ib_owner( i ) ) THEN
           ibl = ib_local( i )
           lambda( ibl, : ) = prod( : )
       END IF
       DEALLOCATE( prod )
       RETURN
     END SUBROUTINE update_rlambda

     SUBROUTINE update_clambda( i, lambda, c0, cdesc, c2 )
       USE electrons_module, ONLY: ib_owner, ib_local
       USE mp_global, ONLY: mpime
       USE mp, ONLY: mp_sum
       USE wave_base, ONLY: hpsi
       USE wave_types, ONLY: wave_descriptor
       IMPLICIT NONE
       COMPLEX(DP) :: lambda(:,:)
       COMPLEX(DP) :: c0(:,:), c2(:)
       TYPE (wave_descriptor), INTENT(IN) :: cdesc
       INTEGER :: i
       !
       COMPLEX(DP), ALLOCATABLE :: prod(:)
       INTEGER :: ibl
       !
       ALLOCATE( prod( SIZE( c0, 2 ) ) )
       prod = hpsi( cdesc%gzero, c0(:,:), c2 )
       CALL mp_sum( prod )
       IF( mpime == ib_owner( i ) ) THEN
           ibl = ib_local( i )
           lambda( ibl, : ) = prod( : )
       END IF
       DEALLOCATE( prod )
       RETURN
     END SUBROUTINE update_clambda



!=----------------------------------------------------------------------------=!
   END MODULE wave_constrains
!=----------------------------------------------------------------------------=!
