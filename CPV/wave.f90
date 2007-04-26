!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"



!=----------------------------------------------------------------------------=!
     SUBROUTINE interpolate_lambda_x( lambdap, lambda, lambdam )
!=----------------------------------------------------------------------------=!
       USE kinds, ONLY: DP
       IMPLICIT NONE
       REAL(DP) :: lambdap(:,:,:), lambda(:,:,:), lambdam(:,:,:) 
       !
       ! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
       !
       lambdap= 2.d0*lambda - lambdam
       lambdam=lambda 
       lambda =lambdap
       RETURN
     END SUBROUTINE interpolate_lambda_x


!=----------------------------------------------------------------------------=!
     SUBROUTINE update_lambda_x( i, lambda, c0, c2, n, noff, tdist )
!=----------------------------------------------------------------------------=!
       USE kinds,              ONLY: DP
       USE electrons_module,   ONLY: ib_owner, ib_local
       USE mp_global,          ONLY: me_image, intra_image_comm
       USE mp,                 ONLY: mp_sum
       USE wave_base,          ONLY: hpsi
       USE reciprocal_vectors, ONLY: gzero
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, noff
       REAL(DP)            :: lambda(:,:)
       COMPLEX(DP)         :: c0(:,:), c2(:)
       INTEGER, INTENT(IN) :: i
       LOGICAL, INTENT(IN) :: tdist   !  if .true. lambda is distributed
       !
       REAL(DP), ALLOCATABLE :: prod(:)
       INTEGER :: ibl
       !
       ALLOCATE( prod( n ) )
       prod = hpsi( gzero, c0, SIZE( c0, 1 ), c2, n, noff )
       CALL mp_sum( prod, intra_image_comm )
       IF( tdist ) THEN
          IF( me_image == ib_owner( i ) ) THEN
             ibl = ib_local( i )
             lambda( ibl, : ) = prod( : )
          END IF
       ELSE
          lambda( i, : ) = prod( : )
       END IF
       DEALLOCATE( prod )
       RETURN
     END SUBROUTINE update_lambda_x




!=----------------------------------------------------------------------------=!
  subroutine elec_fakekine_x( ekincm, ema0bg, emass, c0, cm, ngw, n, noff, delt )
!=----------------------------------------------------------------------------=!
    !
    !  This subroutine computes the CP(fake) wave functions kinetic energy
    
    USE kinds,              only : DP
    use mp,                 only : mp_sum
    use mp_global,          only : intra_image_comm
    use reciprocal_vectors, only : gstart
    use wave_base,          only : wave_speed2
    !
    IMPLICIT NONE
    !
    integer, intent(in)      :: ngw    !  number of plane wave coeff.
    integer, intent(in)      :: n      !  number of bands
    integer, intent(in)      :: noff   !  offset for band index
    real(DP), intent(out)    :: ekincm
    real(DP), intent(in)     :: ema0bg( ngw ), delt, emass
    complex(DP), intent(in)  :: c0( ngw, n ), cm( ngw, n )
    !
    real(DP), allocatable :: emainv(:)
    real(DP) :: ftmp
    integer  :: i

    ALLOCATE( emainv( ngw ) )
    emainv = 1.0d0 / ema0bg
    ftmp = 1.0d0
    if( gstart == 2 ) ftmp = 0.5d0

    ekincm=0.0d0
    do i = noff, n + noff - 1
      ekincm = ekincm + 2.0d0 * wave_speed2( c0(:,i), cm(:,i), emainv, ftmp )
    end do
    ekincm = ekincm * emass / ( delt * delt )

    CALL mp_sum( ekincm, intra_image_comm )
    DEALLOCATE( emainv )

    return
  end subroutine elec_fakekine_x




!=----------------------------------------------------------------------------=!
   SUBROUTINE update_wave_functions_x( cm, c0, cp )
!=----------------------------------------------------------------------------=!

      USE kinds,              ONLY: DP

      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN)    :: cp(:,:)
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:)
      COMPLEX(DP), INTENT(OUT)   :: cm(:,:)

      cm(:,:) = c0(:,:)
      c0(:,:) = cp(:,:)

      RETURN
   END SUBROUTINE update_wave_functions_x



!=----------------------------------------------------------------------------=!
   SUBROUTINE crot_gamma ( c0, ngwl, nx, noff, lambda, nrl, eig )
!=----------------------------------------------------------------------------=!

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
      !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and
      !  so on).
      !  nrl = local number of rows
      !  ----------------------------------------------

      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE mp_global,        ONLY: nproc_image, me_image, intra_image_comm
      USE parallel_toolkit, ONLY: pdspev_drv, dspev_drv

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER, INTENT(IN) :: ngwl, nx, nrl, noff
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:)
      REAL(DP) :: lambda(:,:)
      REAL(DP) :: eig(:)

      ! ... declare other variables
      COMPLEX(DP), ALLOCATABLE :: c0rot(:,:)
      REAL(DP),    ALLOCATABLE :: uu(:,:), vv(:,:), ap(:)
      INTEGER   :: i, j, k, ip
      INTEGER   :: jl, nrl_ip

      IF( nx < 1 ) THEN
        RETURN
      END IF

      ALLOCATE( vv( nrl, nx ) )
      ALLOCATE( c0rot( ngwl, nx ) )

      c0rot = 0.0d0

      IF( nrl /= nx ) THEN

         ! Distributed lambda

         ALLOCATE( uu( nrl, nx ) )

         uu    = lambda

         CALL pdspev_drv( 'V', uu, nrl, eig, vv, nrl, nrl, nx, nproc_image, me_image)

         DEALLOCATE(uu)

         DO ip = 1, nproc_image

            nrl_ip = nx/nproc_image
            IF((ip-1).LT.mod(nx,nproc_image)) THEN
              nrl_ip = nrl_ip + 1
            END IF
 
            ALLOCATE(uu(nrl_ip,nx))
            IF(me_image.EQ.(ip-1)) THEN
              uu = vv
            END IF
            CALL mp_bcast(uu, (ip-1), intra_image_comm)
 
            j      = ip
            DO jl = 1, nrl_ip
              DO i = 1, nx
                CALL DAXPY(2*ngwl,uu(jl,i),c0(1,j+noff-1),1,c0rot(1,i),1)
              END DO
              j = j + nproc_image
            END DO
            DEALLOCATE(uu)
 
         END DO

      ELSE

         ! NON distributed lambda

         ALLOCATE( ap( nx * ( nx + 1 ) / 2 ) )

         K = 0
         DO J = 1, nx
            DO I = J, nx
               K = K + 1
               ap( k ) = lambda( i, j )
            END DO
          END DO

         CALL dspev_drv( 'V', 'L', nx, ap, eig, vv, nx )

         DEALLOCATE( ap )

         DO j = 1, nrl
            DO i = 1, nx
               CALL DAXPY( 2*ngwl, vv(j,i), c0(1,j+noff-1), 1, c0rot(1,i), 1 )
            END DO
         END DO

      END IF

      c0(:,noff:noff+nx-1) = c0rot(:,:)

      DEALLOCATE( vv )
      DEALLOCATE( c0rot )

      RETURN
   END SUBROUTINE crot_gamma


!=----------------------------------------------------------------------------=!
   SUBROUTINE protate_x ( c0, bec, c0rot, becrot, ngwl, nss, noff, lambda, nrl, &
                        na, nsp, ish, nh, np_rot, me_rot, comm_rot  )
!=----------------------------------------------------------------------------=!

      !  this routine rotates the wave functions using the matrix lambda
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
      !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and
      !  so on).
      !  nrl = local number of rows
      !  ----------------------------------------------

      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE mp_global,        ONLY: nproc_image, me_image, intra_image_comm
      USE parallel_toolkit, ONLY: pdspev_drv, dspev_drv

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER, INTENT(IN) :: ngwl, nss, nrl, noff
      INTEGER, INTENT(IN) :: na(:), nsp, ish(:), nh(:)
      INTEGER, INTENT(IN) :: np_rot, me_rot, comm_rot 
      COMPLEX(DP), INTENT(IN) :: c0(:,:)
      COMPLEX(DP), INTENT(OUT) :: c0rot(:,:)
      REAL(DP), INTENT(IN) :: lambda(:,:)
      REAL(DP), INTENT(IN) :: bec(:,:)
      REAL(DP), INTENT(OUT) :: becrot(:,:)

      ! ... declare other variables
      INTEGER   :: i, j, k, ip
      INTEGER   :: jl, nrl_ip, is, ia, jv, jnl, nj
      REAL(DP), ALLOCATABLE :: uu(:,:)

      IF( nss < 1 ) THEN
        RETURN
      END IF

      CALL start_clock('protate')

      becrot = 0.0d0
      c0rot  = 0.0d0

         DO ip = 1, np_rot

            nrl_ip = nss / np_rot
            IF( (ip-1) < mod( nss, np_rot ) ) THEN
              nrl_ip = nrl_ip + 1
            END IF
 
            ALLOCATE( uu( nrl_ip, nss ) )
            IF( me_rot .EQ. (ip-1) ) THEN
              uu = lambda( 1:nrl_ip, 1:nss )
            END IF
            CALL mp_bcast( uu, (ip-1), intra_image_comm)
 
            j      = ip
            DO jl = 1, nrl_ip
              DO i = 1, nss
                CALL DAXPY(2*ngwl,uu(jl,i),c0(1,j+noff-1),1,c0rot(1,i+noff-1),1)
              END DO

              do is=1,nsp
                 do jv=1,nh(is)
                    do ia=1,na(is)
                       jnl=ish(is)+(jv-1)*na(is)+ia
                       do i = 1, nss
                          becrot(jnl,i+noff-1) = becrot(jnl,i+noff-1)+ uu(jl, i) * bec( jnl, j+noff-1 )
                       end do
                    end do
                 end do
              end do

              j = j + np_rot
            END DO

            DEALLOCATE(uu)
 
         END DO

      CALL stop_clock('protate')

      RETURN
   END SUBROUTINE protate_x



!=----------------------------------------------------------------------------=!
   SUBROUTINE crot_gamma2 ( c0rot, c0, ngw, n, noffr, noff, lambda, nx, eig )
!=----------------------------------------------------------------------------=!

      !  this routine rotates the wave functions to the Kohn-Sham base
      !  it works with a block-like distributed matrix
      !  of the Lagrange multipliers ( lambda ).
      !
      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE mp_global,        ONLY: nproc_image, me_image, intra_image_comm
      USE parallel_toolkit, ONLY: dspev_drv

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER,     INTENT(IN)    :: ngw, n, nx, noffr, noff
      COMPLEX(DP), INTENT(INOUT) :: c0rot(:,:)
      COMPLEX(DP), INTENT(IN)    :: c0(:,:)
      REAL(DP),    INTENT(IN)    :: lambda(:,:)
      REAL(DP),    INTENT(OUT)   :: eig(:)

      ! ... declare other variables
      !
      REAL(DP), ALLOCATABLE :: vv(:,:), ap(:)
      INTEGER               :: i, j, k

      IF( nx < 1 ) THEN
        RETURN
      END IF

      ALLOCATE( vv( nx, nx ) )

      ! NON distributed lambda

      ALLOCATE( ap( nx * ( nx + 1 ) / 2 ) )

      K = 0
      DO J = 1, n
         DO I = J, n
            K = K + 1
            ap( k ) = lambda( i, j )
         END DO
      END DO

      CALL dspev_drv( 'V', 'L', n, ap, eig, vv, nx )

      DEALLOCATE( ap )

      DO i = 1, n
         c0rot( :, i+noffr-1 ) = 0.0d0
      END DO

      DO j = 1, n
         DO i = 1, n
            CALL DAXPY( 2*ngw, vv(j,i), c0(1,j+noff-1), 1, c0rot(1,i+noffr-1), 1 )
         END DO
      END DO

      DEALLOCATE( vv )

      RETURN
   END SUBROUTINE crot_gamma2



!=----------------------------------------------------------------------------=!
   SUBROUTINE proj_gamma( a, b, ngw, n, noff, lambda)
!=----------------------------------------------------------------------------=!

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
        !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and so on).
        !  ----------------------------------------------
         
! ...   declare modules
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: nproc_image, me_image, intra_image_comm
        USE wave_base,          ONLY: dotp
        USE reciprocal_vectors, ONLY: gzero

        IMPLICIT NONE

! ...   declare subroutine arguments
        INTEGER,     INTENT( IN )  :: ngw, n, noff
        COMPLEX(DP), INTENT(INOUT) :: a(:,:), b(:,:)
        REAL(DP),    OPTIONAL      :: lambda(:,:)

! ...   declare other variables
        REAL(DP), ALLOCATABLE :: ee(:)
        INTEGER :: i, j, jl
        COMPLEX(DP) :: alp

! ... end of declarations
!  ----------------------------------------------

        IF( n < 1 ) THEN
          RETURN
        END IF

        ALLOCATE( ee( n ) )
        DO i = 1, n
          DO j = 1, n
            ee(j) = -dotp( gzero, ngw, b(:,j+noff-1), a(:,i+noff-1) )
          END DO
          IF( PRESENT(lambda) ) THEN
            IF( MOD( (i-1), nproc_image ) == me_image ) THEN
              DO j = 1, n
                lambda( (i-1) / nproc_image + 1, j ) = ee(j)
              END DO
            END IF
          END IF
          DO j = 1, n
            alp = CMPLX(ee(j),0.0d0)
            CALL ZAXPY( ngw, alp, b(1,j+noff-1), 1, a(1,i+noff-1), 1 )
          END DO
        END DO
        DEALLOCATE(ee)

        RETURN
   END SUBROUTINE proj_gamma




!=----------------------------------------------------------------------------=!
   SUBROUTINE wave_rand_init_x( cm, n, noff )
!=----------------------------------------------------------------------------=!

      !  this routine sets the initial wavefunctions at random

! ... declare modules
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_sum
      USE mp_wave,            ONLY: splitwf
      USE mp_global,          ONLY: me_image, nproc_image, root_image, intra_image_comm
      USE reciprocal_vectors, ONLY: ig_l2g, ngw, ngwt, gzero
      USE io_global,          ONLY: stdout
      USE random_numbers,     ONLY: rranf
      
      IMPLICIT NONE

      ! ... declare subroutine arguments 
      INTEGER,     INTENT(IN)  :: n, noff
      COMPLEX(DP), INTENT(OUT) :: cm(:,:)

      ! ... declare other variables
      INTEGER :: ntest, ig, ib
      REAL(DP) ::  rranf1, rranf2, ampre
      COMPLEX(DP), ALLOCATABLE :: pwt( : )

      ! ... Check array dimensions

      IF( SIZE( cm, 1 ) < ngw ) THEN 
        CALL errore(' wave_rand_init ', ' wrong dimensions ', 3)
      END IF

      ! ... Reset them to zero
 
      cm( :, noff : noff + n - 1 ) = 0.0d0

      ! ... initialize the wave functions in such a way that the values
      ! ... of the components are independent on the number of processors

      ampre = 0.01d0
      ALLOCATE( pwt( ngwt ) )

      ntest = ngwt / 4
      IF( ntest < SIZE( cm, 2 ) ) THEN
         ntest = ngwt
      END IF
      !
      ! ... assign random values to wave functions
      !
      DO ib = noff, noff + n - 1
        pwt( : ) = 0.0d0
        DO ig = 3, ntest
          rranf1 = 0.5d0 - rranf()
          rranf2 = rranf()
          pwt( ig ) = ampre * CMPLX(rranf1, rranf2)
        END DO
        CALL splitwf ( cm( :, ib ), pwt, ngw, ig_l2g, me_image, nproc_image, root_image, intra_image_comm )
      END DO
      IF ( gzero ) THEN
        cm( 1, noff : noff + n - 1 ) = (0.0d0, 0.0d0)
      END IF

      DEALLOCATE( pwt )

      RETURN
    END SUBROUTINE wave_rand_init_x



!=----------------------------------------------------------------------------=!
   SUBROUTINE kohn_sham_x( c, ngw, eforces, n, nl, noff )
!=----------------------------------------------------------------------------=!
        !
        ! ...   declare modules

        USE kinds,          ONLY: DP
        USE cp_interfaces,  ONLY: update_lambda, crot

        IMPLICIT NONE

        ! ...   declare subroutine arguments
        INTEGER, INTENT(IN) :: ngw   ! number of plane waves
        INTEGER, INTENT(IN) :: n     ! number of ks states
        INTEGER, INTENT(IN) :: nl    ! local (to the processor) number of states
        INTEGER, INTENT(IN) :: noff  ! offset of the first state in array "c"
        COMPLEX(DP), INTENT(INOUT) ::  c(:,:)
        COMPLEX(DP) :: eforces(:,:)

        ! ...   declare other variables
        INTEGER ::  ib
        REAL(DP),  ALLOCATABLE :: gam(:,:)
        REAL(DP),  ALLOCATABLE :: eig(:)
        LOGICAL ::  tdist

        ! ...   end of declarations

        IF( n < 1 ) THEN

           eforces = 0.0d0

        ELSE
 
           tdist  =  ( nl /= n )

           ALLOCATE( eig( n ) )
           ALLOCATE( gam( nl, n ) )

           DO ib = 1, n
              CALL update_lambda( ib, gam, c, eforces(:,ib+noff-1), n, noff, tdist )
           END DO
           CALL crot( c, ngw, n, noff, gam, nl, eig )

           DEALLOCATE( gam, eig )

        END IF

        RETURN
        ! ...
   END SUBROUTINE kohn_sham_x

