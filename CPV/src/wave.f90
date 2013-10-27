!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!



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
       USE mp_global,          ONLY: me_bgrp, intra_bgrp_comm
       USE mp,                 ONLY: mp_sum
       USE wave_base,          ONLY: hpsi
       USE gvect, ONLY: gstart
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, noff
       REAL(DP)            :: lambda(:,:)
       COMPLEX(DP)         :: c0(:,:), c2(:)
       INTEGER, INTENT(IN) :: i
       LOGICAL, INTENT(IN) :: tdist   !  if .true. lambda is distributed
       !
       REAL(DP), ALLOCATABLE :: prod(:)
       INTEGER :: ibl
       LOGICAL :: gzero
       !
       gzero = (gstart == 2) 
       ALLOCATE( prod( n ) )
       prod = hpsi( gzero, c0, SIZE( c0, 1 ), c2, n, noff )
       CALL mp_sum( prod, intra_bgrp_comm )
       IF( tdist ) THEN
          IF( me_bgrp == ib_owner( i ) ) THEN
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
    use mp_global,          only : intra_bgrp_comm, nbgrp, inter_bgrp_comm
    use gvect, only : gstart
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

    CALL mp_sum( ekincm, intra_bgrp_comm )
    IF( nbgrp > 1 ) &
       CALL mp_sum( ekincm, inter_bgrp_comm )
    DEALLOCATE( emainv )

    return
  end subroutine elec_fakekine_x

!=----------------------------------------------------------------------------=!
  subroutine bandsum( bsum, c0, ngw, tbgrp )
!=----------------------------------------------------------------------------=!
    !
    !  This subroutine computes the CP(fake) wave functions kinetic energy
    
    USE kinds,              only : DP
    use mp,                 only : mp_sum
    use mp_global,          only : intra_bgrp_comm, nbgrp, inter_bgrp_comm
    USE electrons_base,     ONLY : nbsp, nbsp_bgrp
    !
    IMPLICIT NONE
    !
    integer, intent(in)      :: ngw    !  number of plane wave coeff.
    real(DP), intent(out)    :: bsum
    complex(DP), intent(in)  :: c0( ngw, * )
    logical, intent(in)      :: tbgrp
    !
    integer  :: i, n

    n = nbsp
    IF( tbgrp ) n = nbsp_bgrp

    bsum=0.0d0
    do i = 1, n
      bsum = bsum + SUM( DBLE( CONJG( c0( :, i ) ) * c0( :, i ) ) )
    end do
    CALL mp_sum( bsum, intra_bgrp_comm )
    IF( tbgrp ) &
       CALL mp_sum( bsum, inter_bgrp_comm )

    return
  end subroutine bandsum



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
      !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_bgrp+1 to PE 1 and
      !  so on).
      !  nrl = local number of rows
      !  ----------------------------------------------

      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE mp_global,        ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm
      USE dspev_module,     ONLY: pdspev_drv, dspev_drv

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

      DO i = 1, nss
         c0rot( :, i+noff-1 ) = 0.0d0
         becrot(:,i+noff-1 ) = 0.0d0
      END DO



!      becrot = 0.0d0
!      c0rot  = 0.0d0

         DO ip = 1, np_rot

            nrl_ip = nss / np_rot
            IF( (ip-1) < mod( nss, np_rot ) ) THEN
              nrl_ip = nrl_ip + 1
            END IF
 
            ALLOCATE( uu( nrl_ip, nss ) )
            IF( me_rot .EQ. (ip-1) ) THEN
              uu = lambda( 1:nrl_ip, 1:nss )
            END IF
            CALL mp_bcast( uu, (ip-1), intra_bgrp_comm)
 
            j      = ip
            DO jl = 1, nrl_ip
              DO i = 1, nss
                CALL daxpy(2*ngwl,uu(jl,i),c0(1,j+noff-1),1,c0rot(1,i+noff-1),1)
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
      USE mp_global,        ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm
      USE dspev_module,     ONLY: dspev_drv

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
            CALL daxpy( 2*ngw, vv(j,i), c0(1,j+noff-1), 1, c0rot(1,i+noffr-1), 1 )
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
        !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_bgrp+1 to PE 1 and so on).
        !  ----------------------------------------------
         
! ...   declare modules
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm
        USE wave_base,          ONLY: dotp
        USE gvect, ONLY: gstart

        IMPLICIT NONE

! ...   declare subroutine arguments
        INTEGER,     INTENT( IN )  :: ngw, n, noff
        COMPLEX(DP), INTENT(INOUT) :: a(:,:), b(:,:)
        REAL(DP),    OPTIONAL      :: lambda(:,:)

! ...   declare other variables
        REAL(DP), ALLOCATABLE :: ee(:)
        INTEGER :: i, j, jl
        COMPLEX(DP) :: alp
        LOGICAL :: gzero

! ... end of declarations
!  ----------------------------------------------

        IF( n < 1 ) THEN
          RETURN
        END IF
        gzero = (gstart == 2) 
        ALLOCATE( ee( n ) )
        
        DO i = 1, n
          DO j = 1, n
            ee(j) = -dotp( gzero, ngw, b(:,j+noff-1), a(:,i+noff-1), intra_bgrp_comm )
          END DO
          IF( PRESENT(lambda) ) THEN
            IF( MOD( (i-1), nproc_bgrp ) == me_bgrp ) THEN
              DO j = 1, n
                lambda( (i-1) / nproc_bgrp + 1, j ) = ee(j)
              END DO
            END IF
          END IF
          DO j = 1, n
            alp = CMPLX(ee(j),0.0d0,kind=DP)
            CALL zaxpy( ngw, alp, b(1,j+noff-1), 1, a(1,i+noff-1), 1 )
          END DO
        END DO
        DEALLOCATE(ee)

        RETURN
   END SUBROUTINE proj_gamma




!=----------------------------------------------------------------------------=!
   SUBROUTINE wave_rand_init_x( cm_bgrp, global )
!=----------------------------------------------------------------------------=!

      !  this routine sets the initial wavefunctions at random

! ... declare modules
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_sum, mp_max, mp_min
      USE mp_wave,            ONLY: splitwf
      USE mp_global,          ONLY: me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm
      USE gvect,              ONLY: ig_l2g, gstart
      USE gvect,              ONLY: mill, gg
      USE gvecw,              ONLY: ngw, ngw_g
      USE io_global,          ONLY: stdout
      USE random_numbers,     ONLY: randy
      USE electrons_base,     ONLY: nbsp, ibgrp_g2l
      
      IMPLICIT NONE

      ! ... declare subroutine arguments 
      COMPLEX(DP), INTENT(OUT) :: cm_bgrp(:,:)
      LOGICAL, OPTIONAL, INTENT(IN) :: global

      ! ... declare other variables
      INTEGER :: ntest, ig, ib, ibgrp, lb, ub
      REAL(DP) ::  rranf1, rranf2, ampre, ggx, fac, r1, r2, r3
      COMPLEX(DP), ALLOCATABLE :: pwt( : )
      REAL(DP),    ALLOCATABLE :: RND( : , : )
      INTEGER :: iss, n1, n2, m1, m2, i
      LOGICAL :: local

      ! ... initialize the wave functions in such a way that the values
      ! ... of the components are independent on the number of processors
      ! ... with "local" option the initialization is independend from G sorting too!

      ! ... Check array dimensions

      IF( SIZE( cm_bgrp, 1 ) < ngw ) THEN 
        CALL errore(' wave_rand_init ', ' wrong dimensions ', 3)
      END IF

      local = .TRUE.
      IF( PRESENT( global ) ) THEN
         local = .NOT. global
      END IF

      ! ... Reset them to zero
 
      cm_bgrp = 0.0d0

      ampre   = 0.01d0

      IF( local ) THEN
         ggx = MAXVAL( gg( 1:ngw ) )
         CALL mp_max( ggx, intra_bgrp_comm )
         lb = MINVAL( mill )
         CALL mp_min( lb, intra_bgrp_comm )
         ub = MAXVAL( mill )
         CALL mp_max( ub, intra_bgrp_comm )
         ALLOCATE( RND( 3, lb:ub ) )
      ELSE 
         ALLOCATE( pwt( ngw_g ) )
         ntest = ngw_g / 4
         IF( ntest < SIZE( cm_bgrp, 2 ) ) THEN
            ntest = ngw_g
         END IF
      END IF
      !
      ! ... assign random values to wave functions
      !
      DO ib = 1, nbsp

        IF( local ) THEN
           rnd = 0.0d0
           DO ig = lb, ub
              rnd( 1, ig ) = 0.5d0 - randy()
              rnd( 2, ig ) = 0.5d0 - randy()
              rnd( 3, ig ) = 0.5d0 - randy()
           END DO
        ELSE
           pwt( : ) = 0.0d0
           DO ig = 3, ntest
             rranf1 = 0.5d0 - randy()
             rranf2 = randy()
             pwt( ig ) = ampre * CMPLX(rranf1, rranf2,kind=DP)
           END DO
        END IF
        !
        ibgrp = ibgrp_g2l( ib )
        !
        IF( ibgrp > 0 ) THEN
          DO ig = 1, ngw
            IF( local ) THEN
               IF( gg(ig) < ggx / 2.519d0 ) THEN  ! 2.519 = 4^(2/3), equivalent to keep only (ngw_g/4) values
                  rranf1 = rnd( 1, mill(1,ig) ) * rnd( 2, mill(2,ig) ) * rnd( 3, mill(3,ig) )
                  rranf2 = 0.0d0
                  cm_bgrp( ig, ibgrp ) =  ampre * CMPLX( rranf1, rranf2 ,kind=DP) / ( 1.0d0 + gg(ig) )
               END IF
            ELSE
               cm_bgrp( ig, ibgrp ) = pwt( ig_l2g( ig ) )
            END IF
          END DO
        END IF
        !
      END DO
      IF ( gstart == 2 ) THEN
        cm_bgrp( 1, : ) = (0.0d0, 0.0d0)
      END IF

      IF( ALLOCATED( pwt ) ) DEALLOCATE( pwt )
      IF( ALLOCATED( rnd ) ) DEALLOCATE( rnd )

#ifdef PIPPO_DEBUG
      write(1000+mpime,fmt='(8I5)') nbsp, nbsp_bgrp, nudx, nudx_bgrp, nbsp, nbsp_bgrp, nbspx, nbspx_bgrp
      DO iss = 1, nspin
         write(1000+mpime,fmt='(5I5)') nupdwn(iss), iupdwn(iss), nupdwn_bgrp(iss), iupdwn_bgrp(iss), i2gupdwn_bgrp(iss)
      END DO
      DO ib = 1, nbsp
         ! write(1000+mpime,fmt='(2I5)') ib, ibgrp_g2l(ib)
      END DO

      DO iss = nspin, 1, -1
         write(1000+mpime,*) 'copy'
         n1 = iupdwn_bgrp(iss)
         n2 = n1 + nupdwn_bgrp(iss) - 1
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         DO i = m2, m1, -1
            write(1000+mpime,fmt='(2I5)') i, i-m1+n1
         END DO
      END DO
      DO iss = 1, nspin
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         write(1000+mpime,*) 'zero'
         DO i = iupdwn(iss), m1-1
            write(1000+mpime,fmt='(1I5)') i
         END DO
         write(1000+mpime,*) 'zero'
         DO i = m2+1, iupdwn(iss) + nupdwn(iss) - 1
            write(1000+mpime,fmt='(1I5)') i
         END DO
      END DO
#endif


      RETURN
    END SUBROUTINE wave_rand_init_x


    SUBROUTINE c_bgrp_expand_x( c_bgrp )
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_sum
      USE electrons_base,     ONLY: nspin, i2gupdwn_bgrp, nupdwn, iupdwn_bgrp, iupdwn, nupdwn_bgrp
      USE mp_global,          ONLY: nbgrp, inter_bgrp_comm
      IMPLICIT NONE
      COMPLEX(DP) :: c_bgrp(:,:)
      INTEGER :: iss, n1, n2, m1, m2, i
      IF( nbgrp < 2 ) &
         RETURN
      DO iss = nspin, 1, -1
         n1 = iupdwn_bgrp(iss)
         n2 = n1 + nupdwn_bgrp(iss) - 1
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         DO i = m2, m1, -1
            c_bgrp(:,i) = c_bgrp(:,i-m1+n1)
         END DO
      END DO
      DO iss = 1, nspin
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         DO i = iupdwn(iss), m1-1
            c_bgrp(:,i) = 0.0d0
         END DO
         DO i = m2+1, iupdwn(iss) + nupdwn(iss) - 1
            c_bgrp(:,i) = 0.0d0
         END DO
      END DO
      CALL mp_sum( c_bgrp, inter_bgrp_comm )
      RETURN
    END SUBROUTINE c_bgrp_expand_x


    SUBROUTINE c_bgrp_pack_x( c_bgrp )
      USE kinds,              ONLY: DP
      USE electrons_base,     ONLY: nspin, i2gupdwn_bgrp, nupdwn, iupdwn_bgrp, iupdwn, nupdwn_bgrp
      USE mp_global,          ONLY: nbgrp
      IMPLICIT NONE
      COMPLEX(DP) :: c_bgrp(:,:)
      INTEGER :: iss, n1, n2, m1, m2, i
      IF( nbgrp < 2 ) &
         RETURN
      DO iss = 1, nspin
         n1 = iupdwn_bgrp(iss)
         n2 = n1 + nupdwn_bgrp(iss) - 1
         m1 = iupdwn(iss)+i2gupdwn_bgrp(iss) - 1
         m2 = m1 + nupdwn_bgrp(iss) - 1
         DO i = n1, n2
            c_bgrp(:,i) = c_bgrp(:,i-n1+m1)
         END DO
      END DO
      RETURN
    END SUBROUTINE c_bgrp_pack_x

