!
! Copyright (C) 2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
! FFT base Module.
! Written by Carlo Cavazzoni 
!----------------------------------------------------------------------
!

#if defined __HPM
#  include "/cineca/prod/hpm/include/f_hpm.h"
#endif


!=----------------------------------------------------------------------=!
      MODULE fft_base
!=----------------------------------------------------------------------=!

        USE kinds, ONLY: DP
        USE parallel_include

        USE fft_types, ONLY: fft_dlay_descriptor

        IMPLICIT NONE
        
        TYPE ( fft_dlay_descriptor ) :: dfftp ! descriptor for dense grid
        TYPE ( fft_dlay_descriptor ) :: dffts ! descriptor for smooth grid
        TYPE ( fft_dlay_descriptor ) :: dfftb ! descriptor for box grids

        SAVE

        PRIVATE

#if ! defined __PARA
        INTERFACE fft_itranspose
           MODULE procedure fft_transpose
        END INTERFACE
#endif

        PUBLIC :: fft_itranspose, fft_transpose, fft_scatter
        PUBLIC :: dfftp, dffts, dfftb, fft_dlay_descriptor


!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!

!

#if defined __PARA


        SUBROUTINE fft_transpose( zstick, ldz, r, ldx, ldy, dfft, me, group, nproc, iopt )

          IMPLICIT NONE

          COMPLEX (DP) :: zstick( * )
          COMPLEX (DP) :: r( * )
          TYPE (fft_dlay_descriptor), INTENT(IN) ::  dfft
          INTEGER, INTENT(IN) :: me    ! processor index starting from 1
          INTEGER, INTENT(IN) :: nproc
          INTEGER, INTENT(IN) :: group, iopt, ldz, ldx, ldy

          INTEGER :: i, j, k, ipz, offset, k_start, k_end, is
          INTEGER :: npz, nz_l, ns_l, ns_lp, nz_lx
          INTEGER :: nsx_l, msgsiz
          INTEGER :: mc1, mc2, mc3, mc4, is_offset, ns1, ierr
          COMPLEX (DP) :: bswp( ldz * 4 )
          COMPLEX (DP), ALLOCATABLE :: mp_snd_buffer(:)
          COMPLEX (DP), ALLOCATABLE :: mp_rcv_buffer(:)

          !
          ! ... SUBROUTINE BODY
          !

          npz  = nproc
          nz_lx = MAXVAL( dfft%npp( 1 : nproc ) )

          IF( ABS( iopt ) == 2 ) THEN
            ns_l  = dfft%nsw( me )
            nsx_l = MAXVAL( dfft%nsw( : ) )
          ELSE
            ns_l  = dfft%nsp( me )
            nsx_l = MAXVAL( dfft%nsp( : ) )
          END IF

          msgsiz = nsx_l * nz_lx 
          ALLOCATE( mp_snd_buffer( msgsiz * npz ) )
          ALLOCATE( mp_rcv_buffer( msgsiz * npz ) )

          IF ( iopt < 1 ) THEN

            r( 1 : dfft%nnr ) = 0.0_DP

            k_start = 1
            DO ipz = 1, npz
              nz_l = dfft%npp( ipz )
              k_end   = k_start  + nz_l - 1
              offset  = ( ipz - 1 ) * msgsiz - k_start + 1
              DO is = 1, ns_l
                is_offset = ( is - 1 ) * ldz
                DO k = k_start , k_end
                  mp_snd_buffer(k + offset) = zstick( k + is_offset )
                END DO
                offset = offset + nz_l
              END DO
              k_start = k_start + nz_l
            END DO


            call MPI_ALLTOALL(mp_snd_buffer(1),msgsiz,MPI_DOUBLE_COMPLEX, &
                        mp_rcv_buffer(1),msgsiz,MPI_DOUBLE_COMPLEX, &
                        group, ierr)

            DO ipz = 1, npz
              nz_l = dfft%npp( me )
              offset = ( ipz - 1 ) * msgsiz 
              is_offset = dfft%iss( ipz )
              IF( ABS( iopt ) == 1 ) THEN
                ns_lp  = dfft%nsp(ipz)
              ELSE
                ns_lp  = dfft%nsw(ipz)
              END IF
              ns1 = MOD( ns_lp, 4 )
              IF( ns1 /= 0 ) THEN
                DO is = 1, ns1
                  mc1 = dfft%ismap( is   + is_offset )
                  DO k = 1 , nz_l
                    r( mc1 + (k-1)*dfft%nnp ) = mp_rcv_buffer( k + offset )
                  END DO
                  offset = offset + nz_l
                END DO
              END IF
              IF( ns_lp >= 4 ) THEN
                ns1 = ns1 + 1
                DO is = ns1, ns_lp, 4
                  mc1 = dfft%ismap( is   + is_offset )
                  mc2 = dfft%ismap( is+1 + is_offset )
                  mc3 = dfft%ismap( is+2 + is_offset )
                  mc4 = dfft%ismap( is+3 + is_offset )
                  DO k = 1 , nz_l
                    bswp( k ) = mp_rcv_buffer( k + offset )
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    bswp( k + nz_l ) = mp_rcv_buffer( k + offset )
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    bswp( k + 2*nz_l ) = mp_rcv_buffer( k + offset )
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    bswp( k + 3*nz_l ) = mp_rcv_buffer( k + offset )
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    r( mc1 + (k-1)*dfft%nnp ) = bswp( k          )
                    r( mc2 + (k-1)*dfft%nnp ) = bswp( k +   nz_l )
                    r( mc3 + (k-1)*dfft%nnp ) = bswp( k + 2*nz_l )
                    r( mc4 + (k-1)*dfft%nnp ) = bswp( k + 3*nz_l )
                  END DO
                END DO
              END IF
            END DO

          ELSE IF ( iopt > 0 ) THEN

            DO ipz = 1, npz
              nz_l = dfft%npp( me )
              offset = ( ipz - 1 ) * msgsiz 
              is_offset = dfft%iss( ipz )
              IF( ABS( iopt ) == 1 ) THEN
                ns_lp  = dfft%nsp(ipz)
              ELSE
                ns_lp  = dfft%nsw(ipz)
              END IF
              ns1 = MOD( ns_lp, 4 )
              IF( ns1 /= 0 ) THEN
                DO is = 1, ns1
                  mc1 = dfft%ismap( is   + is_offset )
                  DO k = 1 , nz_l
                    mp_snd_buffer( k + offset ) = r( mc1 + (k-1)*dfft%nnp )
                  END DO
                  offset = offset + nz_l
                END DO
              END IF
              IF( ns_lp >= 4 ) THEN
                ns1 = ns1 + 1
                DO is = ns1, ns_lp, 4
                  mc1 = dfft%ismap( is   + is_offset )
                  mc2 = dfft%ismap( is+1 + is_offset )
                  mc3 = dfft%ismap( is+2 + is_offset )
                  mc4 = dfft%ismap( is+3 + is_offset )
                  DO k = 1, nz_l
                    bswp( k          ) = r( mc1 + (k-1)*dfft%nnp )
                    bswp( k +   nz_l ) = r( mc2 + (k-1)*dfft%nnp )
                    bswp( k + 2*nz_l ) = r( mc3 + (k-1)*dfft%nnp )
                    bswp( k + 3*nz_l ) = r( mc4 + (k-1)*dfft%nnp )
                  END DO
                  DO k = 1 , nz_l
                    mp_snd_buffer( k + offset ) = bswp( k )
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    mp_snd_buffer( k + offset ) = bswp( k + nz_l )
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    mp_snd_buffer( k + offset ) = bswp( k + 2*nz_l )
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    mp_snd_buffer( k + offset ) = bswp( k + 3*nz_l )
                  END DO
                  offset = offset + nz_l
                END DO
              END IF
            END DO

            call MPI_ALLTOALL(mp_snd_buffer(1),msgsiz,MPI_DOUBLE_COMPLEX, &
                        mp_rcv_buffer(1),msgsiz,MPI_DOUBLE_COMPLEX, &
                        group, ierr)

            k_start = 1
            DO IPZ = 1, NPZ
              nz_l = dfft%npp( ipz )
              k_end   = k_start + nz_l - 1
              offset  = (ipz-1) * msgsiz - k_start + 1
              DO is = 1, ns_l
                is_offset = ( is - 1 ) * ldz
                DO k = k_start , k_end
                  zstick( k + is_offset ) = mp_rcv_buffer( k + offset )
                END DO
                offset = offset + nz_l
              END DO
              k_start = k_start + nz_l
            END DO

          END IF

          DEALLOCATE( mp_snd_buffer )
          DEALLOCATE( mp_rcv_buffer )

          RETURN 
        END SUBROUTINE fft_transpose


        !=----------------------------------------------------------------------=!


        SUBROUTINE fft_itranspose( zstick, ldz, r, ldx, ldy, dfft, me, group, nproc, iopt)

          !
          ! Non-blocking version of the transpose routine 
          ! Transpose re-distribute fft data between 1D fft along Z dimension and
          ! 2D fft in the XY plane
          !

          IMPLICIT NONE

          COMPLEX (DP) :: zstick( * )
          COMPLEX (DP) :: r( * )
          TYPE (fft_dlay_descriptor), INTENT(IN) ::  dfft
          INTEGER, INTENT(IN) :: me    ! processor index starting from 1
          INTEGER, INTENT(IN) :: nproc
          INTEGER, INTENT(IN) :: group, iopt
          INTEGER, INTENT(IN) :: ldz, ldx, ldy

          INTEGER :: i, j, k, ipz, offset, k_start, k_end, is, is_start, is_end
          INTEGER :: npz, nz_l, ns_l, ns_lp, nbuf, ierr, itag, mype, nsx_l, nz_lx
          INTEGER, ALLOCATABLE :: ishand( : )
          INTEGER, ALLOCATABLE :: irhand( : )
          INTEGER, ALLOCATABLE :: istatus( :, : )
          LOGICAL, ALLOCATABLE :: rtest( : )
          LOGICAL, ALLOCATABLE :: rdone( : )
          COMPLEX(DP), ALLOCATABLE :: sndbuf(:,:)
          COMPLEX(DP), ALLOCATABLE :: rcvbuf(:,:)
          INTEGER :: i1, i2, j1, j2, mc1, mc2, is_offset

          !
          ! ... SUBROUTINE BODY
          !

          npz  = nproc
          mype = me - 1
          nz_l = dfft%npp( me )
          nz_lx = MAXVAL( dfft%npp( 1 : nproc ) )
          IF( ABS( iopt ) == 2 ) THEN
            ns_l  = dfft%nsw( me )
            nsx_l = MAXVAL( dfft%nsw( : ) )
          ELSE
            ns_l  = dfft%nsp( me )
            nsx_l = MAXVAL( dfft%nsp( : ) )
          END IF

          ALLOCATE( sndbuf( nsx_l * nz_lx, npz ), ishand( npz ) )
          ALLOCATE( rcvbuf( nsx_l * nz_lx, npz ), irhand( npz ) )
          ALLOCATE( rtest( npz ), rdone( npz ) )
          ALLOCATE( istatus( MPI_STATUS_SIZE, npz ) )
          nbuf = nsx_l * nz_lx

          IF ( iopt < 0 ) THEN

            DO ipz = 1, npz
              itag = mype + 1 + npz * ( ipz - 1 )
              call mpi_irecv(rcvbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                group, irhand(ipz), ierr )
            END DO

            k_start = 1
            DO ipz = 1, npz
              nz_l = dfft%npp( ipz )
              k_end   = k_start  + nz_l - 1
              offset  = - k_start + 1
              DO is = 1, ns_l
                DO k = k_start , k_end
                  sndbuf(k + offset, ipz) = zstick( k + (is-1)*ldz )
                END DO
                offset = offset + nz_l
              END DO
              itag = ipz + npz * mype
              CALL mpi_isend( sndbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                group, ishand(ipz), ierr )
              k_start = k_start + nz_l
            END DO

            r( 1 : dfft%nnr ) = 0.0_DP
            rdone = .FALSE.

   111      CONTINUE
            DO IPZ = 1, NPZ
              nz_l = dfft%npp( me )
              call mpi_test(irhand(ipz), rtest(ipz), istatus(1,ipz), ierr)
              IF( rtest(ipz) .AND. .NOT. rdone(ipz) ) THEN
                offset = 0
                is_offset = dfft%iss( ipz )
                IF( ABS( iopt ) == 2 ) THEN
                  ns_lp  = dfft%nsw( ipz )  
                ELSE
                  ns_lp  = dfft%nsp( ipz )  
                END IF
                DO is = 1, ns_lp - 1, 2
                  mc1 = dfft%ismap( is   + is_offset )
                  mc2 = dfft%ismap( is+1 + is_offset )
                  DO k = 1 , nz_l
                    r( mc1 + (k-1)*ldx*ldy ) = rcvbuf(k + offset, ipz)
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    r( mc2 + (k-1)*ldx*ldy ) = rcvbuf(k + offset, ipz)
                  END DO
                  offset = offset + nz_l
                END DO
                IF( MOD( ns_lp, 2 ) /= 0 ) THEN
                  is = ns_lp
                  mc1 = dfft%ismap( is   + is_offset )
                  DO k = 1 , nz_l
                    r( mc1 + (k-1)*ldx*ldy ) = rcvbuf(k + offset, ipz)
                  END DO
                  offset = offset + nz_l
                END IF
                rdone( ipz ) = .TRUE.
              END IF
            END DO
            IF( .NOT. ALL( rtest ) ) GO TO 111


          ELSE IF ( iopt > 0 ) THEN

            DO IPZ = 1, NPZ
              itag = mype + 1 + npz * ( ipz - 1 )
              call mpi_irecv(rcvbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                group, irhand(ipz), ierr )
            END DO

            DO ipz = 1, npz
              nz_l = dfft%npp( me )
              offset = 0
              is_offset = dfft%iss( ipz )
              IF( ABS( iopt ) == 2 ) THEN
                ns_lp  = dfft%nsw( ipz )  
              ELSE
                ns_lp  = dfft%nsp( ipz )  
              END IF
              DO is = 1, ns_lp
                mc1 = dfft%ismap( is   + is_offset )
                DO k = 1 , nz_l
                  sndbuf( k + offset, ipz ) = r( mc1 + (k-1)*ldx*ldy )
                END DO
                offset = offset + nz_l
              END DO
              itag = ipz + npz * mype
              call mpi_isend(sndbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                group, ishand(ipz), ierr )
            END DO

            rdone = .FALSE.

   112      CONTINUE
            k_start = 1
            DO IPZ = 1, NPZ
              nz_l = dfft%npp( ipz )
              call mpi_test(irhand(ipz), rtest(ipz), istatus(1,ipz), ierr)
              IF( rtest(ipz) .AND. .NOT. rdone(ipz) ) THEN
                k_end   = k_start + nz_l - 1
                offset  = - k_start + 1
                DO is = 1, ns_l
                  DO k = k_start , k_end
                    zstick( k + (is-1)*ldz ) = rcvbuf( k + offset, ipz )
                  END DO
                  offset = offset + nz_l
                END DO
                rdone( ipz ) = .TRUE.
              END IF
              k_start = k_start + nz_l
            END DO
            IF( .NOT. ALL( rtest ) ) GO TO 112


          END IF

          DO ipz = 1, npz
            call mpi_wait(ishand(ipz), istatus(1,ipz), ierr)
          END DO

          DEALLOCATE(sndbuf, rcvbuf, ishand, irhand, rtest, rdone, istatus)

          RETURN 
        END SUBROUTINE fft_itranspose


#else 

!     Scalar code

        SUBROUTINE fft_transpose( zstick, ldz, r, ldx, ldy, dfft, me, group, nproc, iopt)

          IMPLICIT NONE

          COMPLEX (DP) :: zstick( * )
          COMPLEX (DP) :: r( * )
          TYPE (fft_dlay_descriptor), INTENT(IN) ::  dfft
          INTEGER, INTENT(IN) :: me    ! processor index starting from 1
          INTEGER, INTENT(IN) :: group, nproc
          INTEGER, INTENT(IN) :: iopt, ldx, ldy, ldz

          INTEGER :: k, is, nz, ns, mc1

          nz = dfft%npp( me )

          IF ( iopt < 0 ) THEN

            r( 1 : dfft%nnr ) = CMPLX(0.0_DP,0.0_DP)

            IF( iopt == 2 ) THEN
              ns = dfft%nsw( me )
            ELSE
              ns = dfft%nsp( me )
            END IF

            DO is = 1, ns
              mc1 = dfft%ismap( is )
              DO k = 1 , nz
                r( mc1 + (k-1)*ldx*ldy ) = zstick( k + (is-1)*ldz )
              END DO
            END DO

          ELSE IF ( iopt > 0 ) THEN

            IF( iopt == -2 ) THEN
              ns = dfft%nsw( me )
            ELSE
              ns = dfft%nsp( me )
            END IF

            DO is = 1, ns
              mc1 = dfft%ismap( is )
              DO k = 1 , nz
                zstick( k + (is-1)*ldz ) = r( mc1 + (k-1)*ldx*ldy )
              END DO
            END DO

          END IF

          RETURN 
        END SUBROUTINE fft_transpose


#endif


!
!

!
!-----------------------------------------------------------------------
subroutine fft_scatter ( f_in, nrx3, nxx_, f_aux, ncp_, npp_, sign, use_tg )
  !-----------------------------------------------------------------------
  !
  ! transpose the fft grid across nodes
  ! a) From columns to planes (sign > 0)
  !
  !    "columns" (or "pencil") representation:
  !    processor "me" has ncp_(me) contiguous columns along z
  !    Each column has nrx3 elements for a fft of order nr3
  !    nrx3 can be =nr3+1 in order to reduce memory conflicts.
  !
  !    The transpose take places in two steps:
  !    1) on each processor the columns are divided into slices along z
  !       that are stored contiguously. On processor "me", slices for
  !       processor "proc" are npp_(proc)*ncp_(me) big
  !    2) all processors communicate to exchange slices
  !       (all columns with z in the slice belonging to "me"
  !        must be received, all the others must be sent to "proc")
  !    Finally one gets the "planes" representation:
  !    processor "me" has npp_(me) complete xy planes
  !
  !  b) From planes to columns (sign < 0)
  !
  !  Quite the same in the opposite direction
  !
  !  The output is overwritten on f_in ; f_aux is used as work space
  !
  !  If optional argument "use_tg" is true the subroutines performs
  !  the trasposition using the Task Groups distribution
  !
#ifdef __PARA
  USE parallel_include
#endif
  use mp_global,   ONLY : nproc_pool, me_pool, intra_pool_comm, nproc, &
                          my_image_id, nogrp, pgrp_comm
  USE kinds,       ONLY : DP
  USE task_groups, ONLY : nplist

  implicit none

  integer, intent(in)           :: nrx3, nxx_, sign, ncp_ (:), npp_ (:)
  complex (DP)                  :: f_in (nxx_), f_aux (nxx_)
  logical, optional, intent(in) :: use_tg

#ifdef __PARA

  integer :: dest, from, k, offset1 (nproc), sendcount (nproc), &
       sdispls (nproc), recvcount (nproc), rdispls (nproc), &
       proc, ierr, me, nprocp, gproc, gcomm
  !
  LOGICAL :: use_tg_

#if defined __HPM
     !       CALL f_hpmstart( 10, 'scatter' )
#endif

  !
  !  Task Groups

  use_tg_ = .FALSE.

  IF( PRESENT( use_tg ) ) use_tg_ = use_tg
     
  me     = me_pool + 1
  !
  IF( use_tg_ ) THEN
    !  This is the number of procs. in the plane-wave group
     nprocp = nproc_pool / nogrp 
  ELSE
     nprocp = nproc_pool
  END IF
  !
  if (nprocp.eq.1) return
  !
  call start_clock ('fft_scatter')
  !
  ! sendcount(proc): amount of data processor "me" must send to processor
  ! recvcount(proc): amount of data processor "me" must receive from
  !
  IF( use_tg_ ) THEN
     do proc = 1, nprocp
        gproc = nplist( proc ) + 1
        sendcount (proc) = npp_ ( gproc ) * ncp_ (me)
        recvcount (proc) = npp_ (me) * ncp_ ( gproc )
     end do 
  ELSE
     do proc = 1, nprocp
        sendcount (proc) = npp_ (proc) * ncp_ (me)
        recvcount (proc) = npp_ (me) * ncp_ (proc)
     end do
  END IF
  !
  ! offset1(proc) is used to locate the slices to be sent to proc
  ! sdispls(proc)+1 is the beginning of data that must be sent to proc
  ! rdispls(proc)+1 is the beginning of data that must be received from pr
  !
  offset1 (1) = 1
  IF( use_tg_ ) THEN
     do proc = 2, nprocp
        gproc = nplist( proc - 1 ) + 1
        offset1 (proc) = offset1 (proc - 1) + npp_ ( gproc )
     enddo
  ELSE
     do proc = 2, nprocp
        offset1 (proc) = offset1 (proc - 1) + npp_ (proc - 1)
     enddo
  END IF

  sdispls (1) = 0
  rdispls (1) = 0
  do proc = 2, nprocp
     sdispls (proc) = sdispls (proc - 1) + sendcount (proc - 1)
     rdispls (proc) = rdispls (proc - 1) + recvcount (proc - 1)
  enddo
  !

  ierr = 0
  if (sign.gt.0) then
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     do proc = 1, nprocp
        from = offset1 (proc)
        dest = 1 + sdispls (proc)
        IF( use_tg_ ) THEN
           gproc = nplist(proc)+1
        ELSE
           gproc = proc
        END IF
        do k = 1, ncp_ (me)
           call DCOPY (2 * npp_ ( gproc ), f_in (from + (k - 1) * nrx3), &
                1, f_aux (dest + (k - 1) * npp_ ( gproc ) ), 1)
        enddo
     enddo
     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     f_in = 0.0_DP
     !
     ! step two: communication
     !
     IF( use_tg_ ) THEN
        gcomm = pgrp_comm
     ELSE
        gcomm = intra_pool_comm
     END IF

     !  call mpi_barrier (gcomm, ierr)  ! why barrier?

     call mpi_alltoallv (f_aux(1), sendcount, sdispls, MPI_DOUBLE_COMPLEX, f_in(1), &
          recvcount, rdispls, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     if( ABS(ierr) /= 0 ) call errore ('fft_scatter', 'info<>0', ABS(ierr) )
     !
  else
     !
     !  "backward" scatter from planes to columns
     !
     !  step two: communication
     !
     IF( use_tg_ ) THEN
        gcomm = pgrp_comm
     ELSE
        gcomm = intra_pool_comm
     END IF

     ! call mpi_barrier (gcomm, ierr)  ! why barrier

     call mpi_alltoallv (f_in(1), recvcount, rdispls, MPI_DOUBLE_COMPLEX, f_aux(1), &
          sendcount, sdispls, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     if( ABS(ierr) /= 0 ) call errore ('fft_scatter', 'info<>0', ABS(ierr) )
     !
     !  step one: store contiguously the columns
     !
     f_in = 0.0_DP
     !
     do proc = 1, nprocp
        from = 1 + sdispls (proc)
        dest = offset1 (proc)
        IF( use_tg_ ) THEN
           gproc = nplist(proc)+1
        ELSE
           gproc = proc
        END IF
        do k = 1, ncp_ (me)
           call DCOPY ( 2 * npp_ ( gproc ), f_aux (from + (k - 1) * npp_ ( gproc ) ), 1, &
                                            f_in  (dest + (k - 1) * nrx3 ), 1 )
        enddo

     enddo

  endif

  call stop_clock ('fft_scatter')

#endif

#if defined __HPM
     !       CALL f_hpmstop( 10 )
#endif

  return

end subroutine fft_scatter

!=----------------------------------------------------------------------=!
      END MODULE fft_base
!=----------------------------------------------------------------------=!
