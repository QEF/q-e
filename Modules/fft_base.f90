!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!----------------------------------------------------------------------
! FFT base Module.
! Written by Carlo Cavazzoni 
!----------------------------------------------------------------------
!

#if defined __HPM
#  include "/cineca/prod/hpm_2_4_2/include/f_hpm.h"
#endif

#define __FFT_BASE_TS1 

!=----------------------------------------------------------------------=!
      MODULE fft_base
!=----------------------------------------------------------------------=!

        USE kinds, ONLY: dbl
        USE mp, ONLY: mp_max, mp_sum, mp_barrier

        USE fft_types, ONLY: fft_dlay_descriptor

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: fft_transpose, fft_scatter


        INTEGER, ALLOCATABLE :: stmask(:,:,:)


!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!

!=----------------------------------------------------------------------=!
!  ...  FFT inizialization
!
        SUBROUTINE transpose_setup( dfft, me, nproc )

          IMPLICIT NONE
          
          TYPE (fft_dlay_descriptor), INTENT(IN) :: dfft
          INTEGER, INTENT(IN) :: me    ! processor index starting from 1
          INTEGER, INTENT(IN) :: nproc ! number of processor

          INTEGER :: i, j, is, ip, mc
          INTEGER :: nspx
          INTEGER :: ierr
!
! ...     Subroutine Body
!

! ...     nx, ny, nz are te sizes of the 3D fft data grid

          nspx = MAXVAL( dfft%nsp )

          ierr = 0
          IF( ALLOCATED(stmask) ) DEALLOCATE(stmask, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' fft_base_setup ' , ' deallocation of stmask failed ', ierr)

          ALLOCATE( stmask ( 2, nspx, nproc ), STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' fft_base_setup ' , ' allocation of stmask failed ', ierr)

! ...     the stick mask is copied locally to increase the efficiency

          DO ip = 1, nproc
            DO is = 1, dfft%nsp( ip )
              mc = dfft%ismap( is + dfft%iss( ip ) )
              j  = ( mc - 1 ) / dfft%nr1x + 1
              i  = MOD( ( mc - 1 ), dfft%nr1x ) + 1
              stmask( 1, is, ip ) = i
              stmask( 2, is, ip ) = j
            END DO
          END DO

          RETURN
        END SUBROUTINE transpose_setup

!
!
!=======================================================================
!

#if defined __PARA

#  if defined __FFT_BASE_TS1

        SUBROUTINE fft_transpose( zstick, r, dfft, me, nproc, iopt)

          USE mp_buffers, ONLY: mp_allocate_buffers, &
            mp_snd_buffer, mp_rcv_buffer, mp_sendrecv_buffers, &
            mp_barrier_buffers, mp_deallocate_buffers, mp_alltoall_buffers, &
            mp_p_snd_buffer, mp_p_rcv_buffer, mp_bufsize_msgmax

          IMPLICIT NONE

          include 'mpif.h'

          COMPLEX (dbl) :: zstick(:,:)
          COMPLEX (dbl) :: r(:,:,:)
          TYPE (fft_dlay_descriptor), INTENT(IN) ::  dfft
          INTEGER, INTENT(IN) :: me    ! processor index starting from 1
          INTEGER, INTENT(IN) :: nproc
          INTEGER, INTENT(IN) :: iopt


          INTEGER :: i, j, k, ipz, offset, k_start, k_end, is
          INTEGER :: npz, nz_l, ns_l, ns_lp
          INTEGER :: nsx_l, msgsiz
          INTEGER :: i1, i2, j1, j2
          COMPLEX (dbl) :: zero

          INTEGER, SAVE :: dfft_id = -1

!
! ... SUBROUTINE BODY
!
          IF( iopt == 0 ) THEN
            CALL transpose_setup( dfft, me, nproc )
            dfft_id = dfft%id
            RETURN
          END IF

          IF( dfft_id /= dfft%id ) THEN
            CALL transpose_setup( dfft, me, nproc )
            dfft_id = dfft%id
          END IF

          npz  = nproc
          nz_l = dfft%npp( me )
          IF( ABS( iopt ) == 2 ) THEN
            ns_l  = dfft%nsw( me )
            nsx_l = MAXVAL( dfft%nsw( : ) )
          ELSE
            ns_l  = dfft%nsp( me )
            nsx_l = MAXVAL( dfft%nsp( : ) )
          END IF

          zero = 0.0d0

          msgsiz = nsx_l * nz_l 
          CALL mp_allocate_buffers( msgsiz * npz )

          IF ( iopt < 1 ) THEN

            r = 0.0d0

#if defined __HPM
            CALL f_hpmstart( 1, 'pack_1' )
#endif

              DO ipz = 1, npz
                k_start = (ipz-1)  * nz_l + 1
                k_end   = k_start  + nz_l - 1
                offset  = (ipz-1) * msgsiz - k_start + 1
                DO is = 1, ns_l
                  DO k = k_start , k_end
                    mp_snd_buffer(k + offset) = zstick(k, is)
                  END DO
                  offset = offset + nz_l
                END DO
              END DO

#if defined __HPM
            CALL f_hpmstop( 1 )
#endif

              CALL mp_alltoall_buffers(mp_snd_buffer, mp_rcv_buffer)

#if defined __HPM
            CALL f_hpmstart( 2, 'pack_1' )
#endif

              DO ipz = 1, npz
                offset = (ipz-1) * msgsiz 
                IF( ABS( iopt ) == 1 ) THEN
                  ns_lp  = dfft%nsp(ipz)
                ELSE
                  ns_lp  = dfft%nsw(ipz)
                END IF
                DO is = 1, ns_lp
                  i = stmask( 1, is, ipz ) 
                  j = stmask( 2, is, ipz ) 
                  DO k = 1 , nz_l
                    r( i, j, k ) = mp_rcv_buffer( k + offset )
                  END DO
                  offset = offset + nz_l
                END DO
              END DO


#if defined __HPM
            CALL f_hpmstop( 2 )
#endif

          ELSE IF ( iopt > 0 ) THEN

              DO ipz = 1, npz
                offset = (ipz-1) * msgsiz 
                IF( ABS( iopt ) == 1 ) THEN
                  ns_lp  = dfft%nsp(ipz)
                ELSE
                  ns_lp  = dfft%nsw(ipz)
                END IF
                DO is = 1, ns_lp
                  i = stmask(1,is,ipz) 
                  j = stmask(2,is,ipz) 
                  DO k = 1 , nz_l
                    mp_snd_buffer( k + offset ) = r(i,j,k)
                  END DO
                  offset = offset + nz_l
                END DO
              END DO

              call mp_alltoall_buffers(mp_snd_buffer, mp_rcv_buffer)

              DO IPZ = 1, NPZ
                k_start = (ipz-1) * nz_l + 1
                k_end   = k_start + nz_l - 1
                offset  = (ipz-1) * msgsiz - k_start + 1
                DO is = 1, ns_l
                  DO k = k_start , k_end
                    zstick(k,is) = mp_rcv_buffer( k + offset )
                  END DO
                  offset = offset + nz_l
                END DO
              END DO

          END IF

          CALL mp_deallocate_buffers()


          RETURN 
        END SUBROUTINE fft_transpose

#  else


        SUBROUTINE fft_transpose( zstick, r, dfft, me, nproc, iopt)

          IMPLICIT NONE

          include 'mpif.h'

          COMPLEX (dbl) :: zstick(:,:)
          COMPLEX (dbl) :: r(:,:,:)
          TYPE (fft_dlay_descriptor), INTENT(IN) ::  dfft
          INTEGER, INTENT(IN) :: me    ! processor index starting from 1
          INTEGER, INTENT(IN) :: nproc
          INTEGER, INTENT(IN) :: iopt

          INTEGER :: i, j, k, ipz, offset, k_start, k_end, is, is_start, is_end
          INTEGER :: npz, nz_l, ns_l, ns_lp, nbuf, ierr, itag, mype, nsx_l
          INTEGER, ALLOCATABLE :: ishand( : )
          INTEGER, ALLOCATABLE :: irhand( : )
          INTEGER, ALLOCATABLE :: istatus( :, : )
          LOGICAL, ALLOCATABLE :: rtest( : )
          LOGICAL, ALLOCATABLE :: rdone( : )
          COMPLEX(dbl), ALLOCATABLE :: sndbuf(:,:)
          COMPLEX(dbl), ALLOCATABLE :: rcvbuf(:,:)
          INTEGER :: i1, i2, j1, j2

          INTEGER, SAVE :: dfft_id = -1

!
! ... SUBROUTINE BODY
!

          IF( iopt == 0 ) THEN
            CALL transpose_setup( dfft, me, nproc )
            dfft_id = dfft%id
            RETURN
          END IF

          IF( dfft_id /= dfft%id ) THEN
            CALL transpose_setup( dfft, me, nproc )
            dfft_id = dfft%id
          END IF


          npz  = nproc
          mype = me - 1
          nz_l = dfft%npp( me )
          IF( ABS( iopt ) == 2 ) THEN
            ns_l  = dfft%nsw( me )
            nsx_l = MAXVAL( dfft%nsw( : ) )
          ELSE
            ns_l  = dfft%nsp( me )
            nsx_l = MAXVAL( dfft%nsp( : ) )
          END IF

          WRITE(6, fmt = "( 'DEBUG fft_transpose: ', 5I5 ) " ) me, nproc, iopt, nz_l, ns_l
          

          ALLOCATE( sndbuf( nsx_l * nz_l, npz ), ishand( npz ) )
          ALLOCATE( rcvbuf( nsx_l * nz_l, npz ), irhand( npz ) )
          ALLOCATE( rtest( npz ), rdone( npz ) )
          ALLOCATE( istatus( MPI_STATUS_SIZE, npz ) )
          nbuf = nsx_l * nz_l

          IF ( iopt < 0 ) THEN

            DO ipz = 1, npz
              itag = mype + 1 + npz * ( ipz - 1 )
              call mpi_irecv(rcvbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                MPI_COMM_WORLD, irhand(ipz), ierr )
            END DO

            DO ipz = 1, npz
              k_start = ( ipz - 1 )  * nz_l + 1
              k_end   = k_start  + nz_l - 1
              offset  = - k_start + 1
              DO is = 1, ns_l
                DO k = k_start , k_end
                  sndbuf(k + offset, ipz) = zstick(k, is)
                END DO
                offset = offset + nz_l
              END DO
              itag = ipz + npz * mype
              CALL mpi_isend( sndbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                MPI_COMM_WORLD, ishand(ipz), ierr )
            END DO

            r = 0.0d0
            rdone = .FALSE.

#if defined __HPM
            CALL f_hpmstart( 1, 'unpack_1' )
#endif

   111      CONTINUE
            DO IPZ = 1, NPZ
              call mpi_test(irhand(ipz), rtest(ipz), istatus(1,ipz), ierr)
              IF( rtest(ipz) .AND. .NOT. rdone(ipz) ) THEN
                offset = 0
                IF( ABS( iopt ) == 2 ) THEN
                  ns_lp  = dfft%nsw( ipz )  
                ELSE
                  ns_lp  = dfft%nsp( ipz )  
                END IF
                DO is = 1, ns_lp - 1, 2
                  i1 = stmask(1,is,ipz) ! descz%mask(1,is,ipz)
                  j1 = stmask(2,is,ipz) ! descz%mask(2,is,ipz)
                  i2 = stmask(1,is+1,ipz) ! descz%mask(1,is,ipz)
                  j2 = stmask(2,is+1,ipz) ! descz%mask(2,is,ipz)
                  DO k = 1 , nz_l
                    r(i1,j1,k) = rcvbuf(k + offset, ipz)
                  END DO
                  offset = offset + nz_l
                  DO k = 1 , nz_l
                    r(i2,j2,k) = rcvbuf(k + offset, ipz)
                  END DO
                  offset = offset + nz_l
                END DO
                IF( MOD( ns_lp, 2 ) /= 0 ) THEN
                  is = ns_lp
                  i = stmask(1,is,ipz) ! descz%mask(1,is,ipz)
                  j = stmask(2,is,ipz) ! descz%mask(2,is,ipz)
                  DO k = 1 , nz_l
                    r(i,j,k) = rcvbuf(k + offset, ipz)
                  END DO
                  offset = offset + nz_l
                END IF
                rdone( ipz ) = .TRUE.
              END IF
            END DO
            IF( .NOT. ALL( rtest ) ) GO TO 111

#if defined __HPM
            CALL f_hpmstop( 1 )
#endif

          ELSE IF ( iopt > 0 ) THEN

            DO IPZ = 1, NPZ
              itag = mype + 1 + npz * ( ipz - 1 )
              call mpi_irecv(rcvbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                MPI_COMM_WORLD, irhand(ipz), ierr )
            END DO

            DO ipz = 1, npz
              offset = 0
              IF( ABS( iopt ) == 2 ) THEN
                ns_lp  = dfft%nsw( ipz )  
              ELSE
                ns_lp  = dfft%nsp( ipz )  
              END IF
              DO is = 1, ns_lp
                i = stmask(1,is,ipz) 
                j = stmask(2,is,ipz) 
                DO k = 1 , nz_l
                  sndbuf( k + offset, ipz ) = r(i,j,k)
                END DO
                offset = offset + nz_l
              END DO
              itag = ipz + npz * mype
              call mpi_isend(sndbuf(1,ipz), nbuf, MPI_DOUBLE_COMPLEX, ipz-1, itag, &
                MPI_COMM_WORLD, ishand(ipz), ierr )
            END DO

#if defined __HPM
            CALL f_hpmstart( 2, 'unpack_2' )
#endif

            rdone = .FALSE.

   112      CONTINUE
            DO IPZ = 1, NPZ
              call mpi_test(irhand(ipz), rtest(ipz), istatus(1,ipz), ierr)
              IF( rtest(ipz) .AND. .NOT. rdone(ipz) ) THEN
                k_start = (ipz-1) * nz_l + 1
                k_end   = k_start + nz_l - 1
                offset  = - k_start + 1
                DO is = 1, ns_l
                  DO k = k_start , k_end
                    zstick(k,is) = rcvbuf( k + offset, ipz )
                  END DO
                  offset = offset + nz_l
                END DO
                rdone( ipz ) = .TRUE.
              END IF
            END DO
            IF( .NOT. ALL( rtest ) ) GO TO 112

#if defined __HPM
            CALL f_hpmstop( 2 )
#endif

          END IF

          DO ipz = 1, npz
            call mpi_wait(ishand(ipz), istatus(1,ipz), ierr)
          END DO

          DEALLOCATE(sndbuf, rcvbuf, ishand, irhand, rtest, rdone, istatus)

          RETURN 
        END SUBROUTINE fft_transpose

#  endif

     
#else 

!     Scalar code

        SUBROUTINE fft_transpose( zstick, r, dfft, me, nproc, iopt)

          IMPLICIT NONE

          COMPLEX (dbl) :: zstick(:,:)
          COMPLEX (dbl) :: r(:,:,:)
          TYPE (fft_dlay_descriptor), INTENT(IN) ::  dfft
          INTEGER, INTENT(IN) :: me    ! processor index starting from 1
          INTEGER, INTENT(IN) :: nproc
          INTEGER, INTENT(IN) :: iopt

          INTEGER :: i, j, k, is, nz, ns
!
! ... SUBROUTINE BODY
!
          nz = dfft%npp( me )

          IF ( iopt > 0 ) THEN

            r = CMPLX(0.0d0,0.0d0)

            IF( iopt == 2 ) THEN
              ns = dfft%nsw( me )
            ELSE
              ns = dfft%nsp( me )
            END IF
            DO is = 1, ns
              i = stmask( 1, is, 1 )
              j = stmask( 2, is, 1 )
              DO k = 1 , nz
                r( i, j, k ) = zstick( k, is )
              END DO
            END DO

          ELSE IF ( iopt < 0 ) THEN

            IF( iopt == -2 ) THEN
              ns = dfft%nsw( me )
            ELSE
              ns = dfft%nsp( me )
            END IF
            DO is = 1, ns
              i = stmask( 1, is, 1 )
              j = stmask( 2, is, 1 )
              DO k = 1 , nz
                zstick( k, is ) = r( i, j, k )
              END DO
            END DO

          END IF

          RETURN 
        END SUBROUTINE fft_transpose

#endif


!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine fft_scatter (f_in, nrx3, nxx_, f_aux, ncp_, npp_, sign)
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
#include "machine.h"
  use mp_global, ONLY: nproc_pool, me_pool, intra_pool_comm, nproc
  use parameters, only : DP
  implicit none

  integer, intent(in) :: nrx3, nxx_, sign, ncp_ (:), npp_ (:)
  complex (kind=DP) :: f_in (nxx_), f_aux (nxx_)

#ifdef __PARA

  include 'mpif.h'

  integer :: dest, from, k, offset1 (nproc), sendcount (nproc), &
       sdispls (nproc), recvcount (nproc), rdispls (nproc), &
       proc, ierr, me, nprocp
  !
  me     = me_pool + 1
  nprocp = nproc_pool
  !
  if (nprocp.eq.1) return
  !
  !call start_clock ('fft_scatter')
  !
  ! sendcount(proc): amount of data processor "me" must send to processor
  ! recvcount(proc): amount of data processor "me" must receive from
  !
  do proc = 1, nprocp
     sendcount (proc) = npp_ (proc) * ncp_ (me)
     recvcount (proc) = npp_ (me) * ncp_ (proc)
  enddo
  !
  ! offset1(proc) is used to locate the slices to be sent to proc
  ! sdispls(proc)+1 is the beginning of data that must be sent to proc
  ! rdispls(proc)+1 is the beginning of data that must be received from pr
  !
  offset1 (1) = 1
  sdispls (1) = 0
  rdispls (1) = 0
  do proc = 2, nprocp
     offset1 (proc) = offset1 (proc - 1) + npp_ (proc - 1)
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
        do k = 1, ncp_ (me)
           call DCOPY (2 * npp_ (proc), f_in (from + (k - 1) * nrx3), &
                1, f_aux (dest + (k - 1) * npp_ (proc) ), 1)
        enddo
     enddo
     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     f_in = 0.0d0
     !
     ! step two: communication
     !
     call mpi_barrier (intra_pool_comm, ierr)
     call mpi_alltoallv (f_aux(1), sendcount, sdispls, MPI_DOUBLE_COMPLEX, f_in(1), &
          recvcount, rdispls, MPI_DOUBLE_COMPLEX, intra_pool_comm, ierr)
     if( ABS(ierr) /= 0 ) call errore ('fft_scatter', 'info<>0', ABS(ierr) )
     !
  else
     !
     !  "backward" scatter from planes to columns
     !
     !  step two: communication
     !
     call mpi_barrier (intra_pool_comm, ierr)
     call mpi_alltoallv (f_in(1), recvcount, rdispls, MPI_DOUBLE_COMPLEX, f_aux(1), &
          sendcount, sdispls, MPI_DOUBLE_COMPLEX, intra_pool_comm, ierr)
     if( ABS(ierr) /= 0 ) call errore ('fft_scatter', 'info<>0', ABS(ierr) )
     !
     !  step one: store contiguously the columns
     !
     f_in = 0.0d0
     !
     do proc = 1, nprocp
        from = 1 + sdispls (proc)
        dest = offset1 (proc)
        do k = 1, ncp_ (me)
           call DCOPY (2 * npp_ (proc), f_aux (from + (k - 1) * npp_ ( &
                proc) ), 1, f_in (dest + (k - 1) * nrx3), 1)
        enddo

     enddo

  endif

  !call stop_clock ('fft_scatter')

#endif

  return

end subroutine fft_scatter


!=----------------------------------------------------------------------=!
      END MODULE fft_base
!=----------------------------------------------------------------------=!
