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

        PUBLIC :: fft_scatter
        PUBLIC :: dfftp, dffts, dfftb, fft_dlay_descriptor


!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
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

     call mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib

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

     call mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib

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
