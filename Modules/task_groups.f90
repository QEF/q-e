!
! Copyright (C) 2002-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------
! Contributed by C. Bekas, October 2005
! Revised by C. Cavazzoni
!--------------------------------------------

MODULE task_groups

   USE kinds,      ONLY: DP

   IMPLICIT NONE
   SAVE


CONTAINS


!========================================================================================
! ADDED SUBROUTINEs FOR TASK GROUP PARALLIZATION
! C. Bekas, IBM Research, Zurich
!        - GROUPS: Define and initialize Task Groups
!        - tg_ivfftw: Inverse FFT driver for Task Groups
!=======================================================================================


!-----------------------------------------------------------------------
!      SUBROUTINE GROUPS (added by C. Bekas)
!      Define groups for task group parallilization
!-----------------------------------------------------------------------

SUBROUTINE task_groups_init( dffts )

   USE parallel_include
   !
   USE mp_global,      ONLY : me_pool, nproc_pool, intra_pool_comm
   USE mp_global,      ONLY : NOGRP, NPGRP, ogrp_comm, pgrp_comm  
   USE mp_global,      ONLY : nolist, nplist
   USE io_global,      only : stdout
   USE fft_types,      only : fft_dlay_descriptor

   ! T.G. 
   ! NPGRP:      Number of processors per group
   ! NOGRP:      Number of group

   IMPLICIT NONE 

   TYPE(fft_dlay_descriptor), INTENT(INOUT) :: dffts

   !----------------------------------
   !Local Variables declaration
   !----------------------------------

   INTEGER  :: I
   INTEGER  :: IERR
   INTEGER  :: num_planes, num_sticks
   INTEGER  :: nnrsx_vec ( nproc_pool )
   INTEGER  :: pgroup( nproc_pool )
   INTEGER  :: strd

   !
   WRITE( stdout, 100 ) nogrp, npgrp

100 FORMAT( /,3X,'Task Groups are in use',/,3X,'groups and procs/group : ',I5,I5 )

   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space

#if defined __MPI
   CALL MPI_Allgather( dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, intra_pool_comm, IERR)
   strd = MAXVAL( nnrsx_vec( 1:nproc_pool ) )
#else
   strd = dffts%nnr 
#endif

   IF( strd /= dffts%nnrx ) CALL errore( ' task_groups_init ', ' inconsistent nnrx ', 1 )

   !-------------------------------------------------------------------------------------
   !C. Bekas...TASK GROUP RELATED. FFT DATA STRUCTURES ARE ALREADY DEFINED ABOVE
   !-------------------------------------------------------------------------------------
   !dfft%nsw(me) holds the number of z-sticks for the current processor per wave-function
   !We can either send these in the group with an mpi_allgather...or put the
   !in the PSIS vector (in special positions) and send them with them.
   !Otherwise we can do this once at the beginning, before the loop.
   !we choose to do the latter one.
   !-------------------------------------------------------------------------------------
   !
   ALLOCATE( dffts%tg_nsw(nproc_pool))
   ALLOCATE( dffts%tg_npp(nproc_pool))

   num_sticks = 0
   num_planes = 0
   DO i = 1, nogrp
      num_sticks = num_sticks + dffts%nsw( nolist(i) + 1 )
      num_planes = num_planes + dffts%npp( nolist(i) + 1 )
   ENDDO

#if defined __MPI
   CALL MPI_ALLGATHER(num_sticks, 1, MPI_INTEGER, dffts%tg_nsw(1), 1, MPI_INTEGER, intra_pool_comm, IERR)
   CALL MPI_ALLGATHER(num_planes, 1, MPI_INTEGER, dffts%tg_npp(1), 1, MPI_INTEGER, intra_pool_comm, IERR)
#else
   dffts%tg_nsw(1) = num_sticks
   dffts%tg_npp(1) = num_planes
#endif

   ALLOCATE( dffts%tg_snd( nogrp ) )
   ALLOCATE( dffts%tg_rcv( nogrp ) )
   ALLOCATE( dffts%tg_psdsp( nogrp ) )
   ALLOCATE( dffts%tg_usdsp( nogrp ) )
   ALLOCATE( dffts%tg_rdsp( nogrp ) )

   dffts%tg_snd(1)   = dffts%nr3x * dffts%nsw( me_pool + 1 )
   IF( dffts%nr3x * dffts%nsw( me_pool + 1 ) > dffts%nnrx ) THEN
      CALL errore( ' task_groups_init ', ' inconsistent dffts%nnrx ', 1 )
   END IF
   dffts%tg_psdsp(1) = 0
   dffts%tg_usdsp(1) = 0
   dffts%tg_rcv(1)  = dffts%nr3x * dffts%nsw( nolist(1) + 1 )
   dffts%tg_rdsp(1) = 0
   DO i = 2, nogrp
      dffts%tg_snd(i)  = dffts%nr3x * dffts%nsw( me_pool + 1 )
      dffts%tg_psdsp(i) = dffts%tg_psdsp(i-1) + dffts%nnrx
      dffts%tg_usdsp(i) = dffts%tg_usdsp(i-1) + dffts%tg_snd(i-1)
      dffts%tg_rcv(i)  = dffts%nr3x * dffts%nsw( nolist(i) + 1 )
      dffts%tg_rdsp(i) = dffts%tg_rdsp(i-1) + dffts%tg_rcv(i-1)
   ENDDO
  
   ! ALLOCATE( dffts%tg_sca_snd( nproc_pool / nogrp )  )
   ! ALLOCATE( dffts%tg_sca_rcv( nproc_pool / nogrp )  )
   ! ALLOCATE( dffts%tg_sca_sdsp( nproc_pool / nogrp )  )
   ! ALLOCATE( dffts%tg_sca_rdsp( nproc_pool / nogrp )  )
   ! ALLOCATE( dffts%tg_sca_off( nproc_pool / nogrp )  )

   ! do i = 1, nproc_pool / nogrp
   !    dffts%tg_sca_snd (i) = dffts%tg_npp ( nplist( i ) + 1 ) * dffts%tg_nsw ( me_pool + 1 )
   !    dffts%tg_sca_rcv (i) = dffts%tg_npp ( me_pool + 1 )     * dffts%tg_nsw ( nplist( i ) + 1 )
   ! end do
   ! dffts%tg_sca_off(1) = 0
   ! do i = 2, nproc_pool / nogrp
   !    dffts%tg_sca_off(i) = dffts%tg_sca_off(i - 1) + dffts%tg_npp ( nplist( i - 1 ) + 1 )
   ! end do
   ! dffts%tg_sca_sdsp (1) = 0
   ! dffts%tg_sca_rdsp (1) = 0
   ! do i = 2, nproc_pool / nogrp
   !    dffts%tg_sca_sdsp (i) = dffts%tg_sca_sdsp (i - 1) + dffts%tg_sca_snd (i - 1)
   !    dffts%tg_sca_rdsp (i) = dffts%tg_sca_rdsp (i - 1) + dffts%tg_sca_rcv (i - 1)
   ! enddo

   dffts%have_task_groups = .TRUE.

   RETURN

END SUBROUTINE task_groups_init

!

SUBROUTINE tg_gather( dffts, v, tg_v )
   !
   USE parallel_include
   !
   USE mp_global,      ONLY : me_pool, nogrp, ogrp_comm, nolist
   USE fft_types,      only : fft_dlay_descriptor

   ! T.G. 
   ! NPGRP:      Number of processors per group
   ! NOGRP:      Number of group

   IMPLICIT NONE 

   TYPE(fft_dlay_descriptor), INTENT(IN) :: dffts

   REAL(DP) :: v(:)
   REAL(DP) :: tg_v(:)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( nogrp ), recv_displ( nogrp )

   nsiz_tg = dffts%nnrx * nogrp

   IF( SIZE( tg_v ) < nsiz_tg ) &
      call errore( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - SIZE( tg_v ) ) )

   nsiz = dffts%npp( me_pool+1 ) * dffts%nr1x * dffts%nr2x

   IF( SIZE( v ) < nsiz ) &
      call errore( ' tg_gather ', ' v too small ',  ( nsiz - SIZE( v ) ) )

   !
   !  The potential in v is distributed accros all processors
   !  We need to redistribute it so that it is completely contained in the
   !  processors of an orbital TASK-GROUP
   !
   recv_cnt(1)   = dffts%npp( nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
   recv_displ(1) = 0
   DO i = 2, NOGRP
      recv_cnt(i) = dffts%npp( nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
      recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
   ENDDO
 
   ! clean only elements that will not be overwritten
   !
   DO i = recv_displ(nogrp) + recv_cnt( nogrp ) + 1, SIZE( tg_v )
      tg_v( i ) = 0.0d0
   END DO

#if defined (__PARA) && defined (__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, ogrp_comm, IERR)
   !
   IF( ierr /= 0 ) &
      call errore( ' tg_gather ', ' MPI_Allgatherv ', ABS( ierr ) )

#endif


   RETURN
END SUBROUTINE


END MODULE task_groups
