!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------=
   MODULE stick_set
!=----------------------------------------------------------------------=

!  ... Distribute G-vectors across processors as sticks and planes,
!  ... initialize FFT descriptors for both dense and smooth grids

!  ... Most important dependencies: next three modules
      USE grid_dimensions
      USE smooth_grid_dimensions
      USE stick_base
!
      USE kinds, ONLY: DP
      USE io_global, ONLY: ionode, stdout
      USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_allocate, &
                           fft_dlay_set, fft_dlay_scalar
      USE mp_global, ONLY: me_pool, nproc_pool, intra_pool_comm, nogrp, use_task_groups
      USE mp_global,      ONLY : me_pool, nproc_pool, intra_pool_comm
      USE mp_global,      ONLY : NOGRP, NPGRP, ogrp_comm, pgrp_comm
      USE mp_global,      ONLY : nolist

      PRIVATE
      SAVE

      PUBLIC :: pstickset

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

      SUBROUTINE pstickset( gamma_only, bg, gcut, gkcut, gcuts, &
          dfftp, dffts, ngw, ngm, ngs )

          LOGICAL, INTENT(in) :: gamma_only
! ...     bg(:,1), bg(:,2), bg(:,3) reciprocal space base vectors.
          REAL(DP), INTENT(in) :: bg(3,3)
          REAL(DP), INTENT(in) :: gcut, gkcut, gcuts
          TYPE(fft_dlay_descriptor), INTENT(inout) :: dfftp, dffts
          INTEGER, INTENT(out) :: ngw, ngm, ngs

          LOGICAL :: tk

          INTEGER :: ub(3), lb(3)
! ...     ub(i), lb(i) upper and lower miller indexes

!
! ...     Plane Waves
!
        INTEGER, ALLOCATABLE :: stw(:,:)
! ...   stick map (wave functions), stw(i,j) = number of G-vector in the
! ...     stick whose x and y miller index are i and j

        INTEGER, ALLOCATABLE :: nstpw(:)
! ...   number of sticks (wave functions), nstpw(ip) = number of stick
! ...     for processor ip

        INTEGER, ALLOCATABLE :: sstpw(:)
! ...   number of G-vectors (wave functions), sstpw(ip) = sum of the
! ...     sticks length for processor ip = number of G-vectors
! ...     owned by the processor ip

        INTEGER :: nstw, nstpwx
! ...   nstw     local number of sticks (wave functions)
! ...   nstpwx   maximum among all processors of nstw

!
! ...     Potentials
!

        INTEGER, ALLOCATABLE :: st(:,:)
! ...   stick map (potentials), st(i,j) = number of G-vector in the
! ...     stick whose x and y miller index are i and j

        INTEGER, ALLOCATABLE :: nstp(:)
! ...   number of sticks (potentials), nstp(ip) = number of stick
! ...     for processor ip

        INTEGER, ALLOCATABLE :: sstp(:)
! ...   number of G-vectors (potentials), sstp(ip) = sum of the
! ...     sticks length for processor ip = number of G-vectors
! ...     owned by the processor ip

        INTEGER :: nst, nstpx
! ...   nst      local number of sticks (potentials)
! ...   nstpx    maximum among all processors of nst

!
! ...     Smooth Mesh
!

        INTEGER, ALLOCATABLE :: sts(:,:)
! ...   stick map (smooth mesh), sts(i,j) = number of G-vector in the
! ...     stick whose x and y miller index are i and j

        INTEGER, ALLOCATABLE :: nstps(:)
! ...   number of sticks (smooth mesh), nstp(ip) = number of stick
! ...     for processor ip

        INTEGER, ALLOCATABLE :: sstps(:)
! ...   number of G-vectors (smooth mesh), sstps(ip) = sum of the
! ...     sticks length for processor ip = number of G-vectors
! ...     owned by the processor ip

        INTEGER :: nsts
! ...   nsts      local number of sticks (smooth mesh)


        INTEGER, ALLOCATABLE :: ist(:,:)    ! sticks indices ordered
          INTEGER :: ip, ngm_ , ngs_
          INTEGER, ALLOCATABLE :: idx(:)

          tk    = .not. gamma_only
          ub(1) = ( nr1 - 1 ) / 2
          ub(2) = ( nr2 - 1 ) / 2
          ub(3) = ( nr3 - 1 ) / 2
          lb    = - ub

          ! ...       Allocate maps

          ALLOCATE( stw ( lb(1):ub(1), lb(2):ub(2) ) )
          ALLOCATE( st  ( lb(1):ub(1), lb(2):ub(2) ) )
          ALLOCATE( sts ( lb(1):ub(1), lb(2):ub(2) ) )

          st  = 0
          stw = 0
          sts = 0

! ...       Fill in the stick maps, for given g-space base and cut-off

          CALL sticks_maps( tk, ub, lb, bg(:,1), bg(:,2), bg(:,3), &
                            gcut, gkcut, gcuts, st, stw, sts, me_pool, &
                            nproc_pool, intra_pool_comm )

! ...       Now count the number of stick nst and nstw

          nst  = count( st  > 0 )
          nstw = count( stw > 0 )
          nsts = count( sts > 0 )

          IF (ionode) THEN
            WRITE( stdout,*)
            WRITE( stdout,10)
 10         FORMAT(3X,'Stick Mesh',/, &
                   3X,'----------')
            WRITE( stdout,15) nst, nstw, nsts
 15         FORMAT( 3X, 'nst =', I6, ',  nstw =', I6, ', nsts =', I6 )
          ENDIF

          ALLOCATE(ist(nst,5))

          ALLOCATE(nstp(nproc_pool))
          ALLOCATE(sstp(nproc_pool))

          ALLOCATE(nstpw(nproc_pool))
          ALLOCATE(sstpw(nproc_pool))

          ALLOCATE(nstps(nproc_pool))
          ALLOCATE(sstps(nproc_pool))

! ...       initialize the sticks indexes array ist

          CALL sticks_countg( tk, ub, lb, st, stw, sts, &
            ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5) )

! ...       Sorts the sticks according to their length

          ALLOCATE( idx( nst ) )

          CALL sticks_sort( ist(:,4), ist(:,3), ist(:,5), nst, idx, nproc_pool )

          ! ... Set as first stick the stick containing the G=0
          !
          !  DO iss = 1, nst
          !    IF( ist( idx( iss ), 1 ) == 0 .AND. ist( idx( iss ), 2 ) == 0 )  EXIT
          !  END DO
          !  itmp         = idx( 1 )
          !  idx( 1 )   = idx( iss )
          !  idx( iss ) = itmp

          CALL sticks_dist( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts, nproc_pool )

          ngw = sstpw( me_pool + 1 )
          ngm = sstp( me_pool + 1 )
          ngs = sstps( me_pool + 1 )

          CALL sticks_pairup( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts, nproc_pool )

          ! ...   Allocate and Set fft data layout descriptors

#if defined __PARA

          CALL fft_dlay_allocate( dfftp, me_pool, nproc_pool, intra_pool_comm, nr1x,  nr2x )
          CALL fft_dlay_allocate( dffts, me_pool, nproc_pool, intra_pool_comm, nr1sx, nr2sx )

          CALL fft_dlay_set( dfftp, tk, nst, nr1, nr2, nr3, nr1x, nr2x, nr3x, (me_pool+1), &
            nproc_pool, intra_pool_comm, nogrp, ub, lb, idx, ist(:,1), ist(:,2), nstp, nstpw, sstp, sstpw, st, stw )
          CALL fft_dlay_set( dffts, tk, nsts, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, (me_pool+1), &
            nproc_pool, intra_pool_comm, nogrp, ub, lb, idx, ist(:,1), ist(:,2), nstps, nstpw, sstps, sstpw, sts, stw )

#else

          DEALLOCATE( stw )
          ALLOCATE( stw( lb(2) : ub(2), lb(3) : ub(3) ) )

          CALL sticks_maps_scalar( (.not.tk), ub, lb, bg(:,1),bg(:,2),bg(:,3),&
                                    gcut, gkcut, gcuts, stw, ngm_ , ngs_ )

          IF( ngm_ /= ngm ) CALL errore( ' pstickset ', ' inconsistent ngm ', abs( ngm - ngm_ ) )
          IF( ngs_ /= ngs ) CALL errore( ' pstickset ', ' inconsistent ngs ', abs( ngs - ngs_ ) )

          CALL fft_dlay_allocate( dfftp, me_pool, nproc_pool, intra_pool_comm, max(nr1x, nr3x),  nr2x  )
          CALL fft_dlay_allocate( dffts, me_pool, nproc_pool, intra_pool_comm, max(nr1sx, nr3sx), nr2sx )

          CALL fft_dlay_scalar( dfftp, ub, lb, nr1, nr2, nr3, nr1x, nr2x, nr3x, stw )
          CALL fft_dlay_scalar( dffts, ub, lb, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, stw )

#endif
      !   set the dimensions of the arrays allocated for the FFT
      nrxx  = dfftp % nnr
      nrxxs = dffts % nnr          

! ...     Maximum number of sticks (potentials)
          nstpx  = maxval( nstp )
! ...     Maximum number of sticks (wave func.)
          nstpwx = maxval( nstpw  )

          IF( use_task_groups ) THEN
            !
            !  Initialize task groups.
            !  Note that this call modify dffts adding task group data.
            !
            CALL task_groups_init( dffts )
            !
          END IF

          IF (ionode) WRITE( stdout,118)
 118      FORMAT(3X,'            n.st   n.stw   n.sts    n.g    n.gw   n.gs')
          WRITE( stdout,121) minval(nstp),  minval(nstpw), minval(nstps), minval(sstp), minval(sstpw), minval(sstps)
          WRITE( stdout,122) maxval(nstp),  maxval(nstpw), maxval(nstps), maxval(sstp), maxval(sstpw), maxval(sstps)
!          DO ip = 1, nproc_pool
!            IF (ionode) THEN
!              WRITE( stdout,120) ip, nstp(ip),  nstpw(ip), nstps(ip), sstp(ip), sstpw(ip), sstps(ip)
!            END IF
!          END DO
          IF (ionode) THEN
            WRITE( stdout,120)  sum(nstp),  sum(nstpw), sum(nstps), sum(sstp), sum(sstpw), sum(sstps)
          ENDIF
 120      FORMAT(3X,7I8)
 121      FORMAT(3X,'min    ',6I8)
 122      FORMAT(3X,'max    ',6I8)


          DEALLOCATE( ist )
          DEALLOCATE( idx )

          DEALLOCATE( st, stw, sts )
          DEALLOCATE( sstp )
          DEALLOCATE( nstp )
          DEALLOCATE( sstpw )
          DEALLOCATE( nstpw )
          DEALLOCATE( sstps )
          DEALLOCATE( nstps )

          IF(ionode) WRITE( stdout,*)

          RETURN
        END SUBROUTINE pstickset


!-----------------------------------------
! Task groups Contributed by C. Bekas, October 2005
! Revised by C. Cavazzoni
!--------------------------------------------

SUBROUTINE task_groups_init( dffts )

   USE parallel_include
   !
   USE io_global,      ONLY : stdout
   USE fft_types,      ONLY : fft_dlay_descriptor

   ! T.G.
   ! NPGRP:      Number of processors per group
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE

   TYPE(fft_dlay_descriptor), INTENT(inout) :: dffts

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
   IF ( nogrp > 1 ) WRITE( stdout, 100 ) nogrp, npgrp

100 FORMAT( /,3X,'Task Groups are in USE',/,3X,'groups and procs/group : ',I5,I5 )

   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space

#if defined __MPI
   CALL MPI_Allgather( dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, intra_pool_comm, IERR)
   strd = maxval( nnrsx_vec( 1:nproc_pool ) )
#else
   strd = dffts%nnr
#endif

   IF( strd /= dffts%tg_nnr ) CALL errore( ' task_groups_init ', ' inconsistent nnr ', 1 )

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
   dffts%nogrp = nogrp
   dffts%npgrp = npgrp
   dffts%ogrp_comm = ogrp_comm
   dffts%pgrp_comm = pgrp_comm
   !
   ALLOCATE( dffts%tg_nsw(nproc_pool))
   ALLOCATE( dffts%tg_npp(nproc_pool))
   ALLOCATE( dffts%nolist(nogrp))

   num_sticks = 0
   num_planes = 0
   DO i = 1, nogrp
      dffts%nolist( i ) = nolist( i )
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
   IF( dffts%nr3x * dffts%nsw( me_pool + 1 ) > dffts%tg_nnr ) THEN
      CALL errore( ' task_groups_init ', ' inconsistent dffts%tg_nnr ', 1 )
   ENDIF
   dffts%tg_psdsp(1) = 0
   dffts%tg_usdsp(1) = 0
   dffts%tg_rcv(1)  = dffts%nr3x * dffts%nsw( nolist(1) + 1 )
   dffts%tg_rdsp(1) = 0
   DO i = 2, nogrp
      dffts%tg_snd(i)  = dffts%nr3x * dffts%nsw( me_pool + 1 )
      dffts%tg_psdsp(i) = dffts%tg_psdsp(i-1) + dffts%tg_nnr
      dffts%tg_usdsp(i) = dffts%tg_usdsp(i-1) + dffts%tg_snd(i-1)
      dffts%tg_rcv(i)  = dffts%nr3x * dffts%nsw( nolist(i) + 1 )
      dffts%tg_rdsp(i) = dffts%tg_rdsp(i-1) + dffts%tg_rcv(i-1)
   ENDDO

   dffts%have_task_groups = .true.

   RETURN

END SUBROUTINE task_groups_init


!=----------------------------------------------------------------------=
   END MODULE stick_set
!=----------------------------------------------------------------------=
