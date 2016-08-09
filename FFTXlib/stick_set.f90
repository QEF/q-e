!
! Copyright (C) Quantum ESPRESSO group
!
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

      USE stick_base
      USE fft_types, ONLY: fft_type_descriptor, fft_type_set, fft_type_scalar

      IMPLICIT NONE

      INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

      PRIVATE
      SAVE

      PUBLIC :: pstickset, pstickset_custom

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

      SUBROUTINE pstickset( gamma_only, bg, gcut, gkcut, gcuts, &
          dfftp, dffts, ngw, ngm, ngs, mype, root, nproc, comm, nogrp_ , &
          ionode, stdout, dtgs, dfft3d )

          USE task_groups, ONLY: task_groups_descriptor, task_groups_init

          LOGICAL, INTENT(in) :: gamma_only
! ...     bg(:,1), bg(:,2), bg(:,3) reciprocal space base vectors.
          REAL(DP), INTENT(in) :: bg(3,3)
          REAL(DP), INTENT(in) :: gcut, gkcut, gcuts
          TYPE(fft_type_descriptor), INTENT(inout) :: dfftp, dffts
          INTEGER, INTENT(out) :: ngw, ngm, ngs

          INTEGER, INTENT(IN) :: mype, root, nproc, comm
          INTEGER, INTENT(IN) :: nogrp_
          LOGICAL, INTENT(IN) :: ionode
          INTEGER, INTENT(IN) :: stdout

          TYPE(task_groups_descriptor), OPTIONAL, INTENT(inout) :: dtgs
          TYPE(fft_type_descriptor), OPTIONAL, INTENT(inout) :: dfft3d

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
        INTEGER :: ip, ngm_ , ngs_, ipg
        INTEGER, ALLOCATABLE :: idx(:)

!
          tk    = .not. gamma_only
          ub(1) = ( dfftp%nr1 - 1 ) / 2
          ub(2) = ( dfftp%nr2 - 1 ) / 2
          ub(3) = ( dfftp%nr3 - 1 ) / 2
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
                            gcut, gkcut, gcuts, st, stw, sts, mype, &
                            nproc, comm )

! ...       Now count the number of stick nst and nstw

          nst  = count( st  > 0 )
          nstw = count( stw > 0 )
          nsts = count( sts > 0 )

          ALLOCATE(ist(nst,5))

          ALLOCATE(nstp(nproc))
          ALLOCATE(sstp(nproc))

          ALLOCATE(nstpw(nproc))
          ALLOCATE(sstpw(nproc))

          ALLOCATE(nstps(nproc))
          ALLOCATE(sstps(nproc))

! ...       initialize the sticks indexes array ist

          CALL sticks_countg( tk, ub, lb, st, stw, sts, &
            ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5) )

! ...       Sorts the sticks according to their length

          ALLOCATE( idx( nst ) )

          CALL sticks_sort( ist(:,4), ist(:,3), ist(:,5), nst, idx, nproc )

          ! ... Set as first stick the stick containing the G=0
          !
          !  DO iss = 1, nst
          !    IF( ist( idx( iss ), 1 ) == 0 .AND. ist( idx( iss ), 2 ) == 0 )  EXIT
          !  END DO
          !  itmp         = idx( 1 )
          !  idx( 1 )   = idx( iss )
          !  idx( iss ) = itmp

          CALL sticks_dist( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts, mype, nproc )

          ngw = sstpw( mype + 1 )
          ngm = sstp( mype + 1 )
          ngs = sstps( mype + 1 )

          CALL sticks_pairup( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts, nproc )

          ! ...   Allocate and Set fft data layout descriptors

#if defined(__MPI)

          CALL fft_type_set( dfftp, tk, nst, ub, lb, idx, ist(:,1), ist(:,2), nstp, nstpw, sstp, sstpw, st, stw )

          CALL fft_type_set( dffts, tk, nsts, ub, lb, idx, ist(:,1), ist(:,2), nstps, nstpw, sstps, sstpw, sts, stw )

          IF( PRESENT( dfft3d ) ) THEN
             DEALLOCATE( stw )
             ALLOCATE( stw( lb(1) : ub(1), lb(2) : ub(2) ) )
             CALL sticks_maps_scalar( (.not.tk), ub, lb, bg(:,1),bg(:,2),bg(:,3), gcut, gkcut, gcuts, stw, ngm_ , ngs_ )
             CALL fft_type_scalar( dfft3d, ub, lb, stw )
          END IF

#else

          DEALLOCATE( stw )
          ALLOCATE( stw( lb(1) : ub(1), lb(2) : ub(2) ) )

          CALL sticks_maps_scalar( (.not.tk), ub, lb, bg(:,1),bg(:,2),bg(:,3), gcut, gkcut, gcuts, stw, ngm_ , ngs_ )

          IF( ngm_ /= ngm ) CALL fftx_error__( ' pstickset ', ' inconsistent ngm ', abs( ngm - ngm_ ) )
          IF( ngs_ /= ngs ) CALL fftx_error__( ' pstickset ', ' inconsistent ngs ', abs( ngs - ngs_ ) )

          CALL fft_type_scalar( dfftp, ub, lb, stw )

          CALL fft_type_scalar( dffts, ub, lb, stw )

#endif

! ...     Maximum number of sticks (potentials)
          nstpx  = maxval( nstp )
! ...     Maximum number of sticks (wave func.)
          nstpwx = maxval( nstpw  )
          !
          !  Initialize task groups.
          !  Note that this call modify dffts adding task group data.
          !
          IF( PRESENT( dtgs ) ) THEN
             CALL task_groups_init( dffts, dtgs, nogrp_ )
          ELSE
             IF( nogrp_ > 1 ) THEN
                CALL fftx_error__( ' pstickset ', &
                & ' more than one Task Groups is requested, but the Task Group descriptor is not initialized ', 1 )
             END IF
          END IF
          !
          IF (ionode) THEN
             WRITE( stdout,*)
             IF ( nproc > 1 ) THEN
                WRITE( stdout, '(5X,"Parallelization info")')
             ELSE
                WRITE( stdout, '(5X,"G-vector sticks info")')
             ENDIF
             WRITE( stdout, '(5X,"--------------------")')
             WRITE( stdout, '(5X,"sticks:   dense  smooth     PW", &
                            & 5X,"G-vecs:    dense   smooth      PW")') 
             IF ( nproc > 1 ) THEN
                WRITE( stdout,'(5X,"Min",4X,2I8,I7,12X,2I9,I8)') &
                   minval(nstp), minval(nstps), minval(nstpw), &
                   minval(sstp), minval(sstps), minval(sstpw)
                WRITE( stdout,'(5X,"Max",4X,2I8,I7,12X,2I9,I8)') &
                   maxval(nstp), maxval(nstps), maxval(nstpw), &
                   maxval(sstp), maxval(sstps), maxval(sstpw)
             END IF
             WRITE( stdout,'(5X,"Sum",4X,2I8,I7,12X,2I9,I8)') &
                sum(nstp), sum(nstps), sum(nstpw), &
                sum(sstp), sum(sstps), sum(sstpw)
             ! in the case k=0, the lines above and below differ:
             ! above all sticks, below only those in the half sphere
             IF ( .NOT. tk ) &
                 WRITE( stdout,'(5X,"Tot",4X,2I8,I7)') nst, nsts, nstw
          ENDIF

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

!----------------------------------------------------------------------

      SUBROUTINE pstickset_custom( gamma_only, bg, gcut, gkcut, gcuts, &
          dfftp, dffts, ngw, ngm, ngs, mype, root, nproc, comm, nogrp_ )

          LOGICAL, INTENT(in) :: gamma_only
! ...     bg(:,1), bg(:,2), bg(:,3) reciprocal space base vectors.
          REAL(DP), INTENT(in) :: bg(3,3)
          REAL(DP), INTENT(in) :: gcut, gkcut, gcuts
          TYPE(fft_type_descriptor), INTENT(inout) :: dfftp, dffts
          INTEGER, INTENT(inout) :: ngw, ngm, ngs

          INTEGER, INTENT(IN) :: mype, root, nproc, comm
          INTEGER, INTENT(IN) :: nogrp_


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
          ub(1) = ( dfftp%nr1 - 1 ) / 2
          ub(2) = ( dfftp%nr2 - 1 ) / 2
          ub(3) = ( dfftp%nr3 - 1 ) / 2
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
                            gcut, gkcut, gcuts, st, stw, sts, mype, &
                            nproc, comm )

! ...       Now count the number of stick nst and nstw

          nst  = count( st  > 0 )
          nstw = count( stw > 0 )
          nsts = count( sts > 0 )

          ALLOCATE(ist(nst,5))

          ALLOCATE(nstp(nproc))
          ALLOCATE(sstp(nproc))

          ALLOCATE(nstpw(nproc))
          ALLOCATE(sstpw(nproc))

          ALLOCATE(nstps(nproc))
          ALLOCATE(sstps(nproc))

! ...       initialize the sticks indexes array ist

          CALL sticks_countg( tk, ub, lb, st, stw, sts, &
            ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5) )

! ...       Sorts the sticks according to their length

          ALLOCATE( idx( nst ) )

          CALL sticks_sort( ist(:,4), ist(:,3), ist(:,5), nst, idx, nproc )

! ...       Distribute the sticks as in dfftp

          CALL sticks_ordered_dist( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts, nproc )

          ngw = sstpw( mype + 1 )
          ngm = sstp( mype + 1 )
          ngs = sstps( mype + 1 )

          CALL sticks_pairup( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts, nproc )

          ! ...   Allocate and Set fft data layout descriptors

#if defined(__MPI)

          CALL fft_type_set( dffts, tk, nsts, ub, lb, idx, ist(:,1), ist(:,2), nstps, nstpw, sstps, sstpw, sts, stw )

#else

          DEALLOCATE( stw )
          ALLOCATE( stw( lb(1) : ub(1), lb(2) : ub(2) ) )

          CALL sticks_maps_scalar( (.not.tk), ub, lb, bg(:,1),bg(:,2),bg(:,3), gcut, gkcut, gcuts, stw, ngm_ , ngs_ )

          IF( ngs_ /= ngs ) CALL fftx_error__( ' pstickset_custom ', ' inconsistent ngs ', abs( ngs - ngs_ ) )

          CALL fft_type_scalar( dffts, ub, lb, stw )

#endif

! ...     Maximum number of sticks (potentials)
          nstpx  = maxval( nstp )
! ...     Maximum number of sticks (wave func.)
          nstpwx = maxval( nstpw  )

          DEALLOCATE( ist )
          DEALLOCATE( idx )

          DEALLOCATE( st, stw, sts )
          DEALLOCATE( sstp )
          DEALLOCATE( nstp )
          DEALLOCATE( sstpw )
          DEALLOCATE( nstpw )
          DEALLOCATE( sstps )
          DEALLOCATE( nstps )

          RETURN
        END SUBROUTINE pstickset_custom


!=----------------------------------------------------------------------=
   END MODULE stick_set
!=----------------------------------------------------------------------=
