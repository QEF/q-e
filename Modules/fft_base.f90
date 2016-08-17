!
! Copyright (C) 2006-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! FFT data Module.
! Written by Carlo Cavazzoni
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE fft_base
!=----------------------------------------------------------------------=!

        USE parallel_include

        USE fft_types, ONLY: fft_type_descriptor
        USE task_groups, ONLY: task_groups_descriptor
        USE fft_smallbox_type, ONLY: fft_box_descriptor
        USE stick_base, ONLY: sticks_map, sticks_map_deallocate


        IMPLICIT NONE

        ! ... data structure containing all information
        ! ... about fft data distribution for a given
        ! ... potential grid, and its wave functions sub-grid.

        TYPE ( fft_type_descriptor ) :: dfftp ! descriptor for dense grid
             !  Dimensions of the 3D real and reciprocal space FFT grid
             !  relative to the charge density and potential ("dense" grid)
        TYPE ( fft_type_descriptor ) :: dffts ! descriptor for smooth grid
             !  Dimensions of the 3D real and reciprocal space
             !  FFT grid relative to the smooth part of the charge density
             !  (may differ from the full charge density grid for USPP )
        TYPE ( fft_box_descriptor ) :: dfftb ! descriptor for box grids
             !  Dimensions of the 3D real and reciprocal space
             !  FFT grid relative to the "small box" computation
             !  of the atomic augmentation part of the 
             !  charge density used in USPP (to speed up CPV iterations)
        TYPE ( fft_type_descriptor ) :: dfft3d
             !
        TYPE ( task_groups_descriptor ) :: dtgs
             !  Dimensions of the task groups
        TYPE (sticks_map) :: smap
             !  Stick map descriptor

        SAVE

        PRIVATE

        PUBLIC :: dfftp, dffts, dfft3d, fft_type_descriptor
        PUBLIC :: dtgs, task_groups_descriptor
        PUBLIC :: dfftb, fft_box_descriptor, fft_base_info
        PUBLIC :: smap, pstickdealloc

   CONTAINS


      SUBROUTINE pstickdealloc()
         CALL sticks_map_deallocate( smap )
      END SUBROUTINE pstickdealloc


      SUBROUTINE fft_base_info( ionode, stdout )

          LOGICAL, INTENT(IN) :: ionode
          INTEGER, INTENT(IN) :: stdout
          !
          !  Display fft basic information
          !
          IF (ionode) THEN
             WRITE( stdout,*)
             IF ( dfftp%nproc > 1 ) THEN
                WRITE( stdout, '(5X,"Parallelization info")')
             ELSE
                WRITE( stdout, '(5X,"G-vector sticks info")')
             ENDIF
             WRITE( stdout, '(5X,"--------------------")')
             WRITE( stdout, '(5X,"sticks:   dense  smooth     PW", &
                            & 5X,"G-vecs:    dense   smooth      PW")') 
             IF ( dfftp%nproc > 1 ) THEN
                WRITE( stdout,'(5X,"Min",4X,2I8,I7,12X,2I9,I8)') &
                   minval(dfftp%nsp), minval(dffts%nsp), minval(dffts%nsw), &
                   minval(dfftp%ngl), minval(dffts%ngl), minval(dffts%nwl)
                WRITE( stdout,'(5X,"Max",4X,2I8,I7,12X,2I9,I8)') &
                   maxval(dfftp%nsp), maxval(dffts%nsp), maxval(dffts%nsw), &
                   maxval(dfftp%ngl), maxval(dffts%ngl), maxval(dffts%nwl)
             END IF
             WRITE( stdout,'(5X,"Sum",4X,2I8,I7,12X,2I9,I8)') &
                sum(dfftp%nsp), sum(dffts%nsp), sum(dffts%nsw), &
                sum(dfftp%ngl), sum(dffts%ngl), sum(dffts%nwl)
          ENDIF

          IF(ionode) WRITE( stdout,*)

          RETURN
        END SUBROUTINE fft_base_info

!=----------------------------------------------------------------------=!
   END MODULE fft_base
!=----------------------------------------------------------------------=!
