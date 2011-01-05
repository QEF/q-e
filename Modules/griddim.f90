!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE grid_dimensions
!=----------------------------------------------------------------------------=!

     !  Dimensions of the 3D real and reciprocal space FFT grid
     !  relative to the charge density and potential ("dense" grid)

     IMPLICIT NONE
     SAVE

     !  dimensions of the "dense" 3D grid (global)
     INTEGER :: nr1  = 0, nr2  = 0, nr3  = 0

     !  dimensions of the arrays for the "dense" 3D grid (global)
     !  may differ from nr1 ,nr2 ,nr3 in order to boost performances
     INTEGER :: nr1x = 0, nr2x = 0, nr3x = 0

     ! size of the arrays allocated for the FFT, local to each processor:
     ! in parallel execution may differ from nr1x*nr2x*nr3x
     ! Not to be confused either with nr1*nr2*nr3 
     INTEGER :: nrxx  = 0

     PRIVATE
     PUBLIC :: nr1, nr2,nr3, nr1x,nr2x,nr3x, nrxx

!=----------------------------------------------------------------------------=!
   END MODULE grid_dimensions
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE smooth_grid_dimensions
!=----------------------------------------------------------------------------=!

     !  This module contains the dimensions of the 3D real and reciprocal space
     !  FFT grid relative to the smooth part of the charge density
     !  (may differ from the full charge density grid for USPP )

     IMPLICIT NONE
     SAVE

     !  parameter description: same as above but for smooth grid
     INTEGER :: nr1s = 0, nr2s = 0, nr3s = 0
     INTEGER :: nr1sx= 0, nr2sx= 0, nr3sx= 0
     INTEGER :: nrxxs = 0

     PRIVATE
     PUBLIC :: nr1s, nr2s,nr3s, nr1sx,nr2sx,nr3sx, nrxxs

!=----------------------------------------------------------------------------=!
   END MODULE smooth_grid_dimensions
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE grid_subroutines
!=----------------------------------------------------------------------------=!

     ! This module contains subroutines that are related to grids 
     ! parameters

     USE kinds, ONLY: DP
     USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x
     USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx

     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: realspace_grids_init, realspace_grids_info

   CONTAINS


     SUBROUTINE realspace_grids_init( at, bg, gcutm, gcuts )
       !
       USE fft_scalar, only: good_fft_dimension, good_fft_order
       USE io_global, only: stdout
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
       REAL(DP), INTENT(IN) :: gcutm, gcuts
       !
       IF( nr1 == 0 .OR. nr2 == 0 .OR. nr3 == 0 ) THEN
         !
         ! ... calculate the size of the real-space dense grid for FFT
         ! ... first, an estimate of nr1,nr2,nr3, based on the max values
         ! ... of n_i indices in:   G = i*b_1 + j*b_2 + k*b_3   
         ! ... We use G*a_i = n_i => n_i .le. |Gmax||a_i|
         !
         nr1 = int ( sqrt (gcutm) * &
               sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
         nr2 = int ( sqrt (gcutm) * &
               sqrt (at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2) ) + 1
         nr3 = int ( sqrt (gcutm) * &
               sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1
         !
         CALL grid_set( bg, gcutm, nr1, nr2, nr3 )
         !
       ELSE
         WRITE( stdout, '( /, 3X,"Info: using nr1, nr2, nr3 values from input" )' )
       END IF

       nr1 = good_fft_order( nr1 )
       nr2 = good_fft_order( nr2 )
       nr3 = good_fft_order( nr3 )

       nr1x  = good_fft_dimension( nr1 )
       nr2x  = nr2
       nr3x  = good_fft_dimension( nr3 )

       ! ... As above, for the smooth grid

       IF( nr1s == 0 .OR. nr2s == 0 .OR. nr3s == 0 ) THEN
         !
         IF ( gcuts == gcutm ) THEN
            ! ... No double grid, the two grids are the same
            nr1s = nr1 ; nr2s = nr2 ; nr3s = nr3
            nr1sx= nr1x; nr2sx= nr2x; nr3sx= nr3x
            RETURN
         END IF
         !
         nr1s= int (2 * sqrt (gcuts) * &
               sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
         nr2s= int (2 * sqrt (gcuts) * &
               sqrt (at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2) ) + 1
         nr3s= int (2 * sqrt (gcuts) * &
               sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1
         !
         CALL grid_set( bg, gcuts, nr1s, nr2s, nr3s )
         !
       ELSE
         WRITE( stdout, '( /, 3X,"Info: using nr1s, nr2s, nr3s values from input" )' )
       END IF

       nr1s = good_fft_order( nr1s )
       nr2s = good_fft_order( nr2s )
       nr3s = good_fft_order( nr3s )

       nr1sx = good_fft_dimension(nr1s)
       nr2sx = nr2s
       nr3sx = good_fft_dimension(nr3s)

       IF ( nr1s > nr1 .or. nr2s > nr2 .or. nr3s > nr3 ) THEN
          CALL errore(' realspace_grids_init ', ' smooth grid larger than dense grid?',1)
       END IF

       RETURN

     END SUBROUTINE realspace_grids_init

!=----------------------------------------------------------------------------=!

    SUBROUTINE realspace_grids_info ( dfftp, dffts, nproc_ )

      !  Print info on local and global dimensions for real space grids

      USE io_global, ONLY: ionode, stdout
      USE fft_types, ONLY: fft_dlay_descriptor
      USE grid_dimensions, ONLY: nrxx
      USE smooth_grid_dimensions, ONLY: nrxxs

      IMPLICIT NONE

      TYPE(fft_dlay_descriptor), INTENT(IN) :: dfftp, dffts
      INTEGER, INTENT(IN) :: nproc_

      INTEGER :: i

      IF(ionode) THEN

        WRITE( stdout,*)
        WRITE( stdout,*) '  Real Mesh'
        WRITE( stdout,*) '  ---------'
        WRITE( stdout,1000) nr1, nr2, nr3, nr1, nr2, dfftp%npl, 1, 1, nproc_
        WRITE( stdout,1010) nr1x, nr2x, nr3x
        WRITE( stdout,1020) nrxx
        WRITE( stdout,*) '  Number of x-y planes for each processors: '
        WRITE( stdout, fmt = '( 3X, "nr3l = ", 10I5 )' ) &
           ( dfftp%npp( i ), i = 1, nproc_ )

        WRITE( stdout,*)
        WRITE( stdout,*) '  Smooth Real Mesh'
        WRITE( stdout,*) '  ----------------'
        WRITE( stdout,1000) nr1s, nr2s, nr3s, nr1s, nr2s, dffts%npl,1,1, nproc_
        WRITE( stdout,1010) nr1sx, nr2sx, nr3sx
        WRITE( stdout,1020) nrxxs
        WRITE( stdout,*) '  Number of x-y planes for each processors: '
        WRITE( stdout, fmt = '( 3X, "nr3sl = ", 10I5 )' ) &
           ( dffts%npp( i ), i = 1, nproc_ )

      END IF

1000  FORMAT(3X, &
         'Global Dimensions   Local  Dimensions   Processor Grid',/,3X, &
         '.X.   .Y.   .Z.     .X.   .Y.   .Z.     .X.   .Y.   .Z.',/, &
         3(1X,I5),2X,3(1X,I5),2X,3(1X,I5) )
1010  FORMAT(3X, 'Array leading dimensions ( nr1x, nr2x, nr3x )   = ', 3(1X,I5) )
1020  FORMAT(3X, 'Local number of cell to store the grid ( nrxx ) = ', 1X, I9 )

      RETURN
      END SUBROUTINE realspace_grids_info


   SUBROUTINE grid_set( bg, gcut, nr1, nr2, nr3 )

!  this routine returns in nr1, nr2, nr3 the minimal 3D real-space FFT 
!  grid required to fit the G-vector sphere with G^2 <= gcut
!  On input, nr1,nr2,nr3 must be set to values that match or exceed
!  the largest i,j,k (Miller) indices in G(i,j,k) = i*b1 + j*b2 + k*b3
!  ----------------------------------------------

! ... declare modules
      USE kinds, ONLY: DP
      USE mp, ONLY: mp_max, mp_min, mp_sum
      USE mp_global, ONLY: me_image, nproc_image,  intra_image_comm

      IMPLICIT NONE

! ... declare arguments
      INTEGER, INTENT(INOUT) :: nr1, nr2, nr3
      REAL(DP), INTENT(IN) :: bg(3,3), gcut

! ... declare other variables
      INTEGER :: i, j, k, nr, nb(3)
      REAL(DP) :: gsq, g(3)

!  ----------------------------------------------

      nb     = 0

! ... calculate moduli of G vectors and the range of indices where
! ... |G|^2 < gcut (in parallel whenever possible)

      DO k = -nr3, nr3
        !
        ! ... me_image = processor number, starting from 0
        !
        IF( MOD( k + nr3, nproc_image ) == me_image ) THEN
          DO j = -nr2, nr2
            DO i = -nr1, nr1

              g( 1 ) = DBLE(i)*bg(1,1) + DBLE(j)*bg(1,2) + DBLE(k)*bg(1,3)
              g( 2 ) = DBLE(i)*bg(2,1) + DBLE(j)*bg(2,2) + DBLE(k)*bg(2,3)
              g( 3 ) = DBLE(i)*bg(3,1) + DBLE(j)*bg(3,2) + DBLE(k)*bg(3,3)

! ...         calculate modulus

              gsq =  g( 1 )**2 + g( 2 )**2 + g( 3 )**2 

              IF( gsq < gcut ) THEN

! ...           calculate maximum index
                nb(1) = MAX( nb(1), ABS( i ) )
                nb(2) = MAX( nb(2), ABS( j ) )
                nb(3) = MAX( nb(3), ABS( k ) )
              END IF

            END DO
          END DO
        END IF
      END DO

      CALL mp_max( nb,  intra_image_comm )

! ... the size of the required (3-dimensional) matrix depends on the
! ... maximum indices. Note that the following choice is slightly
! ... "small": 2*nb+2 would be needed in order to guarantee that the
! ...  sphere in G-space never overlaps its periodic image

      nr1 = 2 * nb(1) + 1
      nr2 = 2 * nb(2) + 1
      nr3 = 2 * nb(3) + 1

      RETURN
   
   END SUBROUTINE grid_set

!=----------------------------------------------------------------------------=!
   END MODULE grid_subroutines
!=----------------------------------------------------------------------------=!

