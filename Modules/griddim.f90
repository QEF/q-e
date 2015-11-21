!
! Copyright (C) 2002-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE grid_subroutines
!=----------------------------------------------------------------------------=!

     ! This module contains subroutines that are related to grids 
     ! parameters

     USE kinds,     ONLY: DP
     USE fft_types, ONLY: fft_dlay_descriptor

     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: realspace_grid_init, realspace_grids_info

   CONTAINS

     SUBROUTINE realspace_grid_init( dfft, at, bg, gcutm, fft_fact )
       !
       ! ... Sets optimal values for dfft%nr[123] and dfft%nr[123]x
       ! ... If fft_fact is present, force nr[123] to be multiple of fft_fac([123])
       !
       USE fft_support, only: good_fft_dimension, good_fft_order
       USE io_global, only: stdout
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
       REAL(DP), INTENT(IN) :: gcutm
       INTEGER, INTENT(IN), OPTIONAL :: fft_fact(3)
       TYPE(fft_dlay_descriptor), INTENT(INOUT) :: dfft
       !
       IF( dfft%nr1 == 0 .OR. dfft%nr2 == 0 .OR. dfft%nr3 == 0 ) THEN
         !
         ! ... calculate the size of the real-space dense grid for FFT
         ! ... first, an estimate of nr1,nr2,nr3, based on the max values
         ! ... of n_i indices in:   G = i*b_1 + j*b_2 + k*b_3
         ! ... We use G*a_i = n_i => n_i .le. |Gmax||a_i|
         !
         dfft%nr1 = int ( sqrt (gcutm) * &
               sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
         dfft%nr2 = int ( sqrt (gcutm) * &
               sqrt (at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2) ) + 1
         dfft%nr3 = int ( sqrt (gcutm) * &
               sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1
         !
         CALL grid_set( bg, gcutm, dfft%nr1, dfft%nr2, dfft%nr3 )
         !
       ELSE
          WRITE( stdout, '( /, 3X,"Info: using nr1, nr2, nr3 values from input" )' )
       END IF

       IF (PRESENT(fft_fact)) THEN
          dfft%nr1 = good_fft_order( dfft%nr1, fft_fact(1) )
          dfft%nr2 = good_fft_order( dfft%nr2, fft_fact(2) )
          dfft%nr3 = good_fft_order( dfft%nr3, fft_fact(3) )
       ELSE
          dfft%nr1 = good_fft_order( dfft%nr1 )
          dfft%nr2 = good_fft_order( dfft%nr2 )
          dfft%nr3 = good_fft_order( dfft%nr3 )
       END IF

       dfft%nr1x  = good_fft_dimension( dfft%nr1 )
       dfft%nr2x  = dfft%nr2
       dfft%nr3x  = good_fft_dimension( dfft%nr3 )

     END SUBROUTINE realspace_grid_init

!=----------------------------------------------------------------------------=!

    SUBROUTINE realspace_grids_info ( dfftp, dffts, nproc_ )

      !  Print info on local and global dimensions for real space grids

      USE io_global, ONLY: ionode, stdout
      USE fft_types, ONLY: fft_dlay_descriptor

      IMPLICIT NONE

      TYPE(fft_dlay_descriptor), INTENT(IN) :: dfftp, dffts
      INTEGER, INTENT(IN) :: nproc_

      INTEGER :: i

      IF(ionode) THEN

        WRITE( stdout,*)
        WRITE( stdout,*) '  Real Mesh'
        WRITE( stdout,*) '  ---------'
        WRITE( stdout,1000) dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1, dfftp%nr2, dfftp%npl, 1, 1, nproc_
        WRITE( stdout,1010) dfftp%nr1x, dfftp%nr2x, dfftp%nr3x
        WRITE( stdout,1020) dfftp%nnr
        WRITE( stdout,*) '  Number of x-y planes for each processors: '
        WRITE( stdout, fmt = '( 3X, "nr3l = ", 10I5 )' ) &
           ( dfftp%npp( i ), i = 1, nproc_ )

        WRITE( stdout,*)
        WRITE( stdout,*) '  Smooth Real Mesh'
        WRITE( stdout,*) '  ----------------'
        WRITE( stdout,1000) dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1, dffts%nr2, dffts%npl,1,1, nproc_
        WRITE( stdout,1010) dffts%nr1x, dffts%nr2x, dffts%nr3x
        WRITE( stdout,1020) dffts%nnr
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
      USE mp_images, ONLY: me_image, nproc_image,  intra_image_comm

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

