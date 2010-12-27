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

     !  dimensions of the "dense" 3D grid (local on each processor)
     INTEGER :: nr1l = 0, nr2l = 0, nr3l = 0

     ! size of the arrays allocated for the FFT, local to each processor:
     ! in parallel execution may differ from nr1x*nr2x*nr3x
     ! Not to be confused either with nr1*nr2*nr3 
     INTEGER :: nrxx  = 0

     PRIVATE
     PUBLIC :: nr1, nr2,nr3, nr1x,nr2x,nr3x, nrxx
     PUBLIC :: nr1l, nr2l,nr3l
     !PUBLIC :: grid_init

   CONTAINS

     SUBROUTINE grid_init( b1, b2, b3, gcutm )
       !
       USE fft_scalar, only: good_fft_dimension, good_fft_order
       USE io_global, only: stdout
       USE kinds, ONLY: DP
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: b1(3), b2(3), b3(3)
       REAL(DP), INTENT(IN) :: gcutm

       IF( nr1 == 0 .OR. nr2 == 0 .OR. nr3 == 0 ) THEN
         ! ... calculate the size of the real and reciprocal dense grids
         !!!CALL ngnr_set( b1, b2, b3, gcutm, qk, ng, nr1, nr2, nr3 )
       ELSE
         WRITE( stdout, '(/,3X,"Info: using nr1, nr2, nr3 values from input")')
       END IF

       nr1 = good_fft_order( nr1 )
       nr2 = good_fft_order( nr2 )
       nr3 = good_fft_order( nr3 )

       nr1x  = good_fft_dimension( nr1 )
       nr2x  = nr2
       nr3x  = good_fft_dimension( nr3 )

     END SUBROUTINE grid_init

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
     INTEGER :: nr1sl= 0, nr2sl= 0, nr3sl= 0
     INTEGER :: nrxxs = 0

     PRIVATE
     PUBLIC :: nr1s, nr2s,nr3s, nr1sx,nr2sx,nr3sx, nrxxs
     PUBLIC :: nr1sl, nr2sl,nr3sl
     !PUBLIC :: smooth_grid_init

   CONTAINS

     SUBROUTINE smooth_grid_init( b1, b2, b3, gcutms )
       !
       USE fft_scalar, only: good_fft_dimension, good_fft_order
       USE io_global, only: stdout
       USE kinds, ONLY: DP
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: b1(3), b2(3), b3(3)
       REAL(DP), INTENT(IN) :: gcutms

       IF( nr1s == 0 .OR. nr2s == 0 .OR. nr3s == 0 ) THEN
         ! ... calculate the size of the real and reciprocal dense grids
         !!!CALL ngnr_set( b1, b2, b3, gcutms, qk, ng, nr1s, nr2s, nr3s )
       ELSE
         WRITE( stdout, '(/,3X,"Info: using nr1, nr2, nr3 values from input")')
       END IF

       nr1s = good_fft_order( nr1s )
       nr2s = good_fft_order( nr2s )
       nr3s = good_fft_order( nr3s )

       nr1sx  = good_fft_dimension( nr1s )
       nr2sx  = nr2s
       nr3sx  = good_fft_dimension( nr3s )

     END SUBROUTINE smooth_grid_init

!=----------------------------------------------------------------------------=!
   END MODULE smooth_grid_dimensions
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE grid_subroutines
!=----------------------------------------------------------------------------=!

     ! This module contains subroutines that are related to grids 
     ! parameters

     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: realspace_grids_init, realspace_grids_para

   CONTAINS


     SUBROUTINE realspace_grids_init( b1, b2, b3, gcutd, gcuts, ng, ngs )
       !
       USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x
       USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx
       USE fft_scalar, only: good_fft_dimension, good_fft_order
       USE io_global, only: stdout
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: b1(3), b2(3), b3(3)
       REAL(DP), INTENT(IN) :: gcutd, gcuts
       INTEGER, INTENT(OUT) :: ng, ngs
       !
       REAL(DP) :: qk(3) = 0.0_DP

       IF( nr1 == 0 .OR. nr2 == 0 .OR. nr3 == 0 ) THEN
         ! ... This subroutines calculates the size of the real and reciprocal dense grids
         CALL ngnr_set( b1, b2, b3, gcutd, qk, ng, nr1, nr2, nr3 )
       ELSE
         WRITE( stdout, '( /, 3X,"Info: using nr1, nr2, nr3 values from input" )' )
       END IF

       nr1 = good_fft_order( nr1 )
       nr2 = good_fft_order( nr2 )
       nr3 = good_fft_order( nr3 )

       nr1x  = good_fft_dimension( nr1 )
       nr2x  = nr2
       nr3x  = good_fft_dimension( nr3 )

       IF( nr1s == 0 .OR. nr2s == 0 .OR. nr3s == 0 ) THEN
         ! ... This subroutines calculates the size of the real and reciprocal smoth grids
         CALL ngnr_set( b1, b2, b3, gcuts, qk, ngs, nr1s, nr2s, nr3s )
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

    SUBROUTINE realspace_grids_para( dfftp, dffts )

      !  This subroutines sets local dimensions for real space grids

      USE io_global, ONLY: ionode, stdout
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: nproc_image
      USE fft_types, ONLY: fft_dlay_descriptor
      USE grid_dimensions, ONLY: nr1,  nr2,  nr3, nr1x, nr2x, nr3x
      USE grid_dimensions, ONLY: nr1l, nr2l, nr3l, nrxx
      USE smooth_grid_dimensions, ONLY: nr1s,  nr2s,  nr3s, nr1sx, nr2sx, nr3sx
      USE smooth_grid_dimensions, ONLY: nr1sl, nr2sl, nr3sl, nrxxs

      IMPLICIT NONE

      TYPE(fft_dlay_descriptor), INTENT(IN) :: dfftp, dffts

      INTEGER :: i

      ! ... Subroutine body

      !   set the actual (local) FFT dimensions

      nr1l = dfftp % nr1
      nr2l = dfftp % nr2
      nr3l = dfftp % npl

      nr1sl = dffts % nr1
      nr2sl = dffts % nr2
      nr3sl = dffts % npl

      !   set the dimensions of the array allocated for the FFT
      !   this could in principle be different than the FFT dimensions

      nrxx  = dfftp % nnr
      nrxxs = dffts % nnr

      IF ( nr1s > nr1 .or. nr2s > nr2 .or. nr3s > nr3)                    &
     &   CALL errore(' pmeshset ', ' smooth grid larger than dense grid? ', 1 )

      IF(ionode) THEN

        WRITE( stdout,*)
        WRITE( stdout,*) '  Real Mesh'
        WRITE( stdout,*) '  ---------'
        WRITE( stdout,1000) nr1, nr2, nr3, nr1l, nr2l, nr3l, 1, 1, nproc_image
        WRITE( stdout,1010) nr1x, nr2x, nr3x
        WRITE( stdout,1020) nrxx
        WRITE( stdout,*) '  Number of x-y planes for each processors: '
        WRITE( stdout, fmt = '( 3X, "nr3l = ", 10I5 )' ) ( dfftp%npp( i ), i = 1, nproc_image )

        WRITE( stdout,*)
        WRITE( stdout,*) '  Smooth Real Mesh'
        WRITE( stdout,*) '  ----------------'
        WRITE( stdout,1000) nr1s, nr2s, nr3s, nr1sl, nr2sl, nr3sl, 1, 1, nproc_image
        WRITE( stdout,1010) nr1sx, nr2sx, nr3sx
        WRITE( stdout,1020) nrxxs
        WRITE( stdout,*) '  Number of x-y planes for each processors: '
        WRITE( stdout, fmt = '( 3X, "nr3sl = ", 10I5 )' ) ( dffts%npp( i ), i = 1, nproc_image )

      END IF

1000  FORMAT(3X, &
         'Global Dimensions   Local  Dimensions   Processor Grid',/,3X, &
         '.X.   .Y.   .Z.     .X.   .Y.   .Z.     .X.   .Y.   .Z.',/, &
         3(1X,I5),2X,3(1X,I5),2X,3(1X,I5) )
1010  FORMAT(3X, 'Array leading dimensions ( nr1x, nr2x, nr3x )   = ', 3(1X,I5) )
1020  FORMAT(3X, 'Local number of cell to store the grid ( nrxx ) = ', 1X, I9 )

      RETURN
      END SUBROUTINE realspace_grids_para



   SUBROUTINE ngnr_set( b1, b2, b3, gcut, qk, ng, nr1, nr2, nr3 )

!  this routine calculates the storage required for G vectors arrays
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE kinds, ONLY: DP
      USE mp, ONLY: mp_max, mp_min, mp_sum
      USE mp_global, ONLY: me_image, nproc_image,  intra_image_comm

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: nr1, nr2, nr3, ng
      REAL(DP), INTENT(IN) :: b1(3), b2(3), b3(3), gcut, qk(3)

! ... declare other variables
      INTEGER :: i, j, k
      INTEGER :: nr1tab, nr2tab, nr3tab, nr
      INTEGER :: nb(3)
      REAL(DP) :: gsq, sqgc
      REAL(DP) :: c(3), g(3)
      LOGICAL :: tqk

! ... end of declarations
!  ----------------------------------------------

! ... me_image = processor number, starting from 0

! ... evaluate cutoffs in reciprocal space and the required mesh size
      sqgc  = sqrt(gcut)
      nr     = int(sqgc) + 2      ! nr   = mesh size parameter

! ... verify that, for G<gcut, coordinates never exceed nr
! ... (increase nr if needed)
      CALL vec_prod(c,b1,b2)
      nr3tab=nint(2.0_DP*sqgc/abs(dot_prod(c,b3))*vec_mod(c))
      CALL vec_prod(c,b3,b1)
      nr2tab=nint(2.0_DP*sqgc/abs(dot_prod(c,b2))*vec_mod(c))
      CALL vec_prod(c,b2,b3)
      nr1tab=nint(2.0_DP*sqgc/abs(dot_prod(c,b1))*vec_mod(c))
      nr = max(nr,nr3tab)
      nr = max(nr,nr2tab)
      nr = max(nr,nr1tab)

! ... initialize some variables
      ng     = 0
      nb     = 0

      IF( ALL( qk == 0.0_DP ) ) THEN
        tqk = .FALSE.
      ELSE
        tqk = .TRUE.
      END IF

! *** START OF LOOP ***
! ... calculate moduli of G vectors and the range of indexes where
! ... |G| < gcut 

      DO k = -nr, nr
        IF( MOD( k + nr, nproc_image ) == me_image ) THEN
          DO j = -nr, nr
            DO i = -nr, nr

              g( 1 ) = DBLE(i) * b1(1) + DBLE(j) * b2(1) + DBLE(k) * b3(1)
              g( 2 ) = DBLE(i) * b1(2) + DBLE(j) * b2(2) + DBLE(k) * b3(2)
              g( 3 ) = DBLE(i) * b1(3) + DBLE(j) * b2(3) + DBLE(k) * b3(3)

! ...         calculate modulus
              IF( tqk ) THEN
                gsq = ( g( 1 ) + qk( 1 ) )**2 + &
                    & ( g( 2 ) + qk( 2 ) )**2 + &
                    & ( g( 3 ) + qk( 3 ) )**2 
              ELSE
                gsq =  g( 1 )**2 + g( 2 )**2 + g( 3 )**2 
              END IF

              IF( gsq < gcut ) THEN
! ...           increase counters
                ng  = ng  + 1
! ...           calculate minimum and maximum indexes
                nb(1) = MAX( nb(1), ABS( i ) )
                nb(2) = MAX( nb(2), ABS( j ) )
                nb(3) = MAX( nb(3), ABS( k ) )
              END IF

            END DO
          END DO
        END IF
      END DO

      CALL mp_sum( ng , intra_image_comm )

! ... the size of the required (3-dimensional) matrix depends on the
! ... minimum and maximum indices
      CALL mp_max( nb,  intra_image_comm )

      nr1 = 2 * nb(1) + 1
      nr2 = 2 * nb(2) + 1
      nr3 = 2 * nb(3) + 1

      RETURN

   CONTAINS

!  ----------------------------------------------
!  ----------------------------------------------
      FUNCTION dot_prod(a,b)

!  this function calculates the dot product of two vectors
!  ----------------------------------------------
      
      REAL(DP) :: dot_prod
! ... declare function arguments
      REAL(DP) :: a(3),b(3)     

! ... evaluate dot product
      dot_prod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

      RETURN
      END FUNCTION dot_prod

!  ----------------------------------------------
!  ----------------------------------------------
      FUNCTION vec_mod(a)

!  this function calculates the norm of a vector
!  ----------------------------------------------

      REAL(DP) :: vec_mod
! ... declare function argument
      REAL(DP) :: a(3)

! ... evaluate norm
      vec_mod=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      RETURN
      END FUNCTION vec_mod

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE vec_prod(c,a,b)
   
!  this subroutine calculates the vector (cross) product of vectors
!  a,b and stores the result in vector c
!  ----------------------------------------------

! ... declare subroutine arguments
      REAL(DP) :: a(3),b(3),c(3)

! ... evaluate cross product
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

      RETURN 
      END SUBROUTINE vec_prod

   END  SUBROUTINE ngnr_set

!=----------------------------------------------------------------------------=!
   END MODULE grid_subroutines
!=----------------------------------------------------------------------------=!

