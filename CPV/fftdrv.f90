!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  ----------------------------------------------
!  This Module written by Carlo Cavazzoni 
!  Last modified April 2003
!  ----------------------------------------------

#include "f_defs.h"

!=---------------------------------------------------------------------==!
!
!
!     FFT high level Driver 
!     ( Charge density and Wave Functions )
!
!
!=---------------------------------------------------------------------==!


!==---------------------------------------------------------------------==!
     SUBROUTINE pc3fft_drv(r, c, isign, dfft, mode)
!==---------------------------------------------------------------------==!

! ...
     USE fft_base, ONLY: fft_transpose
     USE fft_scalar, ONLY: cft_1z, cft_2xy
     USE mp_global, ONLY: mpime, nproc
     USE fft, ONLY: FFT_MODE_WAVE, FFT_MODE_POTE, tims
     USE stick, ONLY: dfftp
     USE fft_types, ONLY: fft_dlay_descriptor
     USE kinds, ONLY: dbl

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: isign
     TYPE (fft_dlay_descriptor), INTENT(IN) ::  dfft
     INTEGER, INTENT(IN) :: mode
     COMPLEX (dbl) :: R( dfft%nnr )
     COMPLEX (dbl) :: C( dfft%nr3x * dfft%nst )

!    R( * )   3D real space grid, 3rd dimension is 
!             distributed among processors (nproc).
!             The size of the local block of R is specified
!             by the structure dfft .
!
!    C( * )   3D reciprocal space grid stored as an array of
!             sticks ( "z" columns, one stick for each "x" and "y" 
!             coordinates ). The sticks are distributed among processors. 
!             The distribution and the size of the local block
!             of C are specified in the structure "dfft".
!
!    ISIGN    FFT direction and/or FFT initialization
!             > 0  backward direction  G-space to R-space, 
!                  output = \sum_G f(G)exp(+iG*R)
!             = 0  initialization
!             < 0  forward direction   R-space to G-space, 
!                  output = \int_R f(R)exp(-iG*R)/Omega
!
!    MODE     FFT_MODE_POTE  
!               ( potential mode, use the full G-vec sphere )
!             FFT_MODE_WAVE  
!               ( wave func mode, use the small G-vec sphere )
!
! ...

     INTEGER, SAVE :: nx, ny, nz, nz_l, ns_l, ldx, ldy, ldz
     LOGICAL, SAVE :: reinit
     INTEGER, SAVE :: FFT_MODE = 0
     INTEGER :: ierr
     REAL(dbl) :: s1, s2, s3, s4, s5
     REAL(dbl) :: cclock
     EXTERNAL  :: cclock

!
! ...     Subroutine body 
!

     IF ( ( isign == 0 ) .OR. ( MODE /= FFT_MODE ) ) THEN

       s1 = 0.0d0; s2 = 0.0d0; s3 = 0.0d0; s4 = 0.0d0; s5 = 0.0d0;

       IF( ( MODE <= 0 ) .OR. ( MODE > 2 ) ) THEN
         CALL errore( ' PC3FFT_STICK ', ' WRONG MODE ', MODE )
       ELSE IF( MODE == FFT_MODE ) THEN
         reinit = .FALSE.
       ELSE
         reinit = .TRUE.
       END IF

       FFT_MODE = MODE

       IF( reinit ) THEN

         nx   = dfft%nr1
         ny   = dfft%nr2
         nz   = dfft%nr3

         ldx  = dfft%nr1x
         ldy  = dfft%nr2x
         ldz  = dfft%nr3x

         nz_l = dfft%npp( mpime + 1 )

         IF( FFT_MODE == FFT_MODE_POTE ) THEN
           ns_l = dfft%nsp( mpime + 1 )
         ELSE
           ns_l = dfft%nsw( mpime + 1 )
         END IF

       END IF
            
     END IF

     IF ( isign > 0 ) THEN
!
! ...       BACKWARD FFT
!
       s1 = cclock()

       CALL cft_1z( c, ns_l, nz, ldz, isign, c )

       s2 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, -1)
       ELSE
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, -2)
       END IF

       s3 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplp ) 
       ELSE
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplw ) 
       END IF

       s4 = cclock()

     ELSE IF( isign < 0 ) THEN
!
! ...       FORWARD FFT
!
       s4 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplp ) 
       ELSE
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplw ) 
       END IF

       s3 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, 1)
       ELSE
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, 2)
       END IF

       s2 = cclock()

       CALL cft_1z( c, ns_l, nz, ldz, isign, c )

       s1 = cclock()

     END IF

     tims( 2, FFT_MODE ) = tims( 2, FFT_MODE ) + ABS(s4-s3)
     tims( 3, FFT_MODE ) = tims( 3, FFT_MODE ) + ABS(s2-s1)
     tims( 4, FFT_MODE ) = tims( 4, FFT_MODE ) + ABS(s3-s2)
!
     RETURN

!==---------------------------------------------------------------------==!
   END SUBROUTINE pc3fft_drv
!==---------------------------------------------------------------------==!
