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
     REAL(dbl), EXTERNAL :: cclock

     INTEGER, PARAMETER :: FFT_MODE_WAVE = 2
     INTEGER, PARAMETER :: FFT_MODE_POTE = 1

     !
     ! ...     Subroutine body 
     !

     IF( ( MODE <= 0 ) .OR. ( MODE > 2 ) ) THEN
         CALL errore( ' PC3FFT_STICK ', ' WRONG MODE ', MODE )
     END IF

     FFT_MODE = MODE

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

     IF ( isign > 0 ) THEN

       !
       ! ...       BACKWARD FFT
       !
       !s1 = cclock()

       CALL cft_1z( c, ns_l, nz, ldz, isign, c )

       !s2 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, -1)
       ELSE
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, -2)
       END IF

       !s3 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplp ) 
       ELSE
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplw ) 
       END IF

       !s4 = cclock()

     ELSE IF( isign < 0 ) THEN
       !
       ! ...       FORWARD FFT
       !
       !s4 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplp ) 
       ELSE
         CALL cft_2xy( r, nz_l, nx, ny, ldx, ldy, isign, dfft%iplw ) 
       END IF

       !s3 = cclock()

       IF( FFT_MODE == FFT_MODE_POTE ) THEN
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, 1)
       ELSE
         CALL fft_transpose(c, ldz, r, ldx, ldy, dfft, (mpime+1), nproc, 2)
       END IF

       !s2 = cclock()

       CALL cft_1z( c, ns_l, nz, ldz, isign, c )

       !s1 = cclock()

     END IF
!
     RETURN

!==---------------------------------------------------------------------==!
   END SUBROUTINE pc3fft_drv
!==---------------------------------------------------------------------==!
!

MODULE fft_cp

  USE fft_types, ONLY: fft_dlay_descriptor

  IMPLICIT NONE
  SAVE

CONTAINS

!----------------------------------------------------------------------
      subroutine cfft_cp (f,nr1,nr2,nr3,nr1x,nr2x,nr3x,sign, dfft)
!----------------------------------------------------------------------
!
!   sign = +-1 : parallel 3d fft for rho and for the potential
!   sign = +-2 : parallel 3d fft for wavefunctions
!
!   sign = + : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!              fft along z using pencils        (cft_1)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along y (using planes) and x (cft_2)
!   sign = - : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega 
!              fft along x and y(using planes)  (cft_2)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along z using pencils        (cft_1)
!
!   The array "planes" signals whether a fft is needed along y :
!     planes(i)=0 : column f(i,*,*) empty , don't do fft along y
!     planes(i)=1 : column f(i,*,*) filled, fft along y needed
!   "empty" = no active components are present in f(i,*,*) 
!             after (sign>0) or before (sign<0) the fft on z direction
!
!   Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
!   and all planes(i) are set to 1
!
!   based on code written by Stefano de Gironcoli for PWSCF
!
      use mp_global, only: mpime, nproc
      use fft_base, only: fft_scatter
      use fft_scalar, only: cft_1z, cft_2xy
!
      implicit none
      integer, intent(in) :: nr1,nr2,nr3,nr1x,nr2x,nr3x,sign
      type (fft_dlay_descriptor), intent(in) :: dfft
      complex(kind=8) :: f( dfft%nnr )
      complex(kind=8), allocatable :: aux( : )

      integer  mc, i, j, ii, proc, k, nppx, me
      integer planes(nr1x)
!
! see comments in cfftp for the logic (or lack of it) of the following
!
      if ( nr1  /= dfft%nr1  ) call errore(' cfft ',' wrong dims ', 1)
      if ( nr2  /= dfft%nr2  ) call errore(' cfft ',' wrong dims ', 2)
      if ( nr3  /= dfft%nr3  ) call errore(' cfft ',' wrong dims ', 3)
      if ( nr1x /= dfft%nr1x ) call errore(' cfft ',' wrong dims ', 4)
      if ( nr2x /= dfft%nr2x ) call errore(' cfft ',' wrong dims ', 5)
      if ( nr3x /= dfft%nr3x ) call errore(' cfft ',' wrong dims ', 6)

      me = mpime + 1

      allocate( aux( dfft%nnr ) )

      if ( nproc == 1 ) then
         nppx = dfft%nr3x
      else
         nppx = dfft%npp(me)
      end if

      if ( sign > 0 ) then
         if ( sign /= 2 ) then
            call cft_1z( f, dfft%nsp(me), nr3, nr3x, sign, aux)
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsp, dfft%npp, sign)
            f(:) = (0.d0, 0.d0)
            do i = 1, dfft%nst
               mc = dfft%ismap( i )
               do j = 1, dfft%npp(me)
                  f( mc + (j-1) * dfft%nnp ) = aux( j + (i-1) * nppx)
               end do
            end do
            planes = dfft%iplp
         else
            call cft_1z( f, dfft%nsw(me), nr3, nr3x, sign, aux)
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsw, dfft%npp, sign)
            f(:) = (0.d0, 0.d0)
            ii = 0
            do proc=1,nproc
               do i=1,dfft%nsw(proc)
                  mc = dfft%ismap( i + dfft%iss(proc) )
                  ii = ii + 1 
                  do j=1,dfft%npp(me)
                     f(mc+(j-1)*dfft%nnp) = aux(j + (ii-1)*nppx)
                  end do
               end do
            end do
            planes = dfft%iplw
         end if
!
         call cft_2xy( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
!
      else
!
         if (sign.ne.-2) then
            planes = dfft%iplp
         else
            planes = dfft%iplw
         endif
!
         call cft_2xy( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
!
         if (sign.ne.-2) then
            do i=1,dfft%nst
               mc = dfft%ismap( i )
               do j=1,dfft%npp(me)
                  aux(j + (i-1)*nppx) = f(mc+(j-1)*dfft%nnp)
               end do
            end do
            call fft_scatter(aux,nr3x,dfft%nnr,f,dfft%nsp,dfft%npp,sign)
            call cft_1z( aux, dfft%nsp(me), nr3, nr3x, sign, f)
         else
            ii = 0
            do proc=1,nproc
               do i=1,dfft%nsw(proc)
                  mc = dfft%ismap( i + dfft%iss(proc) )
                  ii = ii + 1 
                  do j=1,dfft%npp(me)
                     aux(j + (ii-1)*nppx) = f(mc+(j-1)*dfft%nnp)
                  end do
               end do
            end do
            call fft_scatter(aux,nr3x,dfft%nnr,f,dfft%nsw,dfft%npp,sign)
            call cft_1z(aux,dfft%nsw(me),nr3,nr3x,sign,f)
         end if
      end if
!
      deallocate( aux )

      return
      end subroutine

!
END MODULE



!----------------------------------------------------------------------
      subroutine cfftpb(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,sign)
!----------------------------------------------------------------------
!   
!   not-so-parallel 3d fft for box grid, implemented only for sign=1
!   G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!   The array f (overwritten on output) is NOT distributed:
!   a copy is present on each processor.
!   The fft along z  is done on the entire grid. 
!   The fft along xy is done only on planes that have components on the
!   dense grid for each processor. Note that the final array will no
!   longer be the same on all processors.
!     
      use para_mod
      use grid_dimensions, only: nr3
      use fft_scalar, only: cft_b
!     
      implicit none
      integer nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,sign
      complex(kind=8) f(nr1bx*nr2bx*nr3bx)
!     
      integer ir3, ibig3, imin3, imax3, np3
!     
      call parabox(nr3b,irb3,nr3,imin3,imax3)
      np3=imax3-imin3+1
! np3 is the number of planes to be transformed
      if (np3.le.0) return
      call cft_b(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,imin3,imax3,sign)
!     
      return
      end

