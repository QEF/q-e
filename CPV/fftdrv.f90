!
! Copyright (C) 2002-2005 FPMD-CPV groups
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
!

MODULE fft_cp

  USE fft_types, ONLY: fft_dlay_descriptor

  IMPLICIT NONE
  SAVE

  INTEGER, PRIVATE :: what_scatter = 1

CONTAINS

!----------------------------------------------------------------------
   SUBROUTINE cfft_cp ( f, nr1, nr2, nr3, nr1x, nr2x, nr3x, sign, dfft )
!----------------------------------------------------------------------
!
!   sign = +-1 : parallel 3d fft for rho and for the potential
!   sign = +-2 : parallel 3d fft for wavefunctions
!
!   sign = + : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!              fft along z using pencils        (cft_1z)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along y (using planes) and x (cft_2xy)
!   sign = - : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega 
!              fft along x and y(using planes)  (cft_2xy)
!              transpose across nodes           (fft_scatter)
!                 and reorder
!              fft along z using pencils        (cft_1z)
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
      use mp_global, only: mpime, nproc, group
      use fft_scalar, only: cft_1z, cft_2xy
!
      implicit none
      !
      integer, intent(in) :: nr1, nr2, nr3, nr1x, nr2x, nr3x, sign
      type (fft_dlay_descriptor), intent(in) :: dfft
      complex(8) :: f( dfft%nnr )
      complex(8), allocatable :: aux( : )

      integer  mc, i, j, ii, proc, k, me
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

      if ( sign > 0 ) then
         !
         if ( sign /= 2 ) then
            call cft_1z( f, dfft%nsp(me), nr3, nr3x, sign, aux )
            CALL fw_scatter( sign ) ! forwart scatter from stick to planes
            planes = dfft%iplp
         else
            call cft_1z( f, dfft%nsw(me), nr3, nr3x, sign, aux )
            CALL fw_scatter( sign ) ! forwart scatter from stick to planes
            planes = dfft%iplw
         end if
         !
         call cft_2xy( f, dfft%npp( me ), nr1, nr2, nr1x, nr2x, sign, planes )
         !
      else
         !
         if ( sign .ne. -2 ) then
            planes = dfft%iplp
         else
            planes = dfft%iplw
         endif
         !
         call cft_2xy( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
         !
         if ( sign /= -2 ) then
            call bw_scatter( sign )
            call cft_1z( aux, dfft%nsp( me ), nr3, nr3x, sign, f )
         else
            call bw_scatter( sign )
            call cft_1z( aux, dfft%nsw( me ), nr3, nr3x, sign, f )
         end if
         !
      end if
!
      deallocate( aux )

      RETURN
      !
      !
   CONTAINS
      !
      !
      SUBROUTINE fw_scatter( iopt )
         !
         use fft_base, only: fft_scatter, fft_transpose, fft_itranspose
         !
         INTEGER, INTENT(IN) :: iopt
         INTEGER :: nppx
         !
         !
         IF( iopt == 2 ) THEN
            !
            IF( what_scatter == 1 ) THEN
               call fft_transpose ( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, -2)
            ELSE IF( what_scatter == 2 ) THEN
               call fft_itranspose( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, -2)
            ELSE 
               if ( nproc == 1 ) then
                  nppx = dfft%nr3x
               else
                  nppx = dfft%npp( me )
               end if
               call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
               f(:) = (0.d0, 0.d0)
               ii = 0
               do proc = 1, nproc
                  do i = 1, dfft%nsw( proc )
                     mc = dfft%ismap( i + dfft%iss( proc ) )
                     ii = ii + 1 
                     do j = 1, dfft%npp( me )
                        f( mc + ( j - 1 ) * dfft%nnp ) = aux( j + ( ii - 1 ) * nppx )
                     end do
                  end do
               end do
            END IF
            !
         ELSE IF( iopt == 1 ) THEN
            !
            if ( nproc == 1 ) then
               nppx = dfft%nr3x
            else
               nppx = dfft%npp( me )
            end if
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
            f(:) = (0.d0, 0.d0)
            do i = 1, dfft%nst
               mc = dfft%ismap( i )
               do j = 1, dfft%npp( me )
                  f( mc + ( j - 1 ) * dfft%nnp ) = aux( j + ( i - 1 ) * nppx )
               end do
            end do
            !
         END IF
         !
         RETURN 
      END SUBROUTINE fw_scatter
      !
      !  
      !
      SUBROUTINE bw_scatter( iopt )
         !
         use fft_base, only: fft_scatter, fft_transpose, fft_itranspose
         !
         INTEGER, INTENT(IN) :: iopt
         INTEGER :: nppx
         !
         !
         IF( iopt == -2 ) THEN
            !
            IF( what_scatter == 1 ) THEN
               call fft_transpose ( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, 2)
            ELSE IF( what_scatter == 2 ) THEN
               call fft_itranspose( aux, nr3, f, nr1x, nr2x, dfft, me, group, nproc, 2)
            ELSE 
               if ( nproc == 1 ) then
                  nppx = dfft%nr3x
               else
                  nppx = dfft%npp( me )
               end if
               ii = 0
               do proc = 1, nproc
                  do i = 1, dfft%nsw( proc )
                     mc = dfft%ismap( i + dfft%iss( proc ) )
                     ii = ii + 1 
                     do j = 1, dfft%npp( me )
                        aux( j + ( ii - 1 ) * nppx ) = f( mc + ( j - 1 ) * dfft%nnp )
                     end do
                  end do
               end do
               call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
            END IF
            !
         ELSE IF( iopt == -1 ) THEN
            !
            if ( nproc == 1 ) then
               nppx = dfft%nr3x
            else
               nppx = dfft%npp( me )
            end if
            do i = 1, dfft%nst
               mc = dfft%ismap( i )
               do j = 1, dfft%npp( me )
                  aux( j + ( i - 1 ) * nppx ) = f( mc + ( j - 1 ) * dfft%nnp )
               end do
            end do
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
            !
         END IF
         !
         RETURN 
      END SUBROUTINE bw_scatter
      !
      !  
      !
   END SUBROUTINE cfft_cp
   !
   !
   !
END MODULE fft_cp
