!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "../include/machine.h"

MODULE fft_cp

  USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_allocate, fft_dlay_deallocate, &
                       fft_dlay_set

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
      use work_fft
      use fft_base, only: fft_scatter
!
      implicit none
      integer, intent(in) :: nr1,nr2,nr3,nr1x,nr2x,nr3x,sign
      type (fft_dlay_descriptor), intent(in) :: dfft
      complex(kind=8) :: f( dfft%nnr )

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

      if ( dfft%tptr < 1 .OR. dfft%tptr > 2 ) call errore(' cfft ',' wrong table pointer ', 6)

      me = mpime + 1

      if ( nproc == 1 ) then
         nppx = dfft%nr3x
      else
         nppx = dfft%npp(me)
      end if

      if ( sign > 0 ) then
         if ( sign /= 2 ) then
            if( dfft%tptr == 1 ) then
              call cft_1( f, dfft%nsp(me), nr3, nr3x, sign, aux)
            else if ( dfft%tptr == 2 ) then
              call cft_1s( f, dfft%nsp(me), nr3, nr3x, sign, aux)
            endif
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsp, dfft%npp, sign)
            call zero(2*dfft%nnr,f)
            do i = 1, dfft%nst
               mc = dfft%ismap( i )
               do j = 1, dfft%npp(me)
                  f( mc + (j-1) * dfft%nnp ) = aux( j + (i-1) * nppx)
               end do
            end do
            planes = dfft%iplp
         else
            if( dfft%tptr == 1 ) then
              call cft_1( f, dfft%nsw(me), nr3, nr3x, sign, aux)
            else if ( dfft%tptr == 2 ) then
              call cft_1s( f, dfft%nsw(me), nr3, nr3x, sign, aux)
            endif
            call fft_scatter( aux, nr3x, dfft%nnr, f, dfft%nsw, dfft%npp, sign)
            call zero( 2*dfft%nnr, f )
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
         if( dfft%tptr == 1 ) then
           call cft_2( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign)
         else if ( dfft%tptr == 2 ) then
           call cft_2s( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
         endif
!
      else
!
         if (sign.ne.-2) then
            planes = dfft%iplp
         else
            planes = dfft%iplw
         endif
!
         if( dfft%tptr == 1 ) then
           call cft_2( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign)
         else if ( dfft%tptr == 2 ) then
           call cft_2s( f, dfft%npp(me), nr1, nr2, nr1x, nr2x, sign, planes)
         endif
!
         if (sign.ne.-2) then
            do i=1,dfft%nst
               mc = dfft%ismap( i )
               do j=1,dfft%npp(me)
                  aux(j + (i-1)*nppx) = f(mc+(j-1)*dfft%nnp)
               end do
            end do
            call fft_scatter(aux,nr3x,dfft%nnr,f,dfft%nsp,dfft%npp,sign)
            if( dfft%tptr == 1 ) then
              call cft_1( aux, dfft%nsp(me), nr3, nr3x, sign, f)
            else if ( dfft%tptr == 2 ) then
              call cft_1s( aux, dfft%nsp(me), nr3, nr3x, sign, f)
            endif
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
            if( dfft%tptr == 1 ) then
              call cft_1(aux,dfft%nsw(me),nr3,nr3x,sign,f)
            else if ( dfft%tptr == 2 ) then
              call cft_1s(aux,dfft%nsw(me),nr3,nr3x,sign,f)
            endif
         end if
      end if
!
      return
      end subroutine
!
END MODULE
