!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifndef PARA
!
!----------------------------------------------------------------------

subroutine cfts_3 (f, nr1, nr2, nr3, nrx1, nrx2, nrx3, igrid, &
 sign)
!
!     driver routine for psi 3d complex "reduced" fft
!     sign > 0 : psi(G) => psi(R)   ; sign < 0 : psi(R) => psi(G)
!
!     The 3D fft are computed only on lines and planes which have
!     non zero elements. These lines and planes are defined by
!     the two vectors do_fft_x and do_fft_y which are passed with the
!     common pencil
!
!     The routine is implemented for:
!
!       IBM: essl library
!       NEC: GPFA library
!       DEC: fortran routine
!
!     This routine does not always improve performances
!     It can be disabled by defining NOPENCILS
!
!----------------------------------------------------------------------
!
#include "machine.h"
use parameters
use pencils
!implicit none

integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3, igrid, sign
                        !
                        !   logical dimensions of the fft
                        !
                        !
                        !   physical dimensions of the f array
                        !
                        !   grid used (always the smooth grid)
                        !   sign of the transformation

integer :: j, k!, do_fft_x (nr2, nr3), do_fft_y (nr3)
                                !
                                ! counters on directions
                                ! pencils on the plane
                                ! the planes to transform

#ifdef AIX
complex(kind=DP) :: f (nrx1, nr2, nr3)
                              ! the fft array
!
! ESSL fft's require a different initialization for sign=-1 and sign=1
! aux1 contains the initialized quantities
! aux2 is work space
!
integer :: naux1, naux2
parameter (naux1 = 20000, naux2 = 20000)
integer :: m, incx1, incx2, isign
real(kind=DP) :: aux2 (naux2), scale, sum

external dscal, dcft
integer :: ibid
real(kind=DP), save :: aux1 (naux1, 3, 2)
logical, save :: first (3, 2)
data first / 6 * .true. /


if (igrid.ne.2) call error ('cft_3', 'wrong grid ?', 1)
if (sign.ne. - 1.and.sign.ne.1) call error ('cft_3', 'which fft ?' &
&, 1)
scale = 1.d0
!
! ESSL sign convention for fft's is the opposite of the "usual" one
!
isign = - sign
!
if (isign.eq. - 1) then
   ibid = 1
!
!  i - direction ...
!
   incx1 = 1
   incx2 = nrx1
   m = 1
   if (first (1, ibid) ) then
      call dcft (1, f, incx1, incx2, f, incx1, incx2, nr1, m, &
       isign, scale, aux1 (1, 1, ibid), naux1, aux2, naux2)
      first (1, ibid) = .false.

   endif
   do k = 1, nr3
   do j = 1, nr2
   if (do_fft_x (j, k) .eq.1) call dcft (0, f (1, j, k), incx1, &
    incx2, f (1, j, k), incx1, incx2, nr1, m, isign, scale, aux1 ( &
    1, 1, ibid), naux1, aux2, naux2)
   enddo
   enddo
!
!  ... j-direction ...
!
   incx1 = nrx1
   incx2 = 1
   m = nr1
   if (first (2, ibid) ) then
      call dcft (1, f, incx1, incx2, f, incx1, incx2, nr2, m, &
       isign, scale, aux1 (1, 2, ibid), naux1, aux2, naux2)
      first (2, ibid) = .false.

   endif
   do k = 1, nr3
   if (do_fft_y (k) .eq.1) call dcft (0, f (1, 1, k), incx1, &
    incx2, f (1, 1, k), incx1, incx2, nr2, m, isign, scale, aux1 ( &
    1, 2, ibid), naux1, aux2, naux2)
   enddo
!
!     ... k-direction
!
   incx1 = nrx1 * nr2
   incx2 = 1
   m = nrx1 * nr2
   if (first (3, ibid) ) then
      call dcft (1, f, incx1, incx2, f, incx1, incx2, nr3, m, &
       isign, scale, aux1 (1, 3, ibid), naux1, aux2, naux2)
      first (3, ibid) = .false.
   endif

   call dcft (0, f, incx1, incx2, f, incx1, incx2, nr3, m, isign, &
    scale, aux1 (1, 3, ibid), naux1, aux2, naux2)
else
   ibid = 2
!
!     ... k-direction
!
   incx1 = nrx1 * nr2
   incx2 = 1
   m = nrx1 * nr2
   if (first (3, ibid) ) then
      call dcft (1, f, incx1, incx2, f, incx1, incx2, nr3, m, &
       isign, scale, aux1 (1, 3, ibid), naux1, aux2, naux2)
      first (3, ibid) = .false.
   endif
   call dcft (0, f, incx1, incx2, f, incx1, incx2, nr3, m, isign, &
    scale, aux1 (1, 3, ibid), naux1, aux2, naux2)
!
!     ... j-direction ...
!
   incx1 = nrx1
   incx2 = 1
   m = nr1
   if (first (2, ibid) ) then
      call dcft (1, f, incx1, incx2, f, incx1, incx2, nr2, m, &
       isign, scale, aux1 (1, 2, ibid), naux1, aux2, naux2)
      first (2, ibid) = .false.
   endif
   do k = 1, nr3
   if (do_fft_y (k) .eq.1) call dcft (0, f (1, 1, k), incx1, &
    incx2, f (1, 1, k), incx1, incx2, nr2, m, isign, scale, aux1 ( &
    1, 2, ibid), naux1, aux2, naux2)
   enddo
!
!     i - direction ...
!
   incx1 = 1
   incx2 = nrx1
   m = 1
   if (first (1, ibid) ) then
      call dcft (1, f, incx1, incx2, f, incx1, incx2, nr1, m, &
       isign, scale, aux1 (1, 1, ibid), naux1, aux2, naux2)
      first (1, ibid) = .false.
   endif
   do k = 1, nr3
   do j = 1, nr2
   if (do_fft_x (j, k) .eq.1) call dcft (0, f (1, j, k), incx1, &
    incx2, f (1, j, k), incx1, incx2, nr1, m, isign, scale, aux1 ( &
    1, 1, ibid), naux1, aux2, naux2)
   enddo

   enddo
   call dscal (2 * nrx1 * nr2 * nr3, 1d0 / (nr1 * nr2 * nr3), &
    f, 1)

endif
#endif
#ifdef DEC
#ifndef DXML


real(kind=DP) :: norm, f (2, nr1, nr2, nr3)
                     ! normalization factor
                              ! the function is divided into real and im
if (nr1.ne.nrx1.or.nr2.ne.nrx2) call error ('cft_3', 'no longer im &
&plemented', 1)
!     +    call mcpack(f,nrx1,nrx2,nrx3,f,nr1,nr2,nr3,1)
!
!    Start putting the correct normalization on f and testing the sign
!
norm = dfloat (nr1 * nr2 * nr3)
if (sign.eq. - 1) then
   call dscal (2 * nrx1 * nrx2 * nrx3, 1.d0 / norm, f, 1)
elseif (sign.ne.1) then
   call error ('cfts_3', 'wrong sign', 1)
endif
!
!   This part goes from reciprocal to real space
!   i direction first
!
if (sign.eq.1) then
   do k = 1, nr3
   do j = 1, nr2
   if (do_fft_x (j, k) .eq.1) call cft (f (1, 1, j, k), f (2, 1, &
    j, k), nr1, nr1, nr1, 2 * sign)
   enddo
   enddo
!
!  j direction
!
   do k = 1, nr3
   if (do_fft_y (k) .eq.1) call cft (f (1, 1, 1, k), f (2, 1, 1, &
    k), nr1 * nr2, nr2, nr1 * nr2, 2 * sign)
   enddo
!
!   k direction
!
   call cft (f (1, 1, 1, 1), f (2, 1, 1, 1), nr1 * nr2 * nr3, nr3, &
    nr1 * nr2 * nr3, 2 * sign)
else
!
!   here we compute the transformation from real to reciprocal space
!
!   first the k direction
!
   call cft (f (1, 1, 1, 1), f (2, 1, 1, 1), nr1 * nr2 * nr3, nr3, &
    nr1 * nr2 * nr3, 2 * sign)
!
!   then the j direction
!
   do k = 1, nr3
   if (do_fft_y (k) .eq.1) call cft (f (1, 1, 1, k), f (2, 1, 1, &
    k), nr1 * nr2, nr2, nr1 * nr2, 2 * sign)
   enddo
!
!  and the i direction.
!
   do k = 1, nr3
   do j = 1, nr2
   if (do_fft_x (j, k) .eq.1) call cft (f (1, 1, j, k), f (2, 1, &
    j, k), nr1, nr1, nr1, 2 * sign)
   enddo
   enddo

endif
!      if( nr1.ne.nrx1 .or. nr2.ne.nrx2 )
!     +    call mcpack(f,nrx1,nrx2,nrx3,f,nr1,nr2,nr3,-1)
#endif
#endif
#ifdef NEC

real(kind=DP) :: f (2, nrx1, nrx2, nrx3)
                                    ! inp/out: the function to transform

integer, save :: nr1s, nr2s, nr3s, inc, jump, lot
                          !
                          ! used to check the initialization values
                          !
                          ! the increment between different values
                          ! the jump between fft arrays
                          ! how many fft

real(kind=DP), allocatable, save :: trig1 (2 * nr1), trig2 (2 * nr2), trig3 (2 * nr3), &
 fact
                            !
                            !    trigonometric factors
                            !
                            !    the multiplication factor

logical, save :: first = .true.
                          ! is true at the first iteration

save first, nr1s, nr2s, nr3s, p_trig1, p_trig2, p_trig3
if (igrid.ne.2) call error ('cft_3', 'wrong grid ?', 1)
!
!    test the sign and put the correct normalization on f
!
if (sign.eq. - 1) then
   fact = 1.d0 / float (nr1 * nr2 * nr3)
   call sscal (2 * nrx1 * nrx2 * nrx3, fact, f, 1)
elseif (sign.ne.1) then
   call error ('cft_3', 'wrong sign', 1)
endif
!
!   At the first iteration initialize
!
if (first) then
   nr1s = nr1
   nr2s = nr2
   nr3s = nr3
   allocate (trig1  2 * nr1))    
   allocate (trig2( 2 * nr2))    
   allocate (trig3( 2 * nr3))    
   call setgpfa (trig1, nr1)
   call setgpfa (trig2, nr2)
   call setgpfa (trig3, nr3)
   first = .false.
endif
!
!    test that the initialization value coincide with actual value
!
if (nr1.ne.nr1s.or.nr2.ne.nr2s.or.nr3.ne.nr3s) call error ('cft3', &
 'nr1,nr2,nr3 .ne. from their initializing values', 1)
!
!    Transformation from reciprocal to real space
!
if (sign.eq.1) then
!
!     i-direction
!
   inc = 2
   jump = 2 * nrx1
   lot = nrx2 * nr3
   call gpfa (f (1, 1, 1, 1), f (2, 1, 1, 1), trig1, inc, jump, &
    nr1, lot, sign)
!
!     ... j-direction ...
!
   inc = 2 * nrx1
   jump = 2
   lot = nr1
   do k = 1, nr3
   if (do_fft_y (k) .eq.1) call gpfa (f (1, 1, 1, k), f (2, 1, 1, &
    k), trig2, inc, jump, nr2, lot, sign)
   enddo
!
!     ... k-direction
!
   inc = 2 * nrx1 * nrx2
   jump = 2
   lot = nrx1 * nr2
   call gpfa (f (1, 1, 1, 1), f (2, 1, 1, 1), trig3, inc, jump, &
    nr3, lot, sign)
else
!
!   This is the fft from real to reciprocal space
!
!   k direction
!
   inc = 2 * nrx1 * nrx2
   jump = 2
   lot = nrx1 * nr2
   call gpfa (f (1, 1, 1, 1), f (2, 1, 1, 1), trig3, inc, jump, &
    nr3, lot, sign)
!
!     ... j-direction ...
!
   inc = 2 * nrx1
   jump = 2
   lot = nr1
   do k = 1, nr3
   if (do_fft_y (k) .eq.1) call gpfa (f (1, 1, 1, k), f (2, 1, 1, &
    k), trig2, inc, jump, nr2, lot, sign)
   enddo
!
!    ...i-direction
!
   inc = 2
   jump = 2 * nrx1
   lot = nrx2 * nr3
   call gpfa (f (1, 1, 1, 1), f (2, 1, 1, 1), trig1, inc, jump, &
    nr1, lot, sign)

endif
#endif
return
end subroutine cfts_3
#else
subroutine sbidons ()
return
end subroutine sbidons

#endif
