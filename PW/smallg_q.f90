!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine smallg_q (xq, iswitch, at, bg, nrot, s, ftau, nr1, nr2, &
 nr3, sym, minus_q)
!-----------------------------------------------------------------------
!
! This routine selects, among the symmetry matrices of the point group
! of a crystal, the symmetry operations which leave q unchanged.
! Furthermore it checks if one of the above matrices send q --> -q+G.
! In this case minus_q is set true.
#include"machine.h"
!
!  input-output variables
!
USE kinds
implicit none

real(kind=DP) :: bg (3, 3), at (3, 3), xq (3)
                        ! input: the reciprocal lattice vectors
                        ! input: the direct lattice vectors
                        ! input: the q point of the crystal

integer :: s (3, 3, 48), nrot, nr1, nr2, nr3, ftau (3, 48), &
 iswitch
                        ! input: the symmetry matrices
                        ! input: number of symmetry operations
                        ! input: fft grid dimension (units for ftau)
                        ! input: fractionary translation of each symmetr
                        ! input: main switch of the program, used for
                        !        q<>0 to restrict the small group of q
                        !        to operation such that Sq=q (exactly,
                        !        without G vectors) when iswitch = -3.
logical :: sym (48), minus_q
                        ! input-output: .true. if symm. op. S q = q + G
                        ! output: .true. if there's an op. sym.: S q = -
!
!  local variables
!

real(kind=DP) :: aq (3), raq (3), zero (3)
                        ! q vector in crystal basis
                        ! the rotated of the q vector
                        ! the zero vector

integer :: irot, ipol, jpol
                       ! counter on symmetry op.
                       ! counter on polarizations
                       ! counter on polarizations

logical :: eqvect
                       ! logical function, check if two vectors are equa
!
! return immediately (with minus_q=.true.) if xq.eq.(0,0,0)
!
minus_q = .true.
if (xq (1) .eq.0.d0.and.xq (2) .eq.0.d0.and.xq (3) .eq.0.d0) &
 return
!
!   Set to zero some variables
!
minus_q = .false.
zero(:) = 0.d0
!
!   Transform xq to the crystal basis
!
call DCOPY (3, xq, 1, aq, 1)
call cryst_to_cart (1, aq, at, - 1)
!
!   Test all symmetries to see if this operation send S q in q+G or in -
!
do irot = 1, nrot
   if (.not.sym (irot) ) goto 100
   raq(:) = 0.d0
   do ipol = 1, 3
      do jpol = 1, 3
         raq(ipol) = raq(ipol) + float( s(ipol,jpol,irot) ) * aq( jpol)
      enddo
   enddo
   sym (irot) = eqvect (raq, aq, zero)
   !
   !  if "iswitch.le.-3" S must be such that Sq=q exactly !
   !
   if (iswitch.le. -3.and.sym (irot) ) then
      do ipol = 1, 3
         sym(irot) = sym(irot) .and. abs(raq(ipol)-aq(ipol)).lt.1.0d-5
      enddo
   endif
   if (sym (irot) .and..not.minus_q) then
   ! l'istruzione "originale" in kreductor era la seguente...
   !         if (.not. minus_q) then
      call DSCAL (3, - 1.d0, raq, 1)
      minus_q = eqvect (raq, aq, zero)
   endif
100 continue
enddo
!
!  if "iswitch.le.-3" time reversal symmetry is not included !
!

if (iswitch.le. -3) minus_q = .false.
return
end subroutine smallg_q

