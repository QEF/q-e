!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine compute_becsum
!-----------------------------------------------------------------------
!
!   This routine computes the becsum term which is used to compute the
!   change of the charge due to the displacement of the augmentation
!   term.
!   It implements Eq.12 of the notes.
!
!
#include "machine.h"



use pwcom
USE kinds, only : DP
use phcom
implicit none

integer :: ik, ikk, ikq, ijkb0, ijh, ikb, jkb, ih, jh, na, nt, ibnd
           ! counter on k points, beta functions, atoms and bands
real(kind=DP) :: wgg1 ! auxiliary weight

if (.not.okvan) return
call setv ( (nhm * (nhm + 1) ) / 2 * nat * nspin, 0.d0, becsum, 1)
do ik = 1, nksq
   if (lgamma) then
      ikk = ik
      ikq = ik
   else
      ikk = 2 * ik - 1
      ikq = ikk + 1
   endif
   if (lsda) current_spin = isk (ikk)
   ijkb0 = 0
   do nt = 1, ntyp
      if (tvanp (nt) ) then
         do na = 1, nat
            if (ityp (na) .eq.nt) then
               ijh = 0
               do ih = 1, nh (nt)
                  ikb = ijkb0 + ih
                  ijh = ijh + 1
                  do ibnd = 1, nbnd_occ (ikk)
                     wgg1 = wg (ibnd, ikk)
                     becsum(ijh,na,current_spin) = &
                       becsum(ijh,na,current_spin) + wgg1 * &
                       DREAL ( conjg(becp1(ikb,ibnd,ik)) * becp1(ikb,ibnd,ik) )
                  enddo
                  do jh = 1, nh (nt)
                     jkb = ijkb0 + jh
                     if (jh.gt.ih) ijh = ijh + 1
                     do ibnd = 1, nbnd
                        if (jh.gt.ih) then
                           wgg1 = wg (ibnd, ikk)
                           becsum(ijh,na,current_spin) = &
                             becsum(ijh,na,current_spin) + wgg1 * 2.d0 * &
                             DREAL ( conjg(becp1(ikb,ibnd,ik)) * &
                                           becp1(jkb,ibnd,ik) )
                        endif
                     enddo
                  enddo
               enddo
               ijkb0 = ijkb0 + nh (nt)
            endif
         enddo
      else
         do na = 1, nat
            if (ityp(na).eq.nt) ijkb0 = ijkb0 + nh (nt)
         enddo
      endif
   enddo
enddo
!      do na=1,nat
!         nt=ityp(na)
!         do ijh=1,nh(nt)*(nh(nt)+1)/2
!            WRITE( stdout,'(2i5,f20.10)') na, ijh, becsum(ijh,na,1)
!         enddo
!      enddo
!      call stop_ph(.true.)
return
end subroutine compute_becsum
