!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine compute_alphasum
  !-----------------------------------------------------------------------
  !
  !   This routine computes the alphasum term which is used to compute the
  !   change of the charge due to the displacement of the augmentation
  !   term. (See Eq. 29)
  !   It implements Eq.13 of the notes.
  !
  !
#include "machine.h"

  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE kinds, only : DP
  USE uspp_param, ONLY: nh, tvanp
  use phcom
  implicit none

  integer :: ik, ikk, ikq, ijkb0, ijh, ikb, jkb, ih, jh, na, nt, &
       ipol, ibnd
  ! counter on k points
  ! counters on beta functions
  ! counters on beta functions
  ! counters for atoms
  ! counter on polarizations
  ! counter on bands

  real(kind=DP) :: wgg1
  ! auxiliary weight


  if (.not.okvan) return

  alphasum (:,:,:,:) = 0.d0
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
              if (ityp (na) == nt) then
                 ijh = 0
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    ijh = ijh + 1
                    do ibnd = 1, nbnd_occ (ikk)
                       wgg1 = wg (ibnd, ikk)
                       do ipol = 1, 3
                          alphasum(ijh,ipol,na,current_spin) = &
                               alphasum(ijh,ipol,na,current_spin) + 2.d0*wgg1*&
                               DREAL (conjg (alphap (ikb,ibnd,ipol,ik) ) * &
                               becp1  (ikb,ibnd,ik) )
                       enddo
                    enddo
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       if (jh > ih) ijh = ijh + 1
                       do ibnd = 1, nbnd
                          if (jh > ih) then
                             wgg1 = wg (ibnd, ikk)
                             do ipol = 1, 3
                                alphasum(ijh,ipol,na,current_spin) = &
                                     alphasum(ijh,ipol,na,current_spin) + &
                                     2.d0 * wgg1 * &
                                     DREAL (conjg (alphap(ikb,ibnd,ipol,ik) )*&
                                     becp1 (jkb,ibnd,ik)         + &
                                     conjg ( becp1 (ikb,ibnd,ik) ) *       &
                                     alphap (jkb,ibnd,ipol,ik) )
                             enddo
                          endif
                       enddo
                    enddo
                 enddo
                 ijkb0 = ijkb0 + nh (nt)
              endif
           enddo
        else
           do na = 1, nat
              if (ityp (na) == nt) ijkb0 = ijkb0 + nh (nt)
           enddo
        endif
     enddo
  enddo
  !      do na=1,nat
  !         nt=ityp(na)
  !         do ijh=1,nh(nt)*(nh(nt)+1)/2
  !            do ipol=1,3
  !               WRITE( stdout,'(3i5,f20.10)') na, ijh, ipol,
  !     +                              alphasum(ijh,ipol,na,1)
  !            enddo
  !         enddo
  !      enddo
  !      call stop_ph(.true.)
  return
end subroutine compute_alphasum
