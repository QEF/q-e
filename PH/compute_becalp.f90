!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine compute_becalp (becq, alpq)
  !---------------------------------------------------------------------
  !
  !     This routine is used only if .not.lgamma and in this case
  !     computes the scalar product of vkb and psi_{k+q}
  !
#include "machine.h"


  use pwcom
  use parameters, only : DP
  use phcom
  implicit none

  complex(kind=DP) :: becq(nkb, nbnd, nksq), alpq(nkb,nbnd,3,nksq)
  ! the becp with psi_{k+q}
  ! the alphap with psi_{k+q}

  integer :: ik, ikq, ipol, ibnd, ig, ios
  ! counter on k points
  ! counter on polarizations, bands and
  ! used for i/o control

  complex(kind=DP) :: fact
  complex(kind=DP), allocatable :: aux (:,:)
  !
  if (lgamma) return

  allocate (aux ( npwx , nbnd))    
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq

     ikq = 2 * ik
     if (nksq.gt.1) then
        read (iunigk, err = 100, iostat = ios) npw, igk

100     call error ('compute_becalp', 'reading igk', abs (ios) )
        read (iunigk, err = 200, iostat = ios) npwq, igkq
200     call error ('compute_becalp', 'reading igkq', abs (ios) )

     endif
     call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     call davcio (evq, lrwfc, iuwfc, ikq, - 1)
     call ccalbec (nkb, npwx, npwq, nbnd, becq (1, 1, ik), vkb, evq)
     do ipol = 1, 3
        do ibnd = 1, nbnd
           do ig = 1, npwq
              aux (ig, ibnd) = evq (ig, ibnd) * &
                   (xk (ipol, ikq) + g (ipol, igkq(ig) ) )
           enddo

        enddo

        call ccalbec (nkb, npwx, npwq, nbnd, alpq(1, 1, ipol, ik),vkb, aux)
     enddo
  enddo
  fact = DCMPLX (0.d0, tpiba)

  call ZSCAL (nkb * nbnd * 3 * nksq, fact, alpq, 1)
  deallocate (aux)
  return
end subroutine compute_becalp
