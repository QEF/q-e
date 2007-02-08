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
  !     computes the scalar product of vkb and psi_{k+q}, and of
  !     the derivative of vkb and psi_{k+q}
  !
#include "f_defs.h"


  use pwcom
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, only : DP
  USE io_files, ONLY: iunigk
  use phcom
  implicit none

  complex(DP) :: becq(nkb, npol, nbnd, nksq), alpq(nkb,npol,nbnd,3,nksq)
  ! the becp with psi_{k+q}
  ! the alphap with psi_{k+q}

  integer :: ik, ikq, ipol, ibnd, ig, ios
  ! counter on k points
  ! counter on polarizations, bands and
  ! used for i/o control

  complex(DP) :: fact
  complex(DP), allocatable :: aux (:,:)
  !
  if (lgamma) return

  allocate (aux ( npwx*npol , nbnd))    
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq

     ikq = 2 * ik
     if (nksq.gt.1) then
        read (iunigk, err = 100, iostat = ios) npw, igk

100     call errore ('compute_becalp', 'reading igk', abs (ios) )
        read (iunigk, err = 200, iostat = ios) npwq, igkq
200     call errore ('compute_becalp', 'reading igkq', abs (ios) )

     endif
     call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     call davcio (evq, lrwfc, iuwfc, ikq, - 1)
     IF (noncolin) THEN
        call ccalbec_nc(nkb,npwx,npwq,npol,nbnd,becq(1,1,1,ik),vkb,evq)
     ELSE
        call ccalbec (nkb, npwx, npwq, nbnd, becq(1, 1, 1,ik), vkb, evq)
     END IF
     do ipol = 1, 3
        aux=(0.d0,0.d0)
        do ibnd = 1, nbnd
           do ig = 1, npwq
              aux (ig, ibnd) = evq (ig, ibnd) * &
                   (xk (ipol, ikq) + g (ipol, igkq(ig) ) )
           enddo
           IF (noncolin) THEN
              do ig = 1, npwq
                 aux (ig+npwx, ibnd) = evq (ig+npwx, ibnd) * &
                   (xk (ipol, ikq) + g (ipol, igkq(ig) ) )
              enddo
           ENDIF
        enddo

        IF (noncolin) THEN
           call ccalbec_nc(nkb,npwx,npwq,npol,nbnd,alpq(1,1,1,ipol,ik),vkb, aux)
        ELSE
           call ccalbec (nkb, npwx, npwq, nbnd, alpq(1,1,1,ipol,ik),vkb, aux)
        END IF
     enddo
  enddo
  fact = CMPLX (0.d0, tpiba)

  call ZSCAL (nkb * nbnd * 3 * nksq * npol, fact, alpq, 1)
  deallocate (aux)
  return
end subroutine compute_becalp
