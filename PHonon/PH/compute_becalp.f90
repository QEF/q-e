!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine compute_becalp (becq, alpq)
  !---------------------------------------------------------------------
  !
  !     This routine is used only at finite q and in this case
  !     computes the scalar product of vkb and psi_{k+q}, and of
  !     the derivative of vkb and psi_{k+q}. Eq. B8 and B10 (at k+q)
  !     of PRB 64 235118 (2001).
  !

  USE kinds, only : DP
  USE cell_base, ONLY : tpiba
  USE klist,     ONLY : xk
  USE gvect,     ONLY : g
  USE becmod, ONLY: calbec, bec_type, becscal
  USE buffers, ONLY: get_buffer
  USE uspp, ONLY: nkb, vkb
  USE noncollin_module, ONLY : noncolin, npol
  USE io_files, ONLY: iunigk
  USE wvfct,    ONLY : nbnd, npw, npwx, igk
  USE paw_variables, ONLY : okpaw

  USE units_ph, ONLY : lrwfc, iuwfc
  USE control_ph, ONLY : rec_code_read
  USE control_lr, ONLY : lgamma
  USE eqv, ONLY : evq
  USE qpoint, ONLY : nksq, npwq, igkq, ikqs

  implicit none

  type (bec_type) :: becq(nksq), alpq(3,nksq)
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
  IF (rec_code_read >= -20.AND..NOT.okpaw) RETURN

  allocate (aux ( npwx*npol , nbnd))
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq

     ikq = ikqs(ik)
     if (nksq.gt.1) then
        read (iunigk, err = 100, iostat = ios) npw, igk

100     call errore ('compute_becalp', 'reading igk', abs (ios) )
        read (iunigk, err = 200, iostat = ios) npwq, igkq
200     call errore ('compute_becalp', 'reading igkq', abs (ios) )

     endif
     call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     call get_buffer (evq, lrwfc, iuwfc, ikq)
     call calbec ( npwq, vkb, evq, becq(ik) )
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

        call calbec ( npwq, vkb, aux, alpq(ipol,ik) )
     enddo
  enddo
  fact = CMPLX(0.d0, tpiba,kind=DP)

  DO ik=1,nksq
     DO ipol=1,3
        CALL becscal(fact,alpq(ipol,ik),nkb,nbnd)
     ENDDO
  ENDDO

  deallocate (aux)
  return
end subroutine compute_becalp
