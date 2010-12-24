!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine cft_wave (evc_g, evc_r, isw)
  !-----------------------------------------------------------------------
  !
  ! Fourier-transformation of a wavefunction
  ! evc_g(npwx):  the wavefunction in G space
  ! evc_r(nrxxs): the wavefunction in R space ("smooth" grid)
  ! isw =+1: input:  evc_g
  !          output: evc_f = Fourier(evc_g)
  !          evc_g is transformed according to igk-indexes
  !          evc_r is set to zero at the beginning
  ! isw =-1: input:  evc_r
  !          output: evc_g = evc_g + Fourier-1(evc_r)
  !          evc_r is transformed according to igkq indexes
  !

  USE kinds, ONLY : DP
  USE wvfct, ONLY : npwx, npw, igk
  USE fft_base,   ONLY: dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvecs, ONLY : nls
  use noncollin_module, ONLY : noncolin, npol
  use qpoint, ONLY : npwq, igkq
  implicit none

  integer :: isw
  complex(DP) :: evc_g (npwx*npol), evc_r (dffts%nnr,npol)

  integer :: ig

  if (isw.eq.1) then
     evc_r = (0.d0, 0.d0)
     do ig = 1, npw
        evc_r (nls (igk (ig) ),1 ) = evc_g (ig)
     enddo
     CALL invfft ('Wave', evc_r(:,1), dffts)
     IF (noncolin) THEN
        DO ig = 1, npw
           evc_r (nls(igk(ig)),2) = evc_g (ig+npwx)
        ENDDO
        CALL invfft ('Wave', evc_r(:,2), dffts)
     ENDIF
  else if(isw.eq.-1) then
     CALL fwfft ('Wave', evc_r(:,1), dffts)
     do ig = 1, npwq
        evc_g (ig) = evc_g (ig) + evc_r (nls (igkq (ig) ), 1 )
     enddo
     IF (noncolin) THEN
        CALL fwfft ('Wave', evc_r(:,2), dffts)
        DO ig = 1, npwq
           evc_g (ig+npwx) = evc_g (ig+npwx) + evc_r (nls(igkq(ig)),2)
        ENDDO
     ENDIF
  else
     call errore (' cft_wave',' Wrong switch',1)
  endif

  return
end subroutine cft_wave
