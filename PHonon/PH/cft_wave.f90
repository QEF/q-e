!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
!
!-----------------------------------------------------------------------
subroutine cft_wave_tg (evc_g, evc_r, isw, v_size, ibnd, nbnd_occ)
  !-----------------------------------------------------------------------
  !
  ! Fourier-transformation of a wavefunction using the task group
  ! features
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
  USE mp_bands, ONLY : me_bgrp
  use noncollin_module, ONLY : noncolin, npol
  use qpoint, ONLY : npwq, igkq

  implicit none

  integer, intent(in) :: v_size
  integer, intent(in) :: isw, ibnd, nbnd_occ
  complex(DP), intent(inout) :: evc_g(npwx*npol,nbnd_occ), evc_r(v_size,npol)

  integer :: ig, ioff, idx

  if (isw.eq.1) then
     evc_r = (0.d0, 0.d0)
     !
     ioff   = 0
     !
     DO idx = 1, dffts%nogrp
        !
        IF( idx + ibnd - 1 <= nbnd_occ ) THEN
           DO ig = 1, npw
              evc_r(nls (igk(ig))+ioff,1) =  evc_g(ig,idx+ibnd-1)
           ENDDO
           IF (noncolin) THEN
              DO ig = 1, npw
                 evc_r(nls (igk(ig))+ioff,2) =  evc_g(npwx+ig,idx+ibnd-1)
              ENDDO
           ENDIF
        ENDIF
        !
        ioff = ioff + dffts%tg_nnr
        !
     ENDDO

     CALL invfft ('Wave', evc_r(:,1), dffts)
     IF (noncolin) CALL invfft ('Wave', evc_r(:,2), dffts)

  else if(isw.eq.-1) then

     CALL fwfft ('Wave', evc_r(:,1), dffts)
     IF (noncolin) CALL fwfft ('Wave', evc_r(:,2), dffts)
     !
     ioff   = 0
     !
     DO idx = 1, dffts%nogrp
        !
        IF( idx + ibnd - 1 <= nbnd_occ ) THEN
           !
           DO ig = 1, npwq
              evc_g(ig, ibnd+idx-1) = evc_g(ig, ibnd+idx-1) +  &
                        evc_r( nls(igkq(ig)) + ioff, 1 )
           ENDDO
           !
           IF (noncolin) THEN
              DO ig = 1, npwq
                 evc_g (ig+npwx, ibnd+idx-1) = evc_g (ig+npwx, ibnd+idx-1) &
                      + evc_r (nls(igkq(ig))+ ioff,2)
              ENDDO
           ENDIF
           !
        ENDIF
        !
        ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
        !
     ENDDO
  else
     call errore (' cft_wave_tg',' Wrong switch',1)
  endif

  return
end subroutine cft_wave_tg
