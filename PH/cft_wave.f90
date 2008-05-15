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

  use pwcom
  use noncollin_module, ONLY : noncolin, npol
  use phcom
  implicit none

  integer :: isw
  complex(DP) :: evc_g (npwx*npol), evc_r (nrxxs,npol)

  integer :: ir, ig

  if (isw.eq.1) then
     evc_r = (0.d0, 0.d0)
     do ig = 1, npw
        evc_r (nls (igk (ig) ),1 ) = evc_g (ig)
     enddo
     call cft3s (evc_r, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
     IF (noncolin) THEN
        DO ig = 1, npw
           evc_r (nls(igk(ig)),2)=evc_g(ig+npwx)
        ENDDO
        CALL cft3s (evc_r(1,2),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,+2)
     ENDIF
  else if(isw.eq.-1) then
     call cft3s (evc_r, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
     do ig = 1, npwq
        evc_g (ig) = evc_g (ig) + evc_r (nls (igkq (ig) ), 1 )
     enddo
     IF (noncolin) THEN
        CALL cft3s(evc_r(1,2),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
        DO ig = 1, npwq
           evc_g(ig+npwx)=evc_g(ig+npwx)+evc_r(nls(igkq(ig)),2)
        ENDDO
     ENDIF
  else
     call errore (' cft_wave',' Wrong switch',1)
  endif

  return
end subroutine cft_wave
