!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
SUBROUTINE cft_wave (evc_g, evc_r, isw)
  !-----------------------------------------------------------------------
  !
  ! Inverse Fourier (isw=+1) or Fourier (isw=-1) transform of a wavefunction
  !
  ! evc_g(npwx*npol):      wavefunction in G space, ordering: see below
  ! evc_r(dffts%nnr,npol): wavefunction in R space, ordering: same as
  !                        real-space points on the "smooth" FFT grid
  !
  ! isw =+1  input:  evc_g (ordered as \psi(k+G), using igk indices)
  !          output: evc_r = Fourier(evc_g), evc_g = unchanged
  ! isw =-1  input:  evc_r (overwritten on output)
  !          output: evc_g = evc_g + InvFourier(evc_r)
  !                  ordered as \psi(k+q+G), using igkq indices
  !
  ! Further required variables from modules (unchanged on output):
  !
  USE kinds, ONLY : DP
  USE wvfct, ONLY : npwx, npw, igk
  USE fft_base,   ONLY: dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvecs, ONLY : nls
  USE noncollin_module, ONLY : noncolin, npol
  USE qpoint, ONLY : npwq, igkq
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: isw
  COMPLEX(DP) :: evc_g (npwx*npol), evc_r (dffts%nnr,npol)
  ! intent depends upon the value of isw

  INTEGER :: ig

  IF (isw == 1) THEN
     evc_r = (0.0_dp, 0.0_dp)
     DO ig = 1, npw
        evc_r (nls (igk (ig) ),1 ) = evc_g (ig)
     ENDDO
     CALL invfft ('Wave', evc_r(:,1), dffts)
     IF (noncolin) THEN
        DO ig = 1, npw
           evc_r (nls(igk(ig)),2) = evc_g (ig+npwx)
        ENDDO
        CALL invfft ('Wave', evc_r(:,2), dffts)
     ENDIF
  ELSE IF (isw == -1) then
     CALL fwfft ('Wave', evc_r(:,1), dffts)
     DO ig = 1, npwq
        evc_g (ig) = evc_g (ig) + evc_r (nls (igkq (ig) ), 1 )
     ENDDO
     IF (noncolin) THEN
        CALL fwfft ('Wave', evc_r(:,2), dffts)
        DO ig = 1, npwq
           evc_g (ig+npwx) = evc_g (ig+npwx) + evc_r (nls(igkq(ig)),2)
        ENDDO
     ENDIF
  ELSE
     CALL  errore (' cft_wave',' Wrong value for isw',1)
  ENDIF

  RETURN
END SUBROUTINE cft_wave
!
!-----------------------------------------------------------------------
SUBROUTINE cft_wave_tg (evc_g, evc_r, isw, v_size, ibnd, nbnd_occ)
  !-----------------------------------------------------------------------
  !
  !
  ! Inverse Fourier (isw=+1) or Fourier (isw=-1) transform of a wavefunction
  ! using task group parallelization
  !
  ! evc_g(npwx*npol,nbnd_occ): wavefunction in G space, ordering: see below
  ! evc_r(v_size,npol):        wavefunction in R space, ordering: same as
  !                            real-space points on the "smooth" FFT grid
  !                            
  ! isw =+1  input:  evc_g (ordered as \psi(k+G), using igk indices)
  !          output: evc_r = Fourier(evc_g), evc_g = unchanged
  ! isw =-1  input:  evc_r (overwritten on output)
  !          output: evc_g = evc_g + InvFourier(evc_r)
  !                  ordered as \psi(k+q+G), using igkq indices
  ! v_size:   dimension of FFT arrays in real space using task groups
  ! ibnd  :   index of first band in the current task group  
  ! nbnd_occ: number of occupied states
  !
  ! Further required variables from modules (unchanged on output):
  !
  USE kinds, ONLY : DP
  USE wvfct, ONLY : npwx, npw, igk
  USE fft_base,   ONLY: dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvecs, ONLY : nls
  USE mp_bands, ONLY : me_bgrp
  USE noncollin_module, ONLY : noncolin, npol
  USE qpoint, ONLY : npwq, igkq

  IMPLICIT NONE

  INTEGER, INTENT(in) :: v_size, isw, ibnd, nbnd_occ
  COMPLEX(DP), INTENT(inout) :: evc_g(npwx*npol,nbnd_occ), evc_r(v_size,npol)

  INTEGER :: ig, ioff, idx

  IF (isw == 1) then
     evc_r = (0.0_dp, 0.0_dp)
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

  ELSE IF(isw == -1) THEN

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
  ELSE
     CALL  errore (' cft_wave_tg',' Wrong value for isw',1)
  ENDIF

  RETURN
END SUBROUTINE cft_wave_tg
