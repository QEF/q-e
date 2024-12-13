!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
SUBROUTINE cft_wave (ik, evc_g, evc_r, isw)
  !-----------------------------------------------------------------------
  !
  ! Inverse Fourier (isw=+1) or Fourier (isw=-1) transform of a wavefunction
  ! Typical usage: compute product Vpsi = (Vq*psi)(k+q+G)
  ! where Vq = perturbation of wave-vector q . Example:
  !    cft_wave (ik, psi, aux,+1)        psi in G space, psi(i) =\psi(k+G_i)
  !    compute  aux(r)= Vq(r)*aux(r)     Vq  in R-space, Vq(i)  =V(r_i)
  !    cft_wave (ik, vpsi, aux,-1)       Vpsi in G-space,Vpsi(i)=Vpsi(k+q+G_i)
  !
  ! Input/Output variables:
  ! ik:                    index of k-point under consideration
  ! evc_g(npwx*npol):      wavefunction in G space, ordering: see below
  ! evc_r(dffts%nnr,npol): wavefunction in R space, ordering: same as
  !                        real-space points on the "smooth" FFT grid
  !
  ! isw =+1  input:  evc_g (ordered as \psi(k+G), using k+G indices)
  !          output: evc_r = Fourier(evc_g), evc_g = unchanged
  ! isw =-1  input:  evc_r (overwritten on output)
  !          output: evc_g = evc_g + InvFourier(evc_r)
  !                  ordered as \psi(k+q+G), using k+q+G indices
  !
  ! ikks(ik) and ikqs(ik) must point in array igk_k(:,i) to the list of k+G
  ! vectors (i=ikks(ik)) and to k+q+G indices (i=ikqs(ik)), respectively 
  !
  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : npwx
  USE qpoint,   ONLY : ikks, ikqs
  USE klist,    ONLY : ngk, igk_k
  USE fft_base, ONLY : dffts
  USE noncollin_module, ONLY : npol
  USE fft_wave,         ONLY : fwfft_wave, invfft_wave

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ik, isw
  COMPLEX(DP) :: evc_g (npwx*npol), evc_r (dffts%nnr,npol)
  ! intent depends upon the value of isw
  !
  INTEGER :: ig, ikk, ikq, npw, npwq

  CALL start_clock ('cft_wave')

  !$acc data present(igk_k)

  IF (isw == 1) THEN
     ikk = ikks(ik) ! points to k+G indices
     npw = ngk(ikk) 
     CALL invfft_wave (npwx, npw, igk_k (1,ikk), evc_g, evc_r )
  ELSE IF (isw == -1) then
     ikq = ikqs(ik) ! points to k+q+G indices
     npwq= ngk(ikq)
     CALL fwfft_wave (npwx, npwq, igk_k (1,ikq), evc_g, evc_r )
  ELSE
     CALL  errore ('cft_wave',' Wrong value for isw',1)
  ENDIF

  !$acc end data
 
  CALL stop_clock ('cft_wave')
 
  RETURN

END SUBROUTINE cft_wave

!-----------------------------------------------------------------------
SUBROUTINE cft_wave_tg (ik, evc_g, evc_r, isw, v_size, ibnd, nbnd_occ)
  !-----------------------------------------------------------------------
  !
  ! Inverse Fourier (isw=+1) or Fourier (isw=-1) transform of a wavefunction
  ! using task group parallelization
  !
  ! ik:                        index of k-point under consideration
  ! evc_g(npwx*npol,nbnd_occ): wavefunction in G space, ordering: see below
  ! evc_r(v_size,npol):        wavefunction in R space, ordering: same as
  !                            real-space points on the "smooth" FFT grid
  ! isw =+1  input:  evc_g (ordered as \psi(k+G), using k+G indices)
  !          output: evc_r = Fourier(evc_g), evc_g = unchanged
  ! isw =-1  input:  evc_r (overwritten on output)
  !          output: evc_g = evc_g + InvFourier(evc_r)
  !                  ordered as \psi(k+q+G), using k+q+G indices
  ! v_size:   dimension of FFT arrays in real space using task groups
  ! ibnd  :   index of first band in the current task group  
  ! nbnd_occ: number of occupied states
  !
  ! Further required variables from modules (must be properly initialized,
  ! unchanged on output):
  !
  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : npwx
  USE fft_base, ONLY : dffts
  USE qpoint,   ONLY : ikks, ikqs
  USE klist,    ONLY : ngk, igk_k
  USE mp_bands, ONLY : me_bgrp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE noncollin_module, ONLY : noncolin, npol
  USE fft_helper_subroutines

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ik, v_size, isw, ibnd, nbnd_occ
  COMPLEX(DP), INTENT(inout) :: evc_g(npwx*npol,nbnd_occ), evc_r(v_size,npol)

  INTEGER :: ig, ikk, ikq, ioff, idx, npw, npwq, ntgrp, right_inc

  CALL start_clock ('cft_wave_tg')

  ntgrp = fftx_ntgrp( dffts )
  CALL tg_get_recip_inc( dffts, right_inc )

  IF (isw == 1) then
     ikk = ikks(ik) ! points to k+G indices
     npw = ngk(ikk)
     evc_r = (0.0_dp, 0.0_dp)
     !
     ioff   = 0
     DO idx = 1, ntgrp
        !
        IF( idx + ibnd - 1 <= nbnd_occ ) THEN
           DO ig = 1, npw
              evc_r(dffts%nl (igk_k(ig,ikk))+ioff,1) =  evc_g(ig,idx+ibnd-1)
           ENDDO
           IF (noncolin) THEN
              DO ig = 1, npw
                 evc_r(dffts%nl (igk_k(ig,ikk))+ioff,2) =  evc_g(npwx+ig,idx+ibnd-1)
              ENDDO
           ENDIF
        ENDIF
        !
        ioff = ioff + right_inc
        !
     ENDDO

     CALL invfft ('tgWave', evc_r(:,1), dffts)
     IF (noncolin) CALL invfft ('tgWave', evc_r(:,2), dffts)

  ELSE IF(isw == -1) THEN
     ikq = ikqs(ik) ! points to k+q+G indices
     npwq= ngk(ikq)
     CALL fwfft ('tgWave', evc_r(:,1), dffts)
     IF (noncolin) CALL fwfft ('tgWave', evc_r(:,2), dffts)
     !
     ioff   = 0
     DO idx = 1, ntgrp
        !
        IF( idx + ibnd - 1 <= nbnd_occ ) THEN
           !
           DO ig = 1, npwq
              evc_g(ig, ibnd+idx-1) = evc_g(ig, ibnd+idx-1) +  &
                        evc_r( dffts%nl(igk_k(ig,ikq)) + ioff, 1 )
           ENDDO
           !
           IF (noncolin) THEN
              DO ig = 1, npwq
                 evc_g (ig+npwx, ibnd+idx-1) = evc_g (ig+npwx, ibnd+idx-1) &
                      + evc_r (dffts%nl(igk_k(ig,ikq))+ ioff,2)
              ENDDO
           ENDIF
           !
        ENDIF
        !
        ioff = ioff + right_inc
        !
     ENDDO
  ELSE
     CALL  errore ('cft_wave_tg',' Wrong value for isw',1)
  ENDIF

  CALL stop_clock ('cft_wave_tg')

  RETURN
END SUBROUTINE cft_wave_tg
