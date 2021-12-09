!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dveqpsi_us (ik)
  !----------------------------------------------------------------------
  !
  ! This routine applies e^iqr to a wavefunction.
  ! In the US and PAW case applies e^iqr K(r,r',r'')
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, ityp
  USE cell_base, ONLY : tpiba
  USE fft_base,  ONLY: dfftp, dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE noncollin_module, ONLY : noncolin
  USE wvfct,     ONLY : nbnd, npwx

  USE wavefunctions,  ONLY: evc
  USE control_lr, ONLY : nbnd_occ
  USE eqv,       ONLY : dvpsi
  USE qpoint,    ONLY : ikks, ikqs
  USE klist,     ONLY : ngk, igk_k 
  IMPLICIT NONE
  !
  !   The dummy variables
  !

  INTEGER :: ik, npw, npwq, ikk, ikq
  ! input: the k point

  INTEGER :: ibnd, ig
  ! counter on bands
  ! counter on G vectors
  ! counter on polarizations
  ! counter on reciprocal mesh

  COMPLEX(DP) , ALLOCATABLE ::  aux2 (:)
  ! work space

  CALL start_clock ('dveqpsi_us')
  ALLOCATE (aux2(dffts%nnr))
  !
  !  The unperturbed wavefunctions must be multiplied by e^{iqr}.
  !  This means that we have to order the coefficients with the mesh
  !  of k+q. We do this by distributing the vectors on the FFT mesh
  !  with the indices of k+G (igk) and then saving them with the mesh 
  !  of k+q+G (igkq)
  !
  dvpsi(:,:) = (0.D0, 0.D0)

  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)


  DO ibnd = 1, nbnd_occ(ikk)
     aux2(:) = (0.D0, 0.D0)
     DO ig = 1, npw
        aux2 (dffts%nl (igk_k (ig,ikk) ) ) = evc (ig, ibnd)
     ENDDO

     DO ig = 1, npwq
        dvpsi (ig, ibnd) = aux2 (dffts%nl (igk_k (ig,ikq) ) )
     ENDDO
     IF (noncolin) THEN
        aux2(:) = (0.d0, 0.d0)
        DO ig = 1, npw
           aux2 (dffts%nl (igk_k (ig,ikk) ) ) = evc (ig+npwx, ibnd)
        ENDDO

        DO ig = 1, npwq
           dvpsi (ig+npwx, ibnd) = aux2 (dffts%nl (igk_k (ig,ikq) ) )
        ENDDO
     END IF
  ENDDO
  !
  DEALLOCATE (aux2)
  !
  !   We add the contribution of the nonlocal potential in the US form
  !
  CALL dveqpsi_us_only (npwq, ik)

  CALL stop_clock ('dveqpsi_us')
  RETURN
END SUBROUTINE dveqpsi_us
