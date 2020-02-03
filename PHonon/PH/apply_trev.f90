!
! Copyright (C) 2018 Andrea Dal Corso and Andrea Urru
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE apply_trev(evc, ikk_evc, ikk_tevc)
!
!  This routine applies the time reversal operator to the wavefunctions
!  evc at the k point ikk_evc and puts the output in evc with the order
!  of G vectors of ikk_tevc
!
USE kinds,     ONLY : DP
USE wvfct,     ONLY : nbnd, npwx
USE klist,     ONLY : ngk, igk_k
USE fft_base,  ONLY: dffts
USE fft_interfaces, ONLY: invfft, fwfft
USE noncollin_module,  ONLY : npol

IMPLICIT NONE
INTEGER :: ikk_evc, ikk_tevc
COMPLEX(DP) :: evc(npwx*npol,nbnd)
COMPLEX(DP), ALLOCATABLE :: aux2(:,:)

INTEGER :: npw, npwt, ibnd, ig


npw=ngk(ikk_evc)
npwt=ngk(ikk_tevc)

ALLOCATE(aux2(dffts%nnr,2))

DO ibnd=1,nbnd
   aux2(:,:) = (0.D0, 0.D0)
   DO ig = 1, npw
      aux2 (dffts%nl (igk_k (ig, ikk_evc) ), 1 ) = evc (ig, ibnd)
      aux2 (dffts%nl (igk_k (ig, ikk_evc) ), 2 ) = evc (npwx+ig, ibnd)
   END DO
   CALL invfft ('Wave', aux2(:,1), dffts)
   CALL invfft ('Wave', aux2(:,2), dffts)
   aux2=CONJG(aux2)
   CALL fwfft ('Wave', aux2(:,1), dffts)
   CALL fwfft ('Wave', aux2(:,2), dffts)
   evc(:,ibnd)=(0.0_DP,0.0_DP)
   DO ig = 1, npwt
      evc(ig,ibnd)=-aux2(dffts%nl(igk_k (ig, ikk_tevc)),2)
      evc(ig+npwx,ibnd)=aux2(dffts%nl(igk_k (ig, ikk_tevc)),1)
   END DO
END DO

DEALLOCATE(aux2)

RETURN
END SUBROUTINE apply_trev
