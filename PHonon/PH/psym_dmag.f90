!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE psym_dmag (nper, irr, dvtosym)
  !-----------------------------------------------------------------------
  !
  ! ...  p-symmetrize the charge density.
  !
  USE kinds,      ONLY : DP
  USE noncollin_module, ONLY : nspin_mag
  USE mp_bands,  ONLY : me_bgrp
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym

  USE lr_symm_base, ONLY : minus_q, nsymq
  !
  IMPLICIT NONE
  !
  INTEGER :: nper, irr
    ! the number of perturbations
    ! the representation under consideration
  COMPLEX(DP) :: dvtosym (dfftp%nnr, nspin_mag, nper)
    ! the potential to symmetrize
    !-local variable
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, npp0

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
  ! the potential to symm


  IF (nsymq.EQ.1.AND. (.NOT.minus_q) ) RETURN
  CALL start_clock ('psym_dmag')

  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin_mag, nper))
  npp0 = 1
  DO i = 1, me_bgrp
     npp0 = npp0 + dfftp%npp (i) * dfftp%nnp

  ENDDO
  DO iper = 1, nper
     DO is = 1, nspin_mag
        CALL cgather_sym (dfftp,dvtosym (:, is, iper), ddvtosym (:, is, iper) )
     ENDDO

  ENDDO

  CALL sym_dmag (nper, irr, ddvtosym)
  DO iper = 1, nper
     DO is = 1, nspin_mag
        CALL zcopy (dfftp%npp (me_bgrp+1) * dfftp%nnp, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     ENDDO

  ENDDO
  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psym_dmag')
#else
  CALL sym_dmag (nper, irr, dvtosym)
#endif

  RETURN

END SUBROUTINE psym_dmag
