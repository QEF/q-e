!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE psym_dmage (dvtosym)
  !-----------------------------------------------------------------------
  !
  ! ...  p-symmetrize the magnetization change due to an electric field.
  !
  USE kinds,     ONLY : DP
  USE lsda_mod,   ONLY : nspin
  USE mp_bands,  ONLY : me_bgrp
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym (dfftp%nnr, nspin, 3)
    ! the potential to symmetrize
    !-local variable
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, npp0

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
  ! the potential to symm

  CALL start_clock ('psym_dmage')

  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin, 3))
  npp0 = 1
  DO i = 1, me_bgrp
     npp0 = npp0 + dfftp%npp (i) * dfftp%nnp
  ENDDO
  DO iper = 1, 3
     DO is = 1, nspin
        CALL cgather_sym (dfftp,dvtosym (:, is, iper), ddvtosym (:, is, iper) )
     ENDDO
  ENDDO

  CALL sym_dmage (ddvtosym)
  DO iper = 1, 3
     DO is = 1, nspin
        CALL zcopy (dfftp%npp (me_bgrp+1) * dfftp%nnp, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     ENDDO
  ENDDO
  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psym_dmage')

#else

  CALL sym_dmage (dvtosym)

#endif

  RETURN

END SUBROUTINE psym_dmage
