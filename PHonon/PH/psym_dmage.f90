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
  INTEGER :: i, is, iper, ir3, ioff, ioff_tg, nxyp

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
  ! the potential to symm

  CALL start_clock ('psym_dmage')

  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin, 3))

  DO iper = 1, 3
     DO is = 1, nspin
        CALL cgather_sym (dfftp,dvtosym (:, is, iper), ddvtosym (:, is, iper) )
     ENDDO
  ENDDO

  CALL sym_dmage (ddvtosym)

  nxyp = dfftp%nr1x * dfftp%my_nr2p
  DO iper = 1, 3
     DO is = 1, nspin
        DO ir3 = 1, dfftp%my_nr3p
           ioff    = dfftp%nr1x * dfftp%my_nr2p * (ir3-1)
           ioff_tg = dfftp%nr1x * dfftp%nr2x    * (dfftp%my_i0r3p+ir3-1) + dfftp%nr1x * dfftp%my_i0r2p
           CALL zcopy (nxyp, ddvtosym (ioff_tg+1, is, iper), 1, dvtosym (ioff+1, is, iper), 1)
        END DO
     ENDDO
  ENDDO

  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psym_dmage')

#else

  CALL sym_dmage (dvtosym)

#endif

  RETURN

END SUBROUTINE psym_dmage
