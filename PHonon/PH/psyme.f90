!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE psyme (dvtosym)
  !-----------------------------------------------------------------------
  !! p-symmetrize the charge density.
  !
  USE kinds,     ONLY : DP
  USE fft_base, ONLY : dfftp
  USE noncollin_module, ONLY : nspin_mag
  USE mp_bands, ONLY : me_bgrp
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym (dfftp%nnr, nspin_mag, 3)
  !! the potential to symmetrize
  !
  ! ... local variables
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, ir3, ioff, ioff_tg, nxyp
  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
    ! the potential to symmet
  !
  !
  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin_mag, 3))

  DO iper = 1, 3
     DO is = 1, nspin_mag
        CALL cgather_sym (dfftp,dvtosym (:, is, iper), ddvtosym (:, is, iper) )
     ENDDO

  ENDDO

  CALL syme (ddvtosym)

  nxyp = dfftp%nr1x * dfftp%my_nr2p
  DO iper = 1, 3
     DO is = 1, nspin_mag
        DO ir3 = 1, dfftp%my_nr3p
           ioff    = dfftp%nr1x * dfftp%my_nr2p * (ir3-1)
           ioff_tg = dfftp%nr1x * dfftp%nr2x    * (dfftp%my_i0r3p+ir3-1) + dfftp%nr1x * dfftp%my_i0r2p
           CALL zcopy (nxyp, ddvtosym (ioff_tg+1, is, iper), 1, dvtosym (ioff+1, is, iper), 1)
        END DO
     ENDDO
  ENDDO

  DEALLOCATE (ddvtosym)

#else
  CALL syme (dvtosym)
#endif

  RETURN

END SUBROUTINE psyme
