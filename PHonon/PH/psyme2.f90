!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine psyme2 (dvtosym)
  !-----------------------------------------------------------------------
  !  p-symmetrize the second derivative of charge density.
  !

  use kinds, only : DP
  USE mp_bands, ONLY: me_bgrp
  USE fft_base,  ONLY: dfftp
  USE scatter_mod,  ONLY: cgather_sym
  implicit none

  complex(DP) :: dvtosym (dfftp%nnr, 6)
  ! the potential to symmetrize
  !-local variable

#if defined(__MPI)
  integer :: i, iper, ir3, ioff, ioff_tg, nxyp

  complex(DP), allocatable :: ddvtosym (:,:)
  ! the potential to symmetrize

  allocate (ddvtosym (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 6))

  do iper = 1, 6
     call cgather_sym (dfftp, dvtosym (:, iper), ddvtosym (:, iper) )
  enddo

  call syme2 (ddvtosym)

  nxyp = dfftp%nr1x * dfftp%my_nr2p
  DO iper = 1, 6
     DO ir3 = 1, dfftp%my_nr3p
        ioff    = dfftp%nr1x * dfftp%my_nr2p * (ir3-1)
        ioff_tg = dfftp%nr1x * dfftp%nr2x    * (dfftp%my_i0r3p+ir3-1) + dfftp%nr1x * dfftp%my_i0r2p
        CALL zcopy (nxyp, ddvtosym (ioff_tg+1, iper), 1, dvtosym (ioff+1, iper), 1)
     END DO
  ENDDO


  deallocate (ddvtosym)
#else
  call syme2 (dvtosym)
#endif
  return
end subroutine psyme2
