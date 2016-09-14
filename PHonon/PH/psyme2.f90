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
  integer :: i, iper, npp0

  complex(DP), allocatable :: ddvtosym (:,:)
  ! the potential to symmetrize

  allocate (ddvtosym (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 6))

  npp0 = 0
  do i = 1, me_bgrp
     npp0 = npp0 + dfftp%npp (i)
  enddo
  npp0 = npp0 * dfftp%nnp + 1
  do iper = 1, 6
     call cgather_sym (dfftp, dvtosym (:, iper), ddvtosym (:, iper) )
  enddo

  call syme2 (ddvtosym)

  do iper = 1, 6
     call zcopy (dfftp%npp (me_bgrp+1) * dfftp%nnp, ddvtosym (npp0, iper), 1, &
                 dvtosym (1, iper), 1)
  enddo

  deallocate (ddvtosym)
#else
  call syme2 (dvtosym)
#endif
  return
end subroutine psyme2
