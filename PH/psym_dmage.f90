!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE psym_dmage (dvtosym)
  !-----------------------------------------------------------------------
  !
  ! ...  p-symmetrize the magnetization change due to an electric field.
  !
  USE pwcom
  USE kinds,     ONLY : DP
  USE phcom
  USE mp_global, ONLY : me_pool
  USE pfft,      ONLY : npp, ncplane
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym (nrxx, nspin, 3)
    ! the potential to symmetrize
    !-local variable
  !
#if defined (__PARA)
  !
  INTEGER :: i, is, iper, npp0

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
  ! the potential to symm

  CALL start_clock ('psym_dmage')

  ALLOCATE (ddvtosym ( nrx1 * nrx2 * nrx3, nspin, 3))    
  npp0 = 1
  DO i = 1, me_pool
     npp0 = npp0 + npp (i) * ncplane
  ENDDO
  DO iper = 1, 3
     DO is = 1, nspin
        CALL cgather_sym (dvtosym (1, is, iper), ddvtosym (1, is, iper) )
     ENDDO
  ENDDO

  CALL sym_dmage (ddvtosym)
  DO iper = 1, 3
     DO is = 1, nspin
        CALL ZCOPY (npp (me_pool+1) * ncplane, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     ENDDO
  ENDDO
  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psym_dmage')

#endif

  RETURN

END SUBROUTINE psym_dmage
