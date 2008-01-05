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
SUBROUTINE psymdvscf (nper, irr, dvtosym)
  !-----------------------------------------------------------------------
  !
  ! ...  p-symmetrize the charge density.
  !
  USE pwcom
  USE kinds,     ONLY : DP
  USE phcom
  USE mp_global, ONLY : me_pool
  USE fft_base,  ONLY : dfftp
  !
  IMPLICIT NONE
  !
  INTEGER :: nper, irr
    ! the number of perturbations
    ! the representation under consideration
  COMPLEX(DP) :: dvtosym (nrxx, nspin, nper)
    ! the potential to symmetrize
    !-local variable
  !
#if defined (__PARA)
  !
  INTEGER :: i, is, iper, npp0

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
  ! the potential to symm


  IF (nsymq.EQ.1.AND. (.NOT.minus_q) ) RETURN
  CALL start_clock ('psymdvscf')

  ALLOCATE (ddvtosym ( nrx1 * nrx2 * nrx3, nspin, nper))    
  npp0 = 1
  DO i = 1, me_pool
     npp0 = npp0 + dfftp%npp (i) * dfftp%nnp

  ENDDO
  DO iper = 1, nper
     DO is = 1, nspin
        CALL cgather_sym (dvtosym (1, is, iper), ddvtosym (1, is, iper) )
     ENDDO

  ENDDO

  CALL symdvscf (nper, irr, ddvtosym)
  DO iper = 1, nper
     DO is = 1, nspin
        CALL ZCOPY (dfftp%npp (me_pool+1) * dfftp%nnp, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     ENDDO

  ENDDO
  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psymdvscf')

#endif

  RETURN

END SUBROUTINE psymdvscf
