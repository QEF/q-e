!
! Copyright (C) 2001 PWSCF group
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
  USE gvect,      ONLY : nrxx, nrx1,nrx2,nrx3
  USE noncollin_module,   ONLY : nspin_mag
  USE modes,     ONLY : minus_q, nsymq
  USE mp_global, ONLY : me_pool
  USE fft_base,  ONLY : dfftp, cgather_sym
  !
  IMPLICIT NONE
  !
  INTEGER :: nper, irr
    ! the number of perturbations
    ! the representation under consideration
  COMPLEX(DP) :: dvtosym (nrxx, nspin_mag, nper)
    ! the potential to symmetrize
    !-local variable
  !
#if defined (__PARA)
  !
  INTEGER :: i, is, iper, npp0

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
  ! the potential to symm


  IF (nsymq.EQ.1.AND. (.NOT.minus_q) ) RETURN
  CALL start_clock ('psym_dmag')

  ALLOCATE (ddvtosym ( nrx1 * nrx2 * nrx3, nspin_mag, nper))    
  npp0 = 1
  DO i = 1, me_pool
     npp0 = npp0 + dfftp%npp (i) * dfftp%nnp

  ENDDO
  DO iper = 1, nper
     DO is = 1, nspin_mag
        CALL cgather_sym (dvtosym (:, is, iper), ddvtosym (:, is, iper) )
     ENDDO

  ENDDO

  CALL sym_dmag (nper, irr, ddvtosym)
  DO iper = 1, nper
     DO is = 1, nspin_mag
        CALL zcopy (dfftp%npp (me_pool+1) * dfftp%nnp, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     ENDDO

  ENDDO
  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psym_dmag')

#endif

  RETURN

END SUBROUTINE psym_dmag
