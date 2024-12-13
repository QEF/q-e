!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!adapted from PH code
!-----------------------------------------------------------------------
SUBROUTINE kcw_psymdvscf (dvtosym)
  !-----------------------------------------------------------------------
  ! Symmetrization of the response charge density.
  !
  USE kinds,     ONLY : DP
  USE noncollin_module, ONLY : nspin_mag
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym
  USE lr_symm_base, ONLY : nsymq, minus_q
  USE control_kcw,  ONLY : irr_bz
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym (dfftp%nnr, nspin_mag)
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, ir3, ioff, ioff_tg, nxyp
  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:)
  !
  if ( .not. irr_bz ) RETURN
  !
  IF (nsymq == 1 .AND. (.NOT.minus_q) ) RETURN
  !
  CALL start_clock ('hp_psymdvscf')
  !
  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin_mag))
  !
  DO is = 1, nspin_mag
     CALL cgather_sym (dfftp, dvtosym(:,is), ddvtosym(:,is))
  ENDDO
  !
  CALL kcw_symdvscf (ddvtosym) !(nper, irr, ddvtosym)
  !
  nxyp = dfftp%nr1x * dfftp%my_nr2p
  DO is = 1, nspin_mag
     DO ir3 = 1, dfftp%my_nr3p
        ioff    = dfftp%nr1x * dfftp%my_nr2p * (ir3-1)
        ioff_tg = dfftp%nr1x * dfftp%nr2x    * (dfftp%my_i0r3p+ir3-1) + dfftp%nr1x * dfftp%my_i0r2p
        CALL zcopy (nxyp, ddvtosym (ioff_tg+1, is), 1, dvtosym (ioff+1, is), 1)
     ENDDO
  ENDDO
  !
  DEALLOCATE (ddvtosym)
  !
  CALL stop_clock ('psymdvscf')
  !
#else
  if ( .not. irr_bz ) RETURN
  CALL kcw_symdvscf (dvtosym)
#endif
  !
  RETURN
  !
END SUBROUTINE kcw_psymdvscf
