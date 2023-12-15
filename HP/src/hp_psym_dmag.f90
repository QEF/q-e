!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_psym_dmag (dvtosym)
    !-----------------------------------------------------------------------
    !
    ! ...  p-symmetrize the charge density.
    !
    USE kinds,      ONLY : DP
    USE noncollin_module, ONLY : nspin_mag
    USE fft_base,  ONLY : dfftp
    USE scatter_mod,  ONLY : cgather_sym
    USE lr_symm_base, ONLY : minus_q, nsymq
    !
    IMPLICIT NONE
    !
    COMPLEX(DP) :: dvtosym (dfftp%nnr, nspin_mag)
      ! the potential to symmetrize
      !-local variable
    !
#if defined (__MPI)
    !
    INTEGER :: i, is, ir3, ioff, ioff_tg, nxyp
  
    COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:)
    ! the potential to symm
  
  
    IF (nsymq.EQ.1.AND. (.NOT.minus_q) ) RETURN
    CALL start_clock ('hp_psym_dmag')
  
    ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin_mag))
  
    DO is = 1, nspin_mag
       CALL cgather_sym (dfftp,dvtosym (:, is), ddvtosym (:, is) )
    ENDDO
  
    CALL hp_sym_dmag (ddvtosym)
  
    nxyp = dfftp%nr1x * dfftp%my_nr2p
    DO is = 1, nspin_mag
      DO ir3 = 1, dfftp%my_nr3p
        ioff    = dfftp%nr1x * dfftp%my_nr2p * (ir3-1)
        ioff_tg = dfftp%nr1x * dfftp%nr2x    * (dfftp%my_i0r3p+ir3-1) + dfftp%nr1x * dfftp%my_i0r2p
        CALL zcopy (nxyp, ddvtosym (ioff_tg+1, is), 1, dvtosym (ioff+1, is), 1)
       END DO
    ENDDO
  
    DEALLOCATE (ddvtosym)
  
    CALL stop_clock ('hp_psym_dmag')
#else
    CALL hp_sym_dmag (dvtosym)
#endif
  
    RETURN
  
  END SUBROUTINE hp_psym_dmag
  