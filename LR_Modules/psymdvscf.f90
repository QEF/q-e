!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE psymdvscf (dvtosym, dfft)
  !-----------------------------------------------------------------------
  !! p-symmetrize the potential.
  !!
  !! The real space points of dv is distributed, but symmetry may map a point in one
  !! core to a point in a different core. Hence, gather dvscf in real space, symmetrize,
  !! and then scatter back.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : fft_type_descriptor
  USE noncollin_module, ONLY : nspin_mag, noncolin, domag
  USE scatter_mod,      ONLY : cgather_sym
  USE lr_symm_base,     ONLY : nsymq, minus_q, lr_npert
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  !! FFT descriptor for dvtosym (dfftp for hard grid, dffts for soft grid)
  COMPLEX(DP), INTENT(INOUT) :: dvtosym(dfft%nnr, nspin_mag, lr_npert)
  !! the potential or density to symmetrize
  !
  ! ... local variable
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, ir3, ioff, ioff_tg, nxyp
  !
  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
  ! the potential to symmetrize, in gathered form
  IF (nsymq == 1 .AND. (.NOT.minus_q) ) RETURN
  CALL start_clock ('psymdvscf')
  !
  ALLOCATE (ddvtosym ( dfft%nr1x * dfft%nr2x * dfft%nr3x, nspin_mag, lr_npert))
  !
  ! Gather real-space points
  !
  DO iper = 1, lr_npert
     DO is = 1, nspin_mag
        CALL cgather_sym (dfft, dvtosym (:, is, iper), ddvtosym (:, is, iper) )
     ENDDO
  ENDDO
  !
  ! Symmetrize
  !
  ! Nonmagnetic : nspin_mag = 1, nspin_lsda = 1.
  !               symdvscf symmetrizes this component.
  ! Collinear magnet (LSDA) : nspin_mag = 2, nspin_lsda = 2.
  !                           symdvscf symmetrizes both components (spin up and down potentials).
  ! Noncollinear magnet : nspin_mag = 4, nspin_lsda = 1.
  !                       symdvscf symmetrizes the first component (potential),
  !                       and sym_dmag symmetrizes the other three components (magnetic field).
  !
  CALL symdvscf (ddvtosym, dfft)
  IF (noncolin .AND. domag) CALL sym_dmag(ddvtosym, dfft)
  !
  ! Scatter back the real-space points
  !
  nxyp = dfft%nr1x * dfft%my_nr2p
  DO iper = 1, lr_npert
     DO is = 1, nspin_mag
        DO ir3 = 1, dfft%my_nr3p
           ioff    = dfft%nr1x * dfft%my_nr2p * (ir3-1)
           ioff_tg = dfft%nr1x * dfft%nr2x    * (dfft%my_i0r3p+ir3-1) + dfft%nr1x * dfft%my_i0r2p
           CALL zcopy (nxyp, ddvtosym (ioff_tg+1, is, iper), 1, dvtosym (ioff+1, is, iper), 1)
        END DO
     ENDDO
  ENDDO
  DEALLOCATE (ddvtosym)
  !
  CALL stop_clock ('psymdvscf')
#else
  !
  CALL symdvscf(dvtosym, dfft)
  IF (noncolin .AND. domag) CALL sym_dmag(dvtosym, dfft)
  !
#endif
  !
  !
END SUBROUTINE psymdvscf
