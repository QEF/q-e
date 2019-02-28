!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------
SUBROUTINE lr_psym_eels (dvtosym)
  !---------------------------------------------------------------
  !
  ! Symmetrize the response charge density (parallel case).
  ! Inspired by PH/psymdvscf.f90
  !
  ! TODO: Try to use the routine PH/psymdvscf.f90 directly.
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : nspin_mag
  USE mp_global,        ONLY : me_bgrp
  USE fft_base,         ONLY : dfftp
  USE scatter_mod,      ONLY : cgather_sym
  USE lr_symm_base,     ONLY : nsymq, minus_q
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym(dfftp%nnr, nspin_mag)  
  ! the charge density response to symmetrize
    
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, ir3, ioff, ioff_tg, nxyp 
  COMPLEX(DP), ALLOCATABLE :: ddvtosym(:,:)
  ! the potential to symm
  !
  IF (nsymq==1) RETURN
  !IF (nsymq==1 .AND. (.NOT.minus_q) ) RETURN
  !
  CALL start_clock ('lr_psym_eels')
  !
  ALLOCATE(ddvtosym(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin_mag))
  !
  ! Gather complex data for the symmetrization.
  !
  DO is = 1, nspin_mag
     CALL cgather_sym ( dfftp, dvtosym(:,is), ddvtosym(:,is) )
  ENDDO
  !
  ! Symmetrization
  !
  CALL lr_sym_eels (ddvtosym)
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
  CALL stop_clock ('lr_psym_eels')
  !
#endif
  !
  RETURN
  !
END SUBROUTINE lr_psym_eels
