!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_int_forces()
!----------------------------------------------------------------------------
!
USE kinds,            ONLY : DP
USE io_global,        ONLY : stdout
USE ions_base,        ONLY : nat, ityp
USE force_mod,        ONLY : force
USE control_flags,    ONLY : iverbosity
USE plugin_flags
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
INTEGER  :: ipol, na
  ! counter on polarization
  ! counter on atoms
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!
END SUBROUTINE plugin_int_forces

!----------------------------------------------------------------------
SUBROUTINE external_wg_corr_force( rhor, force )
!----------------------------------------------------------------------
  !
  USE kinds,             ONLY : DP
  USE cell_base,         ONLY : omega
  USE ions_base,         ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE gvect,             ONLY : ngm, g
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft
  USE martyna_tuckerman, ONLY : wg_corr_force
  USE vlocal,            ONLY : strf
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN) ::  rhor ( dfftp%nnr )
  REAL( DP ), INTENT(OUT) :: force (3, nat)
  !
  ! ... Local variables
  !
  COMPLEX (DP), ALLOCATABLE :: auxg( : ), auxr( : )
  !
  allocate(auxr(dfftp%nnr))
  auxr = cmplx(rhor,0.0_dp,kind=dp)
  call fwfft ("Rho", auxr, dfftp)
  !
  allocate(auxg(ngm))
  auxg = cmplx(0.0_dp,0.0_dp,kind=dp)
  auxg(:)=auxr(dfftp%nl(:))
  deallocate(auxr)
  !
  call wg_corr_force(.false.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                     auxg, force)
  !
  deallocate(auxg)
  !
  RETURN
  !
END SUBROUTINE external_wg_corr_force

!----------------------------------------------------------------------
SUBROUTINE external_force_lc( rhor, force )
!----------------------------------------------------------------------
  !
  !  interface to call force_lc from external programs
  !
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : alat, omega
  USE ions_base,     ONLY : nat, ityp, tau
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm, gstart, ngl, igtongl, g
  USE control_flags, ONLY : gamma_only
  !
  USE vlocal,        ONLY : vloc
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN) ::  rhor ( dfftp%nnr )
  REAL( DP ), INTENT(OUT) :: force ( 3, nat )
  !
  ! ... Local variables
  !
  CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
       g, rhor, gstart, gamma_only, vloc, force )
  !
  RETURN
  !
END SUBROUTINE external_force_lc
