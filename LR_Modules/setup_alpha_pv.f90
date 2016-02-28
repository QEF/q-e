!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE setup_alpha_pv
  !-----------------------------------------------------------------------
  !
  ! This subroutine computes alpha_pv (i.e. the coefficient
  ! in front of the projector on occupied states in the
  ! linear response system, in order to avoid singularity problems)
  ! See Eq. (30) in Baroni et al., Rev. Mod. Phys. 73, 515 (2001)
  ! This coefficient is also used in orthogonalize.f90 when lgauss=.true.
  !
  USE kinds,      ONLY : DP
  USE wvfct,      ONLY : nbnd, et
  USE klist,      ONLY : nks
  USE mp,         ONLY : mp_max, mp_min
  USE mp_pools,   ONLY : inter_pool_comm
  USE klist,      ONLY : lgauss
  USE control_lr, ONLY : alpha_pv, nbnd_occ
  !
  IMPLICIT NONE
  !
  REAL(DP) :: target, emin, emax
  ! auxiliary variable used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy
  INTEGER :: ik, ibnd
  !
  CALL start_clock ('setup_alpha_pv')
  !
  emin = et (1, 1)
  !
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     ENDDO
  ENDDO
  !
#ifdef __MPI
  ! Find the minimum across pools
  CALL mp_min( emin, inter_pool_comm )
#endif
  !
  IF (lgauss) THEN
     emax = target
     alpha_pv = emax - emin
  ELSE
     emax = et (1, 1)
     DO ik = 1, nks
        DO ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        ENDDO
     ENDDO
#ifdef __MPI
     ! Find the maximum across pools
     CALL mp_max( emax, inter_pool_comm )
#endif
     alpha_pv = 2.d0 * (emax - emin)
  ENDIF
  !
  ! Avoid zero value for alpha_pv
  !
  alpha_pv = max (alpha_pv, 1.0d-2)
  !
  CALL stop_clock ('setup_alpha_pv')
  !
  RETURN
  !
END SUBROUTINE setup_alpha_pv
