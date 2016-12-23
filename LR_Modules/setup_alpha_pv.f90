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
  USE constants,  ONLY : pi
  USE wvfct,      ONLY : nbnd, et
  USE klist,      ONLY : nks, lgauss, ngauss, degauss, ltetra
  USE ener,       ONLY : ef
  USE mp,         ONLY : mp_max, mp_min
  USE mp_pools,   ONLY : inter_pool_comm
  USE control_lr, ONLY : alpha_pv, nbnd_occ
  USE dfpt_tetra_mod, ONLY : dfpt_tetra_main
  !
  IMPLICIT NONE
  !
  REAL(DP) :: small, fac, xmax, emin, emax
  ! auxiliary variable used in the metallic case
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
#if defined(__MPI)
  ! Find the minimum across pools
  CALL mp_min( emin, inter_pool_comm )
#endif
  !
  IF (lgauss) THEN
     !
     ! Discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint:
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     IF (ngauss.eq. - 99) THEN
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     ENDIF
     !
     emax = ef + xmax * degauss
     !
     alpha_pv = emax - emin
     !
  ELSE IF(ltetra) then
     !
     CALL dfpt_tetra_main()
     !
  ELSE
     !
     emax = et (1, 1)
     !
     DO ik = 1, nks
        DO ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        ENDDO
     ENDDO
     !
#if defined(__MPI)
     ! Find the maximum across pools
     CALL mp_max( emax, inter_pool_comm )
#endif
     !
     alpha_pv = 2.d0 * (emax - emin)
     !
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
