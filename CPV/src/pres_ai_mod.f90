!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE pres_ai_mod
!=----------------------------------------------------------------------------=!

  use kinds, only: dp
  use parameters, only: nsx
  IMPLICIT NONE
!
  SAVE
!
  logical :: abivol, abisur, pvar, fill_vac, scale_at, t_gauss, jellium
  logical :: cntr(nsx)
  real(dp), allocatable:: rho_gaus(:), posv(:,:), v_vol(:), f_vol(:,:,:)
  real(dp) :: P_ext, P_in, P_fin, rho_thr, step_rad(nsx)
  real(dp) :: Surf_t, dthr, volclu, surfclu, n_ele, nelect
  real(dp) :: R_j, e_j, h_j
  real(dp) :: stress_vol(3,3)
  real(dp) :: delta_eps, delta_sigma
  real(dp) :: xc0(500), weight(500)
  integer :: n_cntr, axis

CONTAINS

     !----------------------------------------------------------------------
  SUBROUTINE pres_ai_init (abivol_, abisur_, pvar_, fill_vac_,     &
                           scale_at_, t_gauss_, jellium_, cntr_,   &
                           P_ext_, P_in_, P_fin_, rho_thr_,        &
                           step_rad_, Surf_t_, dthr_, R_j_, h_j_,  &
                           delta_eps_, delta_sigma_, n_cntr_, axis_)
     !----------------------------------------------------------------------
     !
     USE constants, ONLY : au_gpa
     !
     IMPLICIT NONE
     !
     LOGICAL  :: abivol_, abisur_, pvar_, fill_vac_, scale_at_, &
                 t_gauss_, jellium_, cntr_(nsx)
     REAL(dp) :: P_ext_, P_in_, P_fin_, rho_thr_, step_rad_(nsx), &
                 Surf_t_, dthr_, R_j_, h_j_, delta_eps_, delta_sigma_
     INTEGER  :: n_cntr_, axis_
     !
     !  Copy variables read from input into module
     !
     abivol = abivol_
     abisur = abisur_
     pvar = pvar_
     fill_vac = fill_vac_
     scale_at = scale_at_
     t_gauss  = t_gauss_ 
     cntr_(:) = cntr_(:)
     jellium = .false.  ! provvisorio
     rho_thr = rho_thr_
     step_rad(:) = step_rad_(:)
     Surf_t = Surf_t_
     dthr = dthr_
     R_j = R_j_
     h_j = h_j_
     delta_eps = delta_eps_
     delta_sigma = delta_sigma_
     n_cntr = n_cntr_
     axis = axis_
     !
     !  Correct (a.u.) units to pressure
     !
     P_ext = P_ext_ / au_gpa
     P_in = P_in_ / au_gpa
     P_fin = P_fin_ / au_gpa
     if (pvar) P_ext = P_in
     !
  END SUBROUTINE pres_ai_init
!=----------------------------------------------------------------------------=!
   END MODULE pres_ai_mod
!=----------------------------------------------------------------------------=!
