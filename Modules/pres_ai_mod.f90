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

  use parameters, only: nsx, natx
!
  SAVE
!
  logical :: abivol, abisur, pvar, fill_vac, scale_at, t_gauss, jellium
  logical :: cntr(nsx)
  real(kind=8), allocatable:: rho_gaus(:), posv(:,:), v_vol(:), f_vol(:,:,:)
  real(kind=8) :: P_ext, P_in, P_fin, rho_thr, step_rad(nsx)
  real(kind=8) :: Surf_t, dthr, volclu, surfclu, n_ele, nelect
  real(kind=8) :: R_j, e_j, h_j
  real(kind=8) :: stress_vol(3,3)
  real(kind=8) :: delta_eps, delta_sigma
  real(kind=8) :: xc0(500), weight(500)
  integer n_cntr, axis


!=----------------------------------------------------------------------------=!
   END MODULE pres_ai_mod
!=----------------------------------------------------------------------------=!
