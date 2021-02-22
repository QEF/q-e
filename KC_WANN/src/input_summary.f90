!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine input_summary ( )
  !-----------------------------------------------------------------------
  !
  !! Summary of the input variables
  !
  USE run_info,      ONLY : title
  USE io_global,     ONLY : stdout, ionode
  USE mp,            ONLY : mp_bcast
  USE control_ph,    ONLY : nmix_ph, tr2_ph, niter_ph
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_kc_wann
  USE control_lr,     ONLY : lrpa
  USE input_parameters, ONLY : assume_isolated
  !
  implicit none
  !
  IF (ionode) THEN 
    !
    WRITE( stdout, '(/, 5X, "KCWANN INPUT SUMMARY ")') 
!!!! CONTOL NAMELIST
    WRITE(stdout,'(5X,  42("="))')
    WRITE(stdout, 41)  "# title               =", TRIM(title)
    WRITE(stdout, 41)  "# out_dir             =", TRIM(tmp_dir)
    WRITE(stdout, 42)  "# prefix              =", TRIM(prefix)
    WRITE(stdout, 42)  "# seedname            =", TRIM(seedname)
    WRITE(stdout, 45)  "# kc_iverbosity       =", kc_iverbosity
    WRITE(stdout, 43)  "# kc_at_ks            =", kc_at_ks     
    WRITE(stdout, 43)  "# homo_only           =", homo_only    
    WRITE(stdout, 43)  "# read_unitary_matrix =", read_unitary_matrix  
    WRITE(stdout, 43)  "# check_ks            =", check_ks
!!! WANNIER
    WRITE(stdout, 45)  "# num_wann_occ        =", num_wann_occ
    IF (have_empty) WRITE(stdout, 45)  "# num_wann_emp        =", num_wann_emp
    WRITE(stdout, 43)  "# have_empty          =", have_empty 
    WRITE(stdout, 43)  "# has_disentangle     =", has_disentangle 
    WRITE(stdout, 45)  "# spin_component      =", spin_component
    WRITE(stdout, 43)  "# l_unique_manifold   =", l_unique_manifold 
!! SCREEN
    WRITE(stdout, 43)  "# lrpa                =", lrpa     
    IF(kc_at_ks) WRITE(stdout, 43)  "# fix_orb             =", fix_orb      
    WRITE(stdout, 46)  "# tr2_ph              =", tr2_ph
    WRITE(stdout, 45)  "# niter_ph            =", niter_ph
    WRITE(stdout, 45)  "# nmix_ph             =", nmix_ph
    IF (i_orb /= -1 )   WRITE(stdout, 45)  "# i_orb               =", i_orb
    WRITE(stdout, 43)  "# check_spread        =", check_spread     
!!! HAM
    WRITE(stdout, 43)  "# qp_symm             =", qp_symm      
    WRITE(stdout, 43)  "# kipz_corr           =", kipz_corr    
    WRITE(stdout, 43)  "# compute_hf          =", compute_hf    
    WRITE(stdout, 47)  "# MP grid             =", mp1, mp2, mp3
    WRITE(stdout, 43)  "# do_bands            =", do_bands
    WRITE(stdout, 43)  "# use_ws_distance     =", use_ws_distance
    WRITE(stdout, 43)  "# write_hr            =", write_hr
    WRITE(stdout, 43)  "# l_vcut              =", l_vcut    
    WRITE(stdout, 43)  "# l_alpha_corr        =", l_alpha_corr
    WRITE(stdout, 42)  "# assume_isolated     =", TRIM(assume_isolated)
    WRITE(stdout, 46)  "# eps_inf             =", eps_inf
    WRITE(stdout,'(5X, 42("="),/)')
    !
  ENDIF
  !
41  FORMAT(5X, A23, A18)
42  FORMAT(5X, A23, A18)
43  FORMAT(5X, A23, L18)
44  FORMAT(5X, A23, F12.8, 3x, " [Ry]")
45  FORMAT(5X, A23, I18)
46  FORMAT(5X, A23, E18.2)
47  FORMAT(5X, A23, 3I6)
  !
  RETURN
  !
end subroutine input_summary
