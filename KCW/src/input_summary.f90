!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
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
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_kcw
  USE control_lr,     ONLY : lrpa
  USE input_parameters, ONLY : assume_isolated
  !
  implicit none
  !
  IF (ionode) THEN 
    !
    WRITE( stdout, '(/, 5X, "KCW INPUT SUMMARY ")') 
    WRITE(stdout,'(5X,  44("="))')
!!!! CONTOL NAMELIST
    WRITE( stdout, '(/, 6X, "CONTROL ")') 
    WRITE(stdout, 41)  "# title               =", TRIM(title)
    WRITE(stdout, 41)  "# out_dir             =", TRIM(tmp_dir)
    WRITE(stdout, 42)  "# prefix              =", TRIM(prefix)
    WRITE(stdout, 42)  "# calculation         =", TRIM(calculation)
    WRITE(stdout, 45)  "# kcw_iverbosity      =", kcw_iverbosity
    WRITE(stdout, 43)  "# kcw_at_ks           =", kcw_at_ks     
    WRITE(stdout, 47)  "# MP grid             =", mp1, mp2, mp3
    WRITE(stdout, 45)  "# spin_component      =", spin_component
    WRITE(stdout, 43)  "# homo_only           =", homo_only    
    WRITE(stdout, 43)  "# read_unitary_matrix =", read_unitary_matrix  
    WRITE(stdout, 43)  "# check_ks            =", check_ks
    WRITE(stdout, 43)  "# l_vcut              =", l_vcut    
    WRITE(stdout, 42)  "# assume_isolated     =", TRIM(assume_isolated)
    WRITE(stdout, 43)  "# io_sp               =", io_sp
    WRITE(stdout, 43)  "# io_real_space       =", io_real_space
    WRITE(stdout, 43)  "# irr_bz              =", irr_bz
    WRITE(stdout, 43)  "# use_wct         =", use_wct
    !
    IF ( .NOT. kcw_at_ks .AND. .NOT. calculation=='cc' ) THEN 
!!! WANNIER
      WRITE( stdout, '(/, 6X, "WANNIER ")') 
      WRITE(stdout, 42)  "# seedname            =", TRIM(seedname)
      WRITE(stdout, 45)  "# num_wann_occ        =", num_wann_occ
      IF (have_empty) WRITE(stdout, 45)  "# num_wann_emp        =", num_wann_emp
      WRITE(stdout, 43)  "# have_empty          =", have_empty 
      WRITE(stdout, 43)  "# has_disentangle     =", has_disentangle 
      WRITE(stdout, 43)  "# l_unique_manifold   =", l_unique_manifold 
    ENDIF
    !
    IF (calculation == 'screen') THEN 
!! SCREEN
      WRITE( stdout, '(/, 6X, "SCREEN ")') 
      WRITE(stdout, 43)  "# lrpa                =", lrpa     
      IF(kcw_at_ks) WRITE(stdout, 43)  "# fix_orb             =", fix_orb      
      WRITE(stdout, 46)  "# tr2                 =", tr2
      WRITE(stdout, 45)  "# niter               =", niter
      WRITE(stdout, 45)  "# nmix                =", nmix
      WRITE(stdout, 46)  "# eps_inf             =", eps_inf
      IF (i_orb /= -1 )   WRITE(stdout, 45)  "# i_orb               =", i_orb
      WRITE(stdout, 43)  "# check_spread        =", check_spread     
    ENDIF
    !
    IF ( calculation == 'ham') THEN
!!! HAM
      WRITE( stdout, '(/, 6X, "HAM ")') 
      WRITE(stdout, 43)  "# qp_symm             =", qp_symm      
      WRITE(stdout, 43)  "# kipz_corr           =", kipz_corr    
      WRITE(stdout, 47)  "# MP grid             =", mp1, mp2, mp3
      WRITE(stdout, 43)  "# do_bands            =", do_bands
      WRITE(stdout, 43)  "# use_ws_distance     =", use_ws_distance
      WRITE(stdout, 43)  "# write_hr            =", write_hr
      WRITE(stdout, 43)  "# l_alpha_corr        =", l_alpha_corr
      WRITE(stdout, 43)  "# on_site_only        =", on_site_only
    ENDIF
    !
    WRITE(stdout,'(5X, 44("="),/)')
    !
  ENDIF
  !
41  FORMAT(7X, A23, A18)
42  FORMAT(7X, A23, A18)
43  FORMAT(7X, A23, L18)
44  FORMAT(7X, A23, F12.8, 3x, " [Ry]")
45  FORMAT(7X, A23, I18)
46  FORMAT(7X, A23, E18.4)
47  FORMAT(7X, A23, 3I6)
  !
  RETURN
  !
end subroutine input_summary
