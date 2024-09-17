!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine input_pp_summary ( )
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends the acfdt input to all
  !     the other processors
  !
  !
  USE run_info,      ONLY : title
  USE io_global,     ONLY : stdout, ionode
  USE mp,            ONLY : mp_bcast
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_kcw
  !
  implicit none
  !
  IF (ionode) THEN 
    !
    WRITE( stdout, '(/, 5X, "KC_PP INPUT SUMMARY ")') 
!!!! CONTOL NAMELIST
    WRITE(stdout,'(5X,  42("="))')
    WRITE(stdout, 41)  "# title               =", TRIM(title)
    WRITE(stdout, 41)  "# out_dir             =", TRIM(tmp_dir)
    WRITE(stdout, 42)  "# prefix              =", TRIM(prefix)
    WRITE(stdout, 42)  "# seedname            =", TRIM(seedname)
    WRITE(stdout, 45)  "# kcw_iverbosity       =", kcw_iverbosity
    WRITE(stdout, 45)  "# num_wann            =", num_wann
    WRITE(stdout, 47)  "# MP grid             =", mp1, mp2, mp3
    WRITE(stdout, 43)  "# use_ws_distance     =", use_ws_distance
    WRITE(stdout, 43)  "# have_empty          =", have_empty
    WRITE(stdout, 43)  "# io_sp               =", io_sp
    WRITE(stdout, 43)  "# io_real_space       =", io_real_space
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
end subroutine input_pp_summary
