
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE set_defaults()
  !-----------------------------------------------------------------------------
  !
  USE input_parameters, ONLY : full_phs_path_flag
  !
  USE control_flags, ONLY : &
                            lscf, &
                            lpath
                            

  lscf      = .true.
  lpath       = .true.
  full_phs_path_flag = .true.
  !
END SUBROUTINE set_defaults
