!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!--------------------------------------------------------------------------
!
MODULE acfdt_ener
  !
  USE kinds, ONLY : DP
  !
  LOGICAL  :: acfdt_in_pw = .FALSE.  ! the default value is .false.
  REAL(DP) :: acfdt_eband                    
  !
END MODULE acfdt_ener

