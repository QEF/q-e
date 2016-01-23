!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

#if defined __AIX
#  define  __BSIZ_VALUE  55
#else 
#  define  __BSIZ_VALUE  35
#endif

