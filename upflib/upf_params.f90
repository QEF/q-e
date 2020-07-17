!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE upf_params
  
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: &
       npsx   = 10,     &! max number of different PPs
       lmaxx  = 3,      &! max non local angular momentum (l=0 to lmaxx)      
       lqmax= 2*lmaxx+1  ! max number of angular momenta of Q

END MODULE upf_params

