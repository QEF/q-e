!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE small_box
!------------------------------------------------------------------------------!

      !  This module contains the basis vector of the small sub-cell (small box)
      !  used for charge augmentation process

      USE kinds, ONLY : dbl
!
      IMPLICIT NONE
      SAVE

        !  a1, a2 and a3 are the simulation cell base vector as calculated from celldm

      REAL(dbl) :: a1b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
      REAL(dbl) :: a2b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
      REAL(dbl) :: a3b(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)

      REAL(dbl) :: ainvb(3,3) = 0.0d0

      REAl(dbl) :: omegab = 0.0d0  !  volume of the small boxes 

      REAL(dbl) :: tpibab = 0.0d0
!
!------------------------------------------------------------------------------!
   END MODULE small_box
!------------------------------------------------------------------------------!
