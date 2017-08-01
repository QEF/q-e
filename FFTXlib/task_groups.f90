!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------=
   MODULE task_groups
!=----------------------------------------------------------------------=

!  ... Distribute G-vectors across processors as sticks and planes,
!  ... initialize FFT descriptors for both dense and smooth grids

      IMPLICIT NONE

      INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

      TYPE task_groups_descriptor
        !
        !  task groups logical
        !
        LOGICAL :: have_task_groups = .FALSE.
        !
      END TYPE

      SAVE

!=----------------------------------------------------------------------=
   END MODULE task_groups
!=----------------------------------------------------------------------=
