!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE electrons_base
!------------------------------------------------------------------------------!

      USE kinds, ONLY: dbl
!
      IMPLICIT NONE
      SAVE

      INTEGER :: nbnd       = 0    !  total number of electronic states
      INTEGER :: nbndx      = 0    !  array dimension nbndx >= nbnd
      INTEGER :: nel(2)     = 0    !  number of electrons (up, down)
      INTEGER :: nelt       = 0    !  total number of electrons ( up + down )
      INTEGER :: nupdwn(2)  = 0    !  number of states with spin up (1) and down (2)
      INTEGER :: iupdwn(2)  = 0    !  first state with spin (1) and down (2)
      INTEGER :: nspin      = 0    !  nspin = number of spins (1=no spin, 2=LSDA)
!

!------------------------------------------------------------------------------!
  END MODULE electrons_base
!------------------------------------------------------------------------------!
