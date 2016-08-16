!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------=
   MODULE stick_set
!=----------------------------------------------------------------------=

!  ... Distribute G-vectors across processors as sticks and planes,
!  ... initialize FFT descriptors for both dense and smooth grids

      USE stick_base

      IMPLICIT NONE

      PRIVATE
      SAVE

      TYPE(sticks_map) :: smap

      PUBLIC :: pstickdealloc
      PUBLIC :: smap

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

      SUBROUTINE pstickdealloc()
         CALL sticks_map_deallocate( smap )
      END SUBROUTINE pstickdealloc

!=----------------------------------------------------------------------=
   END MODULE stick_set
!=----------------------------------------------------------------------=
