!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE nmr_mod
  !
  ! Extra variable needed to describe nmr
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  SAVE
  !
  INTEGER, ALLOCATABLE :: &
       igk_q(:)             !  correspondence k+q+G <-> G
  REAL(DP), ALLOCATABLE :: &
       g2kin_q(:)             !  kinetic energy for k+q
  COMPLEX(DP), ALLOCATABLE :: &
       vkb_q(:,:)           ! all beta functions in reciprocal space for k+q

  CONTAINS

    SUBROUTINE init_nmr()
      
      USE wvfct,                 only: npwx
      USE uspp,                 ONLY : nkb
      IMPLICIT NONE

      allocate (igk_q(npwx))
      allocate (g2kin_q(npwx))
      allocate (vkb_q(npwx,  nkb))

    end SUBROUTINE init_nmr
  
  end MODULE nmr_mod
