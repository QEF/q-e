!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE symme
  !
  USE kinds,      ONLY : DP
  !
  ! ... The variables needed to describe the symmetry properties
  !  
  SAVE
  !
  INTEGER :: &
       s(3,3,48),            &! simmetry matrices, in crystal axis
       ftau(3,48),           &! fractional translations, in FFT coordinates
       nrot,                 &! number of bravais lattice symmetries 
       nsym                   ! number of crystal symmetries
  INTEGER :: &
       t_rev(48) = 0          ! time reversal flag, for noncolinear magnetisation
  INTEGER, ALLOCATABLE :: &
       irt(:,:)               ! symmetric atom for each atom and sym.op.
  LOGICAL :: &
       time_reversal=.true., &! if .TRUE. the system has time_reversal symmetry
       invsym                 ! if .TRUE. the system has inversion symmetry
  REAL(DP),TARGET :: &
       d1(3,3,48),           &! matrices for rotating spherical
       d2(5,5,48),           &! harmonics (d1 for l=1, ...)
       d3(7,7,48)             !
  CHARACTER(LEN=45) ::  sname(48)   ! name of the symmetries
  !
END MODULE symme
