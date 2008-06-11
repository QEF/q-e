!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!--------------------------------------------------------------------------
!
!DCC
MODULE gcoarse
  !
  !
  USE kinds, ONLY : DP
  USE fft_types, ONLY: fft_dlay_descriptor
  !
  SAVE
  !
  INTEGER :: &
       ngmc,        &! the number of smooth G vectors
       ngmc_g,      &! the global number of smooth G vectors
                     !  (sum over all processors)
       ngmc_l,      &! the local number of smooth G vectors
                     !  (only present processor)
       nr1c,        &!
       nr2c,        &! the dimension of the smooth grid
       nr3c,        &!
       nrx1c,       &! maximum dimension of the smooth grid
       nrx2c,       &! maximum dimension of the smooth grid
       nrx3c,       &! maximum dimension of the smooth grid
       nrxxc         ! the total dimension of the smooth grid
  TYPE ( fft_dlay_descriptor ) :: dfftc ! coarse grid (dcc method)

  INTEGER, POINTER :: &
       nlc(:),      &! the correspondence  G <-> smooth mesh
       nlcm(:)       ! the same for gamma point calculation
  REAL(DP) :: &
       gcutmc        ! the cut-off of the smooth mesh
  !

END MODULE gcoarse
!
