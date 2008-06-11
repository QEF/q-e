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
!-----------------------------------------------------------------------
subroutine set_mltgrid_dim
  !-----------------------------------------------------------------------
  !     This routine computes the dimensions of the minimum FFT grid
  !     compatible with the input cut-off
  !
  !     NB: The values of nr1c, nr2c, nr3c are computed only if they are not
  !     given as input parameters. Input values are kept otherwise.
  !
  USE cell_base, ONLY: at, tpiba2
  USE gcoarse, ONLY: gcutmc, nr1c, nr2c, nr3c
  USE ee_mod, ONLY: mr1, mr2, mr3, ecutcoarse, nlev
  use fft_scalar, only: allowed
  implicit none

  integer, parameter :: nmax = 5000
  ! an unreasonably big number for a FFT grid
  !
  !
  gcutmc = ecutcoarse / tpiba2
  !
  nr1c = mr1
  nr2c = mr2
  nr3c = mr3
  !
  ! the values of nr1c, nr2c, nr3c are computed only if they are not given
  ! as input parameters
  !
  if (nr1c == 0) then
     !
     ! estimate nr1c and check if it is an allowed value for FFT
     !
     nr1c = int (2 * sqrt (gcutmc) * sqrt (at (1, 1) **2 + at (2, 1) &
          **2 + at (3, 1) **2) ) + 1
     !
     ! ... Selects an appropriate grid resolution for the multigrid algorithm
     ! ... ( nr = 2 ** nlev * i + 1, where i is an integer ) 
     !
     nr1c = ( nr1c - 1 ) / 2 ** ( nlev - 1 ) * 2 ** ( nlev - 1 ) + 1
     !
  end if
  !
  if (nr2c == 0) then
     !
     ! estimate nr1c and check if it is an allowed value for FFT
     !
     nr2c = int (2 * sqrt (gcutmc) * sqrt (at (1, 2) **2 + at (2, 2) &
          **2 + at (3, 2) **2) ) + 1
     !
     ! ... Selects an appropriate grid resolution for the multigrid algorithm
     ! ... ( nr = 2 ** nlev * i + 1, where i is an integer )
     !
     nr2c = ( nr2c - 1 ) / 2 ** ( nlev - 1 ) * 2 ** ( nlev - 1 ) + 1
     !
  endif
  !
  if (nr3c == 0) then
     !
     ! estimate nr3c and check if it is an allowed value for FFT
     !
     nr3c = int (2 * sqrt (gcutmc) * sqrt (at (1, 3) **2 + at (2, 3) &
          **2 + at (3, 3) **2) ) + 1
     !
     ! ... Selects an appropriate grid resolution for the multigrid algorithm
     ! ... ( nr = 2 ** nlev * i + 1, where i is an integer )
     !
     nr3c = ( nr3c - 1 ) / 2 ** ( nlev - 1 ) * 2 ** ( nlev - 1 ) + 1
     !
  endif
  !
  mr1 = nr1c
  mr2 = nr2c
  mr3 = nr3c
  !
  return
end subroutine set_mltgrid_dim

