!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine init_run
  !-----------------------------------------------------------------------
  !
  use pwcom
  implicit none

  call start_clock ('init_run')

  if (gamma_only) then
     write (6, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials, ", &
                 &   "Gamma point")')
  else
     write (6, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")')
  end if
  write (6, 9010) ntypx, npk, lmaxx, nchix, ndm, nbrx, nqfm

  call iosys
  call setup
  !
  !   allocate memory for G- and R-space fft arrays
  !
  call allocate_fft
  !
  !   generate reciprocal-lattice vectors and fft indices
  !
  call ggen
#ifndef __PARA
  !
  !   generates pencils for 3d-fft of psi and related quantities
  !
  call set_pencils (nks, xk, ngms, gg, nls, ecutwfc / tpiba2, &
       nrx1s, nr1s, nr2s, nr3s)
#endif
  call summary
  call allocate_nlpot
  call allocate_locpot
  call allocate_wfc

  call openfil

  call hinit0
  call potinit

  call newd

  call wfcinit
  call stop_clock ('init_run')
  write (6, * )
  call show_memory ()

  return
9010 format (/5x,'Current dimensions of program pwscf are:' &
             /5x,'ntypx =',i2,'   npk =',i5,'  lmax =',i2   &
             /5x,'nchix =',i2,'  ndim =',i5,'  nbrx =',i2,' nqfm =',i2)
end subroutine init_run

