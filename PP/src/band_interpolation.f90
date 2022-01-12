!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Author: Ivan Carnimeo (September 2021)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
program band_interpolation 
use globalmod
use shepardmod
use fouriermod
use fourierdiffmod
implicit none
  !
  Call read_input ()
  !
  ek = 0.0d0
  !
  if(TRIM(method).eq.'shepard') then 
    !
    Call shepard (1, Nb, Nq, q, eq, Nk, k, ek)
    !
  elseif(TRIM(method).eq.'shepard-sphere') then 
    !
    Call shepard (2, Nb, Nq, q, eq, Nk, k, ek)
    !
  elseif(TRIM(method).eq.'fourier') then 
    !
    Call fourier (Nb, Nq, q, eq, Nk, k, ek, Nsym, at, bg, Op)
    !
  elseif(TRIM(method).eq.'fourier-diff') then 
    !
    Call fourierdiff (Nb, Nq, q, eq, Nk, k, ek, Nsym, at, bg, Op)
    !
  else
    !
    write(*,*) 'ERROR: Wrong method ', TRIM(method)
    stop
    !
  end if
  !
  Call print_bands (TRIM(method))
  !
  Call deallocate_global ()
  !
end program
!----------------------------------------------------------------------------
