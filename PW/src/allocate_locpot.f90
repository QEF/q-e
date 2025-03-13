
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_locpot
  !-----------------------------------------------------------------------
  !! Dynamical allocation of arrays:
  !! local potential for each kind of atom, structure factor
  !
  USE ions_base, ONLY : nat, ntyp => nsp
  USE vlocal,    ONLY : vloc, strf
  USE gvect,     ONLY : eigts1, eigts2, eigts3, ngm, ngl
  USE fft_base , ONLY : dfftp
  !
  IMPLICIT NONE
    INTEGER :: nr1, nr2, nr3
  !
  nr1 = dfftp%nr1
  nr2 = dfftp%nr2
  nr3 = dfftp%nr3
  !
  ALLOCATE( vloc( ngl, ntyp) )
  ALLOCATE( strf( ngm, ntyp) )
  !
  ALLOCATE( eigts1(-dfftp%nr1:dfftp%nr1,nat) )
  ALLOCATE( eigts2(-dfftp%nr2:dfftp%nr2,nat) )
  ALLOCATE( eigts3(-dfftp%nr3:dfftp%nr3,nat) )
  !
  !FIXME  these are persistent device data
  !FIXME  use acc declare create() for eigts[123]  
  !$acc enter data create(eigts1(-nr1:nr1,1:nat), eigts2(-nr2:nr2,1:nat), eigts3(-nr3:nr3,1:nat) )
  !
  RETURN
  !
END SUBROUTINE ALLOCATE_locpot

