!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine setlocal
  !----------------------------------------------------------------------
  !
  !    This routine computes the local potential in real space vltot(ir)
  !
#include "machine.h"
  use pwcom
  use gamma
  implicit none
  complex(kind=DP), allocatable :: aux (:)
  ! auxiliary variable used to contain the potential

  integer :: nt, ng, ir
  ! counter on atomtypes
  ! counter on g vectors
  ! counter on r vectors
  !
  !  calculate the fourier coefficients of the local part of the bare
  !  pseudopotential
  !

  allocate (aux( nrxx))    
  !
  !    here we compute the local potential in real space
  !
  aux(:)=(0.d0,0.d0)
  do nt = 1, ntyp
     do ng =1, ngm
        aux (nl(ng) ) = aux (nl(ng) ) + vloc (igtongl(ng), nt) * strf (ng, nt)
     enddo
  enddo
  do ng = 1, ngm
     aux (nlm(ng)) = conjg (aux (nl(ng)))
  enddo
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  ! vltot is the local part of the electron-ion interaction in real space
  !
  do ir = 1, nrxx
     vltot (ir) = DREAL (aux (ir) )
  enddo
  deallocate(aux)
  !
  !  If required add an electric field to the local potential 
  !
  if (tefield.and.(.not.dipfield)) then
     call add_efield(vltot)
  endif
  !
  return
end subroutine setlocal

