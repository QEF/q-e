!
! Copyright (C) 2001 PWSCF group
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
  USE kinds, ONLY: DP
  USE basis,  ONLY: ntyp
  USE extfield, ONLY: tefield, dipfield, etotefield
  USE gvect,  ONLY: igtongl
  USE scf,    ONLY: vltot
  USE vlocal, ONLY : strf, vloc
  USE wvfct,  ONLY: gamma_only
  USE gvect,  ONLY: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, nlm, ngm
  USE scf,    ONLY: rho
  implicit none
  complex(kind=DP), allocatable :: aux (:)
  ! auxiliary variable
  integer :: nt, ng, ir
  ! counter on atom types
  ! counter on g vectors
  ! counter on r vectors
  !
    allocate (aux( nrxx))    
  !
  aux(:)=(0.d0,0.d0)
  do nt = 1, ntyp
     do ng = 1, ngm
        aux (nl(ng))=aux(nl(ng)) + vloc (igtongl (ng), nt) * strf (ng, nt)
     enddo
  enddo
  if (gamma_only) then
     do ng = 1, ngm
        aux (nlm(ng)) = conjg (aux (nl(ng)))
     enddo
  end if
  !
  ! aux = potential in G-space . FFT to real space
  !
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  do ir = 1, nrxx
     vltot (ir) = DREAL (aux (ir) )
  enddo
!
!  If required add an electric field to the local potential 
!

  if (tefield.and.(.not.dipfield)) then
     call add_efield(rho,vltot,etotefield)
! NB rho is not actually used by add_efield in this case ...
!    it should be fixed and removed from this routine
  endif

  deallocate(aux)
  return
end subroutine setlocal

