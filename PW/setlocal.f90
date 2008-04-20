!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine setlocal
  !----------------------------------------------------------------------
  !
  !    This routine computes the local potential in real space vltot(ir)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  USE ions_base, ONLY : ntyp => nsp
  USE extfield,  ONLY : tefield, dipfield, etotefield
  USE gvect,     ONLY : igtongl, gg
  USE scf,       ONLY : rho, v_of_0, vltot
  USE vlocal,    ONLY : strf, vloc
  USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, nlm, ngm
  USE control_flags, ONLY : gamma_only
  USE mp_global, ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum

  !
  implicit none
  complex(DP), allocatable :: aux (:)
  ! auxiliary variable
  integer :: nt, ng
  ! counter on atom types
  ! counter on g vectors
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
        aux (nlm(ng)) = CONJG(aux (nl(ng)))
     enddo
  end if
  !
  ! ... v_of_0 is (Vloc)(G=0)
  !
  v_of_0=0.0_DP
  if (gg(1) < eps8) v_of_0 = DBLE ( aux (nl(1)) )
  call mp_sum( v_of_0, intra_pool_comm )
  !
  ! ... aux = potential in G-space . FFT to real space
  !
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  vltot (:) =  DBLE (aux (:) )
  !
  ! ... If required add an electric field to the local potential 
  !
  if ( tefield .and. ( .not. dipfield ) ) then
     !
     call add_efield(rho%of_r,vltot,etotefield,0)
     !
     ! NB rho is not actually used by add_efield in this case ...
     !    it should be fixed and removed from this routine
  endif
  !
  deallocate(aux)
  !
  return
end subroutine setlocal

