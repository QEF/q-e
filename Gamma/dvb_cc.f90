!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine dvb_cc (nlcc,npseu,ngm,nr1,nr2,nr3,nrx1,nrx2,nrx3,  &
     nl,rho_core,dmuxc,ga,aux,dvb_nlcc)
  !---------------------------------------------------------------------
  ! calculate the core-correction contribution to Delta V bare
  !
#include "f_defs.h"
  implicit none
  integer:: npseu,ngm,nr1,nr2,nr3,nrx1,nrx2,nrx3,np,ng,i
  logical :: nlcc(npseu)
  integer :: nl(ngm)
  real(8) :: rho_core(ngm), dmuxc(nrx1*nrx2*nrx3)
  complex(8) :: ga(ngm), dvb_nlcc(ngm), aux(nrx1*nrx2*nrx3)
  !
  do np=1,npseu
     if(nlcc(np)) go to 10
  end do
  return
10 continue
  !
  aux(:) = (0.d0, 0.d0)
  do ng=1,ngm
     aux(nl(ng)) = ga(ng) * rho_core(ng)
  end do
  call cft3(aux,nr1,nr2,nr3,nrx1,nr2,nr3,1)
  !
  aux(:) = aux(:) * dmuxc(:)
  !
  call cft3(aux,nr1,nr2,nr3,nrx1,nr2,nr3,-1)
  do ng=1,ngm
     dvb_nlcc(ng) = aux(nl(ng))
  end do
  !
  return
end subroutine dvb_cc

