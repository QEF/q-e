!---------------------------------------------------------------------
subroutine dvb_cc (nlcc,npseu,ngm,nr1,nr2,nr3,nrx1,  &
     nl,rho_core,dmuxc,ga,aux,dvb_nlcc)
  !---------------------------------------------------------------------
  ! calcola il contributo core-correction al Delta V bare
  !
#include "machine.h"
  implicit none
  integer:: npseu,ngm,nr1,nr2,nr3,nrx1,nrxx,np,ng,i
  logical :: nlcc(npseu)
  integer :: nl(ngm)
  real(kind=8) :: rho_core(ngm), dmuxc(nrx1*nr2*nr3)
  complex(kind=8) :: ga(ngm), dvb_nlcc(ngm), aux(nrx1*nr2*nr3)
  !
  do np=1,npseu
     if(nlcc(np)) go to 10
  end do
  return
10 continue
  !
  nrxx=nrx1*nr2*nr3
  call setv(2*nrxx,0.d0,aux,1)
  do ng=1,ngm
     aux(nl(ng)) = ga(ng) * rho_core(ng)
  end do
  call cft3(aux,nr1,nr2,nr3,nrx1,nr2,nr3,1)
  !
  do i=1,nrxx
     aux(i) = aux(i) * dmuxc(i)
  end do
  call cft3(aux,nr1,nr2,nr3,nrx1,nr2,nr3,-1)
  do ng=1,ngm
     dvb_nlcc(ng) = aux(nl(ng))
  end do
  !
  return
end subroutine dvb_cc

