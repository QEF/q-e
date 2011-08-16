!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Generalized to spinor wavefunctions and spin-orbit Oct. 2004 (ADC).
!
subroutine allocate_cond
!
! This subroutine allocates some needed variables
!
  USE fft_base,         ONLY : dfftp
  use lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : npol
  use cond

  implicit none

  allocate( newbg(ngper*npol, n2d) )
  allocate( psiperl( n2d, n2d, nrzpl ) )
  allocate( zkl( n2d, nrzpl ) )
  allocate( zkrl( n2d, nrzpl ) )
  allocate( psipers( n2d, n2d, nrzps ) )
  allocate( zks( n2d, nrzps ) )
  allocate( zkrs( n2d, nrzps ) )
  allocate( psiperr( n2d, n2d, nrzpr ) )
  allocate( zkr( n2d, nrzpr ) )
  allocate( zkrr( n2d, nrzpr ) )


  allocate( fun0(n2d, 2*n2d) )
  allocate( fun1(n2d, 2*n2d) )
  allocate( fund0(n2d, 2*n2d) )
  allocate( fund1(n2d, 2*n2d) )

  IF (lorb) THEN
    allocate( funz0(n2d, 2*n2d+norbf*npol, MAX(nrzpl,nrzps)) )
    allocate( korbl(npol*(nocrosl+noinsl), 2*(n2d+npol*nocrosl)) )
    if (ikind.eq.0) then
     allocate( korbr(npol*nocrosl, 2*(n2d+npol*nocrosl)) )
    else
     allocate( korbr(npol*nocrosr, 2*(n2d+npol*nocrosr)) )
    endif
  ENDIF

  IF (lcharge.and..NOT.ALLOCATED(rho_scatt)) THEN
    allocate( rho_scatt(dfftp%nr3,nspin) )
    rho_scatt(:,:) = 0.d0
  ENDIF

  IF (norbf>0) THEN
     allocate( funl0(n2d, norbf*npol) )
     allocate( funl1(n2d, norbf*npol) )
     allocate( fundl0(n2d, norbf*npol) )
     allocate( fundl1(n2d, norbf*npol) )

     allocate( intw1(norbf*npol, 2*n2d) )
     allocate( intw2(norbf*npol, norbf*npol) )
  ENDIF

  allocate( kvall(2*(n2d+npol*nocrosl)) )
  allocate( kfunl(n2d, 2*(n2d+npol*nocrosl)) )
  allocate( kfundl(n2d, 2*(n2d+npol*nocrosl)) )
  IF (nocrosl>0) THEN
     allocate( kintl(nocrosl*npol, 2*(n2d+npol*nocrosl)) )
     allocate( kcoefl(nocrosl*npol, 2*(n2d+npol*nocrosl)) )
  ENDIF

  if(ikind.ne.0) then
    allocate( kvalr(2*(n2d+npol*nocrosr)) )
    allocate( kfunr(n2d, 2*(n2d+npol*nocrosr)) )
    allocate( kfundr(n2d, 2*(n2d+npol*nocrosr)) )
    IF (nocrosr>0) THEN
       allocate( kintr(nocrosr*npol, 2*(n2d+npol*nocrosr)) )
       allocate( kcoefr(nocrosr*npol, 2*(n2d+npol*nocrosr)) )
    ENDIF
  endif

  return
end subroutine allocate_cond
