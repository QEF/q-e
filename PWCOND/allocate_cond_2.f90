
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Generalized to spinor wavefunctions and spin-orbit Oct. 2004 (ADC).
!
subroutine allocate_cond_2
!
! This subroutine allocates the remaining variables 
! after reduction (ngper --> n2d) of XY basis set
!
#include "f_defs.h"
  USE noncollin_module, ONLY : npol
  use cond 
  implicit none

  allocate( newbg(ngper*npol, n2d) )
  allocate( psiper( n2d, n2d, nrzp ) )
  allocate( zk( n2d, nrzp ) )
  allocate( zkr( n2d, nrzp ) )

  allocate( fun0(n2d, 2*n2d) ) 
  allocate( fun1(n2d, 2*n2d) )
  allocate( fund0(n2d, 2*n2d) )
  allocate( fund1(n2d, 2*n2d) ) 

  allocate( funl0(n2d, norbf*npol) )     
  allocate( funl1(n2d, norbf*npol) )
  allocate( fundl0(n2d, norbf*npol) )
  allocate( fundl1(n2d, norbf*npol) )

  allocate( intw1(norbf*npol, 2*n2d) )
  allocate( intw2(norbf*npol, norbf*npol) )

  allocate( kvall(2*(n2d+npol*nocrosl)) )
  allocate( kfunl(n2d, 2*(n2d+npol*nocrosl)) ) 
  allocate( kfundl(n2d, 2*(n2d+npol*nocrosl)) )  
  allocate( kintl(nocrosl*npol, 2*(n2d+npol*nocrosl)) )  
  allocate( kcoefl(nocrosl*npol, 2*(n2d+npol*nocrosl)) )  

  if(ikind.ne.0) then
    allocate( kvalr(2*(n2d+npol*nocrosr)) ) 
    allocate( kfunr(n2d, 2*(n2d+npol*nocrosr)) )
    allocate( kfundr(n2d, 2*(n2d+npol*nocrosr)) )
    allocate( kintr(nocrosr*npol, 2*(n2d+npol*nocrosr)) )
    allocate( kcoefr(nocrosr*npol, 2*(n2d+npol*nocrosr)) )      
  endif

  return
end subroutine allocate_cond_2
