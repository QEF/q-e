!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine allocate_cond_2
!
! This subroutine allocates the remaining variables 
! after reduction (ngper --> n2d) of XY basis set
!
#include "machine.h"
  use cond 

  allocate( newbg(ngper, n2d) )
  allocate( psiper( n2d, n2d, nrzp ) )
  allocate( zk( n2d, nrzp ) )
  allocate( zkr( n2d, nrzp ) )

  allocate( fun0(n2d, 2*n2d) ) 
  allocate( fun1(n2d, 2*n2d) )
  allocate( fund0(n2d, 2*n2d) )
  allocate( fund1(n2d, 2*n2d) ) 

  allocate( funl0(n2d, norbf) )     
  allocate( funl1(n2d, norbf) )
  allocate( fundl0(n2d, norbf) )
  allocate( fundl1(n2d, norbf) )

  allocate( intw1(norbf, 2*n2d) )
  allocate( intw2(norbf, norbf) )

  allocate( kvall(2*(n2d+nocrosl)) )
  allocate( kfunl(n2d, 2*(n2d+nocrosl)) ) 
  allocate( kfundl(n2d, 2*(n2d+nocrosl)) )  
  allocate( kintl(nocrosl, 2*(n2d+nocrosl)) )  
  allocate( kcoefl(nocrosl, 2*(n2d+nocrosl)) )  

  if(ikind.ne.0) then
    allocate( kvalr(2*(n2d+nocrosr)) ) 
    allocate( kfunr(n2d, 2*(n2d+nocrosr)) )
    allocate( kfundr(n2d, 2*(n2d+nocrosr)) )
    allocate( kintr(nocrosr, 2*(n2d+nocrosr)) )
    allocate( kcoefr(nocrosr, 2*(n2d+nocrosr)) )      
  endif

  return
end subroutine allocate_cond_2
