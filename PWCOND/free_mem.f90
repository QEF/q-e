!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine free_mem
!
! Deallocates memory
!
#include "machine.h"
  use cond
!
! From local_2
!
  deallocate(psiper)
  deallocate(zk)
  deallocate(zkr)
!
! From allocate_cond_2
!
  deallocate(newbg)

  deallocate(fun0) 
  deallocate(fun1)
  deallocate(fund0)
  deallocate(fund1)

  deallocate(funl0)
  deallocate(funl1)
  deallocate(fundl0)
  deallocate(fundl1)

  deallocate(intw1)
  deallocate(intw2)

  deallocate(kvall)
  deallocate(kfunl)
  deallocate(kfundl)
  deallocate(kintl)
  deallocate(kcoefl)

  if (ikind.ne.0) then 
    deallocate(kvalr)
    deallocate(kfunr)
    deallocate(kfundr)
    deallocate(kintr)
    deallocate(kcoefr)
  endif
!
! From init_gper
!
  deallocate(gper)
  deallocate(ninsh)
  deallocate(gnsh)

  return
end subroutine free_mem
