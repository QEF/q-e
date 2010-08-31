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
  use cond
  implicit none
!
! From allocate_cond
!
  deallocate(psiperl)
  deallocate(zkl)
  deallocate(zkrl)
  deallocate(psipers)
  deallocate(zks)
  deallocate(zkrs)
  deallocate(psiperr)
  deallocate(zkr)
  deallocate(zkrr)

  deallocate(newbg)

  deallocate(fun0)
  deallocate(fun1)
  deallocate(fund0)
  deallocate(fund1)

  IF (lorb) THEN
    deallocate( funz0 )
    deallocate( korbl )
    deallocate( korbr )
  ENDIF

  IF (norbf>0) THEN
     deallocate(funl0)
     deallocate(funl1)
     deallocate(fundl0)
     deallocate(fundl1)

     deallocate(intw1)
     deallocate(intw2)
  END IF

  deallocate(kvall)
  deallocate(kfunl)
  deallocate(kfundl)
  IF (nocrosl>0) THEN
     deallocate(kintl)
     deallocate(kcoefl)
  END IF

  if (ikind.ne.0) then
    deallocate(kvalr)
    deallocate(kfunr)
    deallocate(kfundr)
    IF (nocrosr>0) THEN
       deallocate(kintr)
       deallocate(kcoefr)
    END IF
  endif
!
! From init_gper
!
  if (lorb) deallocate( nl_2ds )
  if (lorb) deallocate( nl_2d )
  deallocate(gper)
  deallocate(ninsh)
  deallocate(gnsh)

  return
end subroutine free_mem
