!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine read_ef  
  !-----------------------------------------------------------------------
  ! Reads the shift of the Fermi Energy
  !
#include "machine.h"
  use pwcom
  use d3com
#ifdef PARA
  use para
#endif
  implicit none

#ifdef PARA
  include 'mpif.h'  
  integer :: root, errcode, nat_3  
#endif
  integer :: ios  

  if (degauss.eq.0.d0) return  
#ifdef PARA
  if (me.ne.1.or.mypool.ne.1) goto 210  
#endif
  rewind (unit = iuef)  
  read (iuef, err = 100, iostat = ios) ef_sh  

100 call error ('d3_valence', 'reading iuef', abs (ios) )  
#ifdef PARA
210 continue  
  nat_3 = 3 * nat  
  root = 0  
  call MPI_bcast (ef_sh, nat_3, MPI_REAL8, root, MPI_COMM_WORLD, &
       errcode)

  call error ('read_ef', 'at bcast', errcode)  
#endif
  return  
end subroutine read_ef
