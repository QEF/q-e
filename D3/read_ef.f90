!
! Copyright (C) 2001-2003 PWSCF group
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
#include "f_defs.h"
  use pwcom
  use d3com
#ifdef __PARA
  use para
  use mp, only : mp_bcast
#endif
  implicit none
  integer :: root = 0, ios
  !
  if (degauss.eq.0.d0) return
#ifdef __PARA
  if (me.ne.1.or.mypool.ne.1) goto 210
#endif
  rewind (unit = iuef)
  read (iuef, err = 100, iostat = ios) ef_sh

100 call errore ('d3_valence', 'reading iuef', abs (ios) )
#ifdef __PARA
210 continue
  call mp_bcast (ef_sh, root)
#endif
  return
end subroutine read_ef
