!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
!
! Author: L. Martin-Samos
!
!--------------------------------------------------------------------
subroutine stop_pp
  !--------------------------------------------------------------------
  !
  ! Synchronize processes before stopping.
  !
  use io_files, only: iunwfc
  use mp_global, only: mp_global_end
  USE parallel_include
#if defined(__MPI)

  integer :: info
  logical :: op

  inquire ( iunwfc, opened = op )

  if ( op ) close (unit = iunwfc, status = 'delete')

  call mp_global_end ( )

#endif

  stop
end subroutine stop_pp
