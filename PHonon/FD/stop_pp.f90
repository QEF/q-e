!
! Copyright (C) 2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
  use control_flags, only: twfcollect
  use io_files, only: iunwfc
  use mp, only: mp_end, mp_barrier
  USE parallel_include
#ifdef __PARA

  integer :: info
  logical :: op

  inquire ( iunwfc, opened = op )

  if ( op ) then
     if (twfcollect) then
        close (unit = iunwfc, status = 'delete')
     else
        close (unit = iunwfc, status = 'keep')
     end if
  end if 

  call mp_barrier()

  ! call mpi_finalize (info)
#endif
 
  call mp_end()

#ifdef __T3E
  !
  ! set streambuffers off
  !

  call set_d_stream (0)
#endif

  stop
end subroutine stop_pp
