!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program pwscf  
  !-----------------------------------------------------------------------
  !
  !     Plane Wave Self-Consistent Field                                 c
  !
  use pwcom  
  use io
  implicit none
  character :: cdate * 9, ctime * 9, version * 12  
  external date_and_tim  
  !      call sigcatch ( )
  ! use ".false." to disable all clocks except the total cpu time clock
  ! use ".true."  to enable clocks
  !      call init_clocks(.false.)

  call init_clocks (.true.)  
  call start_clock ('PWSCF')  
  version = 'PWSCF 1.2.b'  
#ifdef PARA
  call startup (nd_nmbr, version)  
#else
  nd_nmbr = '   '  
  call date_and_tim (cdate, ctime)  
  write (6, 9000) version, cdate, ctime  
9000 format (/5x,'Program ',a12,' starts ...',/5x, &
       &            'Today is ',a9,' at ',a9)
#endif
  call init_run  
  istep = 0  
  do while (istep.lt.nstep)  
     istep = istep + 1  
     call electrons  
     if (.not.conv_elec) call stop_pw (conv_elec)  
     call ions  
     if (conv_ions) goto 10  
     call hinit1  

  enddo

10 call punch  

  call stop_pw (conv_ions)  
  stop  
end program pwscf


