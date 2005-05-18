!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      subroutine memory
!-----------------------------------------------------------------------
!
! Prints what is hopefully the size of occupied memory
! Implemented only for SGI Origin and AIX SP3.
! Extremely machine- and operating-system dependent
!
      use io_global, only: stdout, ionode
      implicit none
      character(len=256) :: command
      integer :: pid

#ifdef __AIX
      integer getpid_

      pid=getpid_()
      write(command,10) pid
 10   format('ps -lp ',i8,' | grep -v SZ | awk ''{print $10}'' ')
      WRITE( stdout,'(''Estimated size (kB) of each process: '',$)')
      if( ionode ) call system(command)
#endif

#ifdef __ORIGIN
      integer getpid

      pid=getpid()
      write(command,10) pid
 10   format('ps -lp ',i8,'|grep -v SZ|awk ''{print $10}''|cut -f1 -d:')
      WRITE( stdout,'(''Total estimated size (pages) of each process: '',$)')
      if( ionode ) call system( command )
#endif

      return
      end subroutine memory
