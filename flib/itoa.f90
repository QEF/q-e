!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! trasforma un numero intero in una stringa

       subroutine cpitoa(n,str)
       integer, intent(in) :: n
       character(LEN=*) str
       integer i, npow, j, nq, ntmp
       logical :: lzero 

       j     = 1
       ntmp  = n
       lzero = .FALSE.
       do i = 9, 0, -1
         npow = 10**i
         nq   = ntmp / npow
         ntmp = MOD(ntmp,npow)
         if( lzero .or. (nq.ne.0)) then
           lzero = .TRUE.
           select case(nq)
             case (9) 
               str(j:j) = '9' 
             case (8) 
               str(j:j) = '8' 
             case (7) 
               str(j:j) = '7' 
             case (6) 
               str(j:j) = '6' 
             case (5) 
               str(j:j) = '5' 
             case (4) 
               str(j:j) = '4' 
             case (3) 
               str(j:j) = '3' 
             case (2) 
               str(j:j) = '2' 
             case (1) 
               str(j:j) = '1' 
             case default
               str(j:j) = '0' 
           end select
           j = j + 1
         end if
       end do
       str(j:j) = ' '
       return
       end subroutine


       subroutine unitname(n,name)
       character(len=*) :: name
       character(len=5) :: num
       integer :: l
       call cpitoa(n,num)
       l    = index(num,' ') - 1
       name = 'fort.'//num(1:l)
       return
       end

       subroutine myunitname(me,name)
       character(len=*) :: name
       character(len=5) :: num
       integer :: l
       call cpitoa(me,num)
       l    = index(num,' ') - 1
       name = 'fort_6.'//num(1:l)
       return
       end

