!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
module mytime
  integer :: maxclock
  real (8) :: notrunning
  parameter (maxclock = 100, notrunning = - 1.d0)
  real (8) :: myclock (maxclock), t0 (maxclock)!, scnds
  character (len=12) :: clock_label
  dimension clock_label (maxclock)
  integer :: called (maxclock), nclock
  logical :: no
end module mytime

subroutine init_clocks (go)
  !
  ! flag=.true.  : clocks will run
  ! flag=.false. : only clock #1 will run
  !

  use parameters, only : DP
  use mytime

  implicit none

  real(kind=8) :: scnds
  logical :: go
  integer :: n
  no = .not.go
  do n = 1, maxclock
     called (n) = 0
     myclock (n) = 0.0
     t0 (n) = notrunning

  enddo
  return

end subroutine init_clocks
!
subroutine start_clock (label)
  !
  use mytime
  implicit none

  real(kind=8) :: scnds
  character (len=*) :: label
  integer :: n
  if (no.and.nclock.eq.1) return
  do n = 1, nclock
     if (label.eq.clock_label (n) ) then
        !
        ! found previously defined clock: check if not already started,
        ! store in t0 the starting time
        !
        if (t0 (n) .ne.notrunning) then
           write (6, '("start_clock: clock # ",i2," for ",a12, &
                &            " already started")') n, label
        else
           t0 (n) = scnds ()
        endif
        return

     endif
  enddo
  !
  ! clock not found : add new clock for given label
  !
  if (nclock.eq.maxclock) then
     write (6, '("start_clock: Too many clocks! call ignored")')
  else
     nclock = nclock + 1
     clock_label (nclock) = label
     t0 (nclock) = scnds ()

  endif
  return

end subroutine start_clock
!
subroutine stop_clock (label)
  !
  use mytime
  implicit none

  real(kind=8) :: scnds
  character (len=*) :: label
  integer :: n

  if (no) return
  do n = 1, nclock
     if (label.eq.clock_label (n) ) then
        !
        ! found previously defined clock : check if properly initialised,
        ! add elapsed time, increase the counter of calls
        !
        if (t0 (n) .eq.notrunning) then
           write (6, '("stop_clock: clock # ",i2," for ",a12, &
                &            " not running")') n, label
        else
           myclock (n) = myclock (n) + scnds () - t0 (n)
           t0 (n) = notrunning
           called (n) = called (n) + 1
        endif
        return
     endif
  enddo
  !
  ! clock not found
  !

  write (6, '("stop_clock: no clock for ",a12," found !")') label
  return

end subroutine stop_clock
!

subroutine print_clock (label)
  !
  use mytime
  implicit none

  real(kind=8) :: scnds
  character (len=*) :: label
  integer :: n

  if (label.eq.' ') then
     write (6, * )
     do n = 1, nclock
        call print_this_clock (n)
     enddo
  else
     do n = 1, nclock
        if (label.eq.clock_label (n) ) then
           call print_this_clock (n)
           return
        endif

     enddo
     !
     ! clock not found
     !         if (.not.no) write(6,'("print_clock: no clock for ",
     !     +                            a12," found !")')   label

  endif
  return

end subroutine print_clock
!
subroutine print_this_clock (n)
  !
  use parameters
  use mytime
  implicit none

  real(kind=8) :: scnds
  integer :: n
  real(kind=DP) :: elapsed_cpu_time, nsec

  integer :: nhour, nmin
  if (t0 (n) .eq.notrunning) then
     !
     ! clock stopped, print the stored value for the cpu time
     !
     elapsed_cpu_time = myclock (n)
  else
     !
     ! clock not stopped, print the current value of the cpu time
     !
     elapsed_cpu_time = myclock (n) + scnds () - t0 (n)

  endif
#ifdef __PARA
  !
  ! In the parallel case it is far from clear which value to print
  ! The following is the maximum over all nodes and pools. NOTA BENE:
  ! some trouble could arise if a clock is not started on all nodes
  !
  ! by uncommenting the following line the extreme operation is removed
  ! may be useful for testing purpouses
  ! /* #define DEBUG */
  !
#ifndef DEBUG
  call extreme (elapsed_cpu_time, + 1)

  call poolextreme (elapsed_cpu_time, + 1)
#endif
#endif
  if (n.eq.1) then
     ! The first clock is written as hour/min/sec
     nhour = elapsed_cpu_time / 3600
     nmin = (elapsed_cpu_time-3600 * nhour) / 60
     nsec = (elapsed_cpu_time-3600 * nhour) - 60 * nmin
     if (nhour.gt.0) then
        write (6, '(5x,a12," : ",3x,i2,"h",i2,"m CPU time"/)') &
             clock_label (n) , nhour, nmin
     elseif (nmin.gt.0) then
        write (6, '(5x,a12," : ",i2,"m",f5.2,"s CPU time"/)') &
             clock_label (n) , nmin, nsec
     else
        write (6, '(5x,a12," : ",3x,f5.2,"s CPU time"/)') &
             clock_label (n) , nsec
     endif
  elseif (called (n) .eq.1.or.t0 (n) .ne.notrunning) then
     ! For clocks that have been called only once
     write (6, '(5x,a12," :",f9.2,"s CPU")') clock_label (n) , &
          elapsed_cpu_time
  elseif (called (n) .eq.0) then
     ! For clocks that have never been called
     write (6, '("print_this: clock # ",i2," for ",a12, &
          &       " never called !")') n, clock_label(n)
  else
     ! For all other clocks
     write (6, '(5x,a12," :",f9.2,"s CPU (", &
          &       i8," calls,",f8.3," s avg)")') clock_label (n) , &
          elapsed_cpu_time, called (n) , elapsed_cpu_time / called (n)

  endif
  return
end subroutine print_this_clock
!
real(kind=8) function get_clock (label)
  !
  use parameters
  use mytime

  implicit none

  real(kind=8) :: scnds
  character (len=*) :: label
  integer :: n

  if (no) then
     if (label.eq.clock_label (1) ) then
        get_clock = scnds ()
     else
        get_clock = notrunning
     endif
     return
  endif
  do n = 1, nclock
     if (label.eq.clock_label (n) ) then
        if (t0 (n) .eq.notrunning) then
           get_clock = myclock (n)
        else
           get_clock = myclock (n) + scnds () - t0 (n)
        endif
#ifdef __PARA
        !
        ! In the parallel case, use the maximum over all nodes and pools
        !
        call extreme (get_clock, + 1)
        call poolextreme (get_clock, + 1)
#endif
        return
     endif
  enddo
  !
  ! clock not found
  !
  get_clock = notrunning

  write (6, '("get_clock: no clock for ",a12," found !")') label
  return
end function get_clock

