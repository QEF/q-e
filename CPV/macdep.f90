!
! Copyright (C) 2002-2003 CP group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Machine-dependent routines for:
!    cpu-time and ram measurement
!    random number generation
!    fft dimensions
!
! 1) cpu-time and ram measurement
! =======================
!
!-------------------------------------------------------------------------
      subroutine tictac(i,j)
!-------------------------------------------------------------------------
!
!     measure cpu time, elapsed time, and number of calls
!     i   index of subroutine 
!     j   j=0 start, j=1 stop the clock and add 1 to number of calls
!
!     WARNING: the accuracy of time-measuring routines may vary
!              the resolution for elapsed time may be as low as 1 s
!
      use timex_mod

      implicit none
      integer i, j
      real(kind=8) time1(maxclock),time2(maxclock), tcpu, telaps
      integer count, count_rate, count_max
      integer k
      logical first
      save first, time1, time2
      data first/.true./
!
!    The machine dependent declaration of the timing function
!
#ifdef __CRAYY
      real*8                                                            &
     &       second      ! system function, returns the CPU time in sec.
#endif
#ifdef __AIX
      integer                                                           &
     &       mclock      ! system function, returns the time in sec./100.
      real*8                                                            &
     &       timef       ! system function, returns elapsed time in msec.
#endif
#ifdef __T3E
      real*8                                                            &
     &       tsecnd,    &! system function, returns the CPU time in sec.
     &       timef       ! system function, returns elapsed time in msec.
#endif
#if defined(__PGI) || defined(__INTEL)
      real*4                                                            &
     &       etime,    &! system function, returns the CPU time in sec.
     &       tarry(2)    ! user and system times (not used)
#endif
!
!
      if (i.lt.1.or.i.gt.maxclock)                                      &
     &     call errore('tictac','wrong number of clocks',i)
      if (j.ne.0.and.j.ne.1) call errore('tictac','wrong call',j)
!
! initialization
!
      if (first) then
         do k=1,maxclock
            cputime(k)=0.0
            elapsed(k)=0.0
            ntimes(k)=0
         end do
         first=.false.
      end if
      count_max =0
      count_rate=1
!
! call cpu- and elapsed-time machine-specific routines
!
#ifdef __CRAYY
#define __USE_SYSTEM_TIMER
      tcpu = second()
      call system_clock(count, count_rate, count_max) 
      telaps=float(count)/count_rate
#endif
#ifdef __SX4 
#define __USE_SYSTEM_TIMER
      call clock(tcpu)
      call system_clock(count, count_rate, count_max)
      telaps=float(count)/count_rate
#endif
#ifdef __AIX
#define __USE_SYSTEM_TIMER
      tcpu = mclock() / 100.d0
      telaps=timef()/1000.
#endif
#ifdef __T3E
#define __USE_SYSTEM_TIMER
      tcpu = tsecnd()
      telaps=timef()/1000.
#endif
#if defined(__PGI) || defined(__INTEL)
#define __USE_SYSTEM_TIMER
      tcpu = etime( tarry )
      call system_clock(count, count_rate, count_max)
      telaps=float(count)/count_rate
#endif
#ifndef __USE_SYSTEM_TIMER
!
! call intrinsic f90 cpu- and elapsed-time routines
!
      call cpu_time (tcpu)
      call system_clock(count, count_rate, count_max)
      telaps=float(count)/count_rate
#endif
      if (j.eq.0) then
         time1(i)=tcpu
         time2(i)=telaps
      else if (j.eq.1) then
         cputime(i)=cputime(i) + (  tcpu-time1(i))
         elapsed(i)=elapsed(i) + (telaps-time2(i))
!
! The following workaround is needed because system_clock resets "count"
! every time it reaches "count_max" (usually the largest integer), and
! on some compilers (Absoft, Intel) this happens way too frequently. 
! Will not work if elapsed t between two calls > count_max/count_rate .
!
         if (telaps-time2(i).lt.0.d0)                                   &
     &       elapsed(i)=elapsed(i) + (count_max+1.d0)/count_rate
! BEWARE: (count_max+1) may give integer overflow !!!
         ntimes(i) =ntimes(i)+1
      endif
      return
      end
!
!-----------------------------------------------------------------------
      subroutine memory
!-----------------------------------------------------------------------
!
! Prints what is hopefully the size of occupied memory
! Implemented only for SGI Origin and AIX SP3.
! Extremely machine- and operating-system dependent
!
#ifdef __PARA
      use para_mod, only: me
#endif
      implicit none
      character(len=80) command
      integer pid
#ifdef __AIX
      integer getpid_

      pid=getpid_()
      write(command,10) pid
 10   format('ps -lp ',i8,' | grep -v SZ | awk ''{print $10}'' ')
      write(6,'(''Estimated size (kB) of each process: '',$)')
      call system(command)
#endif

#ifdef __ORIGIN
      integer getpid
      pid=getpid()
      write(command,10) pid
 10   format('ps -lp ',i8,'|grep -v SZ|awk ''{print $10}''|cut -f1 -d:')
      write(6,'(''Total estimated size (pages) of each process: '',$)')
#ifdef __PARA
      if(me.eq.1) &
#endif
      call system(command)
#endif
      return
      end
!
! 2) random number generation
! ===========================
!
!-------------------------------------------------------------------------
      real(kind=8) function randy()
!-------------------------------------------------------------------------
!
! Use machine-specific random-number generator when available
!
#ifdef __CRAYY
#define __USE_SYSTEM_RAND
      randy = ranf()
#endif
#ifdef __SX4
#define __USE_SYSTEM_RAND
      randy=random(0)
#endif
#ifdef __AIX
#define __USE_SYSTEM_RAND
      randy=rand()
#endif
!
! Use fortran random-number generator in all other cases
!
#ifndef __USE_SYSTEM_RAND
      integer m, ia, ic, ntab
      real(kind=8) rm
      parameter (ntab=97,m=714025,ia=1366,ic=150889,rm=1.0/m)
      integer ir(ntab), iff, idum, j, iy
      data iff /0/, idum/0/
      save iff, idum, iy, ir
!
!
      if(idum.lt.0.or.iff.eq.0) then
        iff=1
        idum=mod(ic-idum,m)
        do j=1,ntab
           idum=mod(ia*idum+ic,m)
           ir(j)=idum
        end do
        idum=mod(ia*idum+ic,m)
        iy=idum
      endif
      j=1+(ntab*iy)/m
      if(j.gt.ntab.or.j.lt.1) call errore('randy','j out of range',j)
      iy=ir(j)
      randy=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
#endif
      return
      end
!
