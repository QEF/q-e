!
! Copyright (C) 2002 CP90 group
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
!    erf, erfc, freq functions
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
     &     call error('tictac','wrong number of clocks',i)
      if (j.ne.0.and.j.ne.1) call error('tictac','wrong call',j)
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
#ifdef __NEC 
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
#ifdef __NEC
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
      if(j.gt.ntab.or.j.lt.1) call error('randy','j out of range',j)
      iy=ir(j)
      randy=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
#endif
      return
      end
!
! 3) utilities for fft dimensions
! ================
!
      integer function good_fft_dimension(n)
!
! Determines the optimal maximum dimensions of fft arrays
! Useful on some machines to avoid memory conflicts
!
      integer n, nx
!
! this is the default: max dimension = fft dimension
      nx=n
#if defined(__ESSL)
      if ( n.eq. 8 .or. n.eq.16  .or. n.eq.32 .or.                      &
     &     n.eq.64 .or. n.eq.128 .or. n.eq.256       )  nx=n+1
#endif
#if defined(__CRAYY) || defined(__NEC) 
      if ( mod(n,2).eq.0)  nx=n+1
#endif
      good_fft_dimension=nx
      return
      end
!
!-----------------------------------------------------------------------
      integer function good_fft_order(nr)
!-----------------------------------------------------------------------
!
! Input : tentative order n of a fft 
! Output: the same if n is a good number
!         the closest higher number that is good
! an fft order is not good if not implemented (as on IBM with ESSL)
! or implemented but with awful performances (most other cases)
!
      implicit none
      integer nr
!
      integer factors(5), pwr(5), mr, i, fac, p, maxpwr, maxn
      parameter (maxn=1000)
      logical good
      data    factors /2, 3, 5, 7, 11/
!
! find the factors of the fft dimension
!
 10   mr=nr
      do i=1,5
         pwr(i)=0
      end do
      do i=1,5
         fac=factors(i)
         maxpwr = nint(log(float(mr))/log(float(fac)))+1
         do p=1,maxpwr
            if (mr.eq.1) goto 20
            if (mod (mr,fac).eq.0) then
               mr=mr/fac
               pwr(i)=pwr(i)+1
            end if
         end do
      end do
!
 20   if (nr .ne. mr * 2**pwr(1) * 3**pwr(2) * 5**pwr(3) *              &
     &                 7**pwr(4) *11**pwr(5) )                          &
     &   call error('good_fft_order','what ?!?',1)
      if (mr.ne.1) then
! fft dimension contains factors > 11 : no good in any case
         good=.false.
      else
! specific (machine- and library-dependent cases)
#ifdef __ESSL
!
! IBM machines with essl libraries
!
         good=pwr(1).ge.1 .and.                                         &
     &           pwr(2).le.2 .and.                                      &
     &           pwr(3).le.1 .and.                                      &
     &           pwr(4).le.1 .and.                                      &
     &           pwr(5).le.1 .and.                                      &
     &           ((pwr(2).eq.0 .and. pwr(3)+pwr(4)+pwr(5).le.2) .or.    &
     &            (pwr(2).ne.0 .and. pwr(3)+pwr(4)+pwr(5).le.1)     )
#endif
!
#if defined(__CRAYY) || defined(__NEC)
!
! Cray and t3d machines with scilib libraries
!
         good=pwr(4).eq.0 .and. pwr(5).eq.0
#endif
!
#ifdef __FFTW
         good=pwr(5).eq.0
#endif         
      end if
      if (.not.good) then
         nr=nr+1
         if (nr.gt.maxn)                                                &
     &        call error('good_fft_order','too large',maxn)
         go to 10
      else
         good_fft_order=nr
      end if
!
      end
!
! 4) erf, erfc, freq functions
! ================
!     for machines that do not have these routines in the math libraries
!
#if defined __INTEL || defined __PGI
!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
real(kind=8) function erf (x)  
  !---------------------------------------------------------------------
  !
  !     Error function - computed from the rational approximations of
  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
  !
  !     for abs(x) le 0.47 erf is calculated directly
  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
  !
  implicit none  
  real(kind=8) :: x, x2, p1 (4), q1 (4), erfc  
  external erfc  
  data p1 / 2.42667955230532d2, 2.19792616182942d1, &
       6.99638348861914d0, - 3.56098437018154d-2 /
  data q1 / 2.15058875869861d2, 9.11649054045149d1, &
       1.50827976304078d1, 1.00000000000000d0 /
  !
  if (abs (x) .gt.6.d0) then  
     !
     !  erf(6)=1-10^(-17) cannot be distinguished from 1 with 16-byte words
     !
     erf = sign (1.d0, x)  
  else  
     if (abs (x) .le.0.47d0) then  
        x2 = x**2  
        erf = x * (p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 ( &
             4) ) ) ) / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 ( &
             4) ) ) )
     else  
        erf = 1.d0 - erfc (x)  
     endif
  endif
  !
  return  
end function erf
!
!---------------------------------------------------------------------
real(kind=8) function erfc (x)  
  !---------------------------------------------------------------------
  !
  !     erfc(x) = 1-erf(x)  - See comments in erf
  !
  implicit none  
  real(kind=8) :: x, ax, x2, xm2, erf, p2 (8), q2 (8), p3 (5), q3 (5), &
       pim1
  external erf  
  data p2 / 3.00459261020162d2, 4.51918953711873d2, &
       3.39320816734344d2, 1.52989285046940d2, 4.31622272220567d1, &
       7.21175825088309d0, 5.64195517478974d-1, - 1.36864857382717d-7 /
  data q2 / 3.00459260956983d2, 7.90950925327898d2, &
       9.31354094850610d2, 6.38980264465631d2, 2.77585444743988d2, &
       7.70001529352295d1, 1.27827273196294d1, 1.00000000000000d0 /
  data p3 / - 2.99610707703542d-3, - 4.94730910623251d-2, - &
       2.26956593539687d-1, - 2.78661308609648d-1, - 2.23192459734185d-2 &
       /
  data q3 / 1.06209230528468d-2, 1.91308926107830d-1, &
       1.05167510706793d0, 1.98733201817135d0, 1.00000000000000d0 /

  data pim1 / 0.564189583547756d0 /  
  !        ( pim1= sqrt(1/pi) )
  ax = abs (x)  
  if (ax.gt.26.d0) then  
     !
     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
     !
     erfc = 0.d0  
  elseif (ax.gt.4.d0) then  
     x2 = x**2  
     xm2 = (1.d0 / ax) **2  
     erfc = (1.d0 / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
          (q3 (4) + xm2 * q3 (5) ) ) ) ) )
  elseif (ax.gt.0.47d0) then  
     x2 = x**2  
     erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
  else  
     erfc = 1.d0 - erf (ax)  
  endif
  !
  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
  !
  if (x.lt.0.d0) erfc = 2.d0 - erfc  
  !
  return  
end function erfc
!---------------------------------------------------------------------
real(kind=8) function freq (x)  
  !---------------------------------------------------------------------
  !
  !     freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
  !             - See comments in erf
  !
  real(kind=8) :: x, c, erf, erfc  
  external erf  

  data c / 0.707106781186548d0 /  
  !        ( c= sqrt(1/2) )
  freq = 0.5d0 * erfc ( - x * c)  
  !
  return  
end function freq

#endif
