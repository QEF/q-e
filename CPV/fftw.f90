!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "../include/machine.h"

#if defined(__FFTW)
module fftw_mod

  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)     
  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

end module fftw_mod
!
! Driver routines for fft using fftw libraries (PG). Note that:
!
!     isign=-1 => FFTW_FORWARD
!     isign= 1 => FFTW_BACKWARD
!
!     a separate initialization ("plan") is needed for each direction
!
!     the one-dimensional fortran routines need an output array
!     even if the fft is done in-place!
!
!
! The variables "plan" contain C-style pointers. For Origin machines
! compiled with "-64" and for Compaq Alpha machines, 64-bit integers
! must be used
!
#if defined(__SGI64) || defined(__COMPAQ)
#define PLAN_C_POINTER integer(kind=8)
#else
#define PLAN_C_POINTER integer
#endif
!
!-------------------------------------------------------------------------
      subroutine cfft3(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,isign)
!-----------------------------------------------------------------------
!     driver routine for 3d fft using fftw libraries (PG)
!
      use fftw_mod
      implicit none
!
      integer nr1, nr2, nr3, nr1x, nr2x, nr3x, isign
      complex(kind=8) f(nr1x*nr2x*nr3x)
!
      real(kind=8) fac
      integer ibid
!
! initialization variables
!
      PLAN_C_POINTER plan(2)
      save plan
      data plan/0,0/
!
!
      if (nr1.ne.nr1x .or. nr2.ne.nr2x .or. nr3.ne.nr3x )               &
     &     call errore('cfft3','not implemented',1)
      if (isign.eq.1) then
         ibid=1
      else if (isign.eq.-1) then
         ibid=2
      else 
         call errore('cfft3','isign unexpected',isign)
      end if
!
      if (plan(ibid).eq.0)                                              &
     &     call FFTW3D_F77_CREATE_PLAN                                  &
     &     (plan(ibid),nr1,nr2,nr3,isign,FFTW_ESTIMATE+FFTW_IN_PLACE)
!
      call FFTWND_F77_ONE(plan(ibid), f, 0)
!
      if (isign.eq.-1) then
         fac=1.0/float(nr1*nr2*nr3)
         call DSCAL(2*nr1*nr2*nr3, fac, f, 1)
      end if
!
      return
      end
!
!-------------------------------------------------------------------------
      subroutine cfft3s(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,isign)
!-----------------------------------------------------------------------
!     driver routine for 3d fft using fftw libraries (PG)
!
      use fftw_mod
      implicit none
!
      integer nr1, nr2, nr3, nr1x, nr2x, nr3x, isign
      complex(kind=8) f(nr1x*nr2x*nr3x)
!
      real(kind=8) fac
      integer ibid
!
! initialization variables
!
      PLAN_C_POINTER plan(2)
      save plan
      data plan/0,0/
!
!
      if (nr1.ne.nr1x .or. nr2.ne.nr2x .or. nr3.ne.nr3x )               &
     &     call errore('cfft3s','not implemented',1)
      if (isign.eq.1) then
         ibid=1
      else if (isign.eq.-1) then
         ibid=2
      else 
         call errore('cfft3s','isign unexpected',isign)
      end if
!
      if (plan(ibid).eq.0)                                              &
     &     call FFTW3D_F77_CREATE_PLAN                                  &
     &        (plan(ibid),nr1,nr2,nr3,isign,FFTW_ESTIMATE+FFTW_IN_PLACE)
!
      call FFTWND_F77_ONE(plan(ibid), f, 0)
!
      if (isign.eq.-1) then
         fac=1.0/dfloat(nr1*nr2*nr3)
         call DSCAL(2*nr1*nr2*nr3, fac, f, 1)
      end if
!
      return
      end
!
!-------------------------------------------------------------------------
      subroutine cfft3b(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,isign)
!-----------------------------------------------------------------------
!     driver routine for 3d fft using fftw libraries (PG)
!
      use fftw_mod
      implicit none
!
      integer nr1, nr2, nr3, nr1x, nr2x, nr3x, isign
      complex(kind=8) f(nr1x*nr2x*nr3x)
!
      real(kind=8) fac
      integer ibid
!
! initialization variables
!
      PLAN_C_POINTER plan(2)
      save plan
      data plan/0,0/
!
!
      if (nr1.ne.nr1x .or. nr2.ne.nr2x .or. nr3.ne.nr3x )               &
     &     call errore('cfft3b','not implemented',1)
      if (isign.eq.1) then
         ibid=1
      else if (isign.eq.-1) then
         ibid=2
      else 
         call errore('cfft3b','isign unexpected',isign)
      end if
!
      if (plan(ibid).eq.0)                                              &
     &     call FFTW3D_F77_CREATE_PLAN                                  &
     &        (plan(ibid),nr1,nr2,nr3,isign,FFTW_ESTIMATE+FFTW_IN_PLACE)
!     
      call FFTWND_F77_ONE(plan(ibid), f, 0)
!
      if (isign.eq.-1) then
         fac=1.0/dfloat(nr1*nr2*nr3)
         call DSCAL(2*nr1*nr2*nr3, fac, f, 1)
      end if
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_1     (f,m,n,nx,isign,fout)
!     ===============
!     driver routine for m 1d complex fft's (dense grid) - fftw
!     NOTA BENE: not in-place! output in fout
!----------------------------------------------------------------------
!
      use fftw_mod
      implicit none
      integer m, n, nx, isign
      complex(kind=8)  f(nx*m), fout(nx*m)
#ifdef __PARA
      real(kind=8) fac
      integer ibid
!
! initialization variables
!
      PLAN_C_POINTER plan(2)
      save plan
      data plan/0,0/
!
!
      if (isign.eq.1) then
         ibid =1
      else if (isign.eq.-1) then
         ibid =2
      else
         call errore('cft_1','wrong call',isign)
      end if
!
      if (plan(ibid).eq.0)                                              &
     &     call FFTW_F77_CREATE_PLAN                                    &
     &        (plan(ibid),n,isign,FFTW_ESTIMATE)
!
      call FFTW_F77 (plan(ibid),m,f,1,nx,fout,1,nx)
!
      if (isign.eq.-1) then
         fac=1.0/float(n)
         call DSCAL(2*nx*m, fac, fout, 1)
      end if
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_2     (f,mplane,n1,n2,nx1,nx2,isign)
!     ===============
!     driver routine for mplane 2d complex fft's of lengths n1 and n2
!     (dense grid) - fftw
!----------------------------------------------------------------------
!
      use fftw_mod
      implicit none
      integer n1, n2, mplane, nx1, nx2, isign
      complex(kind=8)  f(nx1*nx2*mplane)
#ifdef __PARA
!
      real(kind=8) fac
      integer ibid
!
! initialization variables
!
      PLAN_C_POINTER plan(2)
      save plan
      data plan/0,0/
!
!
      if (n1.ne.nx1 .or. n2.ne.nx2)                                     &
     &     call errore('cft_2','not implemented',1)
      if (isign.eq.1) then
         ibid =1
      else if (isign.eq.-1) then
         ibid =2
      else
         call errore('cft_2','wrong call',isign)
      end if
!
      if (plan(ibid).eq.0)                                              &
     &     call FFTW2D_F77_CREATE_PLAN                                  &
     &        (plan(ibid),n1,n2,isign,FFTW_ESTIMATE+FFTW_IN_PLACE)
!
      call FFTWND_F77 (plan(ibid),mplane,f,1,nx1*nx2,f,1,nx1*nx2)
!
      if (isign.eq.-1) then
         fac=1.0/float(n1*n2)
         call DSCAL(2*nx1*nx2*mplane, fac, f, 1)
      end if
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_1s     (f,m,n,nx,isign,fout)
!     ===============
!     driver routine for m 1d complex fft's (sparse grid) - fftw
!     NOTA BENE: not in-place! output in fout
!----------------------------------------------------------------------
!
      use fftw_mod
      implicit none
      integer m, n, nx, isign
      complex(kind=8)  f(nx*m), fout(nx*m)
#ifdef __PARA
      real(kind=8) fac
      integer ibid, isign1
!
! initialization variables
!
      PLAN_C_POINTER plan(2)
      save plan
      data plan/0,0/
!
!
      if (isign.eq.1.or.isign.eq.2) then
         isign1=1
         ibid  =1
      else if (isign.eq.-1.or.isign.eq.-2) then 
         isign1=-1
         ibid  =2
      else
         call errore('cft_1s','wrong call',isign)
      end if
!
      if (plan(ibid).eq.0)                                              &
     &     call FFTW_F77_CREATE_PLAN                                    &
     &        (plan(ibid),n,isign1,FFTW_ESTIMATE)
!
      call FFTW_F77 (plan(ibid),m,f,1,nx,fout,1,nx)
!
      if (isign.eq.-1.or.isign.eq.-2) then
         fac=1.0/float(n)
         call DSCAL(2*nx*m, fac, fout, 1)
      end if
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_2s     (f,mplane,n1,n2,nx1,nx2,isign,planes)
!     ===============
!     driver routine for mplane 2d complex fft's of lengths n1 and n2
!     (sparse and wavefunction grid) - fftw
!----------------------------------------------------------------------
!
      use fftw_mod
      implicit none
      integer n1, n2, mplane, nx1, nx2, isign, planes(nx1)
      complex(kind=8)  f(nx1*nx2*mplane)
#ifdef __PARA
      integer nmax
      parameter (nmax=256)
      complex(kind=8) fout(nmax)
      real(kind=8) fac
      integer ibid, isign1, i, k, m, istrt
!
! initialization variables
!
      PLAN_C_POINTER plan1(2), plan2(2)
      save plan1, plan2
      data plan1/0,0/, plan2/0,0/
!
!
      if (n1.gt.nmax.or.n2.gt.nmax)                                     &
     &     call errore('cft_2s','increase nmax',max(n1,n2))
      if (n1.ne.nx1 .or. n2.ne.nx2)                                     &
     &     call errore('cft_2s','not implemented',1)
      if (isign.eq.1.or.isign.eq.2) then
         isign1=1
         ibid  =1
      else if (isign.eq.-1.or.isign.eq.-2) then
         isign1=-1
         ibid  =2
      else
         call errore('cft_2s','wrong call',isign)
      end if
! check how many columns along x are nonzero
      m=0
      do i=1,n1
         m=m+planes(i)
      end do
      if (m.gt.n1.or.m.le.0)                                            &
     &     call errore('cft_2s','something wrong with planes',1)
!
      if (isign1.eq.1) then
! j-direction
         if (plan2(ibid).eq.0)                                          &
     &        call FFTW_F77_CREATE_PLAN                                 &
     &        (plan2(ibid),n2,isign1,FFTW_ESTIMATE+FFTW_IN_PLACE)
!
         do i=1,n1
!
! do only ffts on columns (i,*,k) resulting in nonzero components
!
            if (planes(i).eq.1) then
               call FFTW_F77 &
     &              (plan2(ibid),mplane,f(i),nx1,nx1*nx2,fout,0,0)
            end if
         end do
! i-direction
         if (plan1(ibid).eq.0)                                          &
     &        call FFTW_F77_CREATE_PLAN                                 &
     &        (plan1(ibid),n1,isign1,FFTW_ESTIMATE+FFTW_IN_PLACE)
!
         call FFTW_F77 (plan1(ibid),n2*mplane,f,1,nx1,fout,1,nx1)
      else
! i-direction
         if (plan1(ibid).eq.0)                                          &
     &        call FFTW_F77_CREATE_PLAN                                 &
     &        (plan1(ibid),n1,isign1,FFTW_ESTIMATE+FFTW_IN_PLACE)
!
         call FFTW_F77 (plan1(ibid),n2*mplane,f,1,nx1,fout,1,nx1)
! j-direction
         if (plan2(ibid).eq.0)                                          &
     &        call FFTW_F77_CREATE_PLAN                                 &
     &        (plan2(ibid),n2,isign1,FFTW_ESTIMATE+FFTW_IN_PLACE)
!
         do i=1,n1
!
! do only ffts on columns (i,*,k) resulting in nonzero components
!
            if (planes(i).eq.1) then
               call FFTW_F77                                            &
     &              (plan2(ibid),mplane,f(i),nx1,nx1*nx2,fout,0,0)
            end if
         end do
!
         fac=1.0/float(n1*n2)
         call DSCAL(2*nx1*nx2*mplane, fac, f, 1)
!
      end if
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_b (f,n1,n2,n3,n1x,n2x,n3x,imin3,imax3,isign)
!     ===============
!     driver routine for 3d complex fft's on box grid - fftw
!     fft along xy is done only on planes that correspond to
!     dense grid planes on the current processor, i.e. planes
!     with imin3 .le. n3 .le. imax3
!----------------------------------------------------------------------
!
      use fftw_mod
      implicit none
      integer n1,n2,n3,n1x,n2x,n3x,imin3,imax3, isign
      complex(kind=8) f(n1x*n2x*n3x)
#ifdef __PARA
      complex(kind=8), allocatable:: faux(:)
      integer ibid, nplanes, nstart
!
! initialization variables
!
      integer planz(2), planxy(2)
      save planz, planxy
      data planz/0,0/, planxy/0,0/
!
!
      if (isign.eq.1) then
         ibid =1
      else if (isign.eq.-1) then
         call errore('cft_b','not implemented',isign)
      end if
!
      if (planz(ibid).eq.0)                                             &
     &     call FFTW_F77_CREATE_PLAN                                    &
     &        (planz(ibid),n3,isign,FFTW_ESTIMATE)
!
! the fortran fftw routine for 1d fft cannot do fft in-place
!
      if (planxy(ibid).eq.0)                                            &
     &     call FFTW2D_F77_CREATE_PLAN                                  &
     &        (planxy(ibid),n1,n2,isign,FFTW_ESTIMATE)
!
! for 2d fft it could be done, but it is not convenient
!
      allocate(faux(n1x*n2x*n3x))
!
! 1d fft along z, output in faux
!
      call FFTW_F77(planz(ibid),n1x*n2x,f,n1x*n2x,1,faux,n1x*n2x,1)
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes=imax3-imin3+1
      nstart =(imin3-1)*n1x*n2x+1
      call FFTWND_F77 (planxy(ibid),nplanes,faux(nstart),1,n1x*n2x,      &
     &                 f(nstart),1,n1x*n2x)
      deallocate(faux)
#endif
      return
      end
#else
      subroutine bidon_fftw
         stop 'fftw'
      end subroutine bidon_fftw
#endif
