!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!  cfft3, cfft3s, cfft3b are the 3d fft for ibm rs6000 essl library
!  Note that the three routines are equal, and that no initialization
!  is needed (the initialization arrays are recalculated every time).
!
#if defined (__ESSL)
!
!----------------------------------------------------------------------
      subroutine cfft3 (f,nr1,nr2,nr3,nr1x,nr2x,nr3x,sign)
!----------------------------------------------------------------------
!
      implicit none
      integer nr1, nr2, nr3, nr1x, nr2x, nr3x, sign
      complex(kind=8) f(nr1x*nr2x*nr3x)
! work variables
      integer naux, isign
      parameter (naux=60000)
      real(kind=8) aux(naux), scale
      external dscal, dcft
!
!
      if (sign.ne.-1 .and. sign.ne.1)                                   &
     &  call error('cfft3','which fft ?', 1)
!
! IBM uses a different sign convention from the "traditional" one
!
      isign = -sign
!
      if (isign.gt.0) then
         scale=1.0d0/(nr1*nr2*nr3)
      else
         scale=1.0d0
      end if
!
      call dcft3(f,nr1x,nr1x*nr2x,f,nr1x,nr1x*nr2x,nr1,nr2,nr3,         &
     &           isign, scale,aux,naux)
!
      return
      end
!----------------------------------------------------------------------
      subroutine cfft3s(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,sign)
!----------------------------------------------------------------------
!
      implicit none
      integer nr1, nr2, nr3, nr1x, nr2x, nr3x, sign
      complex(kind=8) f(nr1x*nr2x*nr3x)
! work variables
      integer naux, isign
      parameter (naux=60000)
      real(kind=8) aux(naux), scale
      external dscal, dcft
!
!
      if (sign.ne.-1 .and. sign.ne.1)                                   &
     &  call error('cfft3s','which fft ?', 1)
!
! IBM uses a different sign convention from the "traditional" one
!
      isign = -sign
!
      if (isign.gt.0) then
         scale=1.0d0/(nr1*nr2*nr3)
      else
         scale=1.0d0
      end if
!
      call dcft3(f,nr1x,nr1x*nr2x,f,nr1x,nr1x*nr2x,nr1,nr2,nr3,         &
     &           isign, scale,aux,naux)
!
      return
      end
!----------------------------------------------------------------------
      subroutine cfft3b(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,sign)
!----------------------------------------------------------------------
!
      implicit none
      integer nr1, nr2, nr3, nr1x, nr2x, nr3x, sign
      complex(kind=8) f(nr1x*nr2x*nr3x)
! work variables
      integer naux, isign
      parameter (naux=60000)
      real(kind=8) aux(naux), scale
      external dscal, dcft
!
!
      if (sign.ne.-1 .and. sign.ne.1)                                   &
     &  call error('cfft3b','which fft ?', 1)
!
! IBM uses a different sign convention from the "traditional" one
!
      isign = -sign
!
      if (isign.gt.0) then
         scale=1.0d0/(nr1*nr2*nr3)
      else
         scale=1.0d0
      end if
!
      call dcft3(f,nr1x,nr1x*nr2x,f,nr1x,nr1x*nr2x,nr1,nr2,nr3,         &
     &           isign, scale,aux,naux)
!
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_1     (f,m,n,nx,sgn,fout)
!     ===============
!     driver routine for m 1d complex fft's (dense grid)
!     nx=n+1 is allowed (in order to avoid memory conflicts)
!     A separate initialization is stored for sgn=-1 and sgn=+1
!     NOTA BENE: the output in fout !
!----------------------------------------------------------------------
!
      implicit none
      integer m, n, nx, sgn
      complex(kind=8)  f(nx*m), fout(nx*m)
#ifdef __PARA
!
! initialization variables
!
      logical first(2)
      data first /.true., .true./
      integer naux1
      parameter (naux1=20000)
      real(kind=8) aux3(naux1,2)
      save first, aux3
!
! work variables
!
      integer isign, naux, ibid
      parameter (naux=15000)
      real(kind=8) aux(naux), scale
      external DSCAL
!
!
      isign = -sgn
!
      if (isign.eq.1) then
         ibid =1
         scale=1.d0/n
      else if (isign.eq.-1) then
         ibid =2
         scale=1.d0
      else
         call error('cft_1','wrong call',isign)
      end if
!
      if (first(ibid)) then
         call dcft(1,f,1,nx,fout,1,nx,n,m,isign,                        &
     &        scale,aux3(1,ibid),naux1,aux,naux)
         first(ibid)=.false.
      end if
!
      call dcft(0,f,1,nx,fout,1,nx,n,m,isign,                           &
     &     scale,aux3(1,ibid),naux1,aux,naux)
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_2      (f,mplane,n1,n2,nx1,nx2,sgn)
!     ===============
!     driver routine for mplane 2d complex fft's of lenghts n1 and n2
!     (dense grid) nx1 is the actual dimension of f (may differ from n)
!     for compatibility: nx2=n2, nx2 is not used - uses ESSL
!     A separate initialization is stored for sgn=-1 and sgn=+1
!     
!----------------------------------------------------------------------
      implicit none
      integer n1, n2, mplane, nx1, nx2, sgn
      complex(kind=8)  f(nx1*nx2*mplane)
#ifdef __PARA
!
! initialization variables
!
      logical first(2)
      data first /.true., .true./
      integer naux1
      parameter (naux1=20000)
      real(kind=8) aux1(naux1,2), aux2(naux1,2)
      save first, aux1, aux2
!
      integer isign, ibid, m, incx1, incx2, k, istrt, naux
      parameter (naux=15000)
      real(kind=8) aux(naux), scale
      external DSCAL
!
!
      isign = -sign(1,sgn)
!
      if(n2.ne.nx2)                                                     &
     &  call error('cft_2','no longer implemented',1)
!
      if (isign.eq.1) then
         ibid=1
      else if (isign.eq.-1) then
         ibid=2
      else
         call error('cft_2','wrong call',ibid)
      end if
!
      scale=1.d0
!
!  i - direction ...
!
      incx1=1
      incx2=nx1
      m    =n2*mplane
      if (first(ibid)) then
         call dcft(1,f,incx1,incx2,f,incx1,incx2,n1,m,isign,            &
     &        scale,aux1(1,ibid),naux1,aux,naux)
      end if
      call dcft(0,f,incx1,incx2,f,incx1,incx2,n1,m,isign,               &
     &     scale,aux1(1,ibid),naux1,aux,naux)
!
! ... j-direction ...
!
      incx1 = nx1
      incx2 = 1
      m     = n1
      if (first(ibid)) then
         call dcft(1,f,incx1,incx2,f,incx1,incx2,n2,                    &
     &        m,isign,scale,aux2(1,ibid),naux1,aux,naux)
         first(ibid)=.false.
      end if
!
      do k= 1,mplane
         istrt = 1 + (k-1)*nx1*n2
         call dcft(0,f(istrt),incx1,incx2,f(istrt),incx1,incx2,n2,      &
     &        m,isign,scale,aux2(1,ibid),naux1,aux,naux)
      end do
!
      if (isign.eq.1) call DSCAL(2*nx1*n2*mplane,1d0/(n1*n2),f,1)
!
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_1s    (f,m,n,nx,sgn,fout)
!     ===============
!     driver routine for m 1d complex fft's (sparse grid)
!     nx=n+1 is allowed (in order to avoid memory conflicts)
!     On input, sgn=+/-1 for charge density, sgn=+/-2 for wavefunctions
!     A separate initialization is stored for each of the cases:
!     sgn=-1, sgn=-2, sgn=+1, and sgn=+2 !!!
!     NOTA BENE: the output is in fout
!----------------------------------------------------------------------
!
      implicit none
      integer m, n, nx, sgn
      complex(kind=8)  f(nx*m), fout(nx*m)
#ifdef __PARA
!
! initialization variables
!
      logical first(2,2)
      data first /.true., .true., .true., .true./
      integer naux1
      parameter (naux1=20000)
      real(kind=8) aux3(naux1,2,2)
      save   first, aux3
!
! work variables
!
      integer isign, sign, naux, itype, ibid
      parameter (naux=15000)
      real(kind=8) aux(naux), scale
      external DSCAL
!
!
      isign = -sign(1,sgn)
      itype = abs(sgn) 
      if (itype.le.0.or.itype.gt.2)                                     &
     &     call error('cft_1s','wrong call',sgn)
!
      if (isign.eq.1) then
         ibid =1
         scale=1.d0/n
      else if (isign.eq.-1) then
         ibid =2
         scale=1.d0
      end if
!
      if (first(ibid,itype)) then
         call dcft(1,f,1,nx,fout,1,nx,n,m,isign,                        &
     &        scale,aux3(1,ibid,itype),naux1,aux,naux)
         first(ibid,itype)=.false.
      end if
!
      call dcft(0,f,1,nx,fout,1,nx,n,m,isign,                           &
     &     scale,aux3(1,ibid,itype),naux1,aux,naux)
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_2s (f,mplane,n1,n2,nx1,nx2,sgn,planes)
!     ===============
!     driver routine for mplane 2d complex fft's of lengths n1 and n2
!     (sparse grid, both charge and wavefunctions) - uses ESSL
!     on input, sgn=+/-1 for charge density, sgn=+/-2 for wavefunctions
!     nx1 is the actual dimension of f (may differ from n)
!     for compatibility: nx2=n2, nx2 is not used
!     A separate initialization is stored for sgn=-1 and sgn=+1
!
!----------------------------------------------------------------------
!
      implicit none
      integer n1, n2, mplane, nx1, nx2, sgn, planes(nx1)
      complex(kind=8)  f(nx1*nx2*mplane)
#ifdef __PARA
!
! initialization variables
!
      logical first(2)
      data first /.true., .true./
      integer naux1
      parameter (naux1=20000)
      real(kind=8) aux1(naux1,2), aux2(naux1,2)
      save   first, aux1, aux2
!
! work variables
!
      integer isign, itype, ibid, m, incx1, incx2, k, i, istrt, naux
      parameter (naux=15000)
      real(kind=8)   aux(naux), scale
      external DSCAL
!
!
      isign = -sign(1,sgn)
      itype = abs(sgn) 
      if (itype.le.0.or.itype.gt.2) call error('cft_2','wrong call',1)
!
      if(n2.ne.nx2)                                                     &
     &  call error('cft_2','no longer implemented',1)
!
      scale=1.d0
!
! check how many columns along x are nonzero
!
      m=0
      do i=1,n1
         m=m+planes(i)
      end do
      if (m.gt.n1.or.m.le.0)                                            &
     &     call error('cft_2s','something wrong with planes',1)
!
      if (isign.eq.1) then
!
! ... i - direction
!
         incx1=1
         incx2=nx1
         m=n2*mplane
         if (first(1)) then
            call dcft(1,f,incx1,incx2,f,incx1,incx2,n1,m,isign,         &
     &           scale,aux1(1,1),naux1,aux,naux)
         end if
         call dcft(0,f,incx1,incx2,f,incx1,incx2,n1,m,isign,            &
     &             scale,aux1(1,1),naux1,aux,naux)
!
! j-direction ...
!
         incx1 = nx1
         incx2 = nx1*nx2
         m=mplane
         if (first(1)) then
            call dcft                                                   &
     &           (1,f,incx1,incx2,f,incx1,incx2,                        &
     &           n2,m,isign,scale,aux2(1,1),naux1,aux,naux)
            first(1)=.false.
         end if
         do i= 1,n1
!
! do only ffts on columns (i,*,k) resulting in nonzero components
!
            if (planes(i).eq.1) then
               call dcft                                                &
     &              (0,f(i),incx1,incx2,f(i),incx1,incx2,               &
     &              n2,m,isign,scale,aux2(1,1),naux1,aux,naux)
            end if
         end do
!
         call DSCAL(2*nx1*n2*mplane,1d0/(n1*n2),f,1)
!
      else if (isign.eq.-1) then
!
!     ... j-direction
!
         incx1 = nx1
         incx2 = nx1*nx2
         m=mplane
         if (first(2)) then
            call dcft                                                   &
     &           (1,f,incx1,incx2,f,incx1,incx2,                        &
     &           n2,m,isign,scale,aux2(1,2),naux1,aux,naux)
         end if
         do i = 1,n1
!
! do only ffts for columns (i,*,k) having nonzero components
!
            if (planes(i).eq.1.or.itype.eq.1) then
               call dcft                                                &
     &              (0,f(i),incx1,incx2,f(i),incx1,incx2,               &
     &              n2,m,isign,scale,aux2(1,2),naux1,aux,naux)
            end if
         end do
!
!     i - direction ...
!
         incx1= 1
         incx2=nx1
         m=n2*mplane
         if (first(2)) then
            call dcft(1,f,incx1,incx2,f,incx1,incx2,n1,m,isign,scale,   &
     &           aux1(1,2),naux1,aux,naux)
            first(2)=.false.
         end if
         call dcft(0,f,incx1,incx2,f,incx1,incx2,n1,m,isign,scale,      &
     &                  aux1(1,2),naux1,aux,naux)
!
      endif
#endif
      return
      end
!
!----------------------------------------------------------------------
      subroutine cft_b (f,n1,n2,n3,n1x,n2x,n3x,imin3,imax3,sgn)
!     ===============
!     driver routine for 3d complex fft's on box grid - ibm essl
!     fft along xy is done only on planes that correspond to
!     dense grid planes on the current processor, i.e. planes
!     with imin3 .le. n3 .le. imax3
!----------------------------------------------------------------------
!
      implicit none
      integer n1,n2,n3,n1x,n2x,n3x,imin3,imax3,sgn
      complex(kind=8) f(n1x*n2x*n3x)
#ifdef __PARA
!
! initialization variables
!
      logical first(2)
      data first /.true., .true./
      integer naux1
      parameter (naux1=20000)
      real(kind=8) aux3(naux1,2), aux2(naux1,2), aux1(naux1,2)
      save first, aux1, aux2, aux3
!
! work variables
!
      integer isign, naux, ibid, nplanes, nstart, k
      parameter (naux=15000)
      real(kind=8) aux(naux), scale
!
!
      isign = -sgn
      if (isign.eq.-1) then
         ibid =1
         scale=1.d0
      else if (isign.eq.1) then
         call error('cft_b','not implemented',isign)
      end if
!
      if (first(ibid)) then
!
! initialization for the z-direction...
!
         call dcft(1,f,n1x*n2x,1,f,n1x*n2x,1,n3,n1x*n2x,isign,          &
     &        scale,aux3(1,ibid),naux1,aux,naux)
         first(ibid)=.false.
      end if
!
! fft in the z-direction...
!
      call dcft(0,f,n1x*n2x,1,f,n1x*n2x,1,n3,n1x*n2x,isign,             &
     &        scale,aux3(1,ibid),naux1,aux,naux)
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes=imax3-imin3+1
      nstart =(imin3-1)*n1x*n2x+1
!
! x-direction  - Inizialization must be done every time because it depends
! on nplanes that may vary from call to call !!!! sigh
!
      call dcft(1,f,1,n1x,f,1,n1x,n1,n2x*nplanes,isign,                  &
     &        scale,aux1(1,ibid),naux1,aux,naux)
      call dcft(0,f(nstart),1,n1x,f(nstart),1,n1x,n1,n2x*nplanes,isign,  &
     &        scale,aux1(1,ibid),naux1,aux,naux)
!
! y-direction
!
      call dcft(1,f,n1x,1,f,n1x,1,n2,n1x,isign,                          &
     &        scale,aux2(1,ibid),naux1,aux,naux)
!
      do k= imin3,imax3
        nstart=(k-1)*n1x*n2x+1
        call dcft(0,f(nstart),n1x,1,f(nstart),n1x,1,n2,n1x,isign,        &
     &        scale,aux2(1,ibid),naux1,aux,naux)
     end do
#endif
      return
      end
#else
      subroutine bidon_ibm
         stop 'ibmfft'
      end subroutine bidon_ibm
#endif

