!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CRAYY) || defined(__SX4)
!
!-------------------------------------------------------------------------
      subroutine cfft3(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,isign)
!-----------------------------------------------------------------------
!     driver routine for 3d fft using cray's fft routine  (by kl)
!
      implicit none
      integer nr1, nr2, nr3, nr1x, nr2x, nr3x, isign
      real(8)    f(2,nr1x*nr2x*nr3x)
!
! initialization variables
!
      logical first
      data first/.true./
      integer ifax1(19), ifax2(19), ifax3(19)
      real(8), allocatable :: trig1(:), trig2(:), trig3(:)
!
! work variables
!
      real(8), allocatable:: work(:)
      real(8) fac
      integer inc, jump, lot, i, istart
      integer, save :: nr1_old = -1
      integer, save :: nr2_old = -1
      integer, save :: nr3_old = -1
!
      save    first, ifax1, ifax2, ifax3, trig1, trig2, trig3, work
!
      if( ( nr1_old /= nr1 ) .OR. ( nr2_old /= nr2 ) .OR. ( nr3_old /= nr3 ) ) then
        first = .true.
      end if
!
      if (first) then
         IF( ALLOCATED( work ) ) DEALLOCATE( work )
         IF( ALLOCATED( trig1 ) ) DEALLOCATE( trig1 )
         IF( ALLOCATED( trig2 ) ) DEALLOCATE( trig2 )
         IF( ALLOCATED( trig3 ) ) DEALLOCATE( trig3 )
         allocate(work(4*nr1x*nr2x*nr3x))
         allocate(trig1(2*nr1))
         allocate(trig2(2*nr2))
         allocate(trig3(2*nr3))
         call cftfax(nr1,ifax1,trig1)
         call cftfax(nr2,ifax2,trig2)
         call cftfax(nr3,ifax3,trig3)
         first = .false.
      end if
!                           x - direction
      inc=2
      jump=2*nr1x
      lot=nr3x*nr2x
      call cfftmlt(f(1,1),f(2,1),work,trig1,ifax1,inc,                  &
     &             jump,nr1,lot,isign)
!
!                           y - direction
!     inc=2*nr1x
!     jump=2
!     lot=nr1x
!     do i=1,nr3
!       istart=1+(i-1)*nr2x*nr1x
!       call cfftmlt(f(1,istart),f(2,istart),work,trig2,
!    c             ifax2,inc,jump,nr2,lot,isign)
!     end do
!
      inc=2*nr1x
      jump=2*nr1x*nr2x
      lot=nr3x
      do i=1,nr1
         istart=i
         call cfftmlt(f(1,istart),f(2,istart),work,trig2,               &
     &             ifax2,inc,jump,nr2,lot,isign)
      end do
!
!                           z - direction
      inc=2*nr1x*nr2x
      jump=2
      lot=nr1x*nr2x
      call cfftmlt(f(1,1),f(2,1),work,trig3,ifax3,inc,                  &
     &             jump,nr3,lot,isign)
!
      if (isign.eq.-1) then
         fac=1.d0/DBLE(nr1*nr2*nr3)
         call SSCAL(2*nr1x*nr2x*nr3x,fac,f,1)
      end if
!
      return
      end subroutine cfft3
!-------------------------------------------------------------------------
      subroutine cfft3b(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isign)
!-----------------------------------------------------------------------
!     driver routine for 3d fft using cray's fft routine  (by kl)
!
      implicit none
      integer nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isign
      real(8) f(2,nr1bx*nr2bx*nr3bx)
!
! initialization variables
!
      logical first
      data first/.true./
      integer ifax1(19), ifax2(19), ifax3(19)
      real(8), allocatable :: trig1(:), trig2(:), trig3(:)
!
! work variables
!
      real(8), allocatable:: work(:)
      real(8)  fac
      integer inc, jump, lot, i, istart
      integer, save :: nr1b_old = -1
      integer, save :: nr2b_old = -1
      integer, save :: nr3b_old = -1
      save    first, ifax1, ifax2,ifax3, trig1, trig2, trig3, work
!
      if( ( nr1b_old /= nr1b ) .OR. ( nr2b_old /= nr2b ) .OR. ( nr3b_old /= nr3b ) ) then
        first = .true.
      end if
!
      if (first) then
         IF( ALLOCATED( work ) ) DEALLOCATE( work )
         IF( ALLOCATED( trig1 ) ) DEALLOCATE( trig1 )
         IF( ALLOCATED( trig2 ) ) DEALLOCATE( trig2 )
         IF( ALLOCATED( trig3 ) ) DEALLOCATE( trig3 )
         allocate(work(4*nr1bx*nr2bx*nr3bx))
         allocate(trig1(2*nr1b))
         allocate(trig2(2*nr2b))
         allocate(trig3(2*nr3b))
         call cftfax(nr1b,ifax1,trig1)
         call cftfax(nr2b,ifax2,trig2)
         call cftfax(nr3b,ifax3,trig3)
         first=.false.
      end if
!                           x - direction
      inc=2
      jump=2*nr1bx
      lot=nr3bx*nr2bx
      call cfftmlt(f(1,1),f(2,1),work,trig1,ifax1,inc,                  &
     &             jump,nr1b,lot,isign)
!
!                           y - direction
      inc=2*nr1bx
      jump=2
      lot=nr1bx
      do i=1,nr3b
         istart=1+(i-1)*nr2bx*nr1bx
         call cfftmlt(f(1,istart),f(2,istart),work,trig2,               &
     &                ifax2,inc,jump,nr2b,lot,isign)
      end do
!
!                           z - direction
      inc=2*nr1bx*nr2bx
      jump=2
      lot=nr1bx*nr2bx
      call cfftmlt(f(1,1),f(2,1),work,trig3,ifax3,inc,                  &
     &             jump,nr3b,lot,isign)
!
      if (isign.eq.-1) then
         fac=1.d0/DBLE(nr1b*nr2b*nr3b)
         call SSCAL(2*nr1bx*nr2bx*nr3bx,fac,f,1)
      end if
!
      return
      end subroutine cfft3b
!-------------------------------------------------------------------------
      subroutine cfft3s(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,isign)
!-----------------------------------------------------------------------
!     driver routine for 3d fft using cray's fft routine  (by kl)
!
      implicit none
      integer nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, isign
      real(8) f(2,nr1sx*nr2sx*nr3sx)
!
! initialization variables
!
      logical first
      data first/.true./
      integer ifax1(19), ifax2(19), ifax3(19)
      real(8), allocatable :: trig1(:), trig2(:), trig3(:)
!
! work variables
!
      real(8), allocatable:: work(:)
      real(8) fac
      integer inc, jump, lot, i, istart
      integer, save :: nr1s_old = -1
      integer, save :: nr2s_old = -1
      integer, save :: nr3s_old = -1
!
      save    first, ifax1, ifax2, ifax3, trig1, trig2, trig3, work
!
      if( ( nr1s_old /= nr1s ) .OR. ( nr2s_old /= nr2s ) .OR. ( nr3s_old /= nr3s ) ) then
        first = .true.
      end if
!
      if (first) then
         IF( ALLOCATED( work ) ) DEALLOCATE( work )
         IF( ALLOCATED( trig1 ) ) DEALLOCATE( trig1 )
         IF( ALLOCATED( trig2 ) ) DEALLOCATE( trig2 )
         IF( ALLOCATED( trig3 ) ) DEALLOCATE( trig3 )
         allocate(work(4*nr1sx*nr2sx*nr3sx))
         allocate(trig1(2*nr1s))
         allocate(trig2(2*nr2s))
         allocate(trig3(2*nr3s))
         call cftfax(nr1s,ifax1,trig1)
         call cftfax(nr2s,ifax2,trig2)
         call cftfax(nr3s,ifax3,trig3)
         first=.false.
      end if
!
!                           x - direction
      inc=2
      jump=2*nr1sx
      lot=nr3sx*nr2sx
      call cfftmlt(f(1,1),f(2,1),work,trig1,ifax1,inc,                  &
     &             jump,nr1s,lot,isign)
!
!                           y - direction
!     inc=2*nr1sx
!     jump=2
!     lot=nr1sx
!     do i=1,nr3s
!        istart=1+(i-1)*nr2sx*nr1sx
!        call cfftmlt(f(1,istart),f(2,istart),work,trig2,
!    c             ifax2,inc,jump,nr2s,lot,isign)
!     end do
!
      inc=2*nr1sx
      jump=2*nr1sx*nr2sx
      lot=nr3sx
      do i=1,nr1s
         istart=i
         call cfftmlt(f(1,istart),f(2,istart),work,trig2,               &
     &             ifax2,inc,jump,nr2s,lot,isign)
      end do
!
!                           z - direction
      inc=2*nr1sx*nr2sx
      jump=2
      lot=nr1sx*nr2sx
      call cfftmlt(f(1,1),f(2,1),work,trig3,ifax3,inc,                  &
     &             jump,nr3s,lot,isign)
!
      if (isign.eq.-1) then
         fac=1.d0/DBLE(nr1s*nr2s*nr3s)
         call SSCAL(2*nr1sx*nr2sx*nr3sx,fac,f,1)
      end if
!
      return
      end subroutine cfft3s
#else
      subroutine bidon_cray
         stop 'crayfft'
      end subroutine bidon_cray
#endif
