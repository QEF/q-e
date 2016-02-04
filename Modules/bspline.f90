!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! Modified by Davide Ceresoli:
! - use dp from module kinds
! - write to string instead of stdout
! - return error codes

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   VERSION 2.2
!
!   f90 VERSION
!
!   This library contains routines for B-spline interpolation in
!   one, two, and three dimensions. Part of the routines are based
!   on the book by Carl de Boor: A practical guide to Splines (Springer,
!   New-York 1978) and have the same calling sequence and names as
!   the corresponding routines from the IMSL library. For documen-
!   tation see the additional files. NOTE: The results in the demo
!   routines may vary slightly on different architectures.
!
!   by W. Schadow 12/04/99
!   last changed by W. Schadow 07/28/2000
!
!
!   Wolfgang Schadow
!   TRIUMF
!   4004 Wesbrook Mall
!   Vancouver, B.C. V6T 2A3
!   Canada
!
!   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
!
!   www  : http://www.triumf.ca/people/schadow
!
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   Copyright (C) 2000 Wolfgang Schadow
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Library General Public
!   License as published by the Free Software Foundation; either
!   version 2 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public
!   License along with this library; if not, write to the
!   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!   Boston, MA  02111-1307, USA.
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE bspline
  USE kinds, only : dp

  IMPLICIT NONE


  PRIVATE

  character(80)  :: routine
  character(256) :: errmsg

  public dbsnak
  public dbsint, dbsval, dbsder, dbs1gd
  public dbs2in, dbs2dr, dbs2vl, dbs2gd
  public dbs3in, dbs3vl, dbs3dr, dbs3gd
  public get_error_routine
  public get_error_message

CONTAINS

  !==================================================================
  subroutine dbsnak(nx,xvec,kxord,xknot, ierr)
  !==================================================================
  !
  !  Compute the `not-a-knot' spline knot sequence.
  !  (see de Boor p. 167)
  !
  !   nx     - number of data points.  (input)
  !   xvec   - array of length ndata containing the location of the
  !            data points.  (input)
  !   kxord  - order of the spline.  (input)
  !   xknot  - array of length ndata+korder containing the knot
  !            sequence.  (output)
  !
  implicit none
  integer, intent(in) :: nx, kxord
  real(dp), dimension(nx), intent(in)        :: xvec
  real(dp), dimension(nx+kxord), intent(out) :: xknot
  integer, intent(out) :: ierr
  real(dp) :: eps
  integer        :: ix
  logical        :: first = .true.

  save first,eps

  routine = 'dbsnak'
  ierr = 0

  if (first) then
     first=.false.
     eps = epsilon(1.0_dp)
     !!write(6,*) "subroutine dbsnak: "
     !!write(6,*) "eps = ",eps
  endif

  if((kxord .lt. 0) .or. (kxord .gt. nx)) then
     write(errmsg,*) "0 <= kxord <= nx is required: kxord,nx=", kxord, nx
     ierr = 1
     return
  endif

  do ix = 1, kxord
     xknot(ix) = xvec(1)
  end do

  if(mod(kxord,2) .eq. 0) then
     do ix = kxord+1, nx
        xknot(ix) = xvec(ix-kxord/2)
     end do
  else
     do ix = kxord+1, nx
        xknot(ix) = 0.5_dp * (xvec(ix-kxord/2) + xvec(ix-kxord/2-1))
     end do
  endif

  do ix = nx+1, nx+kxord
     xknot(ix) = xvec(nx) * (1.0_dp + eps)
  end do

  end subroutine dbsnak


  !==================================================================
  subroutine dbsint(nx,xvec,xdata,kx,xknot,bcoef, ierr)
  !==================================================================
  !
  !  Computes the spline interpolant, returning the B-spline coefficients.
  !  (see de Boor p. 204)
  !
  !   nx     - number of data points.  (input)
  !   xvec   - array of length nx containing the data point
  !            abscissas.  (input)
  !   xdata  - array of length ndata containing the data point
  !            ordinates.  (input)
  !   kx     - order of the spline.  (input)
  !            korder must be less than or equal to ndata.
  !   xknot  - array of length nx+kx containing the knot
  !            sequence.  (input)
  !            xknot must be nondecreasing.
  !   bscoef - array of length ndata containing the B-spline
  !            coefficients.  (output)
  !
  implicit none

  integer, intent(in)                    :: nx, kx
  real(dp), dimension(nx), intent(in)    :: xdata, xvec
  real(dp), dimension(nx+kx), intent(in) :: xknot
  real(dp), dimension(nx), intent(out)   :: bcoef
  integer, intent(out)                   :: ierr

  integer                                :: nxp1, kxm1, kpkm2, leftx, lenq
  integer                                :: ix, ik,ilp1mx, jj, iflag
  real(dp)                               :: xveci
  real(dp), dimension((2*kx-1)*nx)       :: work


  routine = 'dbsint'
  ierr = 0
  nxp1  = nx + 1
  kxm1  = kx - 1
  kpkm2 = 2 * kxm1
  leftx = kx
  lenq  = nx * (kx + kxm1)

  do ix = 1, lenq
     work(ix) = 0.0_dp
  end do

  do  ix = 1, nx
     xveci  = xvec(ix)
     ilp1mx = min0(ix+kx,nxp1)
     leftx   = max0(leftx,ix)
     if (xveci .lt. xknot(leftx)) goto 998
30   if (xveci .lt. xknot(leftx+1)) go to 40
     leftx = leftx + 1
     if (leftx .lt. ilp1mx) go to 30
     leftx = leftx - 1
     if (xveci .gt. xknot(leftx+1)) goto 998
40   call bsplvb (xknot,nx+kx,kx,1,xveci,leftx,bcoef)
     jj = ix - leftx + 1 + (leftx - kx) * (kx + kxm1)
     do ik = 1, kx
        jj       = jj + kpkm2
        work(jj) = bcoef(ik)
     end do
  end do

  call banfac(work,kx+kxm1,nx,kxm1,kxm1,iflag)

  if (iflag .ne. 1) then
     write(errmsg,*) 'no solution of linear equation system'
     ierr = 1
     return
  end if

  do ix = 1, nx
     bcoef(ix) = xdata(ix)
  end do

  call banslv(work,kx+kxm1,nx,kxm1,kxm1,bcoef)
  return

998 write(errmsg,*) "xknot(ix) <= xknot(ix+1) required: ix,xknot(ix),xknot(ix+1)=", &
                    ix,xknot(ix),xknot(ix+1)
  ierr = 2
  return

  end subroutine dbsint


  !==================================================================
  function dbsval(x,kx,xknot,nx,bcoef, ierr)
  !==================================================================
  !
  !  Evaluates a spline, given its B-spline representation.
  !
  !   x      - point at which the spline is to be evaluated.  (input)
  !   kx     - order of the spline.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence.  (input)
  !            xknot must be nondecreasing.
  !   nx     - number of B-spline coefficients.  (input)
  !   bcoef  - array of length nx containing the B-spline
  !            coefficients.  (input)
  !   dbsval - value of the spline at x.  (output)
  !
  implicit none

  integer, intent(in)                    :: nx, kx
  real(dp)                               :: dbsval
  real(dp)                               :: x
  real(dp), dimension(nx+kx), intent(in) :: xknot
  real(dp), dimension(nx), intent(in)    :: bcoef
  integer, intent(out)                   :: ierr
  integer                 :: il, ik, ix, leftx
  real(dp)                :: save1, save2
  real(dp), dimension(kx) :: work, dl, dr

  ierr = 0
  routine = 'dbsval'
  dbsval = 0.0_dp
  !
  !     check if xknot(i) <= xknot(i+1) and calculation of i so that
  !     xknot(i) <= x < xknot(i+1)
  !
  leftx = 0

  do ix = 1,nx+kx-1
     if (xknot(ix) .gt. xknot(ix+1)) then
        write(errmsg,*) "xknot(ix) <= xknot(ix+1) required: ix,xknot(ix),xknot(ix+1)=", &
                        ix,xknot(ix),xknot(ix+1)
        ierr = 1
        return
     endif
     if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
  end do

  if(leftx .eq. 0) then
     write(errmsg,*) "ix with xknot(ix) <= x < xknot(ix+1) required: x=", x
     ierr = 2
     return
  endif

  do ik = 1, kx-1
     work(ik) = bcoef(leftx+ik-kx)
     dl(ik)   = x - xknot(leftx+ik-kx)
     dr(ik)   = xknot(leftx+ik) - x
  end do

  work(kx)  = bcoef(leftx)
  dl(kx)    = x - xknot(leftx)

  do ik = 1, kx-1
     save2 = work(ik)
     do il = ik+1, kx
        save1 = work(il)
        work(il) = (dl(il) * work(il) + dr(il-ik) * save2)                  &
             &           / (dl(il) + dr(il - ik))
        save2 = save1
     end do
  end do

  dbsval = work(kx)

  end function dbsval


  !==================================================================
  function dbsder(iderx,x,kx,xknot,nx,bcoef, ierr)
  !==================================================================
  !
  !  Evaluates the derivative of a spline, given its B-spline representation.
  !
  !
  !   iderx  - order of the derivative to be evaluated.  (input)
  !            in particular, iderx = 0 returns the value of the
  !            spline.
  !   x      - point at which the spline is to be evaluated.  (input)
  !   kx     - order of the spline.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence.  (input)
  !            xknot must be nondecreasing.
  !   nx     - number of B-spline coefficients.  (input)
  !   bcoef  - array of length nx containing the B-spline
  !            coefficients.  (input)
  !   dbsder - value of the iderx-th derivative of the spline at x.
  !            (output)
  !
  implicit none

  integer, intent(in)                    :: iderx, kx, nx
  real(dp)                               :: dbsder
  real(dp), intent(in)                   :: x
  real(dp), dimension(nx+kx), intent(in) :: xknot
  real(dp), dimension(nx), intent(in)    :: bcoef
  integer, intent(out)                   :: ierr

  integer                       :: ix, ik, il, leftx
  real(dp)                :: save, save1, save2, y, sum, dik
  real(dp), dimension(kx) :: work, dl, dr,bsp

  ierr = 0
  routine = 'dbsder'
  dbsder = 0.0_dp
  !
  !     check if xknot(i) <= xknot(i+1) and calculation of i so that
  !     xknot(i) <= x < xknot(i+1)
  !
  leftx = 0
  do ix = 1,nx+kx-1
     if (xknot(ix) .gt. xknot(ix+1)) then
        write(errmsg,*) "xknot(ix) <= xknot(ix+1) required: ix,xknot(ix),xknot(ix+1)=", &
                        ix,xknot(ix),xknot(ix+1)
        ierr = 1
        return
     endif
     if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
  end do

  if (leftx .eq. 0) then
     write(errmsg,*) "ix with xknot(ix) <= x < xknot(ix+1) required: x=", x
     ierr = 2
     return
  endif

  if (iderx .eq. 0) then

     do ik = 1,kx-1
        work(ik) = bcoef(leftx+ik-kx)
        dl(ik)   = x - xknot(leftx+ik-kx)
        dr(ik)   = xknot(leftx+ik) - x
     end do

     work(kx)  = bcoef(leftx)
     dl(kx)    = x - xknot(leftx)

     do ik = 1,kx-1
        save2 = work(ik)
        do il = ik+1,kx
           save1 = work(il)
           work(il) = (dl(il) * work(il) + dr(il-ik) * save2)               &
                &              / (dl(il) + dr(il - ik))
           save2 = save1
        end do
     end do

     dbsder = work(kx)

  elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

     bsp(1) = 1.0_dp
     do ik = 1,kx-iderx-1
        dr(ik) = xknot(leftx+ik) - x
        dl(ik) = x - xknot(leftx+1-ik)
        save   = bsp(1)
        bsp(1) = 0.0_dp
        do il = 1, ik
           y         = save / (dr(il) + dl(ik+1-il))
           bsp(il)   = bsp(il) + dr(il) * y
           save      = bsp(il+1)
           bsp(il+1) = dl(ik+1-il) * y
        end do
     end do

     do ik = 1, kx
        work(ik) = bcoef(leftx+ik-kx)
        dr(ik)   = xknot(leftx+ik) - x
        dl(ik)   = x - xknot(leftx+ik-kx)
     end do

     do ik = 1, iderx
        dik   = dble(kx - ik)
        save2 = work(ik)
        do il = ik+1, kx
           save1    = work(il)
           work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
           save2    = save1
        end do
     end do

     sum = 0.0_dp

     do ix = 1, kx-iderx
        sum = sum + bsp(ix) * work(iderx+ix)
     end do

     dbsder = sum

  else
     dbsder = 0.0_dp
  endif

  end function dbsder


  !==================================================================
  subroutine dbs1gd(iderx,nxvec,xvec,kx,xknot,nx,bcoef,val, ierr)
  !==================================================================
  !
  !  Evaluates the derivative of a spline on a grid, given its B-spline
  !  representation.
  !
  !   iderx  - order of the derivative to be evaluated.  (input)
  !            in particular, iderx = 0 returns the value of the
  !            spline.
  !   nxvec  - length of vector xvec.  (input)
  !   xvec   - array of length nxvec containing the points at which the
  !            spline is to be evaluated.  (input)
  !            xvec should be strictly increasing.
  !   kx     - order of the spline.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence.  (input)
  !            xknot must be nondecreasing.
  !   nx     - number of B-spline coefficients.  (input)
  !   bcoef  - array of length nx containing the B-spline
  !            coefficients.  (input)
  !   val    - array of length nxvec containing the values of the
  !            iderx-th derivative of the spline at the points in
  !            xvec.  (output)
  !
  implicit none

  integer, intent(in)                     :: iderx, nxvec, kx, nx
  real(dp), dimension(nxvec), intent(in)  :: xvec
  real(dp), dimension(nx), intent(in)     :: bcoef
  real(dp), dimension(nx+kx), intent(in)  :: xknot
  real(dp), dimension(nxvec), intent(out) :: val
  integer, intent(out)                    :: ierr

  integer                             :: i, il, ik, ix
  integer, dimension(nxvec)           :: leftx
  real(dp)                      :: dik
  real(dp), dimension(nxvec,kx) :: dl, dr, biatx, work
  real(dp), dimension(nxvec)    :: save1, save2, term

  logical :: same, next

  routine = 'dbs1gd'
  ierr = 0
  leftx(1) = 0

  call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

  do ix = 2, nxvec
     leftx(ix) = leftx(ix-1)
     same = (xknot(leftx(ix)) .le. xvec(ix))                                &
          &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
     if(.not. same ) then
        leftx(ix) = leftx(ix) + 1
        next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
             &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
        if (.not. next)                                                     &
             &           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
     endif
  end do

  do ix = 1, nx+kx-1
     if (xknot(ix) .gt. xknot(ix+1)) then
        write(errmsg,*) "xknot(ix) <= xknot(ix+1) required: ix,xknot(ix),xknot(ix+1)=", &
                        ix,xknot(ix),xknot(ix+1)
        ierr = 1
        return
     endif
  end do

  do ix = 1, nxvec
     if ((xvec(ix).lt.xknot(1)).or.(xvec(ix).gt.xknot(nx+kx))) then
       write(errmsg,*) "ix with xknot(ix) <= x < xknot(ix+1) required: x=", xvec(ix)
       ierr = 2
       return
     endif
  end do

  if (iderx .eq. 0) then

     do ix = 1,nxvec
        biatx(ix,1) = 1._dp
        val(ix)     = 0._dp
     end do

     do ik = 1, kx-1
        do ix = 1, nxvec
           dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
           dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
           save1(ix) = 0._dp
        end do

        do il = 1, ik
           do ix = 1,nxvec
              term(ix)     = biatx(ix,il)                                   &
                   &                 / (dr(ix,il) + dl(ix,ik+1-il))
              biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
              save1(ix)    = dl(ix,ik+1-il) * term(ix)
           end do
        end do

        do ix = 1, nxvec
           biatx(ix,ik+1) = save1(ix)
        end do
     end do

     do ik = 1, kx
        do ix = 1, nxvec
           val(ix) = val(ix) + biatx(ix,ik) * bcoef(leftx(ix)-kx+ik)
        end do
     end do

  elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

     do ix = 1, nxvec
        biatx(ix,1) = 1._dp
        val(ix)     = 0._dp
     end do

     do ik = 1, kx-iderx-1
        do ix = 1, nxvec
           dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
           dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+1-ik)
           save1(ix)    = biatx(ix,1)
           biatx(ix,1) = 0.0_dp
           do il = 1, ik
              term(ix)       = save1(ix)                                    &
                   &                 / (dr(ix,il) + dl(ix,ik+1-il))
              biatx(ix,il)   = biatx(ix,il) + dr(ix,il) * term(ix)
              save1(ix)      = biatx(ix,il+1)
              biatx(ix,il+1) = dl(ix,ik+1-il) * term(ix)
           end do
        end do
     end do

     do ik = 1, kx
        do ix = 1, nxvec
           work(ix,ik) = bcoef(leftx(ix)+ik-kx)
           dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
           dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+ik-kx)
        end do
     end do

     do ik = 1, iderx
        dik   = dble(kx - ik)
        do ix = 1, nxvec
           save2(ix) = work(ix,ik)
           do il = ik+1, kx
              save1(ix)   = work(ix,il)
              work(ix,il) = dik * (work(ix,il) - save2(ix))                 &
                   &                 /(dl(ix,il) + dr(ix,il-ik))
              save2(ix)   = save1(ix)
           end do
        end do
     end do

     do i = 1, kx-iderx
        do ix = 1, nxvec
           val(ix) = val(ix) + biatx(ix,i) * work(ix,iderx+i)
        end do
     end do

  else

     do ix = 1, nxvec
        val(ix) = 0.0_dp
     end do

  endif

  end subroutine dbs1gd


  !==================================================================
  function dbsdca(iderx,x,kx,xknot,nx,bcoef,leftx)
  !==================================================================
  !
  ! This routine is equivalent to the routine dbsder, but it does not
  ! check the parameters!!!
  !
  ! Evaluates the derivative of a spline, given its B-spline representation.
  !
  !
  !   iderx  - order of the derivative to be evaluated.  (input)
  !            in particular, iderx = 0 returns the value of the
  !            spline.
  !   x      - point at which the spline is to be evaluated.  (input)
  !   kx     - order of the spline.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence.  (input)
  !            xknot must be nondecreasing.
  !   nx     - number of B-spline coefficients.  (input)
  !   bcoef  - array of length nx containing the B-spline
  !            coefficients.  (input)
  !   leftx  - number of the intervall of xknot that includes x
  !   dbsdca - value of the ideriv-th derivative of the spline at x.
  !            (output)
  !
  implicit none

  integer, intent(in)                    :: iderx, kx, nx
  real(dp)                               :: dbsdca
  real(dp), intent(in)                   :: x
  real(dp), dimension(nx+kx), intent(in) :: xknot
  real(dp), dimension(nx), intent(in)    :: bcoef

  integer                 :: i, ik, il, leftx
  real(dp)                :: save, save1, save2, y, sum, dik
  real(dp), dimension(kx) :: work, dl, dr,bsp


  if (iderx .eq. 0) then

     do ik = 1, kx-1
        work(ik) = bcoef(leftx+ik-kx)
        dl(ik)   = x - xknot(leftx+ik-kx)
        dr(ik)   = xknot(leftx+ik) - x
     end do

     work(kx)  = bcoef(leftx)
     dl(kx)    = x - xknot(leftx)

     do ik = 1, kx-1
        save2 = work(ik)
        do il = ik+1, kx
           save1 = work(il)
           work(il) = (dl(il) * work(il) + dr(il-ik) * save2)               &
                &              / (dl(il) + dr(il - ik))
           save2 = save1
        end do
     end do

     dbsdca = work(kx)

  elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
     bsp(1) = 1.0_dp
     do ik = 1,kx-iderx-1
        dr(ik) = xknot(leftx+ik) - x
        dl(ik) = x - xknot(leftx+1-ik)
        save   = bsp(1)
        bsp(1) = 0.0_dp
        do il = 1, ik
           y         = save / (dr(il) + dl(ik+1-il))
           bsp(il)   = bsp(il) + dr(il) * y
           save      = bsp(il+1)
           bsp(il+1) = dl(ik+1-il) * y
        end do
     end do

     do ik = 1, kx
        work(ik) = bcoef(leftx+ik-kx)
        dr(ik)   = xknot(leftx+ik) - x
        dl(ik)   = x - xknot(leftx+ik-kx)
     end do

     do ik = 1, iderx
        dik   = dble(kx - ik)
        save2 = work(ik)
        do il = ik+1, kx
           save1    = work(il)
           work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
           save2    = save1
        end do
     end do

     sum = 0.0_dp

     do i = 1, kx-iderx
        sum = sum + bsp(i) * work(iderx+i)
     end do

     dbsdca = sum

  else
     dbsdca = 0.0_dp
  endif

  end function dbsdca


  !==================================================================
  subroutine dbs2in(nx,xvec,ny,yvec,xydata,ldf,kx,ky,xknot,yknot,bcoef, ierr)
  !==================================================================
  !
  !  Computes a two-dimensional tensor-product spline interpolant,
  !  returning the tensor-product B-spline coefficients.
  !
  !    nx     - number of data points in the x-direction.  (input)
  !    xvec   - array of length nx containing the data points in
  !             the x-direction.  (input)
  !             xdata must be strictly increasing.
  !    ny     - number of data points in the y-direction.  (input)
  !    yvec   - array of length ny containing the data points in
  !             the y-direction.  (input)
  !             ydata must be strictly increasing.
  !    xydata - array of size nx by nydata containing the values to
  !             be interpolated.  (input)
  !             fdata(i,j) is the value at (xdata(i),ydata(j)).
  !    ldf    - the leading dimension of fdata exactly as specified in
  !             the dimension statement of the calling program.
  !             (input)
  !    kx     - order of the spline in the x-direction.  (input)
  !             kxord must be less than or equal to nxdata.
  !    ky     - order of the spline in the y-direction.  (input)
  !             kyord must be less than or equal to nydata.
  !    xknot  - array of length nx+kx containing the knot
  !             sequence in the x-direction.  (input)
  !             xknot must be nondecreasing.
  !    yknot  - array of length ny+ky containing the knot
  !             sequence in the y-direction.  (input)
  !             yknot must be nondecreasing.
  !    bcoef  - array of length nx*ny containing the
  !             tensor-product B-spline coefficients.  (output)
  !             bscoef is treated internally as a matrix of size nxdata
  !             by nydata.
  !
  implicit none

  integer, intent(in)                           :: nx, ny, kx, ky, ldf

  real(dp), dimension(nx), intent(in)     :: xvec
  real(dp), dimension(ny), intent(in)     :: yvec
  real(dp), dimension(nx+kx), intent(in)  :: xknot
  real(dp), dimension(ny+ky), intent(in)  :: yknot
  real(dp), dimension(ldf,*), intent(in)  :: xydata
  real(dp), dimension(nx,ny), intent(out) :: bcoef
  integer, intent(out)                    :: ierr
  real(dp), dimension(max(nx,ny),max(nx,ny))        :: work1
  real(dp), dimension(max(nx,ny))                   :: work2
  real(dp), dimension(max((2*kx-1)*nx,(2*ky-1)*ny)) :: work3


  call spli2d(xvec,ldf,xydata,xknot,nx,kx,ny,work2,work3,work1, ierr)
  if (ierr /= 0) return
  call spli2d(yvec,ny, work1, yknot,ny,ky,nx,work2,work3,bcoef, ierr)

  end subroutine dbs2in


  !==================================================================
  subroutine spli2d(xyvec,ld,xydata,xyknot,n,k,m,work2,work3,bcoef, ierr)
  !==================================================================
  implicit none

  integer, intent(in)                         :: ld, n, k, m
  real(dp), dimension(n), intent(in)    :: xyvec
  real(dp), dimension(n+k), intent(in)  :: xyknot
  real(dp), dimension(ld,m), intent(in) :: xydata
  real(dp), dimension(m,n), intent(out) :: bcoef

  real(dp), dimension(n), intent(out)         :: work2
  real(dp), dimension((2*k-1)*n), intent(out) :: work3
  integer, intent(out)                    :: ierr


  integer        :: np1, km1, kpkm2, left, lenq, i, iflag, ilp1mx, j, jj
  real(dp) :: xyveci

  routine = 'spli2d'
  ierr = 0
  np1   = n + 1
  km1   = k - 1
  kpkm2 = 2 * km1
  left  = k
  lenq  = n * (k + km1)

  do i = 1,lenq
     work3(i) = 0.0_dp
  end do

  do i = 1, n
     xyveci  = xyvec(i)
     ilp1mx = min0(i+k,np1)
     left   = max0(left,i)
     if (xyveci .lt. xyknot(left)) go to 998
30   if (xyveci .lt. xyknot(left+1)) go to 40
     left = left + 1
     if (left .lt. ilp1mx) go to 30
     left = left - 1
     if (xyveci .gt. xyknot(left+1)) go to 998
40   call bsplvb(xyknot,n+k,k,1,xyveci,left,work2)
     jj = i - left + 1 + (left - k) * (k + km1)
     do j = 1, k
        jj        = jj + kpkm2
        work3(jj) = work2(j)
     end do
  end do

  call banfac(work3,k+km1,n,km1,km1,iflag )

  if (iflag .ne. 1) then
     write(errmsg,*) "no solution of linear equation system"
     ierr = 1
     return
  end if

  do j = 1, m
     do i = 1, n
        work2(i) = xydata(i,j)
     end do

     call banslv(work3,k+km1,n,km1,km1,work2)

     do i = 1, n
        bcoef(j,i) = work2(i)
     end do
  end do

  return

998 write(errmsg,*) "i with knot(i) <= x/y < knot(i+1) required: knot(1),knot(n+k)=", &
                    xyknot(1), xyknot(n+k)
    ierr = 2
    return

  end subroutine spli2d


  !==================================================================
  function dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef, ierr)
  !==================================================================
  !
  !  evaluates a two-dimensional tensor-product spline, given its
  !  tensor-product B-spline representation.
  !
  !   x      - x-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   y      - y-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   kx     - order of the spline in the x-direction.  (input)
  !   ky     - order of the spline in the y-direction.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence in the x-direction.  (input)
  !            xknot must be nondecreasing.
  !   yknot  - array of length ny+ky containing the knot
  !            sequence in the y-direction.  (input)
  !            yknot must be nondecreasing.
  !   nx     - number of B-spline coefficients in the x-direction.
  !            (input)
  !   ny     - number of B-spline coefficients in the y-direction.
  !            (input)
  !   bcoef  - array of length nx*ny containing the
  !            tensor-product B-spline coefficients.  (input)
  !            bscoef is treated internally as a matrix of size nx
  !            by ny.
  !   dbs2vl - value of the spline at (x,y).  (output)
  !

    implicit none

  integer, intent(in)                    :: nx, ny, kx, ky
  real(dp), intent(in)                   :: x, y
  real(dp), dimension(nx+kx), intent(in) :: xknot
  real(dp), dimension(ny+ky), intent(in) :: yknot
  real(dp), dimension(nx,ny), intent(in) :: bcoef
  real(dp)                               :: dbs2vl
  integer, intent(out)                   :: ierr
  integer                       :: ix, iy, iky, leftx, lefty
  real(dp), dimension(ky) :: work

  routine = 'dbs2vl'
  ierr = 0
  dbs2vl = 0.0_dp
  !
  !     check if knot(i) <= knot(i+1) and calculation of i so that
  !     knot(i) <= x < knot(i+1)
  !
  leftx = 0

  do ix = 1, nx+kx-1
     if (xknot(ix) .gt. xknot(ix+1)) then
        write(errmsg,*) "xknot(ix) <= xknot(ix+1) required: ix,xknot(ix),xknot(ix+1)=", &
             ix, xknot(ix), xknot(ix+1)
        ierr = 1
        return
     endif
     if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
  end do

  if(leftx .eq. 0) then
     write(errmsg,*) "ix with xknot(ix) <= x < xknot(ix+1) required: x=", x
     ierr = 2
     return
  endif

  lefty = 0

  do iy = 1, ny+ky-1
     if (yknot(iy) .gt. yknot(iy+1)) then
        write(errmsg,*) "yknot(iy) <= yknot(iy+1) required: iy,yknot(iy),yknot(iy+1)=", &
             iy,yknot(iy),yknot(iy+1)
        ierr = 3
        return
     endif
     if((yknot(iy) .le. y) .and. (y .lt. yknot(iy+1))) lefty = iy
  end do

  if(lefty .eq. 0) then
     write(errmsg,*) "iy with yknot(iy) <= y < yknot(iy+1) required: y=", y
     ierr = 4
     return
  endif

  do iky = 1, ky
     work(iky) = dbsdca(0,x,kx,xknot,nx,bcoef(1,lefty-ky+iky),leftx)
  end do

  dbs2vl = dbsval(y,ky,yknot(lefty-ky+1),ky,work, ierr)

  end function dbs2vl


  !==================================================================
  function dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef, ierr)
  !==================================================================
  !
  !  Evaluates the derivative of a two-dimensional tensor-product spline,
  !  given its tensor-product B-spline representation.
  !
  !   iderx  - order of the derivative in the x-direction.  (input)
  !   idery  - order of the derivative in the y-direction.  (input)
  !   x      - x-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   y      - y-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   kx     - order of the spline in the x-direction.  (input)
  !   ky     - order of the spline in the y-direction.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence in the x-direction.  (input)
  !            xknot must be nondecreasing.
  !   yknot  - array of length ny+ky containing the knot
  !            sequence in the y-direction.  (input)
  !            yknot must be nondecreasing.
  !   nx     - number of B-spline coefficients in the x-direction.
  !            (input)
  !   ny     - number of B-spline coefficients in the y-direction.
  !            (input)
  !   bcoef  - array of length nx*ny containing the
  !            tensor-product B-spline coefficients.  (input)
  !            bscoef is treated internally as a matrix of size nx
  !            by ny.
  !   dbs2dr  - value of the (iderx,idery) derivative of the spline at
  !            (x,y).  (output)
  !
  implicit none

  integer, intent(in)                    :: iderx, idery
  integer, intent(in)                    :: kx, nx, ky, ny
  real(dp)                               :: dbs2dr
  real(dp), intent(in)                   :: x, y
  real(dp), dimension(nx+kx), intent(in) :: xknot
  real(dp), dimension(ny+ky), intent(in) :: yknot
  real(dp), dimension(nx,ny), intent(in) :: bcoef
  integer, intent(out)                   :: ierr
  integer                       :: ix, iy, iky, nintx, ninty
  real(dp), dimension(ky) :: work

  routine = 'dbs2dr'
  ierr = 0
  dbs2dr = 0.0_dp
  !
  !     check if knot(i) <= knot(i+1) and calculation of i so that
  !     knot(i) <= x < knot(i+1)
  !

  nintx = 0

  do ix = 1, nx+kx-1
     if (xknot(ix) .gt. xknot(ix+1)) then
        write(errmsg,*) "xknot(ix) <= xknot(ix+1) required: ix,xknot(ix),xknot(ix+1)=", &
             ix, xknot(ix), xknot(ix+1)
        ierr = 1
        return
     endif
     if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) nintx = ix
  end do

  if(nintx .eq. 0) then
     write(errmsg,*) "ix with xknot(ix) <= x < xknot(ix+1) required: x=", x
     ierr = 2
     return
  endif

  ninty = 0

  do iy = 1, ny+ky-1
     if (yknot(iy) .gt. yknot(iy+1)) then
        write(errmsg,*) "yknot(iy) <= yknot(iy+1) required: iy,yknot(iy),yknot(iy+1)=", &
             iy,yknot(iy),yknot(iy+1)
        ierr = 3
        return
     endif
     if ((yknot(iy) .le. y) .and. (y .lt. yknot(iy+1))) ninty = iy
  end do

  if(ninty .eq. 0) then
     write(errmsg,*) "iy with yknot(iy) <= y < yknot(iy+1) required: y=", y
     ierr = 4
     return
  endif

  do iky = 1, ky
     work(iky) =  dbsdca(iderx,x,kx,xknot,nx,bcoef(1,ninty-ky+iky),nintx)
  end do

  dbs2dr = dbsder(idery,y,ky,yknot(ninty-ky+1),ky,work, ierr)

  end function dbs2dr


  !==================================================================
  subroutine dbs2gd(iderx,idery,nxvec,xvec,nyvec,yvec,kx,ky,xknot,yknot,&
       & nx,ny,bcoef,val,ldf, ierr)
  !==================================================================
  !
  !  Evaluates the derivative of a two-dimensional tensor-product spline,
  !  given its tensor-product B-spline representation on a grid.
  !
  !   iderx   - order of the derivative in the x-direction.  (input)
  !   idery   - order of the derivative in the y-direction.  (input)
  !   nxvec   - number of grid points in the x-direction.  (input)
  !   xvec    - array of length nx containing the x-coordinates at
  !             which the spline is to be evaluated.  (input)
  !             the points in xvec should be strictly increasing.
  !   nyvec   - number of grid points in the y-direction.  (input)
  !   yvec    - array of length ny containing the y-coordinates at
  !             which the spline is to be evaluated.  (input)
  !             the points in yvec should be strictly increasing.
  !   kx      - order of the spline in the x-direction.  (input)
  !   ky      - order of the spline in the y-direction.  (input)
  !   xknot   - array of length nx+kx containing the knot
  !             sequence in the x-direction.  (input)
  !             xknot must be nondecreasing.
  !   yknot   - array of length ny+ky containing the knot
  !             sequence in the y-direction.  (input)
  !             yknot must be nondecreasing.
  !   nx      - number of B-spline coefficients in the x-direction.
  !             (input)
  !   ny      - number of B-spline coefficients in the y-direction.
  !             (input)
  !   bcoef   - array of length nx*ny containing the
  !             tensor-product B-spline coefficients.  (input)
  !             bscoef is treated internally as a matrix of size nx
  !             by ny.
  !   val     - value of the (iderx,idery) derivative of the spline on
  !             the nx by ny grid.  (output)
  !             value(i,j) contains the derivative of the spline at the
  !             point (xvec(i),yvec(j)).
  !   ldf     - leading dimension of value exactly as specified in the
  !             dimension statement of the calling program.  (input)
  !
  implicit none

  integer, intent(in)                           :: iderx, idery
  integer, intent(in)                           :: nxvec, nyvec
  integer, intent(in)                           :: kx, nx, ky, ny
  integer, intent(in)                           :: ldf

  real(dp), dimension(nxvec), intent(in)  :: xvec
  real(dp), dimension(nyvec), intent(in)  :: yvec
  real(dp), dimension(nx+kx), intent(in)  :: xknot
  real(dp), dimension(ny+ky), intent(in)  :: yknot
  real(dp), dimension(nx,ny), intent(in)  :: bcoef
  real(dp), dimension(ldf,*), intent(out) :: val
  integer, intent(out) :: ierr

  integer                                     :: i, ik, il, ix, iy, ikx, iky
  integer, dimension(nxvec)                   :: leftx
  integer, dimension(nyvec)                   :: lefty
  real(dp), dimension(nxvec,kx)         :: dl, dr
  real(dp), dimension(max(nxvec,nyvec)) :: save1
  real(dp), dimension(nxvec,kx)         :: biatx
  real(dp), dimension(nyvec,ky)         :: biaty
  real(dp), dimension(max(nxvec,nyvec)) :: term
  real(dp), dimension(ky)               :: work

  logical :: same,next

  routine = 'dbs2gd'
  ierr = 0
  leftx(1) = 0

  call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

  do ix = 2, nxvec
     leftx(ix) = leftx(ix-1)
     same = (xknot(leftx(ix)) .le. xvec(ix))                                &
          &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
     if(.not. same ) then
        leftx(ix) = leftx(ix) + 1
        next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
             &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
        if (.not. next)                                                     &
             &           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
     endif
  end do

  do i = 1, nx+kx-1
     if (xknot(i) .gt. xknot(i+1)) then
        write(errmsg,*) "xknot(i) <= xknot(i+1) required: i,xknot(i),xknot(i+1)=", &
                        i, xknot(i), xknot(i+1)
        ierr = 1
        return
     endif
  end do

  do i = 1, nxvec
     if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
        write(errmsg,*) "ix with xknot(ix) <= x < xknot(ix+1) required: x=", xvec(i)
        ierr = 2
        return
     endif
  end do

  lefty(1) = 0

  call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

  do iy = 2, nyvec
     lefty(iy) = lefty(iy-1)
     same = (yknot(lefty(iy)) .le. yvec(iy))                                &
          &        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
     if(.not. same ) then
        lefty(iy) = lefty(iy) + 1
        next      = (yknot(lefty(iy)) .le. yvec(iy))                        &
             &           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
        if (.not. next) call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
     endif
  end do

  do i = 1, ny+ky-1
     if (yknot(i) .gt. yknot(i+1)) then
        write(errmsg,*) "yknot(i) <= yknot(i+1) required: i,yknot(i),yknot(i+1)=", &
                        i, yknot(i), yknot(i+1)
        ierr = 3
        return
     endif
  end do

  do i = 1, nyvec
     if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
        write(errmsg,*) "iy with yknot(iy) <= y < yknot(iy+1) required: y=", yvec(i)
        ierr = 4
        return
     endif
  end do

  if ((iderx .eq. 0) .and. (idery .eq. 0)) then

     do ix = 1,nxvec
        biatx(ix,1) = 1._dp
     end do

     do ik = 1, kx-1
        do ix = 1,nxvec
           dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
           dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
           save1(ix) = 0._dp
        end do

        do il = 1,ik
           do ix = 1,nxvec
              term(ix)     = biatx(ix,il)                                   &
                   &                 / (dr(ix,il) + dl(ix,ik+1-il))
              biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
              save1(ix)    = dl(ix,ik+1-il) * term(ix)
           end do
        end do

        do ix = 1, nxvec
           biatx(ix,ik+1) = save1(ix)
        end do
     end do

     do iy = 1, nyvec
        biaty(iy,1) = 1._dp
     end do

     do ik = 1, ky-1
        do iy = 1, nyvec
           dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
           dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
           save1(iy) = 0._dp
        end do

        do il = 1, ik
           do iy = 1,nyvec
              term(iy)     = biaty(iy,il)                                   &
                   &                 / (dr(iy,il) + dl(iy,ik+1-il))
              biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
              save1(iy)    = dl(iy,ik+1-il) * term(iy)
           end do
        end do

        do iy = 1, nyvec
           biaty(iy,ik+1) = save1(iy)
        end do
     end do

     do iy = 1, nyvec
        do ix = 1, nxvec
           val(ix,iy) = 0.0_dp
        end do
     end do

     do iky = 1, ky
        do ikx = 1, kx
           do iy = 1, nyvec
              do ix = 1, nxvec
                 val(ix,iy) = val(ix,iy)                                    &
                      & + biatx(ix,ikx) * biaty(iy,iky)                     &
                      & * bcoef(leftx(ix)-kx+ikx,lefty(iy)-ky+iky)
              end do
           end do
        end do
     end do

  elseif (((iderx .ge. 1) .or. (idery .ge. 1))                              &
       &  .and. ( (iderx .lt. kx) .and. (idery .lt. ky))) then

     do iy = 1, nyvec
        do ix = 1, nxvec
           do iky = 1, ky
              work(iky) = dbsdca(iderx,xvec(ix),kx,xknot,nx,                &
                   &             bcoef(1,lefty(iy)-ky+iky),leftx(ix))
           end do
           val(ix,iy) = dbsder(idery,yvec(iy),ky,                           &
                &              yknot(lefty(iy)-ky+1),ky,work, ierr)
        end do
     end do

  else

     do iy = 1, nyvec
        do ix = 1, nxvec
           val(ix,iy) = 0.0_dp
        end do
     end do

  endif

  end subroutine dbs2gd


  !==================================================================
  subroutine dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,kx,ky,kz, &
       & xknot,yknot,zknot,bcoef, ierr)
  !==================================================================
  !
  !  Computes a three-dimensional tensor-product spline interpolant,
  !  returning the tensor-product B-spline coefficients.
  !
  !   nx      - number of data points in the x-direction.  (input)
  !   xvec    - array of length nxdata containing the data points in
  !             the x-direction.  (input)
  !             xdata must be increasing.
  !   ny      - number of data points in the y-direction.  (input)
  !   yvec    - array of length nydata containing the data points in
  !             the y-direction.  (input)
  !             ydata must be increasing.
  !   nz      - number of data points in the z-direction.  (input)
  !   zvec    - array of length nzdata containing the data points in
  !             the z-direction.  (input)
  !             zdata must be increasing.
  !   xyzdata - array of size nx by ny by nz containing the
  !             values to be interpolated.  (input)
  !             xyzdata(i,j,k) contains the value at
  !             (xvec(i),yvec(j),zvec(k)).
  !   ldf     - leading dimension of fdata exactly as specified in the
  !             dimension statement of the calling program.  (input)
  !   mdf     - middle dimension of fdata exactly as specified in the
  !             dimension statement of the calling program.  (input)
  !   kx      - order of the spline in the x-direction.  (input)
  !             kxord must be less than or equal to nxdata.
  !   ky      - order of the spline in the y-direction.  (input)
  !             kyord must be less than or equal to nydata.
  !   kz      - order of the spline in the z-direction.  (input)
  !             kzord must be less than or equal to nzdata.
  !   xknot   - array of length nx+kx containing the knot
  !             sequence in the x-direction.  (input)
  !             xknot must be nondecreasing.
  !   yknot   - array of length ny+ky containing the knot
  !             sequence in the y-direction.  (input)
  !             yknot must be nondecreasing.
  !   zknot   - array of length nz+kz containing the knot
  !             sequence in the z-direction.  (input)
  !             zknot must be nondecreasing.
  !   bcoef   - array of length nx*ny*nz containing the
  !             tensor-product B-spline coefficients.  (output)
  !             bscoef is treated internally as a matrix of size nx
  !             by ny by nz.
  !
  implicit none

  integer, intent(in) :: nx, ny, nz, kx, ky, kz
  integer, intent(in) :: ldf, mdf

  real(dp), dimension(nx), intent(in)         :: xvec
  real(dp), dimension(ny), intent(in)         :: yvec
  real(dp), dimension(nz), intent(in)         :: zvec
  real(dp), dimension(nx+kx), intent(in)      :: xknot
  real(dp), dimension(ny+ky), intent(in)      :: yknot
  real(dp), dimension(nz+kz), intent(in)      :: zknot
  real(dp), dimension(ldf,mdf,nz), intent(in) :: xyzdata
  real(dp), dimension(nx,ny,nz), intent(out)  :: bcoef

  integer                                :: iz
  real(dp), dimension(nx,ny,nz)    :: work1
  real(dp), dimension(nz)          :: work2
  real(dp), dimension((2*kz-1)*nz) :: work3
  integer, intent(out) :: ierr


  call spli3d(zvec,ldf,mdf,xyzdata,zknot,nz,kz,nx,ny,work2,work3,work1,     &
       &     nx,ny,nz, ierr)
  if (ierr /= 0) return

  do iz = 1, nz
     call dbs2in(nx,xvec,ny,yvec,work1(1,1,iz),nx,kx,ky,xknot,yknot,        &
          &        bcoef(1,1,iz), ierr)
     if (ierr /= 0) return
  end do

  end subroutine dbs3in


  !==================================================================
  subroutine spli3d(xyzvec,ldf,mdf,xyzdata,xyzknot,n,k,m,l,work2,work3, &
       & bcoef,nx,ny,nz, ierr)
  !==================================================================
  implicit none

  integer, intent(in)                               :: ldf, mdf, n, k, m, l
  integer, intent(in)                               :: nx, ny, nz
  real(dp), dimension(n), intent(in)          :: xyzvec
  real(dp), dimension(n+k), intent(in)        :: xyzknot
  real(dp), dimension(ldf,mdf,*), intent(in)  :: xyzdata
  real(dp), dimension(nx,ny,nz), intent(out)  :: bcoef
  real(dp), dimension(n), intent(out)         :: work2
  real(dp), dimension((2*k-1)*n), intent(out) :: work3
  integer, intent(out) :: ierr

  integer        :: np1, km1, kpkm2, left, lenq, i, ilp1mx, j, jj, iflag, in
  real(dp) :: xyzveci

  routine = 'spli3d'
  ierr = 0

  np1   = n + 1
  km1   = k - 1
  kpkm2 = 2 * km1
  left  = k
  lenq  = n * (k + km1)

  do i = 1, lenq
     work3(i) = 0._dp
  end do

  do i = 1, n
     xyzveci = xyzvec(i)
     ilp1mx  = min0(i+k,np1)
     left    = max0(left,i)
     if (xyzveci .lt. xyzknot(left)) go to 998
30   if (xyzveci .lt. xyzknot(left+1)) go to 40
     left = left + 1
     if (left .lt. ilp1mx) go to 30
     left = left - 1
     if (xyzveci .gt. xyzknot(left+1)) go to 998
40   call bsplvb(xyzknot,n+k,k,1,xyzveci,left,work2)
     jj = i - left + 1 + (left - k) * (k + km1)
     do j = 1, k
        jj    = jj + kpkm2
        work3(jj) = work2(j)
     end do
  end do

  call banfac(work3,k+km1,n,km1,km1,iflag)

  if (iflag .ne. 1) then
     write(errmsg,*) "no solution of linear equation system"
     ierr = 1
     return
  end if

  do j = 1, l
     do i = 1, m
        do in = 1, n
           work2(in) = xyzdata(i,j,in)
        end do

        call banslv(work3,k+km1,n,km1,km1,work2)

        do in = 1, n
           bcoef(i,j,in) = work2(in)
        end do

     end do
  end do

  return

998 write(errmsg,*) "i with knot(i) <= x/y/z < knot(i+1) required: xyzknot(1),xyzknot(n+k)=", &
         xyzknot(1), xyzknot(n+k)
  ierr = 2
  return

  end subroutine spli3d


  !==================================================================
  function dbs3vl(x,y,z,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef, ierr)
  !==================================================================
  !
  !  Evaluates a three-dimensional tensor-product spline, given its
  !  tensor-product B-spline representation.
  !
  !   x      - x-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   y      - y-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   z      - z-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   kx     - order of the spline in the x-direction.  (input)
  !   ky     - order of the spline in the y-direction.  (input)
  !   kz     - order of the spline in the z-direction.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence in the x-direction.  (input)
  !            xknot must be nondecreasing.
  !   yknot  - array of length ny+ky containing the knot
  !            sequence in the y-direction.  (input)
  !            yknot must be nondecreasing.
  !   zknot  - array of length nz+kz containing the knot
  !            sequence in the z-direction.  (input)
  !            zknot must be nondecreasing.
  !   nx     - number of B-spline coefficients in the x-direction.
  !            (input)
  !   ny     - number of B-spline coefficients in the y-direction.
  !            (input)
  !   nz     - number of B-spline coefficients in the z-direction.
  !            (input)
  !   bcoef  - array of length nx*ny*nz containing the
  !            tensor-product B-spline coefficients.  (input)
  !            bscoef is treated internally as a matrix of size nx
  !            by ny by nz.
  !   dbs3vl - value of the spline at (x,y,z).  (output)
  !
  implicit none

  integer, intent(in)                             :: nx, ny, nz, kx, ky, kz
  real(dp), intent(in)                      :: x, y, z
  real(dp), dimension(nx+kx), intent(in)    :: xknot
  real(dp), dimension(ny+ky), intent(in)    :: yknot
  real(dp), dimension(nz+kz), intent(in)    :: zknot
  real(dp), dimension(nx,ny,nz), intent(in) :: bcoef
  real(dp)                                  :: dbs3vl
  integer, intent(out) :: ierr
  integer                       :: iz, nintz
  real(dp), dimension(kz) :: work

  routine = 'dbs3vl'
  ierr = 0
  dbs3vl = 0.0_dp
  !
  !     check if knot(i) <= knot(i+1) and calculation of i so that
  !     knot(i) <= x < knot(i+1)
  !

  nintz = 0

  do iz = 1, nz+kz-1
     if (zknot(iz) .gt. zknot(iz + 1)) then
        write(errmsg,*) "zknot(iz) <= zknot(iz+1) required: iz,zknot(iz),zknot(iz+1)=", &
                        iz, zknot(iz), zknot(iz+1)
        ierr = 1
        return
     endif
     if((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
  end do

  if(nintz .eq. 0) then
     write(errmsg,*) "iz with zknot(iz) <= z < zknot(iz+1) required: zknot(iz),z,zknot(iz+1)=", &
                     zknot(iz), z, zknot(iz+1)
     ierr = 2
     return
  endif

  do iz = 1, kz
     work(iz) = dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef(1,1,nintz-kz+iz), ierr)
  end do
  if (ierr /= 0) return

  dbs3vl = dbsval(z,kz,zknot(nintz-kz+1),kz,work, ierr)

  end function dbs3vl


  !==================================================================
  function dbs3dr(iderx,idery,iderz,x,y,z,kx,ky,kz,xknot,yknot,zknot,         &
       & nx,ny,nz,bcoef, ierr)
  !==================================================================
  !
  !  Evaluates the derivative of a three-dimensional tensor-product spline,
  !  given its tensor-product B-spline representation.
  !
  !   iderx  - order of the x-derivative.  (input)
  !   idery  - order of the y-derivative.  (input)
  !   iderz  - order of the z-derivative.  (input)
  !   x      - x-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   y      - y-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   z      - z-coordinate of the point at which the spline is to be
  !            evaluated.  (input)
  !   kx     - order of the spline in the x-direction.  (input)
  !   ky     - order of the spline in the y-direction.  (input)
  !   kz     - order of the spline in the z-direction.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence in the x-direction.  (input)
  !            xknot must be nondecreasing.
  !   yknot  - array of length ny+ky containing the knot
  !            sequence in the y-direction.  (input)
  !            yknot must be nondecreasing.
  !   zknot  - array of length nz+kz containing the knot
  !            sequence in the z-direction.  (input)
  !            zknot must be nondecreasing.
  !   nx     - number of B-spline coefficients in the x-direction.
  !            (input)
  !   ny     - number of B-spline coefficients in the y-direction.
  !            (input)
  !   nz     - number of B-spline coefficients in the z-direction.
  !            (input)
  !   bcoef  - array of length nx*ny*nz containing the
  !            tensor-product B-spline coefficients.  (input)
  !            bscoef is treated internally as a matrix of size nx
  !            by ny by nz.
  !   dbs3dr - value of the (iderx,idery,iderz) derivative of the
  !            spline at (x,y,z).  (output)
  !
  implicit none

  integer, intent(in)                              :: iderx, idery, iderz
  integer, intent(in)                              :: nx, ny, nz, kx, ky, kz
  real(dp), intent(in)                       :: x, y, z
  real(dp), dimension(nx+kx), intent(in)     :: xknot
  real(dp), dimension(ny+ky), intent(in)     :: yknot
  real(dp), dimension(nz+kz), intent(in)     :: zknot
  real(dp), dimension(nx,ny,nz), intent(in)  :: bcoef
  real(dp)                                   :: dbs3dr
  integer, intent(out) :: ierr
  integer                       :: iz, nintz
  real(dp), dimension(kz) :: work

  routine = 'dbs3dr'
  ierr = 0
  dbs3dr = 0.0_dp
  !
  !     check if knot(i) <= knot(i+1) and calculation of i so that
  !     knot(i) <= x < knot(i+1)
  !
  nintz = 0

  do iz = 1, nz+kz-1
     if (zknot(iz) .gt. zknot(iz + 1)) then
        write(errmsg,*) "zknot(iz) <= zknot(iz+1) required: iz,zknot(iz),zknot(iz+1)=", &
                        iz, zknot(iz), zknot(iz+1)
        ierr = 1
        return
     endif
     if((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
  end do

  if(nintz .eq. 0) then
     write(errmsg,*) "iz with zknot(iz) <= z < zknot(iz+1) required: zknot(iz),z,zknot(iz+1)=", &
                     zknot(iz), z, zknot(iz+1)
     ierr = 2
     return
  endif

  do iz = 1, kz
     work(iz) = dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,             &
          &        bcoef(1,1,nintz-kz+iz),ierr)
  end do
  if (ierr /= 0) return

  dbs3dr = dbsder(iderz,z,kz,zknot(nintz-kz+1),kz,work,ierr)

  end function dbs3dr


  !==================================================================
  subroutine dbs3gd(iderx,idery,iderz,nxvec,xvec,nyvec,yvec,nzvec,zvec,       &
       & kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef,val,ldf,mdf, ierr)
  !==================================================================
  !
  !  Evaluates the derivative of a three-dimensional tensor-product spline,
  !  given its tensor-product B-spline representation on a grid.
  !
  !   iderx  - order of the x-derivative.  (input)
  !   idery  - order of the y-derivative.  (input)
  !   iderz  - order of the z-derivative.  (input)
  !   nx     - number of grid points in the x-direction.  (input)
  !   xvec   - array of length nx containing the x-coordinates at
  !            which the spline is to be evaluated.  (input)
  !            the points in xvec should be strictly increasing.
  !   ny     - number of grid points in the y-direction.  (input)
  !   yvec   - array of length ny containing the y-coordinates at
  !            which the spline is to be evaluated.  (input)
  !            the points in yvec should be strictly increasing.
  !   nz     - number of grid points in the z-direction.  (input)
  !   zvec   - array of length nz containing the z-coordinates at
  !            which the spline is to be evaluated.  (input)
  !            the points in yvec should be strictly increasing.
  !   kx     - order of the spline in the x-direction.  (input)
  !   ky     - order of the spline in the y-direction.  (input)
  !   kz     - order of the spline in the z-direction.  (input)
  !   xknot  - array of length nx+kx containing the knot
  !            sequence in the x-direction.  (input)
  !            xknot must be nondecreasing.
  !   yknot  - array of length ny+ky containing the knot
  !            sequence in the y-direction.  (input)
  !            yknot must be nondecreasing.
  !   zknot  - array of length nz+kz containing the knot
  !            sequence in the z-direction.  (input)
  !            zknot must be nondecreasing.
  !   nx     - number of B-spline coefficients in the x-direction.
  !            (input)
  !   ny     - number of B-spline coefficients in the y-direction.
  !            (input)
  !   nz     - number of B-spline coefficients in the z-direction.
  !            (input)
  !   bcoef  - array of length nx*ny*nz containing the
  !            tensor-product B-spline coefficients.  (input)
  !            bscoef is treated internally as a matrix of size nx
  !            by ny by nz.
  !   val    - array of size nx by ny by nz containing the values of
  !            the (iderx,idery,iderz) derivative of the spline on the
  !            nx by ny by nz grid.  (output)
  !            value(i,j,k) contains the derivative of the spline at
  !            the point (xvec(i), yvec(j), zvec(k)).
  !   ldf    - leading dimension of value exactly as specified in the
  !            dimension statement of the calling program.  (input)
  !   mdf    - middle dimension of value exactly as specified in the
  !            dimension statement of the calling program.  (input)
  !
  implicit none

  integer, intent(in)                               :: iderx, idery, iderz
  integer, intent(in)                               :: nxvec, nyvec, nzvec
  integer, intent(in)                               :: kx, nx, ky, ny, kz, nz
  integer, intent(in)                               :: ldf,mdf

  real(dp), dimension(nxvec), intent(in)      :: xvec
  real(dp), dimension(nyvec), intent(in)      :: yvec
  real(dp), dimension(nzvec), intent(in)      :: zvec
  real(dp), dimension(nx+kx), intent(in)      :: xknot
  real(dp), dimension(ny+ky), intent(in)      :: yknot
  real(dp), dimension(nz+kz), intent(in)      :: zknot
  real(dp), dimension(nx,ny,nz), intent(in)   :: bcoef
  real(dp), dimension(ldf,mdf,*), intent(out) :: val
  integer, intent(out) :: ierr

  integer                                           :: i, ik, il, ix, iy, iz
  integer                                           :: ikx, iky, ikz
  integer, dimension(nxvec)                         :: leftx
  integer, dimension(nyvec)                         :: lefty
  integer, dimension(nzvec)                         :: leftz
  real(dp), dimension(nxvec,kx)               :: biatx
  real(dp), dimension(nyvec,ky)               :: biaty
  real(dp), dimension(nzvec,kz)               :: biatz
  real(dp), dimension(max(nxvec,nyvec,nzvec)) :: term, save1

  real(dp), dimension(max(nxvec,nyvec,nzvec), max(kx,ky,kz)) :: dl, dr

  logical :: same,next

  routine = 'dbs3gd'
  ierr = 0

  do i = 1, nx+kx-1
     if (xknot(i) .gt. xknot(i+1)) then
        write(errmsg,*) "xknot(i) <= xknot(i+1) required: i,xknot(i),xknot(i+1)=", &
                        i, xknot(i), xknot(i+1)
        ierr = 1
        return
     endif
  end do

  do i = 1, nxvec
     if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
        write(errmsg,*) "ix with xknot(ix) <= x < xknot(ix+1) required: x=", xvec(i)
        ierr = 2
        return
     endif
  end do

  leftx(1) = 0

  call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

  do ix = 2, nxvec
     leftx(ix) = leftx(ix-1)
     same = (xknot(leftx(ix)) .le. xvec(ix))                                &
          &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
     if(.not. same ) then
        leftx(ix) = leftx(ix) + 1
        next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
             &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
        if (.not. next) call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
     endif
  end do

  do i = 1, ny+ky-1
     if (yknot(i) .gt. yknot(i+1)) then
        write(errmsg,*) "yknot(i) <= yknot(i+1) required: i,yknot(i),yknot(i+1)=", &
                        i, yknot(i), yknot(i+1)
        ierr = 3
        return
     endif
  end do

  do i = 1, nyvec
     if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
        write(errmsg,*) "iy with yknot(iy) <= y < yknot(iy+1) required: y=", yvec(i)
        ierr = 4
        return
     endif
  end do

  lefty(1) = 0

  call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

  do iy = 2, nyvec
     lefty(iy) = lefty(iy-1)
     same = (yknot(lefty(iy)) .le. yvec(iy))                                &
          &        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
     if(.not. same ) then
        lefty(iy) = lefty(iy) + 1
        next      = (yknot(lefty(iy)) .le. yvec(iy))                        &
             &           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
        if (.not. next) call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
     endif
  end do

  do i = 1,nz+kz-1
     if (zknot(i) .gt. zknot(i+1)) then
        write(errmsg,*) "zknot(i) <= zknot(i+1) required: i,zknot(i),zknot(i+1)=", &
                        i, zknot(i), zknot(i+1)
        ierr = 5
        return
     endif
  end do

  do i = 1, nzvec
     if ((zvec(i).lt.zknot(1)).or.(zvec(i).gt.zknot(nz+kz))) then
        write(errmsg,*) "iz with zknot(iz) <= z < zknot(iz+1) required: z=", zvec(i)
        ierr = 6
        return
     endif
  end do

  leftz(1) = 0

  call huntn(zknot,nz+kz,kz,zvec(1),leftz(1))

  do iz = 2, nzvec
     leftz(iz) = leftz(iz-1)
     same = (zknot(leftz(iz)) .le. zvec(iz))                                &
          &        .and. (zvec(iz) .le. zknot(leftz(iz)+1))
     if(.not. same ) then
        leftz(iz) = leftz(iz) + 1
        next      = (zknot(leftz(iz)) .le. zvec(iz))                        &
             &           .and. (zvec(iz) .le. zknot(leftz(iz)+1))
        if (.not. next) call huntn(zknot,nz+kz,kz,zvec(iz),leftz(iz))
     endif
  end do

  if ((iderx .eq. 0) .and. (idery .eq. 0) .and. (iderz .eq.0)) then

     do ix = 1, nxvec
        biatx(ix,1) = 1.0_dp
     end do

     do ik = 1, kx-1
        do ix = 1, nxvec
           dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
           dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
           save1(ix) = 0._dp
        end do

        do il = 1, ik
           do ix = 1, nxvec
              term(ix)     = biatx(ix,il) / (dr(ix,il) + dl(ix,ik+1-il))
              biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
              save1(ix)    = dl(ix,ik+1-il) * term(ix)
           end do
        end do

        do ix = 1, nxvec
           biatx(ix,ik+1) = save1(ix)
        end do
     end do

     do iy = 1, nyvec
        biaty(iy,1) = 1.0_dp
     end do

     do ik = 1, ky-1
        do iy = 1, nyvec
           dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
           dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
           save1(iy) = 0._dp
        end do

        do il = 1,ik
           do iy = 1,nyvec
              term(iy)     = biaty(iy,il) / (dr(iy,il) + dl(iy,ik+1-il))
              biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
              save1(iy)    = dl(iy,ik+1-il) * term(iy)
           end do
        end do

        do iy = 1,nyvec
           biaty(iy,ik+1) = save1(iy)
        end do
     end do

     do iz = 1,nzvec
        biatz(iz,1) = 1.0_dp
     end do

     do ik = 1, kz-1
        do iz = 1, nzvec
           dr(iz,ik) = zknot(leftz(iz)+ik) - zvec(iz)
           dl(iz,ik) = zvec(iz) - zknot(leftz(iz)+1-ik)
           save1(iz) = 0._dp
        end do

        do il = 1, ik
           do iz = 1, nzvec
              term(iz)     = biatz(iz,il) / (dr(iz,il) + dl(iz,ik+1-il))
              biatz(iz,il) = save1(iz) + dr(iz,il) * term(iz)
              save1(iz)    = dl(iz,ik+1-il) * term(iz)
           end do
        end do

        do iz = 1, nzvec
           biatz(iz,ik+1) = save1(iz)
        end do
     end do

     do iz = 1,nzvec
        do iy = 1,nyvec
           do ix = 1,nxvec
              val(ix,iy,iz) = 0.0_dp
           end do
        end do
     end do

     do ikz = 1, kz
        do iky = 1, ky
           do ikx = 1, kx
              do iz = 1, nzvec
                 do iy = 1, nyvec
                    do ix = 1, nxvec
                       val(ix,iy,iz) = val(ix,iy,iz)                        &
                            &  + biatx(ix,ikx) * biaty(iy,iky)              &
                            &  * biatz(iz,ikz)                              &
                            &  * bcoef(leftx(ix)-kx+ikx,                    &
                            &          lefty(iy)-ky+iky,leftz(iz)-kz+ikz)
                    end do
                 end do
              end do
           end do
        end do
     end do

  else

     do iz = 1, nzvec
        do iy = 1, nyvec
           do ix = 1, nxvec
              val(ix,iy,iz) = dbs3dr(iderx,idery,iderz,xvec(ix),            &
                   &  yvec(iy),zvec(iz),kx,ky,kz,xknot,yknot,               &
                   &  zknot,nx,ny,nz,bcoef, ierr)
           end do
        end do
     end do

  endif

  end subroutine dbs3gd


  !==================================================================
  ! Internal routines
  !==================================================================
  subroutine bsplvb(t,n,jhigh,idx,x,left,biatx)
    implicit none
    integer, intent(in) :: n, jhigh, idx, left
    real(dp), intent(in)                    :: x
    real(dp), dimension(n), intent(in)      :: t
    real(dp), dimension(jhigh), intent(out) :: biatx
    integer                          :: j = 1
    integer                          :: i, jp1
    real(dp)                   :: saved, term
    real(dp), dimension(jhigh) :: dl, dr


    if (idx .eq. 1) then
       j = 1
       biatx(1) = 1.0_dp
       if (j .ge. jhigh) return
    end if

20  jp1 = j + 1

    dr(j) = t(left+j) - x
    dl(j) = x - t(left+1-j)
    saved = 0._dp

    do i = 1, j
       term     = biatx(i) / (dr(i) + dl(jp1-i))
       biatx(i) = saved + dr(i) * term
       saved    = dl(jp1-i) * term
    end do

    biatx(jp1) = saved
    j          = jp1

    if (j .lt. jhigh) go to 20
  end subroutine bsplvb


  subroutine banfac(w,nroww,nrow,nbandl,nbandu,iflag)
    implicit none
    integer, intent(in)                                  :: nroww,nrow
    integer, intent(in)                                  :: nbandl,nbandu
    integer, intent(out)                                 :: iflag
    real(dp), dimension(nroww,nrow), intent(inout) :: w
    real(dp) :: pivot, factor
    integer  :: middle, nrowm1, jmax, kmax, ipk, midmk, i, j, k

    iflag  = 1
    middle = nbandu + 1
    nrowm1 = nrow - 1

    if (nrowm1 .lt. 0) goto 999
    if (nrowm1 .eq. 0) goto 900
    if (nrowm1 .gt. 0) goto 10

10  if (nbandl .gt. 0) go to 30

    do i = 1, nrowm1
       if (w(middle,i) .eq. 0._dp) go to 999
    end do

    go to 900

30  if (nbandu .gt. 0) go to 60

    do i = 1, nrowm1
       pivot = w(middle,i)
       if(pivot .eq. 0._dp) go to 999
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do
    end do

    return

60  do i = 1, nrowm1
       pivot = w(middle,i)
       if (pivot .eq. 0._dp) go to 999
       jmax = min0(nbandl,nrow - i)
       do j = 1,jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do

       kmax = min0(nbandu,nrow - i)

       do k = 1, kmax
          ipk    = i + k
          midmk  = middle - k
          factor = w(midmk,ipk)
          do j = 1, jmax
             w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)                  &
                  &              * factor
          end do
       end do
    end do

900 if (w(middle,nrow) .ne. 0._dp) return
999 iflag = 2
  end subroutine banfac


  subroutine banslv(w,nroww,nrow,nbandl,nbandu,b)
    implicit none
    integer, intent(in)                               :: nroww,nrow
    integer, intent(in)                               :: nbandl,nbandu
    real(dp), dimension(nroww,nrow), intent(in) :: w
    real(dp), dimension(nrow), intent(inout)    :: b
    integer :: middle, nrowm1, jmax, i, j

    middle = nbandu + 1
    if (nrow .eq. 1) goto 99
    nrowm1 = nrow - 1
    if (nbandl .eq. 0) goto 30

    do i = 1, nrowm1
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          b(i+j) = b(i+j) - b(i) * w(middle+j,i)
       end do
    end do

30  if (nbandu .gt. 0)  goto 50

    do i = 1, nrow
       b(i) = b(i) / w(1,i)
    end do

    return

50  do i = nrow, 2, -1
       b(i) = b(i)/w(middle,i)
       jmax = min0(nbandu,i-1)
       do j = 1, jmax
          b(i-j) = b(i-j) - b(i) * w(middle-j,i)
       end do
    end do

99  b(1) = b(1) / w(middle,1)
  end subroutine banslv


  subroutine huntn(xx,n,kord,x,jlo)
    implicit none
    integer, intent(in)                :: n, kord
    real(dp), intent(in)               :: x
    real(dp), dimension(n), intent(in) :: xx
    integer, intent(inout)             :: jlo
    integer :: max, null, jhi, jm, inc

    max  = n - kord
    null = kord

    if (jlo.le.null.or.jlo.gt.max) then
       jlo = null
       jhi = max+1
       goto 30
    endif

    inc = 1

    if (x .ge. xx(jlo)) then
10     jhi = jlo + inc
       if (jhi .gt. max) then
          jhi = max + 1
       else if (x .ge. xx(jhi)) then
          jlo = jhi
          inc = inc + inc
          goto 10
       endif
    else
       jhi = jlo
20     jlo = jhi - inc
       if (jlo .le. null) then
          jlo = null
       else if (x .lt. xx(jlo)) then
          jhi = jlo
          inc = inc + inc
          goto 20
       endif
    endif

30  if (jhi-jlo.eq.1) return

    jm = (jhi + jlo) / 2
    if (x .gt. xx(jm)) then
       jlo = jm
    else
       jhi = jm
    endif

    goto 30

  end subroutine huntn



  !==================================================================
  ! error reporting routines
  !==================================================================
  function get_error_routine()
    character(80) :: get_error_routine
    get_error_routine = routine
  end function get_error_routine

  function get_error_message()
    character(256) :: get_error_message
    get_error_message = errmsg
  end function get_error_message
 
end module bspline
