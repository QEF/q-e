!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program plotrho
  !-----------------------------------------------------------------------
  !   2D contour plot - logarithmically or linearly spaced levels
  !                   - Postscript printable output
  !   if " cplot" is called:
  !                   - contour lines plus gray levels
  !                   - negative values are shaded
  !   if "psplot" is called:
  !                   - contour lines of various kinds (solid, dashed, etc)
  !
#include "f_defs.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  ! for spline interpolation using essl toutines
  integer, parameter :: nwrk = 10000
  real(DP) :: wrk (nwrk)
  integer, allocatable :: ityp (:)
  integer :: nxi, nyi, nx, ny, i, j, k, nlevels, na, nat, ierr, ilen
  real(DP), allocatable :: rhoi(:,:), xi(:), yi(:)
  real(DP), allocatable :: rhoo(:,:), x (:), y (:)
  real(DP), allocatable ::  z (:)
  real(DP) :: xmin, xmax, ymin, ymax, rhomin, rhomax, rhoomin, rhoomax
  real(DP) :: xdim, ydim, xs, ys, r0 (3), tau1 (3), tau2 (3)
  real(DP), allocatable :: tau (:,:)
  real(DP) :: at (3, 3), a0
  character (len=256) :: filename, fileout, ans * 1
  logical :: logarithmic_scale


  call get_file ( filename )
  open (unit = 1, file = filename, form = 'formatted', status = 'old')

  read (1, * ) nxi, nyi
  allocate ( xi(0:nxi), yi(0:nyi), rhoi(0:nxi,0:nyi) )
  read (1, * ) (xi (i), i = 0, nxi)
  read (1, * ) (yi (j), j = 0, nyi)
  read (1, * ) ( (rhoi (i, j), i = 0, nxi), j = 0, nyi)
  read (1, * ) r0
  read (1, * ) tau1
  read (1, * ) tau2
  read (1, * ) nat
  if (nat < 0 .or. nat > 1000000) then
     WRITE( stdout, '("Error: unlikely number of atoms ",i4)') nat
     stop
  endif
  allocate (tau (3,nat), ityp(nat) )
  read (1, * ) ( (tau (j, na), j = 1, 3), ityp (na), na = 1, nat)
  read (1, * ) a0
  read (1, * ) at
  close (unit = 1)
  !
  WRITE( stdout, '("r0   : ",3f8.4)') r0
  WRITE( stdout, '("tau1 : ",3f8.4)') tau1
  WRITE( stdout, '("tau2 : ",3f8.4)') tau2
  !
  WRITE( stdout, '("read",i4," atomic positions")') nat
  !      WRITE( stdout,'("Atomic positions:")')
  !      WRITE( stdout,'(3f8.4)') ( (tau(j,na),j=1,3),na=1,nat)
  WRITE( stdout, '("output file > ",$)')

  read (5, '(a)') fileout

  WRITE( stdout, '("Read ",i3," *",i3,"  grid")') nxi+1, nyi+1
#ifdef __ESSL
  !
  ! interpolation implemented only for ESSL routines ...
  !
  WRITE( stdout, '("nx, ny (output) > ",$)')
  read (5, * ) nx, ny
#else
  nx = nxi
  ny = nyi
#endif
  allocate ( x(0:nx), y(0:ny), rhoo(0:nx,0:ny) ) 
  xmin = xi (0)
  xmax = xi (nxi)
  do i = 0, nx
     x (i) = (xi (nxi) - xi (0) ) * DBLE (i) / DBLE (nx)
  enddo

  ymin = yi (0)
  ymax = yi (nyi)
  do i = 0, ny
     y (i) = (yi (nyi) - yi (0) ) * DBLE (i) / DBLE (ny)
  enddo
#ifdef __ESSL
  call dcsin2 (xi, yi, rhoi, nxi + 1, nyi + 1, nxi + 1, x, y, &
               nx + 1, ny + 1, rhoo, nx + 1, wrk, nwrk)
#else
  rhoo (0:nx, 0:ny) = rhoi (0:nx, 0:ny)
#endif
  rhomin = minval (rhoo(0:nx, 0:ny))
  rhomax = maxval (rhoo(0:nx, 0:ny))

  if (rhomin > 0.d0) then
     WRITE( stdout,'("Logarithmic scale (y/n)? > ",$)')
     read (5, '(a)') ans
     logarithmic_scale = ans.ne.'n'.and.ans.ne.'N'
  else
     logarithmic_scale = .false.
  end if
10 continue
  WRITE( stdout, '("Bounds: ",2f12.6)') rhomin, rhomax
  WRITE( stdout, '("min, max, # of levels > ",$)')
  read (5, * ) rhoomin, rhoomax, nlevels
  if ( rhoomax <= rhoomin .or. &
       rhoomin >= rhomax .or. rhoomax <= rhomin ) then
     WRITE( stdout, '("Out of Bounds! try again")')
     go to 10
  end if
  if (nlevels > 1000) then
     WRITE( stdout, '("Are you sure you really need ",i8," levels?")') &
          nlevels
  else if (nlevels < 1) then
     WRITE( stdout, '("Too few levels! assuming 1 level")') 
     nlevels = 1
  end if
  allocate (z(0:nlevels))
  if (logarithmic_scale) then
     do k = 0, nlevels - 1
        z (k) = exp (log (rhoomin) + (log (rhoomax) - log (rhoomin) ) &
             * DBLE (k) / (nlevels - 1) )
     enddo
  else
     do k = 0, nlevels - 1
        z (k) = rhoomin + (rhoomax - rhoomin) * DBLE (k) / (nlevels - 1)
     enddo
  endif

  z (nlevels) = z (nlevels - 1)
  xdim = 15.0d0 * (xmax - xmin) / sqrt ( (xmax - xmin) **2 + (ymax - ymin) **2)
  ydim = 15.0d0 * (ymax - ymin) / sqrt ( (xmax - xmin) **2 + (ymax - ymin) **2)

  xs = 4.0d0
  ys = 3.0d0

  ! uncomment the call to "cplot" if you want contour lines,
  ! plus gray levels and shading for negative values

  call cplot (rhoo, nx, ny, x, xmin, xmax, nx, y, ymin, ymax, &
       ny, nlevels, z, xdim, ydim, xs, ys, filename, fileout)

  ! uncomment the call to "psplot" if you want contour lines
  ! of various kinds: solid, dashed, etc

  ! call psplot ( rhoo, nx, x, nx, y, ny, nlevels, z, xdim, ydim, &
  !      xs, ys, fileout)

  call atomi (nat, tau, ityp, at, a0, r0, tau1, tau2, xdim, ydim)
20 stop
   deallocate (z)
   deallocate (rhoo, x , y )
   deallocate (tau, ityp, z)
   deallocate (rhoi, xi, yi)
end program plotrho
!
!-----------------------------------------------------------------------
subroutine cplot (d, imax, jmax, x, xmin, xmax, iub, y, ymin, &
     ymax, jub, nc, z, xdim, ydim, xs, ys, str, filename)
  !-----------------------------------------------------------------------
  !
  !   draws a contour plot of d(i,j). PostScript output on unit 1
  !   F.Gygi Dec.15 1987 - P. Giannozzi Oct.6 1989 and later
  !   Algorithm by Paul D. Bourke, Byte magazine, june 1987, p. 143
  !   d(0:imax,0:jmax) contains the function to plot
  !   x(0:imax) workspace
  !   xmin,xmax determines the range of the variable x
  !   iub is the number of intervals along the x axis (<=imax)
  !   y(0:jmax) workspace
  !   idem for ymin,ymax and jub (<=jmax)
  !   nc is the number of levels wanted ( <=ncmax),
  !   z(0:nc) are the levels,
  !   xdim and ydim are the physical dimensions of the figure in cm
  !   xs and ys determine a shift of the origin in cm
  USE kinds, only: DP
  implicit none
  integer :: imax, jmax, iub, jub, nc
  real(DP) :: d (0:imax, 0:jmax), x (0:imax), y (0:jmax), z (0:nc)
  real(DP) :: xmin, xmax, ymin, ymax, xdim, ydim, xs, ys

  character (len=*) :: filename, str
  integer, parameter :: ncmax = 19
  integer :: i, j, k
  real(DP) :: gray (0:ncmax), dim
  real(DP) :: cm = 28.453d0, width=0.5d0, gray0=1.0d0, deltagray=0.7d0
  !  cm   : number of points per cm
  !  width: linewidth of the contour plot for PostScript printer

  open (unit = 1, file = filename, status = 'unknown', form = &
       'formatted')
  if (nc > ncmax .or. nc < 1) stop '  nc too big or wrong'
  if (iub > imax .or.iub < 1) stop ' iub too big or wrong'
  if (jub > jmax .or.jub < 1) stop ' jub too big or wrong'
  if (xdim < 3.0d0 .or. ydim < 3.0d0) stop ' really too small!'
  if (xdim > 20.0.or. ydim > 30.0) stop ' really too big!'
  if (abs (xs) > 20.or. abs (ys) > 30) stop ' xs or ys are weird'
  !  initializations for PostScript output
  write (1, '(a)') '%! PS-Adobe-1.0'
  write (1, '("%%BoundingBox:",4f6.1)') xs * cm, ys * cm, (xs + &
       xdim) * cm, (ys + ydim) * cm
  write (1, '(a)') '/localdict 100 dict def'
  write (1, '(a)') 'localdict begin'
  write (1, '(a)') '/cm {28.453 mul} def'
  write (1, '(a)') '/title {('//str//')} def'
  write (1, '(a)') '/Times-Roman findfont 12 scalefont setfont'
  write (1, '(a)') '% cshow prints a centered string at current position'
  write (1, '(a)') '/cshow {gsave dup stringwidth pop 2 div neg 0'
  write (1, '(a)') ' rmoveto show grestore} def'
  write (1, '(a)') '% x1 y1 x2 y2 p :  draws a segment from point 1 to point 2'
  write (1, '(a)') '/p {0 setgray newpath moveto lineto stroke} def'
  write (1, '(a)') '% x1 y1 x2 y2 x3 y3 x4 y4 sn :'
  write (1, '(a)') '% fills the region bounded by points 1 to 4'
  write (1, '(a)') '% with greyscale n'
  !  type of gray for shaded areas
  do k = 0, nc
     gray (k) = gray0 - k * deltagray / nc
     if (k.lt.10) then
        write (1, '("/s",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto lineto fill} def")') k, gray (k)
        write (1, '("/t",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto fill} def")') k, gray (k)
     else
        write (1, '("/u",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto lineto fill} def")') mod (k, 10) , &
             &gray (k)
        write (1, '("/v",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto fill} def")') mod (k, 10) , gray (k &
             &)
     endif

  enddo
  write (1, '(a)') '%%EndPreamble'
  write (1, '(a)') 'gsave'
  write (1, '(1x,f6.2," cm ",f6.2," cm translate")') xs, ys
  write (1, '(a)') '% Uncomment next line if you want a big picture'
  write (1, '(a)') '% 1.8 1.8 scale'
  write (1, '(f7.3," setlinewidth")') width
  write (1, '(a)') '% Comment the next line to remove the title'

  write (1, '(1x,f6.2," cm ",f6.2," cm moveto title cshow")') &
       xdim / 2, ydim + 1.5d0

  call hatch (0.d0, xdim, 0.d0, ydim)
  do i = 0, iub
     x (i) = xdim * DBLE (i) / iub
  enddo
  do j = 0, jub
     y (j) = ydim * DBLE (j) / jub
  enddo

  call conrec (imax, iub, jmax, jub, x, y, d, nc, z)
  !  draw frame of size xdim by ydim
  write (1, '(a)') '1 setlinewidth 0 setgray newpath'
  write (1, '(2f6.1," moveto")') 0.0d0, 0.0d0
  write (1, '(2f6.1," lineto")') xdim * cm, 0.0d0
  write (1, '(2f6.1," lineto")') xdim * cm, ydim * cm
  write (1, '(2f6.1," lineto")') 0.0d0, ydim * cm

  write (1, '(a)') 'closepath stroke'
  !      write (1,'(a)') 'grestore'
  !      write (1,'(a)') '%%Trailer'
  !      write (1,'(a)') 'showpage'
  !      close(1)
  return


end subroutine cplot

subroutine conrec (imax, iub, jmax, jub, x, y, d, nc, z)
  USE kinds, only: DP
  implicit none
  integer :: imax, iub, jmax, jub, nc

  real(DP) :: d (0:imax, 0:jmax), x (0:imax), y (0:jmax), z (0:nc)
  integer, parameter :: ncmax = 19
  character (len=4) :: triangle (0:ncmax), trapez (0:ncmax)
  real(DP) :: h (0:4), xh (0:4), yh (0:4)
  real(DP) :: x1, y1, x2, y2, x3, y3, x4, y4, dx, dy, xx, yy, dmin, dmax
  real(DP) :: cm = 28.453d0
  !   cm : conversion factor from cm to points for PostScript
  integer :: ish (0:4), im (0:3), jm (0:3), castab (0:2, 0:2, 0:2)
  integer :: i, j, k, m, m1, m2, m3, npoint, icase, levelin, nolevel
  data (im (i), i = 0, 3) / 0, 1, 1, 0 /
  data (jm (i), i = 0, 3) / 0, 0, 1, 1 /
  data ( ( (castab (i, j, k), k = 0, 2), j = 0, 2), i = 0, 2) &
       / 0, 0, 8, 0, 2, 5, 7, 6, 9, 0, 3, 4, 1, 3, 1, 4, 3, 0, 9, 6, 7, &
       5, 2, 0, 8, 0, 0 /
 
  dy = (y (jub) - y (0) ) / (nc + 1)
  xx = x (iub) + 1.0d0
  dx = 0.5d0

  write (1, '(a)') '% Start of Color Code'

  call hatch (xx, xx + dx, y (0), y (jub) )
  do k = 0, nc
     yy = y (jub) - k * dy
     write (1, '(8f6.1,$)') xx * cm, yy * cm, (xx + dx) * cm, yy * cm, &
          (xx + dx) * cm, (yy - dy) * cm, xx * cm, (yy - dy) * cm
     if (k < 10) then
        write (triangle (k) , '(" t",i1,1x)') k
        write (trapez (k) , '(" s",i1,1x)') k
     else
        write (triangle (k) , '(" v",i1)') mod (k, 10)
        write (trapez (k) , '(" u",i1)') mod (k, 10)
     endif
     write (1, '(a4)') trapez (k)
     write (1, * ) '0 setgray newpath'
     write (1, '(2f6.1," moveto")') xx * cm, yy * cm
     write (1, '(2f6.1," lineto")') (xx + dx) * cm, yy * cm
     write (1, '(2f6.1," lineto")') (xx + dx) * cm, (yy - dy) &
          * cm
     write (1, '(2f6.1," lineto")') xx * cm, (yy - dy) * cm
     write (1, * ) 'closepath stroke'
     write (1, '(2f6.1," moveto")') (x (iub) + 2.0d0) * cm, (yy - dy / &
          2) * cm
     if (k == 0) then
        write (1, '("(z<",f7.5,") show")') z (0)
     elseif (k.eq.nc) then
        write (1, '("(z>",f7.5,") show")') z (nc - 1)
     else
        write (1, '("(",f7.5,"<z<",f7.5,") show")') z (k - 1) , &
             z (k)
     endif
  enddo
  write (1, '(a)') '% End of Color Code'
  do k = 1, nc - 1
     if (z (k) .le.z (k - 1) ) stop 'zk order'
  enddo
  !   scan the array, top down, left to right, to paint shaded areas
  do j = jub - 1, 0, - 1
     do i = 0, iub - 1
        !  find lowest and highest vertex
        dmin = min (d (i, j), d (i, j + 1), d (i + 1, j), d (i + 1, j + 1) &
             )
        dmax = max (d (i, j), d (i, j + 1), d (i + 1, j), d (i + 1, j + 1) &
             )
        !  search for levels in this box
        nolevel = 0
        do k = 0, nc - 1
           if (z (k) < dmin) nolevel = k + 1
           if (z (k) >= dmin .and. z (k) <= dmax) then
              levelin = k
              goto 10
           endif
        enddo
        !  no level in this box: paint the whole box and pass to another box
        write (1, '(8f6.1,a4)') x (i) * cm, y (j) * cm, x (i + 1) * cm, y &
             (j) * cm, x (i + 1) * cm, y (j + 1) * cm, x (i) * cm, y (j + 1) &
             * cm, trapez (nolevel)
        goto 100
        !  there is at least a level in this box: paint the whole box
10      continue
        write (1, '(8f6.1,a4)') x (i) * cm, y (j) * cm, x (i + 1) * cm, y &
             (j) * cm, x (i + 1) * cm, y (j + 1) * cm, x (i) * cm, y (j + 1) &
             * cm, trapez (levelin)
        do k = levelin, nc - 1
           !  if no more levels in this box, move to another box
           if (z (k) .gt.dmax) goto 100
           !  find contour of zero levels in this box
           do m = 1, 4
              h (m) = d (i + im (m - 1), j + jm (m - 1) ) - z (k)
              xh (m) = x (i + im (m - 1) )
              yh (m) = y (j + jm (m - 1) )
           enddo
           h (0) = (h (1) + h (2) + h (3) + h (4) ) / 4
           xh (0) = (x (i) + x (i + 1) ) / 2
           yh (0) = (y (j) + y (j + 1) ) / 2
           do m = 0, 4
              if (h (m) .gt.0) then
                 ish (m) = 2
              elseif (h (m) .lt.0) then
                 ish (m) = 0
              else
                 ish (m) = 1
              endif
           enddo
           !   scan each triangle in the box to paint shaded areas
           do m = 1, 4
              m1 = m
              m2 = 0
              m3 = mod (m, 4) + 1
              npoint = 0
              icase = castab (ish (m1), ish (m2), ish (m3) )
              if (icase == 0) then
                 if (ish (m1) .eq.2) then
                    !     paint this triangle if positive
                    x1 = xh (m1)
                    y1 = yh (m1)
                    x2 = xh (m2)
                    y2 = yh (m2)
                    x3 = xh (m3)
                    y3 = yh (m3)
                    npoint = 3
                 endif
              elseif (icase == 1) then
                 !     line between vertices m1 and m2
                 x1 = xh (m1)
                 y1 = yh (m1)
                 x2 = xh (m2)
                 y2 = yh (m2)
                 if (ish (m3) .eq.2) then
                    x3 = xh (m3)
                    y3 = yh (m3)
                    npoint = 3
                 endif
              elseif (icase == 2) then
                 !     line between vertices m2 and m3
                 x1 = xh (m2)
                 y1 = yh (m2)
                 x2 = xh (m3)
                 y2 = yh (m3)
                 if (ish (m1)  == 2) then
                    x3 = xh (m1)
                    y3 = yh (m1)
                    npoint = 3
                 endif
              elseif (icase == 3) then
                 !     line between vertices m3 and m1
                 x1 = xh (m3)
                 y1 = yh (m3)
                 x2 = xh (m1)
                 y2 = yh (m1)
                 if (ish (m2)  == 2) then
                    x3 = xh (m2)
                    y3 = yh (m2)
                    npoint = 3
                 endif
              elseif (icase == 4) then
                 !     line between vertex m1 and side m2-m3
                 x1 = xh (m1)
                 y1 = yh (m1)
                 x2 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y2 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 if (ish (m3)  == 2) then
                    x3 = xh (m3)
                    y3 = yh (m3)
                 else
                    x3 = xh (m2)
                    y3 = yh (m2)
                 endif
                 npoint = 3
              elseif (icase == 5) then
                 !   line between vertex m2 and side m3-m1
                 x1 = xh (m2)
                 y1 = yh (m2)
                 x2 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y2 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 if (ish (m1)  == 2) then
                    x3 = xh (m1)
                    y3 = yh (m1)
                 else
                    x3 = xh (m3)
                    y3 = yh (m3)
                 endif
                 npoint = 3
              elseif (icase == 6) then
                 !     line between vertex m3 and line m1-m2
                 x1 = xh (m3)
                 y1 = yh (m3)
                 x2 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y2 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 if (ish (m2)  == 2) then
                    x3 = xh (m2)
                    y3 = yh (m2)
                 else
                    x3 = xh (m1)
                    y3 = yh (m1)
                 endif
                 npoint = 3
              elseif (icase == 7) then
                 !     line between sides m1-m2 and m2-m3
                 x1 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y1 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 x2 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y2 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 if (ish (m2)  == 2) then
                    x3 = xh (m2)
                    y3 = yh (m2)
                    npoint = 3
                 else
                    x3 = xh (m3)
                    y3 = yh (m3)
                    x4 = xh (m1)
                    y4 = yh (m1)
                    npoint = 4
                 endif
              elseif (icase == 8) then
                 !   line between sides m2-m3 and m3-m1
                 x1 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y1 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 x2 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y2 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 if (ish (m3)  == 2) then
                    x3 = xh (m3)
                    y3 = yh (m3)
                    npoint = 3
                 else
                    x3 = xh (m1)
                    y3 = yh (m1)
                    x4 = xh (m2)
                    y4 = yh (m2)
                    npoint = 4
                 endif
              elseif (icase == 9) then
                 !     line between sides m3-m1 and m1-m2
                 x1 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y1 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 x2 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y2 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 if (ish (m1)  == 2) then
                    x3 = xh (m1)
                    y3 = yh (m1)
                    npoint = 3
                 else
                    x3 = xh (m2)
                    y3 = yh (m2)
                    x4 = xh (m3)
                    y4 = yh (m3)
                    npoint = 4
                 endif
              endif
              if (npoint == 3) then
                 write (1, '(6f6.1,a4)') x1 * cm, y1 * cm, x2 * cm, y2 * cm, x3 &
                      * cm, y3 * cm, triangle (k + 1)
              elseif (npoint == 4) then
                 write (1, '(8f6.1,a4)') x1 * cm, y1 * cm, x2 * cm, y2 * cm, x3 &
                      * cm, y3 * cm, x4 * cm, y4 * cm, trapez (k + 1)
              endif
              if (icase.ne.0) write (1, '(4f6.1," p")') x1 * cm, y1 * cm, x2 * &
                   cm, y2 * cm
           enddo
        enddo
100     continue
     enddo
  enddo
  !
  return
end subroutine conrec

!
!-----------------------------------------------------------------------
subroutine atomi (nat, tau, ityp, at, a0, r0, tau1, tau2, xdim, ydim)
  !-----------------------------------------------------------------------
  !
  USE kinds, only: DP
  implicit none
  integer :: nat, ityp (nat)
  real(DP) :: tau (3, nat), at (3, 3), r0 (3), tau1 (3), tau2 (3), a0
  real(DP) :: xdim, ydim
  integer :: n1, n2, n3, i, n
  real(DP) :: r (3), ri (3), tau1n, tau2n, delta0, r1, r2, r3
  real(DP) :: delta = 1.0d0, cm = 28.453d0
  !
  delta0 = delta / a0
  tau1n = sqrt (tau1 (1) **2 + tau1 (2) **2 + tau1 (3) **2)
  tau2n = sqrt (tau2 (1) **2 + tau2 (2) **2 + tau2 (3) **2)
  !  clip the rectangle
  write (1, '(a)') 'gsave newpath'
  write (1, '(2f6.1," moveto")') 0.0d0, 0.0d0
  write (1, '(2f6.1," lineto")') xdim * cm, 0.0d0
  write (1, '(2f6.1," lineto")') xdim * cm, ydim * cm
  write (1, '(2f6.1," lineto")') 0.0d0, ydim * cm
  write (1, '(a)') 'closepath clip stroke'
  ! really bad algorithm to generate (hopefully) all atoms
  ! that are inside the contour plane or at |z| < delta (a.u.)
  do n1 = - 3, + 3
     do n2 = - 3, + 3
        do n3 = - 3, + 3
           do i = 1, 3
              r (i) = n1 * at (i, 1) + n2 * at (i, 2) + n3 * at (i, 3)
           enddo
           do n = 1, nat
              do i = 1, 3
                 ri (i) = tau (i, n) + r (i) - r0 (i)
              enddo
              !  this is the component along direction 1 ...
              r1 = (ri (1) * tau1 (1) + ri (2) * tau1 (2) + ri (3) * tau1 (3))&
                   / tau1n
              if (r1 > - delta0 .and. r1 < tau1n + delta0) then
                 !  and this is along direction 2 ...
                 r2 = (ri (1) * tau2(1) + ri (2) * tau2(2) + ri (3) * tau2(3))&
                      / tau2n
                 if (r2 > - delta0 .and. r2 < tau2n + delta0) then
                    do i = 1, 3
                       ri (i) = ri (i) - r1 * tau1 (i) / tau1n &
                                       - r2 * tau2 (i) / tau2n
                    enddo
                    r3 = sqrt (ri (1) **2 + ri (2) **2 + ri (3) **2)
                    !  and this is along the direction orthogonal to plane
                    if (abs (r3) < delta0) then
                       write (1, '(3f6.1," 0 360 arc gsave ",f4.2, &
                            &      " setgray fill grestore stroke")') &
                            r1 / tau1n * xdim * cm, r2 / tau2n * ydim * cm,&
                            delta0 / tau1n * xdim * cm, abs (r3) / delta0
                    endif
                 endif
              endif
           enddo
        enddo
     enddo
  enddo
  !
  write (1, '(a)') 'grestore'
  write (1, '(a)') 'grestore'
  write (1, '(a)') '%%Trailer'
  write (1, '(a)') 'showpage'

  close (1)
  return


end subroutine atomi
 
subroutine hatch (x1, x2, y1, y2)
  USE kinds, only: DP
  implicit none
  real(DP) :: x1, x2, y1, y2
  real(DP) :: cm =  28.453d0, delta = 0.2d0, dim
  integer :: nhach, n

  write (1, '(a)') '% Beginning of hatching'
  write (1, '(a)') 'gsave newpath'
  write (1, '(2f6.1," moveto")') x1 * cm, y1 * cm
  write (1, '(2f6.1," lineto")') x2 * cm, y1 * cm
  write (1, '(2f6.1," lineto")') x2 * cm, y2 * cm
  write (1, '(2f6.1," lineto")') x1 * cm, y2 * cm
  write (1, '(a)') 'closepath clip'

  dim = max (x2 - x1, y2 - y1)
  nhach = dim / delta
  !      delta=dim/nhach
  do n = 1, nhach
     write (1, '(4f6.1," p")') (x1 + (n - 1) * delta) * cm, y1 * cm, &
          (x1 + dim) * cm, (y1 + dim - (n - 1) * delta) * cm
     write (1, '(4f6.1," p")') x1 * cm, (y1 + (n - 1) * delta) &
          * cm, (x1 + dim - (n - 1) * delta) * cm, (y1 + dim) * cm
  enddo
  write (1, '(a)') 'grestore'
  write (1, '(a)') '% End of hatching'

  return
end subroutine hatch
!
!-----------------------------------------------------------------------
subroutine psplot ( d, imax, x, iub, y, jub, nlevels, z, &
     xdim, ydim, xs, ys, filename)
  !-----------------------------------------------------------------------
  !
  !   draws a contour plot of d(i,j). PostScript output on unit 1
  !
  !   d(0:imax,0:jub) contains the function to plot on a uniform 2d grid
  !   x(0:iub) workspace
  !   iub  is the number of intervals along the x axis (<=imax)
  !   y(0:jub) workspace
  !   jub  is the number of intervals along the y axis
  !   nlevels is the number of levels wanted
  !   z(0:nlevels) are the levels
  !   xdim and ydim are the physical dimensions of the figure in cm
  !   xs and ys determine a shift of the origin in cm
  USE kinds, only : DP
  implicit none
  integer, parameter :: ncontourmax=100, npointmax=500, nmaxtypes=8
  integer :: imax, iub, jub, nlevels
  real(DP) :: d(0:imax,0:jub), x(0:imax), y(0:jub), z(0:nlevels), &
       xdim, ydim, xs, ys
  !
  real(DP) :: line(2,npointmax,ncontourmax), work(2,npointmax), &
       segments(2,2,npointmax*ncontourmax), yy, dy
  integer :: npoints(ncontourmax), ncontours, i, j, k, l, npts
  real(DP), parameter :: cm = 28.453d0
  logical ::  lwork(npointmax*ncontourmax)
  character (len=256) :: filename, linetype(nmaxtypes)
  data linetype &
       / '[1 2    ] 0 setdash 0.5 setlinewidth', &
       '[1 2 4 2] 0 setdash', &
       '[4 2    ] 0 setdash', &
       '[       ] 0 setdash', &
       '[1 2   ] 0 setdash 1.0 setlinewidth', &
       '[1 2 4 2] 0 setdash', &
       '[4 2    ] 0 setdash', &
       '[       ] 0 setdash'/

  open (unit = 1, file = filename, status = 'unknown', form = &
       'formatted')
  if(nlevels < 1)  stop '  nunmber of levels is wrong'
  if(iub > imax .or. iub < 1) stop ' iub too big or wrong'
  if(jub < 1) stop ' jub or wrong'
  if(xdim < 3.0d0 .or. ydim < 3.0d0 ) stop ' really too small!'
  if(xdim > 20.0.or. ydim > 30.0) stop ' really too big!'
  if(abs(xs) > 20 .or. abs(ys) > 30) stop ' xs or ys are weird'

  !  initializations for PostScript output
  write (1,'(a)') '%! PS-Adobe-1.0'
  write (1,'(''%%BoundingBox:'',4f6.1)') &
       xs*cm, ys*cm,(xs+xdim)*cm,(ys+ydim)*cm
  write (1,'(a)') '/localdict 100 dict def'
  write (1,'(a)') 'localdict begin'
  write (1,'(a)') '/cm {28.453 mul} def'
  write (1,'(a)') '/title {('//filename//')} def'
  write (1,'(a)') '/Times-Roman findfont 12 scalefont setfont'
  write (1,'(a)') '% cshow prints a centered string at current position'
  write (1,'(a)') '/cshow {gsave dup stringwidth pop 2 div neg 0'
  write (1,'(a)') ' rmoveto show grestore} def'
  write (1,'(a)') '% x1 y1 nm x2 y2 ...xn yn n-1 lrs :'
  write (1,'(a)') '% draws a line from point 1 to point n'
  write (1,'(a)') '/nm  {newpath moveto} def'
  write (1,'(a)') '/lrs {{lineto} repeat stroke} def'
  write (1,'(a)') '%%EndPreamble'
  write (1,'(a)') 'gsave'
  write (1,'(a)') '0 setgray'
  write (1,'(1x,f6.2,'' cm '',f6.2,'' cm translate'')') xs, ys
  write (1,'(1x,f6.2,'' cm '',f6.2,'' cm moveto title cshow'')') &
       xdim/2, ydim+1.5d0
  !
  do i=0,iub
     x(i)=xdim*DBLE(i)/iub
  end do
  do j=0,jub
     y(j)=ydim*DBLE(j)/jub
  end do

  do k=1,nlevels-1
     if (z(k).le.z(k-1)) stop ' order of levels'
  end do

  yy=y(jub)*0.75d0
  dy=(y(jub)-y(0))/(nlevels-1)/2

  do k=0,nlevels-1
     ! type of line
     write (1,'(''%% This is level '',i3)') k
     write (1,'(50a)') linetype(mod(k,nmaxtypes)+1)
     write (1,*) ' newpath '
     write (1,'(2f6.1,'' moveto'')') (x(iub)+1.0d0)*cm, yy*cm
     write (1,'(2f6.1,'' lineto'')') (x(iub)+1.8d0)*cm, yy*cm
     write (1,*) 'closepath stroke'
     write (1,'(2f6.1,'' moveto'')') (x(iub)+2.0d0)*cm, yy*cm
     write (1,'(''( z = '',f6.3,'') show'')') z(k)
     yy=yy-dy

     call findcontours ( d, imax, iub, jub, x, y, z(k), segments, work, &
          lwork, npointmax, ncontourmax, line, npoints, ncontours)

     do l=1,ncontours
        write (1,'(12f6.1)') &
             ( line(1,npts,l)*cm, line(2,npts,l)*cm, npts=npoints(l),1,-1 )
        write (1,'(''nm'',i4,'' lrs'')') npoints(l)-1
     end do
     
  end do

  !  draw frame of size xdim by ydim
  write (1,*) '[] 0 setdash 1 setlinewidth 0 setgray newpath'
  write (1,'(2f6.1,'' moveto'')')     0.0d0,0.0d0
  write (1,'(2f6.1,'' lineto'')') xdim*cm,0.0d0
  write (1,'(2f6.1,'' lineto'')') xdim*cm,ydim*cm
  write (1,'(2f6.1,'' lineto'')')     0.0d0,ydim*cm
  write (1,*) 'closepath stroke'
  !      write (1,*) 'grestore'
  !      write (1,'(a)') '%%Trailer'
  !      write (1,*) 'showpage'
  return
end subroutine psplot

!
subroutine findcontours (d, imax, iub, jub, x, y, z, segments, work, &
     done, npointmax, ncontourmax, line, npoints, ncontours)
  !
  USE kinds, only : DP
  implicit none
  integer :: imax, iub, jub, npointmax, ncontourmax,  ncontours, &
       npoints(ncontourmax)
  real(DP) :: d(0:imax,0:jub),x(0:iub), y(0:jub), z, &
       line(2,npointmax,ncontourmax), segments(2,2,npointmax*ncontourmax)
 logical :: done(npointmax*ncontourmax)
  !
  integer :: i, j, m, i0, j0, m0, nsegments, nsegmax, nseg, nnext, npts, &
       isame, inext, nleft
  real(DP) :: epsx, epsy, work(2,npointmax)
  logical :: found, start_from_boundary, hit_boundary, loop_closed

  nsegmax=npointmax*ncontourmax
  nsegments=0
  do j0=jub-1,0,-1
     do i0=0,iub-1
        do m0=1,4
           !     scan boxes and triangles until a contour is found
           call triangle ( d, x, y, z, imax, iub, jub, &
                i0, j0, m0, found, segments(1,1,nsegments+1) )
           if (found) then
              nsegments=nsegments+1
              if (nsegments > nsegmax) then
                 print '('' Erore: nsegments, nsegmax '',2i5)', &
                      nsegments, nsegmax
                 stop
              end if
           end if
        end do
     end do
  end do
  
  do nseg=1,nsegments
     done(nseg)=.false.
  end do

  ncontours=0
  nleft=nsegments
  epsx = 0.0001d0 * (x(1)-x(0))
  epsy = 0.0001d0 * (y(1)-y(0))
  do nseg=1,nsegments
     if (.not.done(nseg)) then
        ncontours=ncontours+1
        if (ncontours.gt.ncontourmax) then
           print '('' too many ('',i4,'') contours'')', ncontours
           stop
        end if
        line(1,1,ncontours)=segments(1,1,nseg)
        line(2,1,ncontours)=segments(2,1,nseg)
        start_from_boundary = &
             abs(segments(1,1,nseg)-x(  0)) < epsx .or. &
             abs(segments(1,1,nseg)-x(iub)) < epsx .or. &
             abs(segments(2,1,nseg)-y(  0)) < epsy .or. &
             abs(segments(2,1,nseg)-y(jub)) < epsy
        ! true if the starting point lies on a boundary
        nnext=nseg
        inext=2
        npts=1
10      npts=npts+1
        if (npts > npointmax) then
           print '('' too many ('',i5,'') points'')', npts
           stop
        end if
        line(1,npts,ncontours)=segments(1,inext,nnext)
        line(2,npts,ncontours)=segments(2,inext,nnext)
        nleft=nleft-1
        done(nnext)=.true.
        ! check if the loop is closed
        loop_closed = ( abs(segments(1,inext,nnext) - &
                        line(1,1,ncontours)) < epsx ) .and. &
                      ( abs(segments(2,inext,nnext) - &
                        line(2,1,ncontours)) < epsy )
        ! check if the border has been reached
        hit_boundary = &
             abs(segments(1,inext,nnext)-x(  0)) < epsx .or. &
             abs(segments(1,inext,nnext)-x(iub)) < epsx .or. &
             abs(segments(2,inext,nnext)-y(  0)) < epsy .or. &
             abs(segments(2,inext,nnext)-y(jub)) < epsy
        if (nleft == 0) then
           if (.not.hit_boundary .and. .not.loop_closed) &
                print '('' Warning: no more points'')'
           npoints(ncontours)=npts
           return
        end if
        if (start_from_boundary) then
           if (hit_boundary) then
              npoints(ncontours)=npts
              go to 20
           else if (loop_closed) then
              print '('' Warning: loop closed on the boundary ?'')'
              go to 20
           end if
        else
           if (loop_closed) then
              npoints(ncontours)=npts
              go to 20
           else if (hit_boundary) then
              do nnext=1,npts
                 work(1,nnext)=line(1,npts-nnext+1,ncontours)
                 work(2,nnext)=line(2,npts-nnext+1,ncontours)
              end do
              do nnext=1,npts
                 line(1,nnext,ncontours)=work(1,nnext)
                 line(2,nnext,ncontours)=work(2,nnext)
              end do
              start_from_boundary=.true.
           end if
        end if

        do nnext=1,nsegments
           if (.not.done(nnext)) then
              do isame=1,2
                 found= ( abs(segments(1,isame,nnext) - &
                          line(1,npts,ncontours)) < epsx ) .and. &
                        ( abs(segments(2,isame,nnext) - &
                          line(2,npts,ncontours)) .lt. epsy )
                 if (found) then
                    inext=mod(isame,2)+1
                    !inext =1 if isame=2 and vice versa
                    go to 10
                 endif
              end do
           end if
        end do
        npoints(ncontours)=npts-1
        print '('' Warning: next point not found'',i5)', nseg
     end if
20   continue
  end do

  return
end subroutine findcontours
!
subroutine triangle ( d, x, y, z, imax, iub, jub, i0, j0, m0, found, line )
  !
  USE kinds, only : DP
  implicit none
  integer :: imax, iub, jub, i0, j0, m0
  real(DP) :: d(0:imax,0:jub), x(0:iub), y(0:jub), z, line(2,2)
  logical :: found
  !
  real(DP) :: h(0:4), xh(0:4), yh(0:4), dmin, dmax
  integer :: m, icase, m1, m2, m3, i, j, k, ish(0:4), castab(-1:1,-1:1,-1:1)
  data (((castab(i,j,k),k=-1,1),j=-1,1),i=-1,1) / &
       0 , 0 , 8 , 0 , 2 , 5 , 7 , 6 , 9 , 0 , 3 , 4 , &
       1 , 3 , 1 , 4 , 3 , 0 , 9 , 6 , 7 , 5 , 2 , 0 , 8 , 0 , 0 /

  !     values of the function in the four corners...
  h (1) = d(i0  ,j0  ) - z
  h (2) = d(i0+1,j0  ) - z
  h (3) = d(i0+1,j0+1) - z
  h (4) = d(i0  ,j0+1) - z
  dmin=min(h(1),h(2),h(3),h(4))
  dmax=max(h(1),h(2),h(3),h(4))
  !     return if no contour is found  in this box
  if ( dmin > 0.0d0 .or. dmax < 0.0d0 ) then
     found=.false.
     return
  end if
  !     xy coordinates of the four corners...
  xh(1) = x(i0  )
  xh(2) = x(i0+1)
  xh(3) = x(i0+1)
  xh(4) = x(i0  )
  yh(1) = y(j0  )
  yh(2) = y(j0  )
  yh(3) = y(j0+1)
  yh(4) = y(j0+1)
  !     and this is the center of the box
  h(0) =(h(1)+h(2)+h(3)+h(4))/4.0d0
  xh(0)=(x(i0)+x(i0+1))/2.0d0
  yh(0)=(y(j0)+y(j0+1))/2.0d0
  !     ish=+1 if h>0 ; ish=0 if h=0 ; ish=-1 if h<0 
  do m=0,4
     ish(m) = nint( sign(1.d0, h(m)) )
     if (h(m) == 0.0d0) ish(m) = 0
  end do
  !     starting triangle
  m1=m0
  m2=0
  m3=mod(m0,4)+1
  icase=castab(ish(m1),ish(m2),ish(m3))
  if (icase == 0) then
     !     no contour in this triangle
     found=.false.
     return
     ! NOTA BENE: if there is a point on the boundary,
     ! this is always chosen to be the first point
  else if (icase == 1) then
     !     line between vertices m1 and m2
     line(1,1)=xh(m1)
     line(2,1)=yh(m1)
     line(1,2)=xh(m2)
     line(2,2)=yh(m2)
  else if (icase == 2) then
     !     line between vertices m2 and m3
     line(1,1)=xh(m3)
     line(2,1)=yh(m3)
     line(1,2)=xh(m2)
     line(2,2)=yh(m2)
  else if (icase == 3) then
     !     line between vertices m3 and m1
     line(1,1)=xh(m3)
     line(2,1)=yh(m3)
     line(1,2)=xh(m1)
     line(2,2)=yh(m1)
  else if (icase == 4) then
     !     line between vertex m1 and side m2-m3
     line(1,1)=xh(m1)
     line(2,1)=yh(m1)
     line(1,2)=(h(m3)*xh(m2)-h(m2)*xh(m3))/(h(m3)-h(m2))
     line(2,2)=(h(m3)*yh(m2)-h(m2)*yh(m3))/(h(m3)-h(m2))
  else if (icase == 5) then
     !     line between vertex m2 and side m3-m1
     line(1,1)=(h(m1)*xh(m3)-h(m3)*xh(m1))/(h(m1)-h(m3))
     line(2,1)=(h(m1)*yh(m3)-h(m3)*yh(m1))/(h(m1)-h(m3))
     line(1,2)=xh(m2)
     line(2,2)=yh(m2)
  else if (icase == 6) then
     !     line between vertex m3 and line m1-m2
     line(1,1)=xh(m3)
     line(2,1)=yh(m3)
     line(1,2)=(h(m2)*xh(m1)-h(m1)*xh(m2))/(h(m2)-h(m1))
     line(2,2)=(h(m2)*yh(m1)-h(m1)*yh(m2))/(h(m2)-h(m1))
  else if (icase == 7) then
     !     line between sides m1-m2 and m2-m3
     line(1,1)=(h(m2)*xh(m1)-h(m1)*xh(m2))/(h(m2)-h(m1))
     line(2,1)=(h(m2)*yh(m1)-h(m1)*yh(m2))/(h(m2)-h(m1))
     line(1,2)=(h(m3)*xh(m2)-h(m2)*xh(m3))/(h(m3)-h(m2))
     line(2,2)=(h(m3)*yh(m2)-h(m2)*yh(m3))/(h(m3)-h(m2))
  else if (icase == 8) then
     !     line between sides m2-m3 and m3-m1
     line(1,1)=(h(m1)*xh(m3)-h(m3)*xh(m1))/(h(m1)-h(m3))
     line(2,1)=(h(m1)*yh(m3)-h(m3)*yh(m1))/(h(m1)-h(m3))
     line(1,2)=(h(m3)*xh(m2)-h(m2)*xh(m3))/(h(m3)-h(m2))
     line(2,2)=(h(m3)*yh(m2)-h(m2)*yh(m3))/(h(m3)-h(m2))
  else if (icase == 9) then
     !     line between sides m3-m1 and m1-m2
     line(1,1)=(h(m1)*xh(m3)-h(m3)*xh(m1))/(h(m1)-h(m3))
     line(2,1)=(h(m1)*yh(m3)-h(m3)*yh(m1))/(h(m1)-h(m3))
     line(1,2)=(h(m2)*xh(m1)-h(m1)*xh(m2))/(h(m2)-h(m1))
     line(2,2)=(h(m2)*yh(m1)-h(m1)*yh(m2))/(h(m2)-h(m1))
  endif

  found=.true.

  return
end subroutine triangle
