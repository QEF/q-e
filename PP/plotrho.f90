!
! Copyright (C) 2001 PWSCF group
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
  !                   - negative values are shaded
  !                   - contours plus gray levels
  !                   - Postscript printable output
  !
#include "machine.h"
  use parameters, only: DP
  implicit none
  integer :: nwrk, nximax, nyimax, nxmax, nymax, nlevelx, nax
  parameter (nwrk = 10000, nximax = 64, nyimax = 64, nxmax = 128, &
       nymax = 128, nlevelx = 19, nax = 130)
  integer :: ityp (nax), nxi, nyi, nx, ny, i, j, k, nlevels, na, &
       nat, iargc
  real(kind=DP) :: rhoi (0:nximax, 0:nyimax), xi (0:nximax), yi (0: &
       nyimax), rhoo (0:nxmax, 0:nymax), x (0:nxmax), y (0:nymax), &
       z (0:nlevelx), wrk (nwrk), xmin, xmax, ymin, ymax, rhoomin, &
       rhoomax
  real(kind=DP) :: xdim, ydim, xs, ys
  real(kind=DP) :: r0 (3), tau1 (3), tau2 (3), tau (3, nax)
  real(kind=DP) :: at (3, 3), a0
  character (len=20) :: filename, fileout, ans * 1


  logical :: logarithmic_scale
  i = iargc ()
  if (i.eq.0) then
     print '("  input file > ",$)'
     read (5, '(a)', end = 20, err = 20) filename
  elseif (i.eq.1) then
#ifdef T3D
     call pxfgetarg (1, filename, i, j)
#else
     call getarg (1, filename)
#endif
  else
     print '("   usage: plotrho  [input file] ")'
  endif
  open (unit = 1, file = filename, form = 'formatted', status = &
       'old', iostat = i)
  if (i.ne.0) then
     print '("   file not found ")'
     stop

  endif
  read (1, * ) nxi, nyi
  if (nxi.gt.nximax.or.nyi.gt.nyimax) then
     print '(" nx or ny too big ")'
     stop

  endif
  read (1, * ) (xi (i), i = 0, nxi)
  read (1, * ) (yi (j), j = 0, nyi)
  read (1, * ) ( (rhoi (i, j), i = 0, nxi), j = 0, nyi)
  read (1, * ) r0
  read (1, * ) tau1
  read (1, * ) tau2
  read (1, * ) nat
  if (nat.gt.nax) then
     print '(" too many atoms ")'
     stop
  endif
  read (1, * ) ( (tau (j, na), j = 1, 3), ityp (na), na = 1, nat)
  read (1, * ) a0
  read (1, * ) at

  close (unit = 1)
  write (6, '(" r0   : ",3f8.4)') r0
  write (6, '(" tau1 : ",3f8.4)') tau1
  write (6, '(" tau2 : ",3f8.4)') tau2

  write (6, '(i4," atomic positions read")') nat
  !      write(6,'(" Atomic positions:")')
  !      write(6,'(3f8.4)') ( (tau(j,na),j=1,3),na=1,nat)
  write (6, '("   output file > ",$)')

  read (5, '(a)') fileout
  write (6, '("   nx, ny (output) > ",$)')
  read (5, * ) nx, ny
  if (nx.gt.nxmax.or.ny.gt.nymax) then
     write (6, * ) 'Error: nx or ny too big '
     stop

  endif
  xmin = xi (0)
  xmax = xi (nxi)
  do i = 0, nx
     x (i) = (xi (nxi) - xi (0) ) * float (i) / float (nx)
  enddo

  ymin = yi (0)
  ymax = yi (nyi)
  do i = 0, ny
     y (i) = (yi (nyi) - yi (0) ) * float (i) / float (ny)
  enddo
#ifdef AIX
  call dcsin2 (xi, yi, rhoi, nxi + 1, nyi + 1, nximax + 1, x, y, nx &
       + 1, ny + 1, rhoo, nxmax + 1, wrk, nwrk)
#else
  if (nx.ne.nxi.or.ny.ne.nyi) then
     print '("Sorry, interpolation not (yet) supported")'
     nx = nxi
     ny = nyi
  else
     do j = 0, ny
        do i = 0, nx
           rhoo (i, j) = rhoi (i, j)
        enddo
     enddo
  endif
#endif
  rhoomin = 1.0e+20
  rhoomax = - 1.0e+20
  do j = 0, ny
     do i = 0, nx
        rhoomin = min (rhoomin, rhoo (i, j) )
        rhoomax = max (rhoomax, rhoo (i, j) )
     enddo
  enddo

  print '("   Logarithmic scale (y/n)? > ",$)'
  read (5, '(a)') ans
  logarithmic_scale = ans.ne.'n'.and.ans.ne.'N'
  print '("   Bounds: ",2f12.6)', rhoomin, rhoomax
  print '("   min, max, # of levels > ",$)'

  read (5, * ) rhoomin, rhoomax, nlevels
  if (nlevels.gt.nlevelx) stop ' too many levels '
  if (logarithmic_scale) then
     do k = 0, nlevels - 1
        z (k) = exp (log (rhoomin) + (log (rhoomax) - log (rhoomin) ) &
             * float (k) / (nlevels - 1) )
     enddo
  else
     do k = 0, nlevels - 1
        z (k) = rhoomin + (rhoomax - rhoomin) * float (k) / (nlevels - 1)
     enddo
  endif

  z (nlevels) = z (nlevels - 1)
  xdim = 15.0 * (xmax - xmin) / sqrt ( (xmax - xmin) **2 + (ymax - ymin) **2)
  ydim = 15.0 * (ymax - ymin) / sqrt ( (xmax - xmin) **2 + (ymax - ymin) **2)

  xs = 4.0
  ys = 3.0

  call cplot (rhoo, nxmax, nymax, x, xmin, xmax, nx, y, ymin, ymax, &
       ny, nlevels, z, xdim, ydim, xs, ys, filename, fileout)

  call atomi (nat, tau, ityp, at, a0, r0, tau1, tau2, xdim, ydim)
20 stop
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
  use parameters, only: DP
  implicit none
  integer :: imax, jmax, iub, jub, nc
  real(kind=DP) :: d (0:imax, 0:jmax), x (0:imax), y (0:jmax), z (0:nc)
  real(kind=DP) :: xmin, xmax, ymin, ymax, xdim, ydim, xs, ys

  character (len=20) :: filename, str
  integer :: ncmax
  parameter (ncmax = 19)
  integer :: i, j, k
  real(kind=DP) :: gray (0:ncmax), dim, cm, width, gray0, deltagray

  data cm / 28.453 /, width / 0.5 /, gray0 / 1.0 /, deltagray / 0.7 /
  !  cm   : number of points per cm
  !  width: linewidth of the contour plot for PostScript printer

  open (unit = 1, file = filename, status = 'unknown', form = &
       'formatted')
  if (nc.gt.ncmax.or.nc.lt.1) stop '  nc too big or wrong'
  if (iub.gt.imax.or.iub.lt.1) stop ' iub too big or wrong'
  if (jub.gt.jmax.or.jub.lt.1) stop ' jub too big or wrong'
  if (xdim.lt.3.0.or.ydim.lt.3.0) stop ' really too small!'
  if (xdim.gt.20.0.or.ydim.gt.30.0) stop ' really too big!'
  if (abs (xs) .gt.20.or.abs (ys) .gt.30) stop ' xs or ys are weird'
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
       xdim / 2, ydim + 1.5

  call hatch (0.d0, xdim, 0.d0, ydim)
  do i = 0, iub
     x (i) = xdim * float (i) / iub
  enddo
  do j = 0, jub
     y (j) = ydim * float (j) / jub

  enddo

  call conrec (imax, iub, jmax, jub, x, y, d, nc, z)
  !  draw frame of size xdim by ydim
  write (1, '(a)') '1 setlinewidth 0 setgray newpath'
  write (1, '(2f6.1," moveto")') 0.0, 0.0
  write (1, '(2f6.1," lineto")') xdim * cm, 0.0
  write (1, '(2f6.1," lineto")') xdim * cm, ydim * cm
  write (1, '(2f6.1," lineto")') 0.0, ydim * cm

  write (1, '(a)') 'closepath stroke'
  !      write (1,'(a)') 'grestore'
  !      write (1,'(a)') '%%Trailer'
  !      write (1,'(a)') 'showpage'
  !      close(1)
  return


end subroutine cplot

subroutine conrec (imax, iub, jmax, jub, x, y, d, nc, z)
  use parameters, only: DP
  implicit none
  integer :: imax, iub, jmax, jub, nc

  real(kind=DP) :: d (0:imax, 0:jmax), x (0:imax), y (0:jmax), z (0:nc)
  integer :: ncmax
  parameter (ncmax = 19)
  character (len=4) :: triangle (0:ncmax), trapez (0:ncmax)
  real(kind=DP) :: h (0:4), xh (0:4), yh (0:4)
  real(kind=DP) :: x1, y1, x2, y2, x3, y3, x4, y4, dx, dy, xx, yy, cm, &
       dmin, dmax
  !   cm : conversion factor from cm to points for PostScript
  integer :: ish (0:4), im (0:3), jm (0:3), castab (0:2, 0:2, 0:2)
  integer :: i, j, k, m, m1, m2, m3, npoint, icase, levelin, &
       nolevel
  data cm / 28.453 /
  data (im (i), i = 0, 3) / 0, 1, 1, 0 /
  data (jm (i), i = 0, 3) / 0, 0, 1, 1 /


  data ( ( (castab (i, j, k), k = 0, 2), j = 0, 2), i = 0, 2) &
       / 0, 0, 8, 0, 2, 5, 7, 6, 9, 0, 3, 4, 1, 3, 1, 4, 3, 0, 9, 6, 7, &
       5, 2, 0, 8, 0, 0 /
  dy = (y (jub) - y (0) ) / (nc + 1)
  xx = x (iub) + 1.0
  dx = 0.5

  write (1, '(a)') '% Start of Color Code'

  call hatch (xx, xx + dx, y (0), y (jub) )
  do k = 0, nc
     yy = y (jub) - k * dy
     write (1, '(8f6.1,$)') xx * cm, yy * cm, (xx + dx) * cm, yy * cm, &
          (xx + dx) * cm, (yy - dy) * cm, xx * cm, (yy - dy) * cm
     if (k.lt.10) then
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
     write (1, '(2f6.1," moveto")') (x (iub) + 2.0) * cm, (yy - dy / &
          2) * cm
     if (k.eq.0) then
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
           if (z (k) .lt.dmin) nolevel = k + 1
           if (z (k) .ge.dmin.and.z (k) .le.dmax) then
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
           !  if no more levels in this box, pas to another box
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
              if (icase.eq.0) then
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
              elseif (icase.eq.1) then
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
              elseif (icase.eq.2) then
                 !     line between vertices m2 and m3
                 x1 = xh (m2)
                 y1 = yh (m2)
                 x2 = xh (m3)
                 y2 = yh (m3)
                 if (ish (m1) .eq.2) then
                    x3 = xh (m1)
                    y3 = yh (m1)
                    npoint = 3
                 endif
              elseif (icase.eq.3) then
                 !     line between vertices m3 and m1
                 x1 = xh (m3)
                 y1 = yh (m3)
                 x2 = xh (m1)
                 y2 = yh (m1)
                 if (ish (m2) .eq.2) then
                    x3 = xh (m2)
                    y3 = yh (m2)
                    npoint = 3
                 endif
              elseif (icase.eq.4) then
                 !     line between vertex m1 and side m2-m3
                 x1 = xh (m1)
                 y1 = yh (m1)
                 x2 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y2 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 if (ish (m3) .eq.2) then
                    x3 = xh (m3)
                    y3 = yh (m3)
                 else
                    x3 = xh (m2)
                    y3 = yh (m2)
                 endif
                 npoint = 3
              elseif (icase.eq.5) then
                 !   line between vertex m2 and side m3-m1
                 x1 = xh (m2)
                 y1 = yh (m2)
                 x2 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y2 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 if (ish (m1) .eq.2) then
                    x3 = xh (m1)
                    y3 = yh (m1)
                 else
                    x3 = xh (m3)
                    y3 = yh (m3)
                 endif
                 npoint = 3
              elseif (icase.eq.6) then
                 !     line between vertex m3 and line m1-m2
                 x1 = xh (m3)
                 y1 = yh (m3)
                 x2 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y2 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 if (ish (m2) .eq.2) then
                    x3 = xh (m2)
                    y3 = yh (m2)
                 else
                    x3 = xh (m1)
                    y3 = yh (m1)
                 endif
                 npoint = 3
              elseif (icase.eq.7) then
                 !     line between sides m1-m2 and m2-m3
                 x1 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y1 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 x2 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y2 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 if (ish (m2) .eq.2) then
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
              elseif (icase.eq.8) then
                 !   line between sides m2-m3 and m3-m1
                 x1 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y1 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 x2 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y2 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 if (ish (m3) .eq.2) then
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
              elseif (icase.eq.9) then
                 !     line between sides m3-m1 and m1-m2
                 x1 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y1 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 x2 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y2 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 if (ish (m1) .eq.2) then
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
              if (npoint.eq.3) then
                 write (1, '(6f6.1,a4)') x1 * cm, y1 * cm, x2 * cm, y2 * cm, x3 &
                      * cm, y3 * cm, triangle (k + 1)
              elseif (npoint.eq.4) then
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
  use parameters, only: DP
  implicit none
  integer :: nat, ityp (nat)
  real(kind=DP) :: tau (3, nat), at (3, 3), r0 (3), tau1 (3), tau2 (3), a0
  real(kind=DP) :: xdim, ydim
  integer :: n1, n2, n3, i, n
  real(kind=DP) :: r (3), ri (3), tau1n, tau2n, delta, delta0, cm, r1, r2, r3

  parameter (delta = 1.0, cm = 28.453)
  !  Soluzione piu' che provvisoria con algoritmo mongolissimo:
  !  genera tutti (si spera) gli atomi che stanno dentro al rettangolo
  delta0 = delta / a0
  tau1n = sqrt (tau1 (1) **2 + tau1 (2) **2 + tau1 (3) **2)
  tau2n = sqrt (tau2 (1) **2 + tau2 (2) **2 + tau2 (3) **2)
  !  clippa il rettangolo
  write (1, '(a)') 'gsave newpath'
  write (1, '(2f6.1," moveto")') 0.0, 0.0
  write (1, '(2f6.1," lineto")') xdim * cm, 0.0
  write (1, '(2f6.1," lineto")') xdim * cm, ydim * cm
  write (1, '(2f6.1," lineto")') 0.0, ydim * cm

  write (1, '(a)') 'closepath clip stroke'
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
              !  questa e' la componente lungo la direzione 1 ...
              r1 = (ri (1) * tau1 (1) + ri (2) * tau1 (2) + ri (3) * tau1 (3) ) &
                   / tau1n
              if (r1.gt. - delta0.and.r1.lt.tau1n + delta0) then
                 !  e questa lungo la direzione 2 ...
                 r2 = (ri (1) * tau2 (1) + ri (2) * tau2 (2) + ri (3) * tau2 (3) &
                      ) / tau2n
                 if (r2.gt. - delta0.and.r2.lt.tau2n + delta0) then
                    do i = 1, 3
                       ri (i) = ri (i) - r1 * tau1 (i) / tau1n - r2 * tau2 (i) &
                            / tau2n
                    enddo
                    r3 = sqrt (ri (1) **2 + ri (2) **2 + ri (3) **2)
                    !  e questa lungo la direzione ortogonale al piano
                    if (abs (r3) .lt.delta0) then
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
  write (1, '(a)') 'grestore'
  !
  write (1, '(a)') 'grestore'
  write (1, '(a)') '%%Trailer'
  write (1, '(a)') 'showpage'

  close (1)
  return


end subroutine atomi

subroutine hatch (x1, x2, y1, y2)
  use parameters, only: DP
  implicit none
  real(kind=DP) :: x1, x2, y1, y2
  real(kind=DP) :: cm, delta, dim
  integer :: nhach, n


  data cm / 28.453 /, delta / 0.2 /
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
