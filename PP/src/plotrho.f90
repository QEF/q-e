!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM plotrho
  !-----------------------------------------------------------------------
  !   2D contour plot - logarithmically or linearly spaced levels
  !                   - Postscript printable output
  !   if " cplot" is called:
  !                   - contour lines plus gray levels
  !                   - negative values are shaded
  !   if "psplot" is called:
  !                   - contour lines of various kinds (solid, dashed, etc)
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: stdout=6
  ! for spline interpolation using essl toutines
  INTEGER, PARAMETER :: nwrk = 10000
  real(DP) :: wrk (nwrk)
  INTEGER, ALLOCATABLE :: ityp (:)
  INTEGER :: nxi, nyi, nx, ny, i, j, k, nlevels, na, nat
  real(DP), ALLOCATABLE :: rhoi(:,:), xi(:), yi(:)
  real(DP), ALLOCATABLE :: rhoo(:,:), x (:), y (:)
  real(DP), ALLOCATABLE ::  z (:)
  real(DP) :: xmin, xmax, ymin, ymax, rhomin, rhomax, rhoomin, rhoomax
  real(DP) :: xdim, ydim, xs, ys, r0 (3), tau1 (3), tau2 (3)
  real(DP), ALLOCATABLE :: tau (:,:)
  real(DP) :: at (3, 3), a0
  CHARACTER (len=256) :: filename, fileout, ans * 1
  LOGICAL :: logarithmic_scale


  CALL get_file ( filename )
  OPEN (unit = 1, file = filename, form = 'formatted', status = 'old')

  READ (1, * ) nxi, nyi
  ALLOCATE ( xi(0:nxi), yi(0:nyi), rhoi(0:nxi,0:nyi) )
  READ (1, * ) (xi (i), i = 0, nxi)
  READ (1, * ) (yi (j), j = 0, nyi)
  READ (1, * ) ( (rhoi (i, j), i = 0, nxi), j = 0, nyi)
  READ (1, * ) r0
  READ (1, * ) tau1
  READ (1, * ) tau2
  READ (1, * ) nat
  IF (nat < 0 .or. nat > 1000000) THEN
     WRITE( stdout, '("Error: unlikely number of atoms ",i4)') nat
     STOP
  ENDIF
  ALLOCATE (tau (3,nat), ityp(nat) )
  READ (1, * ) ( (tau (j, na), j = 1, 3), ityp (na), na = 1, nat)
  READ (1, * ) a0
  READ (1, * ) at
  CLOSE (unit = 1)
  !
  WRITE( stdout, '("r0   : ",3f8.4)') r0
  WRITE( stdout, '("tau1 : ",3f8.4)') tau1
  WRITE( stdout, '("tau2 : ",3f8.4)') tau2
  !
  WRITE( stdout, '("read",i4," atomic positions")') nat
  !      WRITE( stdout,'("Atomic positions:")')
  !      WRITE( stdout,'(3f8.4)') ( (tau(j,na),j=1,3),na=1,nat)
  WRITE( stdout, '("output file > ")', advance="NO")

  READ (5, '(a)') fileout

  WRITE( stdout, '("Read ",i3," *",i3,"  grid")') nxi+1, nyi+1
  ! interpolation no longer implemented
  nx = nxi
  ny = nyi
  ALLOCATE ( x(0:nx), y(0:ny), rhoo(0:nx,0:ny) )
  xmin = xi (0)
  xmax = xi (nxi)
  DO i = 0, nx
     x (i) = (xi (nxi) - xi (0) ) * dble (i) / dble (nx)
  ENDDO

  ymin = yi (0)
  ymax = yi (nyi)
  DO i = 0, ny
     y (i) = (yi (nyi) - yi (0) ) * dble (i) / dble (ny)
  ENDDO
  ! interpolation no longer implemented
  rhoo (0:nx, 0:ny) = rhoi (0:nx, 0:ny)
  rhomin = minval (rhoo(0:nx, 0:ny))
  rhomax = maxval (rhoo(0:nx, 0:ny))

  IF (rhomin > 0.d0) THEN
     WRITE( stdout,'("Logarithmic scale (y/n)? > ")', advance="NO")
     READ (5, '(a)') ans
     logarithmic_scale = ans/='n'.and.ans/='N'
  ELSE
     logarithmic_scale = .false.
  ENDIF
10 CONTINUE
  WRITE( stdout, '("Bounds: ",2f12.6)') rhomin, rhomax
  WRITE( stdout, '("min, max, # of levels > ")', advance="NO")
  READ (5, * ) rhoomin, rhoomax, nlevels
  IF ( rhoomax <= rhoomin .or. &
       rhoomin >= rhomax .or. rhoomax <= rhomin ) THEN
     WRITE( stdout, '("Out of Bounds! try again")')
     GOTO 10
  ENDIF
  IF (nlevels > 1000) THEN
     WRITE( stdout, '("Are you sure you really need ",i8," levels?")') &
          nlevels
  ELSEIF (nlevels < 1) THEN
     WRITE( stdout, '("Too few levels! assuming 1 level")')
     nlevels = 1
  ENDIF
  ALLOCATE (z(0:nlevels))
  IF (logarithmic_scale) THEN
     DO k = 0, nlevels - 1
        z (k) = exp (log (rhoomin) + (log (rhoomax) - log (rhoomin) ) &
             * dble (k) / (nlevels - 1) )
     ENDDO
  ELSE
     DO k = 0, nlevels - 1
        z (k) = rhoomin + (rhoomax - rhoomin) * dble (k) / (nlevels - 1)
     ENDDO
  ENDIF

  z (nlevels) = z (nlevels - 1)
  xdim = 15.0d0 * (xmax - xmin) / sqrt ( (xmax - xmin) **2 + (ymax - ymin) **2)
  ydim = 15.0d0 * (ymax - ymin) / sqrt ( (xmax - xmin) **2 + (ymax - ymin) **2)

  xs = 4.0d0
  ys = 3.0d0

  ! uncomment the call to "cplot" if you want contour lines,
  ! plus gray levels and shading for negative values

  CALL cplot (rhoo, nx, ny, x, xmin, xmax, nx, y, ymin, ymax, &
       ny, nlevels, z, xdim, ydim, xs, ys, filename, fileout)

  ! uncomment the call to "psplot" if you want contour lines
  ! of various kinds: solid, dashed, etc

  ! call psplot ( rhoo, nx, x, nx, y, ny, nlevels, z, xdim, ydim, &
  !      xs, ys, fileout)

  CALL atomi (nat, tau, ityp, at, a0, r0, tau1, tau2, xdim, ydim)
20 STOP
   DEALLOCATE (z)
   DEALLOCATE (rhoo, x , y )
   DEALLOCATE (tau, ityp, z)
   DEALLOCATE (rhoi, xi, yi)
END PROGRAM plotrho
!
!-----------------------------------------------------------------------
SUBROUTINE cplot (d, imax, jmax, x, xmin, xmax, iub, y, ymin, &
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
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER :: imax, jmax, iub, jub, nc
  real(DP) :: d (0:imax, 0:jmax), x (0:imax), y (0:jmax), z (0:nc)
  real(DP) :: xmin, xmax, ymin, ymax, xdim, ydim, xs, ys

  CHARACTER (len=*) :: filename, str
  INTEGER, PARAMETER :: ncmax = 19
  INTEGER :: i, j, k
  real(DP) :: gray (0:ncmax), dim
  real(DP) :: cm = 28.453d0, width=0.5d0, gray0=1.0d0, deltagray=0.7d0
  !  cm   : number of points per cm
  !  width: linewidth of the contour plot for PostScript printer

  OPEN (unit = 1, file = filename, status = 'unknown', form = &
       'formatted')
  IF (nc > ncmax .or. nc < 1) STOP '  nc too big or wrong'
  IF (iub > imax .or.iub < 1) STOP ' iub too big or wrong'
  IF (jub > jmax .or.jub < 1) STOP ' jub too big or wrong'
  IF (xdim < 3.0d0 .or. ydim < 3.0d0) STOP ' really too small!'
  IF (xdim > 20.0.or. ydim > 30.0) STOP ' really too big!'
  IF (abs (xs) > 20.or. abs (ys) > 30) STOP ' xs or ys are weird'
  !  initializations for PostScript output
  WRITE (1, '(a)') '%! PS-Adobe-1.0'
  WRITE (1, '("%%BoundingBox:",4f6.1)') xs * cm, ys * cm, (xs + &
       xdim) * cm, (ys + ydim) * cm
  WRITE (1, '(a)') '/localdict 100 dict def'
  WRITE (1, '(a)') 'localdict begin'
  WRITE (1, '(a)') '/cm {28.453 mul} def'
  WRITE (1, '(a)') '/title {('//str//')} def'
  WRITE (1, '(a)') '/Times-Roman findfont 12 scalefont setfont'
  WRITE (1, '(a)') '% cshow prints a centered string at current position'
  WRITE (1, '(a)') '/cshow {gsave dup stringwidth pop 2 div neg 0'
  WRITE (1, '(a)') ' rmoveto show grestore} def'
  WRITE (1, '(a)') '% x1 y1 x2 y2 p :  draws a segment from point 1 to point 2'
  WRITE (1, '(a)') '/p {0 setgray newpath moveto lineto stroke} def'
  WRITE (1, '(a)') '% x1 y1 x2 y2 x3 y3 x4 y4 sn :'
  WRITE (1, '(a)') '% fills the region bounded by points 1 to 4'
  WRITE (1, '(a)') '% with greyscale n'
  !  type of gray for shaded areas
  DO k = 0, nc
     gray (k) = gray0 - k * deltagray / nc
     IF (k<10) THEN
        WRITE (1, '("/s",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto lineto fill} def")') k, gray (k)
        WRITE (1, '("/t",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto fill} def")') k, gray (k)
     ELSE
        WRITE (1, '("/u",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto lineto fill} def")') mod (k, 10) , &
             &gray (k)
        WRITE (1, '("/v",i1," {",f4.2," setgray newpath ", &
             &       "moveto lineto lineto fill} def")') mod (k, 10) , gray (k &
             &)
     ENDIF

  ENDDO
  WRITE (1, '(a)') '%%EndPreamble'
  WRITE (1, '(a)') 'gsave'
  WRITE (1, '(1x,f6.2," cm ",f6.2," cm translate")') xs, ys
  WRITE (1, '(a)') '% Uncomment next line if you want a big picture'
  WRITE (1, '(a)') '% 1.8 1.8 scale'
  WRITE (1, '(f7.3," setlinewidth")') width
  WRITE (1, '(a)') '% Comment the next line to remove the title'

  WRITE (1, '(1x,f6.2," cm ",f6.2," cm moveto title cshow")') &
       xdim / 2, ydim + 1.5d0

  CALL hatch (0.d0, xdim, 0.d0, ydim)
  DO i = 0, iub
     x (i) = xdim * dble (i) / iub
  ENDDO
  DO j = 0, jub
     y (j) = ydim * dble (j) / jub
  ENDDO

  CALL conrec (imax, iub, jmax, jub, x, y, d, nc, z)
  !  draw frame of size xdim by ydim
  WRITE (1, '(a)') '1 setlinewidth 0 setgray newpath'
  WRITE (1, '(2f6.1," moveto")') 0.0d0, 0.0d0
  WRITE (1, '(2f6.1," lineto")') xdim * cm, 0.0d0
  WRITE (1, '(2f6.1," lineto")') xdim * cm, ydim * cm
  WRITE (1, '(2f6.1," lineto")') 0.0d0, ydim * cm

  WRITE (1, '(a)') 'closepath stroke'
  !      write (1,'(a)') 'grestore'
  !      write (1,'(a)') '%%Trailer'
  !      write (1,'(a)') 'showpage'
  !      close(1)
  RETURN


END SUBROUTINE cplot

SUBROUTINE conrec (imax, iub, jmax, jub, x, y, d, nc, z)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER :: imax, iub, jmax, jub, nc

  real(DP) :: d (0:imax, 0:jmax), x (0:imax), y (0:jmax), z (0:nc)
  INTEGER, PARAMETER :: ncmax = 19
  CHARACTER (len=4) :: triangle (0:ncmax), trapez (0:ncmax)
  real(DP) :: h (0:4), xh (0:4), yh (0:4)
  real(DP) :: x1, y1, x2, y2, x3, y3, x4, y4, dx, dy, xx, yy, dmin, dmax
  real(DP) :: cm = 28.453d0
  !   cm : conversion factor from cm to points for PostScript
  INTEGER :: ish (0:4), im (0:3), jm (0:3), castab (0:2, 0:2, 0:2)
  INTEGER :: i, j, k, m, m1, m2, m3, npoint, icase, levelin, nolevel
  data (im (i), i = 0, 3) / 0, 1, 1, 0 /
  data (jm (i), i = 0, 3) / 0, 0, 1, 1 /
  data ( ( (castab (i, j, k), k = 0, 2), j = 0, 2), i = 0, 2) &
       / 0, 0, 8, 0, 2, 5, 7, 6, 9, 0, 3, 4, 1, 3, 1, 4, 3, 0, 9, 6, 7, &
       5, 2, 0, 8, 0, 0 /

  dy = (y (jub) - y (0) ) / (nc + 1)
  xx = x (iub) + 1.0d0
  dx = 0.5d0

  WRITE (1, '(a)') '% Start of Color Code'

  CALL hatch (xx, xx + dx, y (0), y (jub) )
  DO k = 0, nc
     yy = y (jub) - k * dy
     WRITE (1, '(8f6.1)', advance="NO") xx * cm, yy * cm, (xx + dx) * cm, yy * cm, &
          (xx + dx) * cm, (yy - dy) * cm, xx * cm, (yy - dy) * cm
     IF (k < 10) THEN
        WRITE (triangle (k) , '(" t",i1,1x)') k
        WRITE (trapez (k) , '(" s",i1,1x)') k
     ELSE
        WRITE (triangle (k) , '(" v",i1)') mod (k, 10)
        WRITE (trapez (k) , '(" u",i1)') mod (k, 10)
     ENDIF
     WRITE (1, '(a4)') trapez (k)
     WRITE (1, * ) '0 setgray newpath'
     WRITE (1, '(2f6.1," moveto")') xx * cm, yy * cm
     WRITE (1, '(2f6.1," lineto")') (xx + dx) * cm, yy * cm
     WRITE (1, '(2f6.1," lineto")') (xx + dx) * cm, (yy - dy) &
          * cm
     WRITE (1, '(2f6.1," lineto")') xx * cm, (yy - dy) * cm
     WRITE (1, * ) 'closepath stroke'
     WRITE (1, '(2f6.1," moveto")') (x (iub) + 2.0d0) * cm, (yy - dy / &
          2) * cm
     IF (k == 0) THEN
        WRITE (1, '("(z<",f7.5,") show")') z (0)
     ELSEIF (k==nc) THEN
        WRITE (1, '("(z>",f7.5,") show")') z (nc - 1)
     ELSE
        WRITE (1, '("(",f7.5,"<z<",f7.5,") show")') z (k - 1) , &
             z (k)
     ENDIF
  ENDDO
  WRITE (1, '(a)') '% End of Color Code'
  DO k = 1, nc - 1
     IF (z (k) <=z (k - 1) ) STOP 'zk order'
  ENDDO
  !   scan the array, top down, left to right, to paint shaded areas
  DO j = jub - 1, 0, - 1
     DO i = 0, iub - 1
        !  find lowest and highest vertex
        dmin = min (d (i, j), d (i, j + 1), d (i + 1, j), d (i + 1, j + 1) &
             )
        dmax = max (d (i, j), d (i, j + 1), d (i + 1, j), d (i + 1, j + 1) &
             )
        !  search for levels in this box
        nolevel = 0
        DO k = 0, nc - 1
           IF (z (k) < dmin) nolevel = k + 1
           IF (z (k) >= dmin .and. z (k) <= dmax) THEN
              levelin = k
              GOTO 10
           ENDIF
        ENDDO
        !  no level in this box: paint the whole box and pass to another box
        WRITE (1, '(8f6.1,a4)') x (i) * cm, y (j) * cm, x (i + 1) * cm, y &
             (j) * cm, x (i + 1) * cm, y (j + 1) * cm, x (i) * cm, y (j + 1) &
             * cm, trapez (nolevel)
        GOTO 100
        !  there is at least a level in this box: paint the whole box
10      CONTINUE
        WRITE (1, '(8f6.1,a4)') x (i) * cm, y (j) * cm, x (i + 1) * cm, y &
             (j) * cm, x (i + 1) * cm, y (j + 1) * cm, x (i) * cm, y (j + 1) &
             * cm, trapez (levelin)
        DO k = levelin, nc - 1
           !  if no more levels in this box, move to another box
           IF (z (k) >dmax) GOTO 100
           !  find contour of zero levels in this box
           DO m = 1, 4
              h (m) = d (i + im (m - 1), j + jm (m - 1) ) - z (k)
              xh (m) = x (i + im (m - 1) )
              yh (m) = y (j + jm (m - 1) )
           ENDDO
           h (0) = (h (1) + h (2) + h (3) + h (4) ) / 4
           xh (0) = (x (i) + x (i + 1) ) / 2
           yh (0) = (y (j) + y (j + 1) ) / 2
           DO m = 0, 4
              IF (h (m) >0) THEN
                 ish (m) = 2
              ELSEIF (h (m) <0) THEN
                 ish (m) = 0
              ELSE
                 ish (m) = 1
              ENDIF
           ENDDO
           !   scan each triangle in the box to paint shaded areas
           DO m = 1, 4
              m1 = m
              m2 = 0
              m3 = mod (m, 4) + 1
              npoint = 0
              icase = castab (ish (m1), ish (m2), ish (m3) )
              IF (icase == 0) THEN
                 IF (ish (m1) ==2) THEN
                    !     paint this triangle if positive
                    x1 = xh (m1)
                    y1 = yh (m1)
                    x2 = xh (m2)
                    y2 = yh (m2)
                    x3 = xh (m3)
                    y3 = yh (m3)
                    npoint = 3
                 ENDIF
              ELSEIF (icase == 1) THEN
                 !     line between vertices m1 and m2
                 x1 = xh (m1)
                 y1 = yh (m1)
                 x2 = xh (m2)
                 y2 = yh (m2)
                 IF (ish (m3) ==2) THEN
                    x3 = xh (m3)
                    y3 = yh (m3)
                    npoint = 3
                 ENDIF
              ELSEIF (icase == 2) THEN
                 !     line between vertices m2 and m3
                 x1 = xh (m2)
                 y1 = yh (m2)
                 x2 = xh (m3)
                 y2 = yh (m3)
                 IF (ish (m1)  == 2) THEN
                    x3 = xh (m1)
                    y3 = yh (m1)
                    npoint = 3
                 ENDIF
              ELSEIF (icase == 3) THEN
                 !     line between vertices m3 and m1
                 x1 = xh (m3)
                 y1 = yh (m3)
                 x2 = xh (m1)
                 y2 = yh (m1)
                 IF (ish (m2)  == 2) THEN
                    x3 = xh (m2)
                    y3 = yh (m2)
                    npoint = 3
                 ENDIF
              ELSEIF (icase == 4) THEN
                 !     line between vertex m1 and side m2-m3
                 x1 = xh (m1)
                 y1 = yh (m1)
                 x2 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y2 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 IF (ish (m3)  == 2) THEN
                    x3 = xh (m3)
                    y3 = yh (m3)
                 ELSE
                    x3 = xh (m2)
                    y3 = yh (m2)
                 ENDIF
                 npoint = 3
              ELSEIF (icase == 5) THEN
                 !   line between vertex m2 and side m3-m1
                 x1 = xh (m2)
                 y1 = yh (m2)
                 x2 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y2 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 IF (ish (m1)  == 2) THEN
                    x3 = xh (m1)
                    y3 = yh (m1)
                 ELSE
                    x3 = xh (m3)
                    y3 = yh (m3)
                 ENDIF
                 npoint = 3
              ELSEIF (icase == 6) THEN
                 !     line between vertex m3 and line m1-m2
                 x1 = xh (m3)
                 y1 = yh (m3)
                 x2 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y2 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 IF (ish (m2)  == 2) THEN
                    x3 = xh (m2)
                    y3 = yh (m2)
                 ELSE
                    x3 = xh (m1)
                    y3 = yh (m1)
                 ENDIF
                 npoint = 3
              ELSEIF (icase == 7) THEN
                 !     line between sides m1-m2 and m2-m3
                 x1 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y1 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 x2 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y2 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 IF (ish (m2)  == 2) THEN
                    x3 = xh (m2)
                    y3 = yh (m2)
                    npoint = 3
                 ELSE
                    x3 = xh (m3)
                    y3 = yh (m3)
                    x4 = xh (m1)
                    y4 = yh (m1)
                    npoint = 4
                 ENDIF
              ELSEIF (icase == 8) THEN
                 !   line between sides m2-m3 and m3-m1
                 x1 = (h (m3) * xh (m2) - h (m2) * xh (m3) ) / (h (m3) - h (m2) &
                      )
                 y1 = (h (m3) * yh (m2) - h (m2) * yh (m3) ) / (h (m3) - h (m2) &
                      )
                 x2 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y2 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 IF (ish (m3)  == 2) THEN
                    x3 = xh (m3)
                    y3 = yh (m3)
                    npoint = 3
                 ELSE
                    x3 = xh (m1)
                    y3 = yh (m1)
                    x4 = xh (m2)
                    y4 = yh (m2)
                    npoint = 4
                 ENDIF
              ELSEIF (icase == 9) THEN
                 !     line between sides m3-m1 and m1-m2
                 x1 = (h (m1) * xh (m3) - h (m3) * xh (m1) ) / (h (m1) - h (m3) &
                      )
                 y1 = (h (m1) * yh (m3) - h (m3) * yh (m1) ) / (h (m1) - h (m3) &
                      )
                 x2 = (h (m2) * xh (m1) - h (m1) * xh (m2) ) / (h (m2) - h (m1) &
                      )
                 y2 = (h (m2) * yh (m1) - h (m1) * yh (m2) ) / (h (m2) - h (m1) &
                      )
                 IF (ish (m1)  == 2) THEN
                    x3 = xh (m1)
                    y3 = yh (m1)
                    npoint = 3
                 ELSE
                    x3 = xh (m2)
                    y3 = yh (m2)
                    x4 = xh (m3)
                    y4 = yh (m3)
                    npoint = 4
                 ENDIF
              ENDIF
              IF (npoint == 3) THEN
                 WRITE (1, '(6f6.1,a4)') x1 * cm, y1 * cm, x2 * cm, y2 * cm, x3 &
                      * cm, y3 * cm, triangle (k + 1)
              ELSEIF (npoint == 4) THEN
                 WRITE (1, '(8f6.1,a4)') x1 * cm, y1 * cm, x2 * cm, y2 * cm, x3 &
                      * cm, y3 * cm, x4 * cm, y4 * cm, trapez (k + 1)
              ENDIF
              IF (icase/=0) WRITE (1, '(4f6.1," p")') x1 * cm, y1 * cm, x2 * &
                   cm, y2 * cm
           ENDDO
        ENDDO
100     CONTINUE
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE conrec

!
!-----------------------------------------------------------------------
SUBROUTINE atomi (nat, tau, ityp, at, a0, r0, tau1, tau2, xdim, ydim)
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER :: nat, ityp (nat)
  real(DP) :: tau (3, nat), at (3, 3), r0 (3), tau1 (3), tau2 (3), a0
  real(DP) :: xdim, ydim
  INTEGER :: n1, n2, n3, i, n
  real(DP) :: r (3), ri (3), tau1n, tau2n, delta0, r1, r2, r3
  real(DP) :: delta = 1.0d0, cm = 28.453d0
  !
  delta0 = delta / a0
  tau1n = sqrt (tau1 (1) **2 + tau1 (2) **2 + tau1 (3) **2)
  tau2n = sqrt (tau2 (1) **2 + tau2 (2) **2 + tau2 (3) **2)
  !  clip the rectangle
  WRITE (1, '(a)') 'gsave newpath'
  WRITE (1, '(2f6.1," moveto")') 0.0d0, 0.0d0
  WRITE (1, '(2f6.1," lineto")') xdim * cm, 0.0d0
  WRITE (1, '(2f6.1," lineto")') xdim * cm, ydim * cm
  WRITE (1, '(2f6.1," lineto")') 0.0d0, ydim * cm
  WRITE (1, '(a)') 'closepath clip stroke'
  ! really bad algorithm to generate (hopefully) all atoms
  ! that are inside the contour plane or at |z| < delta (a.u.)
  DO n1 = - 3, + 3
     DO n2 = - 3, + 3
        DO n3 = - 3, + 3
           DO i = 1, 3
              r (i) = n1 * at (i, 1) + n2 * at (i, 2) + n3 * at (i, 3)
           ENDDO
           DO n = 1, nat
              DO i = 1, 3
                 ri (i) = tau (i, n) + r (i) - r0 (i)
              ENDDO
              !  this is the component along direction 1 ...
              r1 = (ri (1) * tau1 (1) + ri (2) * tau1 (2) + ri (3) * tau1 (3))&
                   / tau1n
              IF (r1 > - delta0 .and. r1 < tau1n + delta0) THEN
                 !  and this is along direction 2 ...
                 r2 = (ri (1) * tau2(1) + ri (2) * tau2(2) + ri (3) * tau2(3))&
                      / tau2n
                 IF (r2 > - delta0 .and. r2 < tau2n + delta0) THEN
                    DO i = 1, 3
                       ri (i) = ri (i) - r1 * tau1 (i) / tau1n &
                                       - r2 * tau2 (i) / tau2n
                    ENDDO
                    r3 = sqrt (ri (1) **2 + ri (2) **2 + ri (3) **2)
                    !  and this is along the direction orthogonal to plane
                    IF (abs (r3) < delta0) THEN
                       WRITE (1, '(3f6.1," 0 360 arc gsave ",f4.2, &
                            &      " setgray fill grestore stroke")') &
                            r1 / tau1n * xdim * cm, r2 / tau2n * ydim * cm,&
                            delta0 / tau1n * xdim * cm, abs (r3) / delta0
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  WRITE (1, '(a)') 'grestore'
  WRITE (1, '(a)') 'grestore'
  WRITE (1, '(a)') '%%Trailer'
  WRITE (1, '(a)') 'showpage'

  CLOSE (1)
  RETURN


END SUBROUTINE atomi

SUBROUTINE hatch (x1, x2, y1, y2)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  real(DP) :: x1, x2, y1, y2
  real(DP) :: cm =  28.453d0, delta = 0.2d0, dim
  INTEGER :: nhach, n

  WRITE (1, '(a)') '% Beginning of hatching'
  WRITE (1, '(a)') 'gsave newpath'
  WRITE (1, '(2f6.1," moveto")') x1 * cm, y1 * cm
  WRITE (1, '(2f6.1," lineto")') x2 * cm, y1 * cm
  WRITE (1, '(2f6.1," lineto")') x2 * cm, y2 * cm
  WRITE (1, '(2f6.1," lineto")') x1 * cm, y2 * cm
  WRITE (1, '(a)') 'closepath clip'

  dim = max (x2 - x1, y2 - y1)
  nhach = dim / delta
  !      delta=dim/nhach
  DO n = 1, nhach
     WRITE (1, '(4f6.1," p")') (x1 + (n - 1) * delta) * cm, y1 * cm, &
          (x1 + dim) * cm, (y1 + dim - (n - 1) * delta) * cm
     WRITE (1, '(4f6.1," p")') x1 * cm, (y1 + (n - 1) * delta) &
          * cm, (x1 + dim - (n - 1) * delta) * cm, (y1 + dim) * cm
  ENDDO
  WRITE (1, '(a)') 'grestore'
  WRITE (1, '(a)') '% End of hatching'

  RETURN
END SUBROUTINE hatch
!
!-----------------------------------------------------------------------
SUBROUTINE psplot ( d, imax, x, iub, y, jub, nlevels, z, &
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
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: ncontourmax=100, npointmax=500, nmaxtypes=8
  INTEGER :: imax, iub, jub, nlevels
  real(DP) :: d(0:imax,0:jub), x(0:imax), y(0:jub), z(0:nlevels), &
       xdim, ydim, xs, ys
  !
  real(DP) :: line(2,npointmax,ncontourmax), work(2,npointmax), &
       segments(2,2,npointmax*ncontourmax), yy, dy
  INTEGER :: npoints(ncontourmax), ncontours, i, j, k, l, npts
  real(DP), PARAMETER :: cm = 28.453d0
  LOGICAL ::  lwork(npointmax*ncontourmax)
  CHARACTER (len=256) :: filename, linetype(nmaxtypes)
  data linetype &
       / '[1 2    ] 0 setdash 0.5 setlinewidth', &
       '[1 2 4 2] 0 setdash', &
       '[4 2    ] 0 setdash', &
       '[       ] 0 setdash', &
       '[1 2   ] 0 setdash 1.0 setlinewidth', &
       '[1 2 4 2] 0 setdash', &
       '[4 2    ] 0 setdash', &
       '[       ] 0 setdash'/

  OPEN (unit = 1, file = filename, status = 'unknown', form = &
       'formatted')
  IF(nlevels < 1)  STOP '  nunmber of levels is wrong'
  IF(iub > imax .or. iub < 1) STOP ' iub too big or wrong'
  IF(jub < 1) STOP ' jub or wrong'
  IF(xdim < 3.0d0 .or. ydim < 3.0d0 ) STOP ' really too small!'
  IF(xdim > 20.0.or. ydim > 30.0) STOP ' really too big!'
  IF(abs(xs) > 20 .or. abs(ys) > 30) STOP ' xs or ys are weird'

  !  initializations for PostScript output
  WRITE (1,'(a)') '%! PS-Adobe-1.0'
  WRITE (1,'(''%%BoundingBox:'',4f6.1)') &
       xs*cm, ys*cm,(xs+xdim)*cm,(ys+ydim)*cm
  WRITE (1,'(a)') '/localdict 100 dict def'
  WRITE (1,'(a)') 'localdict begin'
  WRITE (1,'(a)') '/cm {28.453 mul} def'
  WRITE (1,'(a)') '/title {('//filename//')} def'
  WRITE (1,'(a)') '/Times-Roman findfont 12 scalefont setfont'
  WRITE (1,'(a)') '% cshow prints a centered string at current position'
  WRITE (1,'(a)') '/cshow {gsave dup stringwidth pop 2 div neg 0'
  WRITE (1,'(a)') ' rmoveto show grestore} def'
  WRITE (1,'(a)') '% x1 y1 nm x2 y2 ...xn yn n-1 lrs :'
  WRITE (1,'(a)') '% draws a line from point 1 to point n'
  WRITE (1,'(a)') '/nm  {newpath moveto} def'
  WRITE (1,'(a)') '/lrs {{lineto} repeat stroke} def'
  WRITE (1,'(a)') '%%EndPreamble'
  WRITE (1,'(a)') 'gsave'
  WRITE (1,'(a)') '0 setgray'
  WRITE (1,'(1x,f6.2,'' cm '',f6.2,'' cm translate'')') xs, ys
  WRITE (1,'(1x,f6.2,'' cm '',f6.2,'' cm moveto title cshow'')') &
       xdim/2, ydim+1.5d0
  !
  DO i=0,iub
     x(i)=xdim*dble(i)/iub
  ENDDO
  DO j=0,jub
     y(j)=ydim*dble(j)/jub
  ENDDO

  DO k=1,nlevels-1
     IF (z(k)<=z(k-1)) STOP ' order of levels'
  ENDDO

  yy=y(jub)*0.75d0
  dy=(y(jub)-y(0))/(nlevels-1)/2

  DO k=0,nlevels-1
     ! type of line
     WRITE (1,'(''%% This is level '',i3)') k
     WRITE (1,'(50a)') linetype(mod(k,nmaxtypes)+1)
     WRITE (1,*) ' newpath '
     WRITE (1,'(2f6.1,'' moveto'')') (x(iub)+1.0d0)*cm, yy*cm
     WRITE (1,'(2f6.1,'' lineto'')') (x(iub)+1.8d0)*cm, yy*cm
     WRITE (1,*) 'closepath stroke'
     WRITE (1,'(2f6.1,'' moveto'')') (x(iub)+2.0d0)*cm, yy*cm
     WRITE (1,'(''( z = '',f6.3,'') show'')') z(k)
     yy=yy-dy

     CALL findcontours ( d, imax, iub, jub, x, y, z(k), segments, work, &
          lwork, npointmax, ncontourmax, line, npoints, ncontours)

     DO l=1,ncontours
        WRITE (1,'(12f6.1)') &
             ( line(1,npts,l)*cm, line(2,npts,l)*cm, npts=npoints(l),1,-1 )
        WRITE (1,'(''nm'',i4,'' lrs'')') npoints(l)-1
     ENDDO

  ENDDO

  !  draw frame of size xdim by ydim
  WRITE (1,*) '[] 0 setdash 1 setlinewidth 0 setgray newpath'
  WRITE (1,'(2f6.1,'' moveto'')')     0.0d0,0.0d0
  WRITE (1,'(2f6.1,'' lineto'')') xdim*cm,0.0d0
  WRITE (1,'(2f6.1,'' lineto'')') xdim*cm,ydim*cm
  WRITE (1,'(2f6.1,'' lineto'')')     0.0d0,ydim*cm
  WRITE (1,*) 'closepath stroke'
  !      write (1,*) 'grestore'
  !      write (1,'(a)') '%%Trailer'
  !      write (1,*) 'showpage'
  RETURN
END SUBROUTINE psplot

!
SUBROUTINE findcontours (d, imax, iub, jub, x, y, z, segments, work, &
     done, npointmax, ncontourmax, line, npoints, ncontours)
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER :: imax, iub, jub, npointmax, ncontourmax,  ncontours, &
       npoints(ncontourmax)
  real(DP) :: d(0:imax,0:jub),x(0:iub), y(0:jub), z, &
       line(2,npointmax,ncontourmax), segments(2,2,npointmax*ncontourmax)
 LOGICAL :: done(npointmax*ncontourmax)
  !
  INTEGER :: i, j, m, i0, j0, m0, nsegments, nsegmax, nseg, nnext, npts, &
       isame, inext, nleft
  real(DP) :: epsx, epsy, work(2,npointmax)
  LOGICAL :: found, start_from_boundary, hit_boundary, loop_closed

  nsegmax=npointmax*ncontourmax
  nsegments=0
  DO j0=jub-1,0,-1
     DO i0=0,iub-1
        DO m0=1,4
           !     scan boxes and triangles until a contour is found
           CALL triangle ( d, x, y, z, imax, iub, jub, &
                i0, j0, m0, found, segments(1,1,nsegments+1) )
           IF (found) THEN
              nsegments=nsegments+1
              IF (nsegments > nsegmax) THEN
                 PRINT '('' Erore: nsegments, nsegmax '',2i5)', &
                      nsegments, nsegmax
                 STOP
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  DO nseg=1,nsegments
     done(nseg)=.false.
  ENDDO

  ncontours=0
  nleft=nsegments
  epsx = 0.0001d0 * (x(1)-x(0))
  epsy = 0.0001d0 * (y(1)-y(0))
  DO nseg=1,nsegments
     IF (.not.done(nseg)) THEN
        ncontours=ncontours+1
        IF (ncontours>ncontourmax) THEN
           PRINT '('' too many ('',i4,'') contours'')', ncontours
           STOP
        ENDIF
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
        IF (npts > npointmax) THEN
           PRINT '('' too many ('',i5,'') points'')', npts
           STOP
        ENDIF
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
        IF (nleft == 0) THEN
           IF (.not.hit_boundary .and. .not.loop_closed) &
                PRINT '('' Warning: no more points'')'
           npoints(ncontours)=npts
           RETURN
        ENDIF
        IF (start_from_boundary) THEN
           IF (hit_boundary) THEN
              npoints(ncontours)=npts
              GOTO 20
           ELSEIF (loop_closed) THEN
              PRINT '('' Warning: loop closed on the boundary ?'')'
              GOTO 20
           ENDIF
        ELSE
           IF (loop_closed) THEN
              npoints(ncontours)=npts
              GOTO 20
           ELSEIF (hit_boundary) THEN
              DO nnext=1,npts
                 work(1,nnext)=line(1,npts-nnext+1,ncontours)
                 work(2,nnext)=line(2,npts-nnext+1,ncontours)
              ENDDO
              DO nnext=1,npts
                 line(1,nnext,ncontours)=work(1,nnext)
                 line(2,nnext,ncontours)=work(2,nnext)
              ENDDO
              start_from_boundary=.true.
           ENDIF
        ENDIF

        DO nnext=1,nsegments
           IF (.not.done(nnext)) THEN
              DO isame=1,2
                 found= ( abs(segments(1,isame,nnext) - &
                          line(1,npts,ncontours)) < epsx ) .and. &
                        ( abs(segments(2,isame,nnext) - &
                          line(2,npts,ncontours)) < epsy )
                 IF (found) THEN
                    inext=mod(isame,2)+1
                    !inext =1 if isame=2 and vice versa
                    GOTO 10
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
        npoints(ncontours)=npts-1
        PRINT '('' Warning: next point not found'',i5)', nseg
     ENDIF
20   CONTINUE
  ENDDO

  RETURN
END SUBROUTINE findcontours
!
SUBROUTINE triangle ( d, x, y, z, imax, iub, jub, i0, j0, m0, found, line )
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER :: imax, iub, jub, i0, j0, m0
  real(DP) :: d(0:imax,0:jub), x(0:iub), y(0:jub), z, line(2,2)
  LOGICAL :: found
  !
  real(DP) :: h(0:4), xh(0:4), yh(0:4), dmin, dmax
  INTEGER :: m, icase, m1, m2, m3, i, j, k, ish(0:4), castab(-1:1,-1:1,-1:1)
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
  IF ( dmin > 0.0d0 .or. dmax < 0.0d0 ) THEN
     found=.false.
     RETURN
  ENDIF
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
  DO m=0,4
     ish(m) = nint( sign(1.d0, h(m)) )
     IF (h(m) == 0.0d0) ish(m) = 0
  ENDDO
  !     starting triangle
  m1=m0
  m2=0
  m3=mod(m0,4)+1
  icase=castab(ish(m1),ish(m2),ish(m3))
  IF (icase == 0) THEN
     !     no contour in this triangle
     found=.false.
     RETURN
     ! NOTA BENE: if there is a point on the boundary,
     ! this is always chosen to be the first point
  ELSEIF (icase == 1) THEN
     !     line between vertices m1 and m2
     line(1,1)=xh(m1)
     line(2,1)=yh(m1)
     line(1,2)=xh(m2)
     line(2,2)=yh(m2)
  ELSEIF (icase == 2) THEN
     !     line between vertices m2 and m3
     line(1,1)=xh(m3)
     line(2,1)=yh(m3)
     line(1,2)=xh(m2)
     line(2,2)=yh(m2)
  ELSEIF (icase == 3) THEN
     !     line between vertices m3 and m1
     line(1,1)=xh(m3)
     line(2,1)=yh(m3)
     line(1,2)=xh(m1)
     line(2,2)=yh(m1)
  ELSEIF (icase == 4) THEN
     !     line between vertex m1 and side m2-m3
     line(1,1)=xh(m1)
     line(2,1)=yh(m1)
     line(1,2)=(h(m3)*xh(m2)-h(m2)*xh(m3))/(h(m3)-h(m2))
     line(2,2)=(h(m3)*yh(m2)-h(m2)*yh(m3))/(h(m3)-h(m2))
  ELSEIF (icase == 5) THEN
     !     line between vertex m2 and side m3-m1
     line(1,1)=(h(m1)*xh(m3)-h(m3)*xh(m1))/(h(m1)-h(m3))
     line(2,1)=(h(m1)*yh(m3)-h(m3)*yh(m1))/(h(m1)-h(m3))
     line(1,2)=xh(m2)
     line(2,2)=yh(m2)
  ELSEIF (icase == 6) THEN
     !     line between vertex m3 and line m1-m2
     line(1,1)=xh(m3)
     line(2,1)=yh(m3)
     line(1,2)=(h(m2)*xh(m1)-h(m1)*xh(m2))/(h(m2)-h(m1))
     line(2,2)=(h(m2)*yh(m1)-h(m1)*yh(m2))/(h(m2)-h(m1))
  ELSEIF (icase == 7) THEN
     !     line between sides m1-m2 and m2-m3
     line(1,1)=(h(m2)*xh(m1)-h(m1)*xh(m2))/(h(m2)-h(m1))
     line(2,1)=(h(m2)*yh(m1)-h(m1)*yh(m2))/(h(m2)-h(m1))
     line(1,2)=(h(m3)*xh(m2)-h(m2)*xh(m3))/(h(m3)-h(m2))
     line(2,2)=(h(m3)*yh(m2)-h(m2)*yh(m3))/(h(m3)-h(m2))
  ELSEIF (icase == 8) THEN
     !     line between sides m2-m3 and m3-m1
     line(1,1)=(h(m1)*xh(m3)-h(m3)*xh(m1))/(h(m1)-h(m3))
     line(2,1)=(h(m1)*yh(m3)-h(m3)*yh(m1))/(h(m1)-h(m3))
     line(1,2)=(h(m3)*xh(m2)-h(m2)*xh(m3))/(h(m3)-h(m2))
     line(2,2)=(h(m3)*yh(m2)-h(m2)*yh(m3))/(h(m3)-h(m2))
  ELSEIF (icase == 9) THEN
     !     line between sides m3-m1 and m1-m2
     line(1,1)=(h(m1)*xh(m3)-h(m3)*xh(m1))/(h(m1)-h(m3))
     line(2,1)=(h(m1)*yh(m3)-h(m3)*yh(m1))/(h(m1)-h(m3))
     line(1,2)=(h(m2)*xh(m1)-h(m1)*xh(m2))/(h(m2)-h(m1))
     line(2,2)=(h(m2)*yh(m1)-h(m1)*yh(m2))/(h(m2)-h(m1))
  ENDIF

  found=.true.

  RETURN
END SUBROUTINE triangle
