!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------!
!   This module is basad on a similar module from CP2K
!-----------------------------------------------------------------------------!

    MODULE splines

! routines for handling splines
!   allocate_spline: allocates x and y vectors for splines
!   init_spline:   generate table for spline (allocate spl%y2)
!   spline:        return value of spline for given abscissa (optional:also y1)
!   spline_1:      return value of 1. derivative of spline for given abscissa
!   spline_int:    return value of integral on given interval of spline
!   kill_spline:   destructor ( spl%x,y und/oder spl%y2)

! NB: splines are always "natural splines", i.e. values of first
!     derivative at the end-points cannot be specified
!-----------------------------------------------------------------------------!
      USE kinds, ONLY : DP

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: spline_data, allocate_spline, init_spline, spline, spline_1, &
        spline_int, kill_spline, splineh, splinedx, splintdx, nullify_spline

      TYPE spline_data
        REAL (DP), POINTER :: x(:)  ! array containing x values
        REAL (DP), POINTER :: y(:)  ! array containing y values 
                                     ! y(i) is the function value corresponding
                                     ! to x(i) in the interpolation table
        REAL (DP), POINTER :: y2(:) ! second derivative of interpolating function
        INTEGER :: n       ! number of element in the interpolation table
        INTEGER :: pos
        REAL (DP) :: h, invh, h26, h16 
        REAL (DP) :: xmin, xmax  ! ... added by Carlo Cavazzoni
      END TYPE spline_data

!-----------------------------------------------------------------------------!

    CONTAINS

!-----------------------------------------------------------------------------!

      SUBROUTINE nullify_spline( spl )
        TYPE (spline_data), INTENT (INOUT) :: spl
        NULLIFY( spl%x )
        NULLIFY( spl%y )
        NULLIFY( spl%y2 )
        spl%n = 0
        spl%pos = 0
        spl%h = 0.0d0
        spl%invh = 0.0d0
        spl%h26 = 0.0d0
        spl%h16 = 0.0d0
        spl%xmin = 0.0d0
        spl%xmax = 0.0d0
        RETURN
      END SUBROUTINE nullify_spline

      SUBROUTINE allocate_spline( spl, nn, xmin, xmax )

        IMPLICIT NONE

        TYPE (spline_data), INTENT (INOUT) :: spl
        INTEGER, INTENT (IN) :: nn
        REAL(DP), INTENT (IN), OPTIONAL :: xmin, xmax

        INTEGER err

        IF( PRESENT( xmin ) .AND. .NOT. PRESENT( xmax ) ) &
          CALL errore(' allocate_spline ', ' wrong number of arguments ', 1 )

        spl%n = nn

        IF ( associated(spl%x) ) THEN
          DEALLOCATE (spl%x,STAT=err)
          IF (err/=0) CALL errore(' allocate_spline ','could not deallocate spl%x',1)
          NULLIFY (spl%x)
        END IF

        ! note that spl%x is not allocated if we use a regular x grid

        IF( PRESENT( xmin ) .AND. PRESENT( xmax ) ) THEN
          IF( xmin >= xmax ) &
            CALL errore(' allocate_spline ', ' wrong interval ', 1)
          spl%xmin = xmin
          spl%xmax = xmax
          spl%h    = ( xmax - xmin ) / DBLE( nn - 1 )
          spl%invh = 1.0d0 / spl%h
        ELSE
          spl%xmin = 0
          spl%xmax = 0
          ALLOCATE (spl%x(1:nn),STAT=err)
          IF (err/=0) CALL errore(' allocate_spline ','could not allocate spl%x',1)
        END IF

        IF (associated(spl%y)) THEN
          DEALLOCATE (spl%y,STAT=err)
          IF (err/=0) CALL errore(' allocate_spline ','could not deallocate spl%y',1)
          NULLIFY (spl%y)
        END IF

        ALLOCATE (spl%y(1:nn),STAT=err)
        IF (err/=0) CALL errore(' allocate_spline ','could not allocate spl%y',1)

        IF (associated(spl%y2)) THEN
          DEALLOCATE (spl%y2,STAT=err)
          IF (err/=0) CALL errore(' allocate_spline ','could not deallocate spl%y2',1)
          NULLIFY (spl%y2)
        END IF

        ALLOCATE (spl%y2(1:nn),STAT=err)
        IF (err/=0) CALL errore(' allocate_spline ','could not allocate spl%y2',1)

        RETURN
      END SUBROUTINE allocate_spline


!-----------------------------------------------------------------------

      SUBROUTINE init_spline( spl, endpt, y1a, y1b )

!   endpt: 's': regular spacing
!          'l': left; 'r': right; 'b': both = specify 1-deriv for each endpoints

        IMPLICIT NONE
        TYPE (spline_data), INTENT (INOUT) :: spl
        CHARACTER (len=*), INTENT (IN), OPTIONAL :: endpt
        REAL (DP), INTENT (IN), OPTIONAL :: y1a, y1b  ! end point derivative
        INTEGER :: err, i, k, n
        REAL (DP) :: p, qn, sig, un, y1l, y1r, dyp, dym, dxp, dxm, dxpm
        REAL (DP), POINTER :: ww(:)
        CHARACTER (len=8) :: ep
        LOGICAL :: reg, lep, rep

        !  shortcat for regular mesh without table of x values

        IF( .NOT. ASSOCIATED( spl%x ) ) THEN
          CALL splinedx( spl%xmin, spl%xmax, spl%y(:), spl%n, 0.0d0, 0.0d0, spl%y2(:) )
          RETURN
        END IF

        !  Find out if y first derivative is given at endpoints

        IF ( .NOT. present(endpt) ) THEN
          ep = ' '
        ELSE
          ep = endpt
        END IF
        reg = ( scan(ep,'sS') > 0 )
        lep = ( scan(ep,'lL') > 0 ) .OR. ( scan(ep,'bB') > 0 )
        rep = ( scan(ep,'rR') > 0 ) .OR. ( scan(ep,'bB') > 0 )

        !  check input parameter consistency

        IF ( ( lep .OR. rep ) .AND. .NOT. present(y1a) ) &
          CALL errore( 'init_spline', 'first deriv. at end-point missing', 1 )
        IF ( lep .AND. rep .AND. .NOT. present(y1b) ) &
          CALL errore( 'init_spline', 'first deriv. at end-point missing', 1 )

        !  define endpoints derivative

        IF ( lep ) y1l = y1a
        IF ( rep .AND. .NOT. lep ) y1r = y1a
        IF ( rep .AND. lep ) y1r = y1b

        spl%pos = 1
        ALLOCATE ( ww( 1 : spl%n ), STAT = err )
        IF (err/=0) CALL errore('init_spline','could not allocate ww',1)

        n = spl % n

        IF ( lep ) THEN
          spl%y2(1) = -0.5d0
          dxp  = spl%x(2) - spl%x(1)
          dyp  = spl%y(2) - spl%y(1)
          ww(1) = ( 3.0d0 / dxp ) * ( dyp / dxp - y1l )
        ELSE
          spl%y2(1) = 0
          ww(1) = 0.d0
        END IF

        DO i = 2, n - 1

          dxp  = spl%x(i+1) - spl%x(i)
          dxm  = spl%x(i)   - spl%x(i-1)
          dxpm = spl%x(i+1) - spl%x(i-1)

          sig = dxm / dxpm 
          p = sig * spl%y2(i-1) + 2.0d0
          spl%y2(i) = ( sig - 1.0d0 ) / p

          dyp  = spl%y(i+1) - spl%y(i)
          dym  = spl%y(i)   - spl%y(i-1)

          ww(i) = ( 6.0d0 * ( dyp / dxp  - dym / dxm ) / dxpm  - sig * ww(i-1) ) / p

        END DO

        IF ( rep ) THEN
          qn = 0.5d0
          dxm = spl%x(n) - spl%x(n-1)
          dym = spl%y(n) - spl%y(n-1)
          un = ( 3.0d0 / dxm ) * ( y1r - dym / dxm )
        ELSE
          qn = 0
          un = 0
        END IF

        spl % y2(n) = ( un - qn * ww(n-1) ) / ( qn * spl%y2(n-1) + 1.0d0 )

        DO k = n - 1, 1, -1
          spl % y2(k) = spl%y2(k) * spl%y2(k+1) + ww(k)
        END DO

        DEALLOCATE ( ww, STAT = err )
        IF (err/=0) CALL errore('init_spline','could not deallocate ww',1)

        IF ( reg ) THEN
          spl%h = ( spl%x(n) - spl%x(1) ) / ( n - 1.0d0 )
          spl%h16 = spl%h / 6.0d0
          spl%h26 = spl%h**2 / 6.0d0
          spl%invh = 1 / spl%h
        ELSE
          spl%h = 0.0d0
          spl%invh = 0.0d0
        END IF

        RETURN
      END SUBROUTINE init_spline

!-----------------------------------------------------------------------

      FUNCTION interv( spl, xx )

        IMPLICIT NONE

        TYPE (spline_data), INTENT (IN) :: spl
        REAL (DP), INTENT (IN) :: xx
        INTEGER :: interv
        INTEGER :: khi, klo, i, p, n, k

        !   if we have a regular mesh use a quick position search

        IF ( spl%h /= 0 ) THEN
          i = ( xx - spl%x(1) ) * spl%invh + 1
          IF ( i < 1 .OR. i > spl%n ) &
            CALL errore('interv', 'illegal x-value passed to interv',1)
          interv = i
          RETURN
        END IF

        p = spl%pos
        IF ( p >= spl%n .OR. p <= 1 ) p = spl%n / 2
        i = 0
        n = spl%n

        !   check if interval is close to previous interval

        IF ( xx < spl%x(p+1) ) THEN
          IF ( xx >= spl%x(p) ) THEN
            i = spl%pos
          ELSE IF ( p > 1 .AND. xx >= spl%x(p-1) ) THEN
            i = p - 1
          ELSE
            klo = 1
            khi = p + 1
          END IF
        ELSE IF ( (p + 2) <= n .AND. xx < spl%x(p+2) ) THEN
          i = p + 1
        ELSE
          klo = p + 1
          khi = n
        END IF

        !   perform binary search

        IF ( i == 0 ) THEN 
          IF ( xx < spl%x(1) .OR. xx > spl%x(n) ) &
            CALL errore('interv', 'xx value out of spline-range',1)
          DO WHILE ( (khi - klo) > 1 )
            k = ( khi + klo ) / 2
            IF ( spl%x(k) > xx ) THEN
              khi = k
            ELSE
              klo = k
            END IF
          END DO
          i = klo
        END IF

        interv = i
        RETURN
      END FUNCTION interv


!-----------------------------------------------------------------------
      FUNCTION spline( spl, xx, y1 )

        IMPLICIT NONE

        TYPE (spline_data), INTENT (INOUT) :: spl
        REAL (DP), INTENT (IN) :: xx
        REAL (DP), INTENT (OUT), OPTIONAL :: y1
        REAL (DP) :: spline

        INTEGER :: khi, klo
        REAL (DP) :: a, b, h, invh, ylo, yhi, y2lo, y2hi, a3ma, b3mb

        !  shortcat for regular mesh without table of x values

        IF( .NOT. ASSOCIATED( spl%x ) ) THEN
          IF( PRESENT( y1 ) ) &
            CALL errore(' spline ', ' y1 without x table not implemented ', 1 )
          CALL splintdx( spl%xmin, spl%xmax, spl%y, spl%y2, spl%n, xx, a )
          spline = a
          RETURN
        END IF

        spl%pos = interv( spl, xx )
        klo = spl%pos
        khi = spl%pos + 1

        IF ( spl%h /= 0 ) THEN
          h    = spl%h
          invh = spl%invh
        ELSE
          h    = spl%x( khi ) - spl%x( klo )
          invh = spl%invh
          IF ( h == 0.0d0 ) &
            CALL errore('spline','bad spl%x input',1)
        END IF

        a = ( spl%x( khi ) - xx ) * invh
        b = 1 - a
        a3ma = a**3 - a
        b3mb = b**3 - b
        ylo  = spl%y( klo )
        yhi  = spl%y( khi )
        y2lo = spl%y2( klo )
        y2hi = spl%y2( khi )
        spline = a * ylo + b * yhi + ( a3ma * y2lo + b3mb * y2hi ) * ( h**2 ) / 6.0d0

        IF ( present( y1 ) ) then
          y1 = ( yhi - ylo ) * invh + &
          ( ( 1.0d0 - 3 * a**2 ) * y2lo + ( 3 * b**2 - 1.0d0 ) * y2hi ) * h / 6.0d0
        END IF

        RETURN
      END FUNCTION spline

!-----------------------------------------------------------------------

      FUNCTION splineh(spl,xx,y1)
        IMPLICIT NONE
        TYPE (spline_data), INTENT (IN) :: spl
        REAL (DP), INTENT (IN) :: xx
        REAL (DP), INTENT (OUT) :: y1
        REAL (DP) :: splineh

        INTEGER :: khi, klo, i
        REAL (DP) :: a, b, h, invh, t, ylo, yhi, y2lo, y2hi, d, d0

! fast spline for pair potentials without checks
        h = spl%h
        invh = spl%invh
        d=xx-spl%x(1); i=INT(d*spl%invh); d0=DBLE(i)*h; i=i+1
        i = (xx-spl%x(1))*invh + 1

        a = (spl%x(i+1)-xx)*invh
        b = 1 - a
        t = -a*b
!    b=(d-d0)*invh; a=1-b; t=-a*b
        ylo = spl%y(i)
        yhi = spl%y(i+1)
        y2lo = spl%y2(i)
        y2hi = spl%y2(i+1)
        splineh = a*ylo + b*yhi + ((a+1)*y2lo+(b+1)*y2hi)*t*spl%h26
        y1 = (yhi-ylo)*invh + ((1.d0-3*a*a)*y2lo+(3*b*b-1.d0)*y2hi)*spl%h16

      END FUNCTION splineh
!-----------------------------------------------------------------------
      FUNCTION spline_1(spl,xx)
        IMPLICIT NONE
        TYPE (spline_data), INTENT (INOUT) :: spl
        REAL (DP), INTENT (IN) :: xx
        REAL (DP) :: spline_1

        INTEGER :: khi, klo
        REAL (DP) :: a, b, h

        spl%pos = interv(spl,xx)
        klo = spl%pos
        khi = spl%pos + 1

        h = spl%x(khi) - spl%x(klo)
        IF (h==0.d0) CALL errore('spline','bad spl%x input',1)
        a = (spl%x(khi)-xx)/h
        b = 1 - a
        spline_1 = (spl%y(khi)-spl%y(klo))/h + ((1.d0-3*a**2)*spl%y2(klo)+(3*b** &
          2-1.d0)*spl%y2(khi))*h/6.d0

        RETURN
      END FUNCTION spline_1



!-----------------------------------------------------------------------
      FUNCTION stamm(spl,p,x)
        IMPLICIT NONE
        TYPE (spline_data), INTENT (IN) :: spl
        INTEGER, INTENT (IN) :: p
        REAL (DP), INTENT (IN) :: x
        REAL (DP) :: stamm
        REAL (DP) :: a, b, aa, bb, h

        h = spl%x(p+1) - spl%x(p)
        b = (x-spl%x(p))/h
        a = 1 - b
        aa = a**2
        bb = b**2
        stamm = 0.5d0*h*(bb*spl%y(p+1)-aa*spl%y(p)) + h**3/12.d0*(aa*(1-0.5d0*aa)* &
          spl%y2(p)-bb*(1-0.5d0*bb)*spl%y2(p+1))

        RETURN
      END FUNCTION stamm



!-----------------------------------------------------------------------
      FUNCTION spline_int(spl,x0,x1)
        IMPLICIT NONE
        TYPE (spline_data), INTENT (INOUT) :: spl
        REAL (DP), INTENT (IN) :: x0, x1
        REAL (DP) :: spline_int

        INTEGER :: j, pa, pb
        REAL (DP) :: h, vorz, xa, xb, i1, i2

        vorz = 1
        xa = min(x0,x1)
        xb = max(x0,x1)
        IF (x0>x1) vorz = -1
        IF (xa<spl%x(1) .OR. xb>spl%x(spl%n)) CALL errore('spline_int', &
          'illegal integration range',1)

        pa = interv(spl,xa)
        pb = interv(spl,xb)

        IF (pa==pb) THEN
          spline_int = vorz*(stamm(spl,pa,xb)-stamm(spl,pa,xa))
          RETURN
        END IF

        i1 = 0
        i2 = 0
        DO j = pa + 1, pb - 1
          h = spl%x(j+1) - spl%x(j)
          i1 = i1 + h*(spl%y(j)+spl%y(j+1))
          i2 = i2 + h**3*(spl%y2(j)+spl%y2(j+1))
        END DO
        h = spl%x(pa+1) - spl%x(pa)
        i1 = i1 + h*spl%y(pa+1)
        i2 = i2 + h**3*spl%y2(pa+1)
        h = spl%x(pb+1) - spl%x(pb)
        i1 = i1 + h*spl%y(pb)
        i2 = i2 + h**3*spl%y2(pb)

        spline_int = vorz*(i1/2.-i2/24.d0+stamm(spl,pb,xb)-stamm(spl,pa,xa))

        RETURN
      END FUNCTION spline_int

!-----------------------------------------------------------------------
      SUBROUTINE kill_spline(spl,what)
!   deallocate splines
!   what=='a' or not present: deallocate all (spl%x, spl%y, spl%y2)
!   what=='d': deallocate only data (spl%x, spl%y)
!   what=='2': deallocate only table of 2. derivatives (spl%y2)
        IMPLICIT NONE
        TYPE (spline_data), INTENT (INOUT) :: spl
        CHARACTER, INTENT (IN), OPTIONAL :: what
        CHARACTER :: w
        INTEGER :: err

        w = 'a'
        IF (present(what)) w = what
        SELECT CASE (w)
        CASE ('d','D')
          IF (associated(spl%x)) THEN
            DEALLOCATE (spl%x,STAT=err)
            IF (err/=0) CALL errore('kill_spline', 'could not deallocate spl%x',1)
            NULLIFY (spl%x)
          END IF
          IF (associated(spl%y)) THEN
            DEALLOCATE (spl%y,STAT=err)
            IF (err/=0) CALL errore('kill_spline', 'could not deallocate spl%y',1)
            NULLIFY (spl%y)
          END IF
        CASE ('2')
          IF (associated(spl%y2)) THEN
            DEALLOCATE (spl%y2,STAT=err)
            IF (err/=0) CALL errore('kill_spline', 'could not deallocate spl%y2',1)
            NULLIFY (spl%y2)
          END IF
        CASE ('a','A')
          IF (associated(spl%x)) THEN
            DEALLOCATE (spl%x,STAT=err)
            IF (err/=0) CALL errore('kill_spline', 'could not deallocate spl%x',1)
            NULLIFY (spl%x)
          END IF
          IF (associated(spl%y)) THEN
            DEALLOCATE (spl%y,STAT=err)
            IF (err/=0) CALL errore('kill_spline', 'could not deallocate spl%y',1)
            NULLIFY (spl%y)
          END IF
          IF (associated(spl%y2)) THEN
            DEALLOCATE (spl%y2,STAT=err)
            IF (err/=0) CALL errore('kill_spline', 'could not deallocate spl%y2',1)
            NULLIFY (spl%y2)
          END IF
        END SELECT

        RETURN
      END SUBROUTINE kill_spline

!=-----------------------------------------------------------------------=!
! Subroutines: splinedx, splintdx
! added for compatibility with SISSA code 
! Carlo Cavazzoni 15-03-2000
!=-----------------------------------------------------------------------=!

      SUBROUTINE splinedx(xmin,xmax,y,n,yp1,ypn,y2)
        USE kinds
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: n
        REAL(DP),  INTENT(IN)  :: yp1,ypn,xmin,xmax,y(:)
        REAL(DP),  INTENT(OUT) :: y2(:)
        INTEGER :: i, k
        REAL(DP)  :: p, qn, sig, un, dx
        REAL(DP)  :: u(n)

        dx = (xmax-xmin)/DBLE(n-1)
        if ( yp1 .gt. 0.99d30 ) then
          y2(1)=0.d0
          u(1)=0.0d0
        else
          y2(1)=-0.5d0
          u(1)=(3.d0/dx) * ( (y(2)-y(1)) / dx - yp1 )
        endif
        do i=2,n-1
          sig=0.5d0
          p=sig*y2(i-1)+2.d0
          y2(i)=(sig-1.d0)/p
          u(i) = (6.0d0 * ( (y(i+1)-y(i))/ dx - (y(i)-y(i-1))/ dx ) &
                        / (2.0d0*dx) - sig * u(i-1) ) / p
        end do 
        if ( ypn .gt. 0.99d30 ) then
          qn=0.d0
          un=0.d0
        else
          qn=0.5d0
          un= ( 3.d0 / dx ) * ( ypn - (y(n)-y(n-1)) / dx )
        endif
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
        do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
        end do
        return
      END SUBROUTINE splinedx


      SUBROUTINE splintdx(xmin,xmax,ya,y2a,n,x,y)
        USE kinds
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: n
        REAL(DP),  INTENT(IN)  :: x,xmin,xmax,ya(:),y2a(:)
        REAL(DP),  INTENT(OUT) :: y
        INTEGER :: khi,klo
        REAL(DP)  :: a,b,h,dx,xhi,xlo
        dx  = (xmax-xmin)/DBLE(n-1)
        klo = INT(x/dx)
        khi = klo+1
        IF(klo.LT.1) THEN
          CALL errore(' splintdx ',' klo less than one ',klo)
        END IF
        IF(khi.GT.n) THEN
          CALL errore(' splintdx ',' khi grether than N ',khi)
        END IF
        xlo  = xmin + DBLE(klo-1) * dx
        xhi  = xmin + DBLE(khi-1) * dx
       
        a   = (xhi-x)/dx
        b   = (x-xlo)/dx
        y   = a*ya(klo) + b*ya(khi) + &
              ( (a*a*a-a)*y2a(klo) + (b*b*b-b)*y2a(khi) ) * (dx*dx)/6.0d0
        RETURN
      END SUBROUTINE splintdx

!-----------------------------------------------------------------------

      SUBROUTINE nr_spline( x, y, n, yp1, ypn, y2 )
      INTEGER :: n
      REAL(DP) :: yp1, ypn, x(n), y(n), y2(n)
      INTEGER :: i, k
      REAL(DP) :: p, qn, sig, un
      REAL(DP) :: u( n )
      if ( yp1 .gt. 0.99d30 ) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) / & 
             (x(i)-x(i-1))) / (x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if ( ypn .gt. 0.99d30 ) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      END SUBROUTINE nr_spline


      SUBROUTINE nr_splint( xa, ya, y2a, n, x, y )
      INTEGER :: n
      REAL(DP) :: x,y,xa(n),y2a(n),ya(n)
      INTEGER :: k,khi,klo
      REAL(DP) :: a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) &
        CALL errore(' nr_splint ', 'bad xa input in splint' , 1 )
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* &
        (h**2)/6.d0
      return
      END SUBROUTINE nr_splint


      SUBROUTINE nr_splie2( x1a, x2a, ya, m, n, y2a )
      INTEGER :: m, n
      REAL(DP) :: x1a(m), x2a(n), y2a(m,n), ya(m,n)
      INTEGER :: j,k
      REAL(DP) :: y2tmp(n), ytmp(n)
      do j = 1, m
        do k = 1, n
          ytmp(k) = ya(j,k)
        end do
        call nr_spline( x2a, ytmp, n, 1.d30, 1.d30, y2tmp )
        do k = 1, n
          y2a(j,k) = y2tmp(k)
        end do
      end do
      return
      END SUBROUTINE nr_splie2


      SUBROUTINE nr_splin2( x1a, x2a, ya, y2a, m, n, x1, x2, y )
      INTEGER :: m, n
      REAL(DP) :: x1, x2, y, x1a(m), x2a(n), y2a(m,n), ya(m,n)
      INTEGER :: j, k
      REAL(DP) :: y2tmp( MAX(n,m) ), ytmp( n ), yytmp( MAX(n,m) )
      do j = 1, m
        do k = 1, n
          ytmp(k)  = ya(j,k)
          y2tmp(k) = y2a(j,k)
        end do
        call nr_splint( x2a, ytmp, y2tmp, n, x2, yytmp(j) )
      end do
      call nr_spline( x1a, yytmp, m, 1.d30, 1.d30, y2tmp )
      call nr_splint( x1a, yytmp, y2tmp, m, x1, y )
      return
      END SUBROUTINE nr_splin2

!-----------------------------------------------------------------------

    END MODULE splines
