SUBROUTINE getvofr( np_in_sp_me, np_in_sp, &
        hcub, rho, v, vnew, vold, vold2,tSelf, d_cur, d_n, d_o, d_o2, cgstep)
    !=======================================================================================
    ! Code Version 1.0 (Princeton University, September 2014)
    !=======================================================================================
    !===============================================================
    ! Given charge density, get the potential using multipole expansion
    ! for boundary region and solving poisson equation for inside box
    ! Adapted from PARSEC by Lingzhu Kong,  http://parsec.ices.utexas.edu/
    !===============================================================
    !
    USE kinds,                   ONLY  :  DP
    USE io_global,               ONLY  : stdout
    USE exx_module,              ONLY  :  coeke, nord2, lmax, n_exx
    USE exx_module,              ONLY  :  thdtood_in_sp 
    USE mp_global,               ONLY  :  me_image
    USE parallel_include
    !
    IMPLICIT NONE
    !
    ! test if we are in the case of self interaction
    LOGICAL tSelf
    ! pair distances and extrapolation coefficients
    REAL(DP) d_cur, d_n, d_o, d_o2, Cex(3)
    INTEGER  np_in_sp_me, np_in_sp
    REAL(DP) hcub
    REAL(DP) rho(np_in_sp)
    !OUTPUT
    REAL(DP) v(np_in_sp_me)
    !INOUT
    REAL(DP) vnew(np_in_sp), vold(np_in_sp), vold2(np_in_sp)
    !
    REAL(DP), ALLOCATABLE  :: v_in_sp(:)
    COMPLEX(DP) qlm(0:lmax, 0:lmax)
    INTEGER  i, j, k, ii, jj, kk, ir, ip, irg, ifg, cgstep
    !========================================================================
    !-----------------------------------------------------------------------
    !First, we calculate the potential outside the inner sphere, which will
    !also give the boundary values for potential inside the sphere.
    !-----------------------------------------------------------------------
    !
    CALL start_clock('getvofr_qlm')
    CALL getqlm(np_in_sp, hcub, rho, qlm)
    CALL stop_clock('getvofr_qlm')
    ALLOCATE( v_in_sp(1:np_in_sp_me) )
    !
    v_in_sp = 0.d0
    Cex = 0.0D0
    CALL start_clock('getvofr_bound')
    CALL exx_boundaryv( np_in_sp_me, np_in_sp, v_in_sp, qlm)
    CALL stop_clock('getvofr_bound')
    !
    !$omp parallel do 
    DO ir = 1+np_in_sp, np_in_sp_me
      v(ir) = v_in_sp(ir)
    END DO
    !$omp end parallel do 
    !-------------------------------------------------------------------------------
    !                   NEXT, THE POTENTIAL INSIDE THE SPHERE
    !-------------------------------------------------------------------------------
    !
    IF (tSelf) THEN
      ! Selfv calculation
      IF (n_exx .EQ. 2) THEN
        !
        !$omp parallel do 
        DO ir = 1, np_in_sp
          !
          v_in_sp(ir) = 2.D0*vnew(ir) - vold(ir)
          !
        ENDDO
        !$omp end parallel do 
        !
      ELSE IF (n_exx .GT. 2) THEN
        !
        !$omp parallel do 
        DO ir = 1, np_in_sp
          !
          v_in_sp(ir) = 3.D0*vnew(ir) - 3.D0*vold(ir) + vold2(ir)
          !
        END DO
        !$omp end parallel do 
        !
      END IF
      !
    ELSE
      ! Pair v calculation
      IF (n_exx .EQ. 2) THEN
        ! pair distances extrapolation 
        IF (ABS(d_n - d_o).GT.1.D-8) THEN
          !
          Cex (1) = (d_cur - d_o)/(d_n - d_o)
          Cex (2) = (d_cur - d_n)/(d_o - d_n)
          IF ((ABS(Cex(1)).GT.3.D0).OR.(ABS(Cex(2)).GT.2.0D0)) THEN 
            !
            Cex (1) = 2.0D0
            Cex (2) = -1.D0
            !
          END IF
          !
        ELSE 
          !
          Cex(1) = 2.0D0
          Cex(2) = -1.D0
          !
        END IF
        !$omp parallel do 
        DO ir = 1, np_in_sp
          ! Hybrid time and pair-dist extrapolation
          v_in_sp(ir) = 0.5D0*((2.0D0+Cex(1))*vnew(ir) &
              + (Cex(2)-1.0D0)*vold(ir))
          !
        END DO
        !$omp end parallel do 
        !
      ELSE IF(n_exx .Gt. 2)THEN
        ! pair distances extrapolation 
        IF ((ABS(d_n - d_o).GT.1.D-8).OR.(ABS(d_n - d_o2).GT.1.D-8).OR.(abs(d_o - d_o2).GT.1.D-8)) THEN
          !
          Cex(1) = (d_cur-d_o2)*(d_cur-d_o)/(d_n-d_o2)/(d_n-d_o)
          Cex(2) = (d_cur-d_o2)*(d_cur-d_n)/(d_o-d_o2)/(d_o-d_n)
          Cex(3) = (d_cur-d_o)*(d_cur-d_n)/(d_o2-d_o)/(d_o2-d_n)
          !
          IF ((ABS(Cex(1)).GT.5.D0).OR.(ABS(Cex(2)).GT.5.0D0).OR.(ABS(Cex(3)).GT.3.D0)) THEN 
            !
            Cex (1) =  3.0D0
            Cex (2) = -3.0D0
            Cex (3) =  1.0D0
            !
          END IF
          !
        ELSE 
          !
          Cex (1) =  3.0D0
          Cex (2) = -3.0D0
          Cex (3) =  1.0D0
          !
        END IF
        !
        !$omp parallel do 
        DO ir = 1, np_in_sp
          ! Hybrid time and pair-dist extrapolation
          v_in_sp(ir) = 0.5D0*((3.D0+Cex(1))*vnew(ir) &
              +(-3.D0+Cex(2))*vold(ir)+(1.D0+Cex(3))*vold2(ir))
          !
        END DO
        !$omp end parallel do 
        !
      END IF
      !
    END IF
    !
    !!!!! This part is original !!!!!!!!!!
    CALL start_clock('getvofr_geterho')
    CALL geterho(np_in_sp_me, np_in_sp,rho, v_in_sp)
    CALL stop_clock('getvofr_geterho')
    CALL start_clock('getvofr_hpotcg')
    CALL hpotcg(np_in_sp_me, np_in_sp, rho, v_in_sp,.TRUE.,cgstep)
    CALL stop_clock('getvofr_hpotcg')
    !
    IF ( n_exx .EQ. 1) THEN
      !
      !$omp parallel do 
      DO ir = 1, np_in_sp
        !
        vnew(ir) = v_in_sp(ir)
        vold(ir) = vnew(ir)
        !
      END DO
      !$omp end parallel do 
      !
      IF (.NOT.tSelf) THEN
        !
        d_n = d_cur
        d_o = d_n
        !
      END IF
      !
    ELSE IF ( n_exx .EQ. 2) THEN
      !
      !$omp parallel do 
      DO ir = 1, np_in_sp
        !
        vold(ir) = vnew(ir)
        vnew(ir) = v_in_sp(ir)
        !
      END DO
      !$omp end parallel do 
      !
      IF (.NOT.tSelf) THEN
        !
        d_o = d_n
        d_n = d_cur
        !
      END IF
      !
    ELSE 
      !
      !$omp parallel do 
      DO ir = 1, np_in_sp
        !
        vold2(ir) = vold(ir)
        vold(ir) = vnew(ir)
        vnew(ir) = v_in_sp(ir)
        !
      END DO
      !$omp end parallel do 
      !
      IF (.NOT.tSelf) THEN
        !
        d_o2 = d_o
        d_o  = d_n
        d_n  = d_cur
        !
      END IF
      !
    END IF
    !
    !$omp parallel do 
    DO ir = 1, np_in_sp
      !
      v(ir) = v_in_sp(ir)
      !
    END DO
    !$omp end parallel do 
    !
    DEALLOCATE( v_in_sp)
    !
    RETURN
    !
END SUBROUTINE getvofr
!===============================================================================

!===============================================================================
SUBROUTINE getqlm(np_in_sp, hcub, rho, qlm)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  lpole=>lmax
    USE exx_module,       ONLY  :  xx_in_sp, yy_in_sp,  zz_in_sp, clm
    !
    IMPLICIT NONE
    !
    INTEGER    np_in_sp
    REAL(DP)   hcub, rho(1:np_in_sp)
    COMPLEX(DP) qlm(0:lpole, 0:lpole)
    !
    REAL(DP)   rrho,xx,yy,zz,xx2,yy2,zz2,xy,r2,x,y
    REAL(DP)   rinv(0:lpole),r(0:lpole)
    !
    REAL(DP)    plm(0:lpole, 0:lpole)  !temporary storage of the associated Legendre polynom
    COMPLEX(DP) cxy(1:lpole)      !coefficient array: e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
    !
    REAL(DP)   hcub2, zero, one, two , tmp1, tmp2
    PARAMETER (zero = 1.0E-10, one = 1.d0, two = 2.d0 )
    INTEGER    i,j,l,m
    !------------------------------------------------------------------------------
    !
    hcub2 = two*hcub
    qlm(:,:) = 0.d0
    !
    !$omp parallel do private(r,plm,cxy,rrho,tmp1,tmp2,xx,yy,zz,xx2,yy2,zz2,r2,rinv,xy,x,y), reduction(+:qlm)
    DO j = 1,np_in_sp
      rrho = rho(j)
      !
      xx   = xx_in_sp(j)
      yy   = yy_in_sp(j)
      zz   = zz_in_sp(j)
      !
      xx2  = xx*xx
      yy2  = yy*yy
      zz2  = zz*zz
      r2   = xx2 + yy2 + zz2
      !
      r(1) = dsqrt(r2)
      qlm(0,0) = qlm(0,0) + rrho
      !
      IF (r(1) .GT. zero) THEN
        DO l = 2, lpole
          r(l) = r(l-1)*r(1)           !r(l) = r^l (higher powers of r)
        END DO
        !
        rinv(1) = one/r(1)
        xy = DSQRT(xx2+yy2)
        !
        x = zz*rinv(1)                  ! x = cos(theta_j)
        y = xy*rinv(1)                  ! y = sin(theta_j)
        !
        CALL setplm(x, y, lpole, plm)   ! associate Legendre polynomials in plm(l,m)
        !
        DO l = 1, lpole
          qlm(l,0) = qlm(l,0) + rrho*plm(l,0)*r(l)   !qlm(l,m=0)
        END DO
        !
        IF (xy .GT. zero) THEN
          !
          cxy(1) = CMPLX(xx,yy,KIND=dp)
          cxy(1) = cxy(1)/xy
          DO m = 2, lpole
            cxy(m) = cxy(m-1)*cxy(1)  !cxy(m) = exp(i*m*phi_j) = (cos(phi_j) + i*sin(phi_j))^m
          END DO                         !       = ((xx + i*yy)/(sqrt(xx^2 + yy^2))^m
          !
          DO l = 1, lpole
            tmp1 = rrho*r(l)
            DO m = 1, l
              tmp2     = plm(l,m)*tmp1
              qlm(l,m) = qlm(l,m) + cxy(m)*tmp2      !qlm(l, m>0)
            END DO
          END DO
          !
        END IF
      END IF
    END DO
    !$omp end parallel do
    !
    DO l = 0, lpole
      qlm(l,0) = qlm(l,0)*clm(l,0)*hcub
    END DO
    !
    DO l = 1, lpole
      DO m = 1, l, 1
        qlm(l,m) = qlm(l,m)*clm(l,m)*hcub2
      END DO
    END DO
    !
    RETURN
END SUBROUTINE getqlm
!===============================================================================

!===============================================================================
SUBROUTINE exx_boundaryv(np_in_sp_me, np_in_sp,v_in_sp,qlm)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  lpole=>lmax 
    USE exx_module,       ONLY  :  xx_in_sp, yy_in_sp,  zz_in_sp, clm
    !
    IMPLICIT  NONE
    !
    INTEGER   np_in_sp_me, np_in_sp
    REAL(DP)  v_in_sp( 1:np_in_sp_me )
    !
    ! For each grid point:
    REAL(dp)   plm(0:lpole, 0:lpole)
    COMPLEX*16 qlm(0:lpole, 0:lpole)
    COMPLEX*16 cxy(1:lpole)   ! e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
    COMPLEX*16 cpole, qlmc
    !
    REAL(dp)   rinv(0:lpole), r(0:lpole)
    REAL(dp)   xh,yh,zh, xx2, yy2, zz2, xy, x, y, r2, zero, one, plmr
    INTEGER    i, l, m
    !
    zero = 0.d0
    one  = 1.d0
    !
    !$omp parallel do private(r,plm,cxy,xh,yh,zh,xx2,yy2,zz2,r2,rinv,xy,x,y,plmr,cpole,qlmc)
    DO i = 1 + np_in_sp, np_in_sp_me
      !
      xh = xx_in_sp(i)
      yh = yy_in_sp(i)
      zh = zz_in_sp(i)
      !
      xx2 = xh*xh
      yy2 = yh*yh
      zz2 = zh*zh
      !
      r2 = xx2 + yy2 + zz2
      r(1) = dsqrt(r2)
      IF (r(1) .LT. 0.000001) PRINT *, 'i =',i, r(1)
      rinv(0) = one/r(1)
      DO l = 1,lpole
        rinv(l) = rinv(l-1)*rinv(0)  !rinv(l) = 1/r^(l+1).
      END DO
      !
      xy= DSQRT(xx2+yy2)
      x = zh*rinv(0)           ! x=cos(theta)
      y = xy*rinv(0)           ! y=sin(theta)
      !
      ! evaluate associate Legendre polynomials for x and y with l from 0 to 
      ! lpole and m from 0 to l. store each in plm(l,m)
      CALL setplm(x, y, lpole, plm)
      !
      cpole = zero
      !
      DO l = 0,lpole
        plmr = plm(l,0)*rinv(l)
        cpole = cpole + qlm(l,0)*plmr    ! m=0 terms
      END DO
      !
      IF (xy .GT. zero) THEN
        cxy(1) = CMPLX(xh,-yh,KIND=dp)
        cxy(1) = cxy(1)/xy
        DO m=2, lpole, 1
          cxy(m) = cxy(m-1)*cxy(1)
        END DO
        !
        DO l = 1, lpole
          DO m = 1, l, 1                          ! m>0 terms
            plmr = plm(l,m)*rinv(l)
            qlmc = qlm(l,m)*cxy(m)
            cpole = cpole + plmr * qlmc
          END DO
        END DO
      END IF
      !
      v_in_sp(i) = REAL(cpole)
      !
    END DO 
    !$omp end parallel do
    !
    RETURN
    !
END SUBROUTINE exx_boundaryv
! ===========================================================================

! ===========================================================================
SUBROUTINE geterho(np_in_sp_me, np_in_sp, rho, v_in_sp)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  nord2,coeke
    USE exx_module,       ONLY  :  odtothd_in_sp, thdtood_in_sp
    USE constants,        ONLY  :  eps12
    !
    IMPLICIT NONE
    !
    INTEGER np_in_sp, np_in_sp_me
    REAL(DP) rho(np_in_sp), v_in_sp( 1:np_in_sp_me )
    !
    INTEGER :: i, j, k, ip, jp, ish, ipp, ipm, jpp, jpm, kpp, kpm
    INTEGER :: ipjp, ipjm, imjp, imjm
    !
    !$omp parallel do private(i,j,k,ipp,ipm,jpp,jpm,kpp,kpm) 
    DO ip = 1, np_in_sp
      !
      !(i,j,k) is within the first sphere
      i = odtothd_in_sp(1,ip)
      j = odtothd_in_sp(2,ip)
      k = odtothd_in_sp(3,ip)
      !
      DO ish  = 1, nord2
        !
        ipp = thdtood_in_sp( i+ish, j,     k     )
        ipm = thdtood_in_sp( i-ish, j,     k     )
        jpp = thdtood_in_sp( i,     j+ish, k     )
        jpm = thdtood_in_sp( i,     j-ish, k     )
        kpp = thdtood_in_sp( i,     j,     k+ish )
        kpm = thdtood_in_sp( i,     j,     k-ish )
        !
        IF(ipp .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1,1)*v_in_sp(ipp) ! apply finite difference stencil
        IF(ipm .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1,1)*v_in_sp(ipm) ! apply finite difference stencil
        IF(jpp .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2,2)*v_in_sp(jpp) ! apply finite difference stencil
        IF(jpm .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2,2)*v_in_sp(jpm) ! apply finite difference stencil
        IF(kpp .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3,3)*v_in_sp(kpp) ! apply finite difference stencil
        IF(kpm .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3,3)*v_in_sp(kpm) ! apply finite difference stencil
        !
      END DO
      !
    END DO
    !$omp end parallel do 
    !
    !! cross derivatives
    !
    IF (ABS(coeke(1,1,2)).GT.eps12) THEN
      !$omp parallel do private(i,j,k,ipjp,ipjm,imjp,imjm)
      DO ip = 1, np_in_sp
        !(i,j,k) is within the first sphere
        i = odtothd_in_sp(1,ip)
        j = odtothd_in_sp(2,ip)
        k = odtothd_in_sp(3,ip)
        !
        DO ish  = 1, nord2
          !
          ipjp = thdtood_in_sp( i+ish, j+ish, k )
          ipjm = thdtood_in_sp( i+ish, j-ish, k )
          imjp = thdtood_in_sp( i-ish, j+ish, k )
          imjm = thdtood_in_sp( i-ish, j-ish, k )
          !
          IF(ipjp .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1,2)*v_in_sp(ipjp)
          IF(ipjm .GT. np_in_sp) rho(ip) = rho(ip) + coeke(ish,1,2)*v_in_sp(ipjm)
          IF(imjp .GT. np_in_sp) rho(ip) = rho(ip) + coeke(ish,1,2)*v_in_sp(imjp)
          IF(imjm .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1,2)*v_in_sp(imjm)
          !
        END DO
        !
      END DO
      !$omp end parallel do
    END IF
    !
    IF (ABS(coeke(1,1,3)).GT.eps12) THEN
      !$omp parallel do private(i,j,k,ipjp,ipjm,imjp,imjm)
      DO ip = 1, np_in_sp
        !(i,j,k) is within the first sphere
        i = odtothd_in_sp(1,ip)
        j = odtothd_in_sp(2,ip)
        k = odtothd_in_sp(3,ip)
        !
        DO ish  = 1, nord2
          !
          ipjp = thdtood_in_sp( i+ish, j, k+ish )
          ipjm = thdtood_in_sp( i+ish, j, k-ish )
          imjp = thdtood_in_sp( i-ish, j, k+ish )
          imjm = thdtood_in_sp( i-ish, j, k-ish )
          !
          IF(ipjp .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1,3)*v_in_sp(ipjp)
          IF(ipjm .GT. np_in_sp) rho(ip) = rho(ip) + coeke(ish,1,3)*v_in_sp(ipjm)
          IF(imjp .GT. np_in_sp) rho(ip) = rho(ip) + coeke(ish,1,3)*v_in_sp(imjp)
          IF(imjm .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1,3)*v_in_sp(imjm)
          !
        END DO
        !
      END DO
      !$omp end parallel do
    END IF
    !
    IF (ABS(coeke(1,2,3)).GT.eps12) THEN
      !$omp parallel do private(i,j,k,ipjp,ipjm,imjp,imjm)
      DO ip = 1, np_in_sp
        !(i,j,k) is within the first sphere
        i = odtothd_in_sp(1,ip)
        j = odtothd_in_sp(2,ip)
        k = odtothd_in_sp(3,ip)
        !
        DO ish  = 1, nord2
          !
          ipjp = thdtood_in_sp( i, j+ish, k+ish )
          ipjm = thdtood_in_sp( i, j+ish, k-ish )
          imjp = thdtood_in_sp( i, j-ish, k+ish )
          imjm = thdtood_in_sp( i, j-ish, k-ish )
          !
          IF(ipjp .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2,3)*v_in_sp(ipjp)
          IF(ipjm .GT. np_in_sp) rho(ip) = rho(ip) + coeke(ish,2,3)*v_in_sp(ipjm)
          IF(imjp .GT. np_in_sp) rho(ip) = rho(ip) + coeke(ish,2,3)*v_in_sp(imjp)
          IF(imjm .GT. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2,3)*v_in_sp(imjm)
          !
        END DO
        !
      END DO
      !$omp end parallel do
    END IF
    !
    RETURN
END SUBROUTINE geterho
! ==================================================================

! ==================================================================
SUBROUTINE setplm(x, y, lpole, plm)
    !
    ! this subroutine calculates all of the associated Legendre polynomials up 
    ! to a supplied lpole (maximum lpole is 9). 
    !
    IMPLICIT NONE
    !
    ! INPUT VARIABLES
    ! cos(theta),sin(theta), for a given grid point
    REAL*8 x,y
    ! order of multipole expansion
    INTEGER lpole
    ! OUTPUT VARIABLES
    ! array containing P_{lm}, the associated Legendre polynomials
    REAL*8 plm(0:lpole, 0:lpole)
    ! WORK VARIABLES:
    ! powers of x, y: xn=x^n, yn=y^n
    REAL*8 x2, x3, x4, x5, x6, x7, x8, x9
    REAL*8 y2, y3, y4, y5, y6, y7, y8, y9
    !
    plm(0,0) = 1.0
    IF (lpole .GE. 1) THEN
      plm(1,0) = x
      plm(1,1) = -y
    END IF
    IF (lpole .GE. 2) THEN
      x2 = x*x
      y2 = y*y
      plm(2,0) = 1.5*x2 - 0.5
      plm(2,1) = -3.0*x*y
      plm(2,2) = 3.0*y2
    END IF
    IF (lpole .GE. 3) THEN
      x3 = x2*x
      y3 = y2*y
      plm(3,0) = 2.5*x3 - 1.5*x
      plm(3,1) = (-7.5*x2 + 1.5)*y
      plm(3,2) = 15.0*x*y2
      plm(3,3) = -15.0*y3
    END IF
    IF (lpole .GE. 4) THEN
      x4 = x2*x2
      y4 = y2*y2
      plm(4,0) =  4.375*x4 - 3.75*x2 + 0.375
      plm(4,1) = (-17.5*x3 + 7.5*x)*y
      plm(4,2) = ( 52.5*x2 - 7.5  )*y2
      plm(4,3) = -105.0*x*y3
      plm(4,4) =  105.0*y4
    END IF
    IF (lpole .GE. 5) THEN
      x5 = x3*x2
      y5 = y3*y2
      plm(5,0) =  7.875*x5 - 8.75*x3 + 1.875*x
      plm(5,1) = (-39.375*x4 + 26.25*x2 - 1.875)*y
      plm(5,2) = ( 157.5*x3 - 52.5*x)*y2
      plm(5,3) = (-472.5*x2 + 52.5)  *y3
      plm(5,4) =  945.0*x*y4
      plm(5,5) = -945.0*y5
    END IF
    IF (lpole .GE. 6) THEN
      x6 = x3*x3
      y6 = y3*y3
      plm(6,0) = 14.4375*x6 - 19.6875*x4 + 6.5625*x2 - 0.3125
      plm(6,1) = (-86.625*x5 + 78.75*x3 - 13.125*x)*y
      plm(6,2) = ( 433.125*x4 - 236.25*x2 + 13.125)*y2
      plm(6,3) = (-1732.5*x3 + 472.5*x            )*y3
      plm(6,4) = ( 5197.5*x2 - 472.5              )*y4
      plm(6,5) = -10395.0*x*y5
      plm(6,6) =  10395.0*y6
    END IF
    IF (lpole .GE. 7) THEN
      x7 = x4*x3
      y7 = y4*y3
      plm(7,0) = 26.8125*x7 - 43.3125*x5 + 19.6875*x3 - 2.1875*x
      plm(7,1) = -187.6875*x6*y + 216.5625*x4*y - 59.0625*x2*y + 2.1875*y
      plm(7,2) = 1126.125*x5*y2 - 866.25*x3*y2 + 118.125*x*y2
      plm(7,3) = -5630.625*x4*y3 + 2598.75*x2*y3 - 118.125*y3
      plm(7,4) = 22522.5*x3*y4 - 5197.5*x*y4
      plm(7,5) = -67567.5*x2*y5 + 5197.5*y5
      plm(7,6) = 135135.0*x*y6
      plm(7,7) = -135135.0*y7
    END IF
    IF (lpole .GE. 8) THEN
      x8 = x4*x4
      y8 = y4*y4
      plm(8,0) = 50.2734375*x8 - 93.84375*x6 + 54.140625*x4 - 9.84375*x2 + 0.2734375
      plm(8,1) = -402.1875*x7*y + 563.0625*x5*y - 216.5625*x3*y + 19.6875*x*y
      plm(8,2) = 2815.3125*x6*y2 - 2815.3125*x4*y2 + 649.6875*x2*y2 - 19.6875*y2
      plm(8,3) = -16891.875*x5*y3 + 11261.25*x3*y3 - 1299.375*x*y3
      plm(8,4) = 84459.375*x4*y4 - 33783.75*x2*y4 + 1299.375*y4
      plm(8,5) = -337837.5*x3*y5 + 67567.5*x*y5
      plm(8,6) = 1013512.5*x2*y6 - 67567.5*y6
      plm(8,7) = -2027025.0*x*y7
      plm(8,8) = 2027025.0*y8
    END IF
    IF (lpole .GE. 9) THEN
      y9 = y5*y4
      x9 = x5*x4
      plm(9,0) = 94.9609375*x9 - 201.09375*x7 + 140.765625*x5 - 36.09375*x3 + 2.4609375*x
      plm(9,1) = -854.6484375*x8*y + 1407.65625*x6*y - 703.828125*x4*y + 108.28125*x2*y - 2.4609375*y
      plm(9,2) = 6837.1875*x7*y2 - 8445.9375*x5*y2 + 2815.3125*x3*y2 - 216.5625*x*y2
      plm(9,3) = -47860.3125*x6*y3 + 42229.6875*x4*y3 - 8445.9375*x2*y3 + 216.5625*y3
      plm(9,4) = 287161.875*x5*y4 - 168918.75*x3*y4 + 16891.875*x*y4
      plm(9,5) = -1435809.375*x4*y5 + 506756.25*x2*y5 - 16891.875*y5
      plm(9,6) = 5743237.5*x3*y6 - 1013512.5*x*y6
      plm(9,7) = -17229712.5*x2*y7 + 1013512.5*y7
      plm(9,8) = 34459425.0*x*y8
      plm(9,9) = -34459425.0*y9
    END IF
    !
    RETURN
    !
END SUBROUTINE setplm 
