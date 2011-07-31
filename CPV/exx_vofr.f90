
      SUBROUTINE getvofr( nnrtot, hcub, n_exx, rho, v, vnew, vold, tran )

!===============================================================
! Given charge density, get the potential using multipole expansion
! for boundary region and solving poisson equation for inside box
!===============================================================

      USE kinds,                   ONLY  :  DP
      USE fft_base,                ONLY  :  dffts
      USE cp_main_variables,       ONLY  :  odtothd_in_sp, thdtood_in_sp
      USE cp_main_variables,       ONLY  :  thdtood, np_in_sp, np_in_sp2
      USE cp_main_variables,       ONLY  :  coeke, nord2, lmax
      USE mp_global,               ONLY  :  me_image

      implicit none

      INTEGER  nnrtot, n_exx, tran(3)
      REAL(DP) hcub
      REAL(DP) rho(np_in_sp)
!OUTPUT
      REAL(DP) v(nnrtot)
!INOUT
      REAL(DP) vnew(np_in_sp), vold(np_in_sp)

      REAL(DP), ALLOCATABLE  :: v_in_sp(:)
      COMPLEX(DP) qlm(0:lmax, 0:lmax)
      INTEGER  i, j, k, ii, jj, kk, ir, ip
!========================================================================

!-----------------------------------------------------------------------
!First, we calculate the potential outside the inner sphere, which will
!also give the boundary values for potential inside the sphere.
!-----------------------------------------------------------------------

       call start_clock('getqlm')
       CALL getqlm( hcub, rho, qlm )
       call stop_clock('getqlm')
       ALLOCATE( v_in_sp(1:np_in_sp+np_in_sp2) )

       call start_clock('exx_bound')
       CALL exx_boundaryv( v_in_sp, qlm)
       call stop_clock('exx_bound')
       DO ir = 1+np_in_sp, np_in_sp2+np_in_sp
          i = odtothd_in_sp(1,ir)
          j = odtothd_in_sp(2,ir)
          k = odtothd_in_sp(3,ir)

          ii = i - tran(1)
          jj = j - tran(2)
          kk = k - tran(3)
   
          if( ii .gt. dffts%nr1) ii = ii - dffts%nr1
          if( jj .gt. dffts%nr2) jj = jj - dffts%nr2
          if( kk .gt. dffts%nr3) kk = kk - dffts%nr3
   
          if( ii .lt. 1) ii = ii + dffts%nr1
          if( jj .lt. 1) jj = jj + dffts%nr2
          if( kk .lt. 1) kk = kk + dffts%nr3
   
          ip = thdtood(ii,jj,kk)
          v(ip) = v_in_sp(ir)
       ENDDO

!-------------------------------------------------------------------------------
!                   NEXT, THE POTENTIAL INSIDE THE SPHERE
!-------------------------------------------------------------------------------

       IF(n_exx .gt. 1)THEN
          DO ir = 1, np_in_sp
             v_in_sp(ir) = vnew(ir) + vnew(ir) - vold(ir)
          ENDDO
       ENDIF

       CALL geterho(rho, v_in_sp)
       CALL hpotcg( rho, v_in_sp)

       IF( n_exx .eq. 1)THEN
          DO ir = 1, np_in_sp
             vnew(ir) = v_in_sp(ir)
          ENDDO
          vold(:) = vnew(:)
       ELSE
          vold(:) = vnew(:)
          DO ir = 1, np_in_sp
             vnew(ir) = v_in_sp(ir)
          ENDDO
       ENDIF

       DO ir = 1, np_in_sp
          i = odtothd_in_sp(1,ir)
          j = odtothd_in_sp(2,ir)
          k = odtothd_in_sp(3,ir)
         
          ii = i - tran(1)
          jj = j - tran(2)
          kk = k - tran(3)
         
          if( ii .gt. dffts%nr1) ii = ii - dffts%nr1
          if( jj .gt. dffts%nr2) jj = jj - dffts%nr2
          if( kk .gt. dffts%nr3) kk = kk - dffts%nr3
         
          if( ii .lt. 1) ii = ii + dffts%nr1
          if( jj .lt. 1) jj = jj + dffts%nr2
          if( kk .lt. 1) kk = kk + dffts%nr3
         
           ip = thdtood(ii,jj,kk)
           v(ip) = v_in_sp(ir)
        ENDDO

        DEALLOCATE( v_in_sp)

!===============================================================================
      RETURN
      END SUBROUTINE getvofr


!===============================================================================
!===============================================================================
       subroutine getqlm(hcub, rho, qlm)

       USE kinds,                   ONLY  :  DP
       USE cp_main_variables,       ONLY  :  np_in_sp, np_in_sp2, lpole=>lmax
       USE cp_main_variables,       ONLY  :  xx=>xx_in_sp, yy=>yy_in_sp, zz=>zz_in_sp, clm

       IMPLICIT NONE

       REAL(DP)   hcub, rho(1:np_in_sp)
       REAL(DP)   rrho,xx2,yy2,zz2,xy,r2,x,y
       REAL(DP)   rinv(0:lpole),r(0:lpole)
 
       REAL(DP)    plm(0:lpole, 0:lpole)  !temporary storage of the associated Legendre polynom
       complex(DP) qlm(0:lpole, 0:lpole)
       complex(DP) cxy(1:lpole)      !coefficient array: e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)

       REAL(DP)   hcub2, zero, one, two , tmp1, tmp2
       parameter (zero = 1.0E-10, one = 1.d0, two = 2.d0 )
       INTEGER    i,j,l,m

!------------------------------------------------------------------------------
       hcub2 = two*hcub
       qlm(:,:) = 0.d0

       do j = 1,np_in_sp
          rrho = rho(j)
          xx2  = xx(j)*xx(j)
          yy2  = yy(j)*yy(j)
          zz2  = zz(j)*zz(j)
          r2   = xx2 + yy2 + zz2

          r(1) = dsqrt(r2)
          qlm(0,0) = qlm(0,0) + rrho

          if (r(1) .gt. zero) then
            do l = 2, lpole
               r(l) = r(l-1)*r(1)           !r(l) = r^l (higher powers of r)
            enddo
 
            rinv(1) = one/r(1)
            xy = dsqrt(xx2+yy2)
 
            x = zz(j)*rinv(1)               ! x = cos(theta_j)
            y = xy*rinv(1)                  ! y = sin(theta_j)

            call setplm(x, y, lpole, plm)   ! associate Legendre polynomials in plm(l,m)

            do l = 1, lpole
               qlm(l,0) = qlm(l,0) + rrho*plm(l,0)*r(l)   !qlm(l,m=0)
            enddo

            if (xy .gt. zero) then
!
               cxy(1) = CMPLX(xx(j),yy(j))
               cxy(1) = cxy(1)/xy
               do m = 2, lpole
                  cxy(m) = cxy(m-1)*cxy(1)  !cxy(m) = exp(i*m*phi_j) = (cos(phi_j) + i*sin(phi_j))^m
               enddo                        !       = ((xx(j) + i*yy(j))/(sqrt(xx(j)^2 + yy(j)^2))^m

               do l = 1, lpole
                  tmp1 = rrho*r(l)
                  do m = 1, l
                     tmp2     = plm(l,m)*tmp1
                     qlm(l,m) = qlm(l,m) + cxy(m)*tmp2      !qlm(l, m>0)
                  enddo
                enddo

             endif
          endif
       enddo

       do l = 0, lpole
          qlm(l,0) = qlm(l,0)*clm(l,0)*hcub
       enddo
!
       do l = 1, lpole
          do m = 1, l, 1
             qlm(l,m) = qlm(l,m)*clm(l,m)*hcub2
          enddo
       enddo

! ==================================================================
       RETURN
       END

!

!===============================================================================
       subroutine exx_boundaryv(v_in_sp,qlm)

       USE kinds,                   ONLY  :  DP
       USE cp_main_variables,       ONLY  :  np_in_sp, np_in_sp2, lpole=>lmax 
       USE cp_main_variables,       ONLY  :  xx=>xx_in_sp, yy=>yy_in_sp, zz=>zz_in_sp
       USE mp_global,               ONLY  : me_image

       implicit   none

       real(DP)   v_in_sp( 1:np_in_sp+np_in_sp2 )

! For each grid point:
       real(DP)   plm(0:lpole, 0:lpole)
       complex*16 qlm(0:lpole, 0:lpole)
       complex*16 cxy(1:lpole)   ! e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
       complex*16 cpole, qlmc

       real(DP)   rinv(0:lpole), r(0:lpole)
       real(DP)   xh,yh,zh, xx2, yy2, zz2, xy, x, y, r2, zero, one, plmr
       integer    i, l, m

       zero = 0.d0
       one  = 1.d0
       do 999 i = 1 + np_in_sp, np_in_sp2 + np_in_sp

          xh = xx(i)
          yh = yy(i)
          zh = zz(i)
   
          xx2 = xh*xh
          yy2 = yh*yh
          zz2 = zh*zh

          r2 = xx2 + yy2 + zz2
          r(1) = dsqrt(r2)
          if(r(1) .lt. 0.000001)print *, 'i =',i, r(1)
          rinv(0) = one/r(1)
          do l = 1,lpole
             rinv(l) = rinv(l-1)*rinv(0)  !rinv(l) = 1/r^(l+1).
          enddo
 
          xy= dsqrt(xx2+yy2)
          x = zh*rinv(0)           ! x=cos(theta)
          y = xy*rinv(0)           ! y=sin(theta)
  
! evaluate associate Legendre polynomials for x and y with l from 0 to 
! lpole and m from 0 to l. store each in plm(l,m)
          call setplm(x, y, lpole, plm)

          cpole = zero

          do l = 0,lpole
             plmr = plm(l,0)*rinv(l)
             cpole = cpole + qlm(l,0)*plmr    ! m=0 terms
          enddo

          if (xy .gt. zero) then
             cxy(1) = CMPLX(xh,-yh)
             cxy(1) = cxy(1)/xy
             do m=2, lpole, 1
                cxy(m) = cxy(m-1)*cxy(1)
             enddo

             do l = 1, lpole
                do m = 1, l, 1                          ! m>0 terms
                   plmr = plm(l,m)*rinv(l)
                   qlmc = qlm(l,m)*cxy(m)
                   cpole = cpole + plmr * qlmc
                enddo
             enddo
          endif

          v_in_sp(i) = REAL(cpole)
  999    continue

! ===========================================================================
       return
       end

! ===========================================================================
! ===========================================================================
       subroutine geterho(rho, v_in_sp)

       USE kinds,                   ONLY  :  DP
       USE cp_main_variables,       ONLY  :  np_in_sp, np_in_sp2, nord2,coeke, lap_dir_num, lap_neig
       USE cp_main_variables,       ONLY  :  odtothd_in_sp, thdtood_in_sp

       implicit none

       REAL(DP) rho(np_in_sp), v_in_sp( 1:np_in_sp+np_in_sp2 )
       INTEGER  i, j, k, ip, jp, ish, ipp, ipm, jpp, jpm, kpp, kpm


       select case (lap_dir_num)
       case (0)

       DO ip = 1, np_in_sp

!(i,j,k) is within the first sphere
          i = odtothd_in_sp(1,ip)
          j = odtothd_in_sp(2,ip)
          k = odtothd_in_sp(3,ip)

          DO ish  = 1, nord2

             ipp = thdtood_in_sp( i+ish, j,     k     )
             ipm = thdtood_in_sp( i-ish, j,     k     )
             jpp = thdtood_in_sp( i,     j+ish, k     )
             jpm = thdtood_in_sp( i,     j-ish, k     )
             kpp = thdtood_in_sp( i,     j,     k+ish )
             kpm = thdtood_in_sp( i,     j,     k-ish )

             IF(ipp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipp)
             IF(ipm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipm)
             IF(jpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpp)
             IF(jpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpm)
             IF(kpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpp)
             IF(kpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpm)

          ENDDO
       ENDDO

      case (1)
       DO ip = 1, np_in_sp

!(i,j,k) is within the first sphere
          i = odtothd_in_sp(1,ip)
          j = odtothd_in_sp(2,ip)
          k = odtothd_in_sp(3,ip)

          DO ish  = 1, nord2

             ipp = thdtood_in_sp( i+ish, j,     k     )
             ipm = thdtood_in_sp( i-ish, j,     k     )
             jpp = thdtood_in_sp( i,     j+ish, k     )
             jpm = thdtood_in_sp( i,     j-ish, k     )
             kpp = thdtood_in_sp( i,     j,     k+ish )
             kpm = thdtood_in_sp( i,     j,     k-ish )

             IF(ipp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipp)
             IF(ipm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipm)
             IF(jpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpp)
             IF(jpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpm)
             IF(kpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpp)
             IF(kpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpm)

             ipp = thdtood_in_sp( i+ish*lap_neig(1,1),j+ish*lap_neig(2,1),k+ish*lap_neig(3,1) )
             ipm = thdtood_in_sp( i-ish*lap_neig(1,1),j-ish*lap_neig(2,1),k-ish*lap_neig(3,1) )

             IF(ipp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,4)*v_in_sp(ipp)
             IF(ipm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,4)*v_in_sp(ipm)

          ENDDO
       ENDDO

      case (2)
       DO ip = 1, np_in_sp

!(i,j,k) is within the first sphere
          i = odtothd_in_sp(1,ip)
          j = odtothd_in_sp(2,ip)
          k = odtothd_in_sp(3,ip)

          DO ish  = 1, nord2

             ipp = thdtood_in_sp( i+ish, j,     k     )
             ipm = thdtood_in_sp( i-ish, j,     k     )
             jpp = thdtood_in_sp( i,     j+ish, k     )
             jpm = thdtood_in_sp( i,     j-ish, k     )
             kpp = thdtood_in_sp( i,     j,     k+ish )
             kpm = thdtood_in_sp( i,     j,     k-ish )

             IF(ipp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipp)
             IF(ipm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipm)
             IF(jpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpp)
             IF(jpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpm)
             IF(kpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpp)
             IF(kpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpm)

             ipp = thdtood_in_sp( i+ish*lap_neig(1,1),j+ish*lap_neig(2,1),k+ish*lap_neig(3,1) )
             ipm = thdtood_in_sp( i-ish*lap_neig(1,1),j-ish*lap_neig(2,1),k-ish*lap_neig(3,1) )

             jpp = thdtood_in_sp( i+ish*lap_neig(1,2),j+ish*lap_neig(2,2),k+ish*lap_neig(3,2) )
             jpm = thdtood_in_sp( i-ish*lap_neig(1,2),j-ish*lap_neig(2,2),k-ish*lap_neig(3,2) )

             IF(ipp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,4)*v_in_sp(ipp)
             IF(ipm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,4)*v_in_sp(ipm)
             IF(jpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,5)*v_in_sp(jpp)
             IF(jpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,5)*v_in_sp(jpm)

          ENDDO
       ENDDO

      case (3)
       DO ip = 1, np_in_sp

!(i,j,k) is within the first sphere
          i = odtothd_in_sp(1,ip)
          j = odtothd_in_sp(2,ip)
          k = odtothd_in_sp(3,ip)

          DO ish  = 1, nord2

             ipp = thdtood_in_sp( i+ish, j,     k     )
             ipm = thdtood_in_sp( i-ish, j,     k     )
             jpp = thdtood_in_sp( i,     j+ish, k     )
             jpm = thdtood_in_sp( i,     j-ish, k     )
             kpp = thdtood_in_sp( i,     j,     k+ish )
             kpm = thdtood_in_sp( i,     j,     k-ish )

             IF(ipp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipp)
             IF(ipm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,1)*v_in_sp(ipm)
             IF(jpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpp)
             IF(jpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,2)*v_in_sp(jpm)
             IF(kpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpp)
             IF(kpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,3)*v_in_sp(kpm)

             ipp = thdtood_in_sp( i+ish*lap_neig(1,1),j+ish*lap_neig(1,2),k+ish*lap_neig(1,3) )
             ipm = thdtood_in_sp( i-ish*lap_neig(1,1),j-ish*lap_neig(1,2),k-ish*lap_neig(1,3) )

             jpp = thdtood_in_sp( i+ish*lap_neig(2,1),j+ish*lap_neig(2,2),k+ish*lap_neig(2,3) )
             jpm = thdtood_in_sp( i-ish*lap_neig(2,1),j-ish*lap_neig(2,2),k-ish*lap_neig(2,3) )

             kpp = thdtood_in_sp( i+ish*lap_neig(3,1),j+ish*lap_neig(3,2),k+ish*lap_neig(3,3) )
             kpm = thdtood_in_sp( i-ish*lap_neig(3,1),j-ish*lap_neig(3,2),k-ish*lap_neig(3,3) )

             IF(ipp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,4)*v_in_sp(ipp)
             IF(ipm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,4)*v_in_sp(ipm)
             IF(jpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,5)*v_in_sp(jpp)
             IF(jpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,5)*v_in_sp(jpm)
             IF(kpp .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,6)*v_in_sp(kpp)
             IF(kpm .gt. np_in_sp) rho(ip) = rho(ip) - coeke(ish,6)*v_in_sp(kpm)

          ENDDO
       ENDDO

      endselect
! ==================================================================
       return
       end


! ==================================================================
       subroutine setplm(x, y, lpole, plm)

! this subroutine calculates all of the associated Legendre polynomials up 
! to a supplied lpole (maximum lpole is 9). 

       implicit none

! INPUT VARIABLES
! cos(theta),sin(theta), for a given grid point
       real*8 x,y
! order of multipole expansion
       integer lpole
! OUTPUT VARIABLES
! array containing P_{lm}, the associated Legendre polynomials
       real*8 plm(0:lpole, 0:lpole)
! WORK VARIABLES:
! powers of x, y: xn=x^n, yn=y^n
       real*8 x2, x3, x4, x5, x6, x7, x8, x9
       real*8 y2, y3, y4, y5, y6, y7, y8, y9

       plm(0,0) = 1.0
       if (lpole .ge. 1) then
         plm(1,0) = x
         plm(1,1) = -y
       endif
       if (lpole .ge. 2) then
         x2 = x*x
         y2 = y*y
         plm(2,0) = 1.5*x2 - 0.5
         plm(2,1) = -3.0*x*y
         plm(2,2) = 3.0*y2
       endif
       if (lpole .ge. 3) then
         x3 = x2*x
         y3 = y2*y
         plm(3,0) = 2.5*x3 - 1.5*x
         plm(3,1) = (-7.5*x2 + 1.5)*y
         plm(3,2) = 15.0*x*y2
         plm(3,3) = -15.0*y3
       endif
       if (lpole .ge. 4) then
         x4 = x2*x2
         y4 = y2*y2
         plm(4,0) =  4.375*x4 - 3.75*x2 + 0.375
         plm(4,1) = (-17.5*x3 + 7.5*x)*y
         plm(4,2) = ( 52.5*x2 - 7.5  )*y2
         plm(4,3) = -105.0*x*y3
         plm(4,4) =  105.0*y4
       endif
       if (lpole .ge. 5) then
         x5 = x3*x2
         y5 = y3*y2
         plm(5,0) =  7.875*x5 - 8.75*x3 + 1.875*x
         plm(5,1) = (-39.375*x4 + 26.25*x2 - 1.875)*y
         plm(5,2) = ( 157.5*x3 - 52.5*x)*y2
         plm(5,3) = (-472.5*x2 + 52.5)  *y3
         plm(5,4) =  945.0*x*y4
         plm(5,5) = -945.0*y5
       endif
       if (lpole .ge. 6) then
         x6 = x3*x3
         y6 = y3*y3
         plm(6,0) = 14.4375*x6 - 19.6875*x4 + 6.5625*x2 - 0.3125
         plm(6,1) = (-86.625*x5 + 78.75*x3 - 13.125*x)*y
         plm(6,2) = ( 433.125*x4 - 236.25*x2 + 13.125)*y2
         plm(6,3) = (-1732.5*x3 + 472.5*x            )*y3
         plm(6,4) = ( 5197.5*x2 - 472.5              )*y4
         plm(6,5) = -10395.0*x*y5
         plm(6,6) =  10395.0*y6
       endif
       if (lpole .ge. 7) then
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
       endif
       if (lpole .ge. 8) then
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
       endif
       if (lpole .ge. 9) then
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
       endif

       end 
