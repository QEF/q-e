! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .


SUBROUTINE mk_u(l,mmax,uvalue,jvalue,u2)
   IMPLICIT NONE
! Passed variables
   DOUBLE PRECISION uvalue,jvalue
   INTEGER          mmax,l
!Local variables
   DOUBLE PRECISION rcl(4),xu,xj,u(7,7,7,7),u2(10,10),tmp
   INTEGER          i,j
   LOGICAL          sw2

!  calculate F^0 .... F^3
   CALL rcl_init(l,uvalue,jvalue,rcl)
   xu = rcl(1)
   xj = 0.d0
   IF(l==1) xj = rcl(2)/5.d0
   IF(l==2) xj = (rcl(2)+rcl(3))/14.d0
   IF(l==3) xj = (4.d0*rcl(2)/15.d0+2.d0*rcl(3)/11.d0+100.d0*rcl(4)/429.d0 )/6.d0
! Produce 4index Coulomb interaction matrix
   CALL u4ind(u,rcl,l)
   DO i = 1, mmax
      DO j = 1, mmax
         u2(i,j) = u(i,j,i,j) -u(i,j,j,i)
         u2(i+mmax,j+mmax) = u(i,j,i,j) -u(i,j,j,i)
         u2(i,j+mmax) = u(i,j,i,j)
         u2(i+mmax,j) = u(j,i,j,i)
      ENDDO
   ENDDO
END SUBROUTINE mk_u


SUBROUTINE u4ind(u,rcl,l)
!----> calculation of <m1m2|1/r12|m3m4> coulomb integrals
   IMPLICIT NONE
! Passed variables
   INTEGER l
   DOUBLE PRECISION u(7,7,7,7),rcl(4)
! Local variables
   INTEGER mmax,k,k2p1,ms1,ms2,ms3,ms4,j,ms5,ms6,ms7,ms8
   DOUBLE PRECISION cgk,cgk0,cgk1,cgk2,uc(7,7,7,7),                  &
                    xk,xm1,xm2,xm3,xm,xm4,xl,                        &
                    yor(7,7),yoi(7,7)
   DOUBLE COMPLEX am1,am2,am3,am4,amz
   DATA amz/(0.d0,0.d0)/
   LOGICAL sw3
   EXTERNAL cgk

   mmax=2*l+1
   xl=dble(l)

   CALL dinit(uc,7*7*7*7)

   DO k = 0, 2*l, 2
      k2p1 = k/2 + 1
      xk = dble(k)
      cgk0 =  cgk(xl,0.d0,xk,0.d0,xl,0.d0)
      DO ms1 = 1,mmax
         xm1 = dble(ms1-l-1)
         DO ms2 = 1,mmax
            xm2 = dble(ms2-l-1)
            DO ms3 = 1,mmax
               xm3 = dble(ms3-l-1)
               xm  = xm1 - xm3
               DO ms4 = 1,mmax
                  IF ((ms1+ms2-ms3-ms4)/=0) CYCLE
                  xm4 = dble(ms4-l-1)
                  cgk1 =  cgk(xl,xm3,xk,xm,xl,xm1)
                  cgk2 =  cgk(xl,xm2,xk,xm,xl,xm4)
                  uc(ms1,ms2,ms3,ms4) = uc(ms1,ms2,ms3,ms4) +         &
                                        rcl(k2p1)*cgk0*cgk0*cgk1*cgk2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   CALL dinit(u,7*7*7*7)

   DO ms1 = 1,mmax
      DO ms2 = 1,mmax
         DO ms3 = 1,mmax
            DO ms4 = 1,mmax
               u(ms1,ms2,ms3,ms4) = uc(ms1,ms2,ms3,ms4)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   CALL ctormt(yor,yoi,l)
   CALL dinit(u,7*7*7*7)

   DO ms1=1,mmax
      DO ms2=1,mmax
         DO ms3=1,mmax
            DO ms4=1,mmax
               DO ms5=1,mmax
                  am1 = dcmplx(yor(ms1,ms5),-yoi(ms1,ms5))
                  IF (am1==amz) CYCLE
                  DO ms6=1,mmax
                     am2 = dcmplx(yor(ms2,ms6),-yoi(ms2,ms6))
                     IF (am2==amz) CYCLE
                     DO ms7=1,mmax
                        am3 = dcmplx(yor(ms3,ms7),yoi(ms3,ms7))
                        IF (am3==amz) CYCLE
                        DO ms8=1,mmax
                           am4 = dcmplx(yor(ms4,ms8),yoi(ms4,ms8))
                           IF (am4==amz) CYCLE
                           u(ms1,ms2,ms3,ms4) = u(ms1,ms2,ms3,ms4)   &
                                        + am1*am2*am3*am4*uc(ms5,ms6,ms7,ms8)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

333   FORMAT(30f8.4)

END SUBROUTINE u4ind


DOUBLE PRECISION FUNCTION cgk(a,al,b,be,c,ga)
   IMPLICIT NONE

! Passed variables
   DOUBLE PRECISION a,al,b,be,c,ga

! Local variables
   INTEGER z,zmin,zmax,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13
   DOUBLE PRECISION fa(0:20)

   DATA fa/1.d0, 1.d0, 2.d0, 6.d0, 24.d0, 12.d1, 72.d1, 504.d1,      &
           4032.d1, 36288.d1, 36288.d2, 399168.d2, 4790016.d2,       &
           62270208.d2, 871782912.d2, 1307674368.d3, 20922789888.d3, &
           355687428096.d3, 6402373705728.d3, 121645100408832.d3,    &
           243290200817664.d4/

   INTRINSIC idint,max0,min0,dsqrt

   i1=0
   i2=idint(a+b-c)
   i3=idint(a-al)
   i4=idint(b+be)
   i5=idint(c-b+al)
   i6=idint(c-a-be)
   zmin=max0(i1,-i5,-i6)
   zmax=min0(i2, i3, i4)
   cgk=0.d0
   IF (dabs(al)>a)   RETURN
   IF (dabs(be)>b)   RETURN
   IF (dabs(ga)>c)   RETURN
   IF ( zmin>zmax )  RETURN
   IF ( (al+be)/=ga ) RETURN
   i7=idint(a-b+c)
   i8=idint(c+b-a)
   i9=idint(c+b+a)
   i10=idint(a+al)
   i11=idint(b-be)
   i12=idint(c+ga)
   i13=idint(c-ga)
   DO z=zmin,zmax
      cgk=cgk+(-1.d0)**z/(fa(z)*fa(i2-z)*fa(i3-z)*fa(i4-z)*fa(i5+z)* &
                          fa(i6+z))
   ENDDO
   cgk=cgk*dsqrt(fa(i2)*fa(i7)*fa(i8)*fa(i10)*fa(i3)*                &
                 fa(i4)*fa(i11)*fa(i12)*fa(i13)*(2.d0*c+1.d0)/       &
                 fa(i9+1))

END FUNCTION cgk


SUBROUTINE ctormt(yor,yoi,l)
!.................................................................ctormt
!
!---->    transformation from (ms) to real harmonics basis set
!
   IMPLICIT NONE
   INTEGER l
   DOUBLE PRECISION yor(7,7),yoi(7,7),sqtwo
   INTRINSIC dsqrt
   CALL dinit(yor,7*7)
   CALL dinit(yoi,7*7)
   sqtwo=1.d0/dsqrt(2.d0)
   IF (l==0) THEN
      yor(1,1)=1.d0
   ELSEIF (l==1) THEN
      yoi(1,1)= sqtwo
      yoi(1,3)= sqtwo
      yor(2,2)=1.d0
      yor(3,1)= sqtwo
      yor(3,3)=-sqtwo
   ELSEIF (l==2) THEN
      yoi(1,1)= sqtwo
      yoi(1,5)=-sqtwo
      yoi(2,2)= sqtwo
      yoi(2,4)= sqtwo
      yor(3,3)=1.d0
      yor(4,2)= sqtwo
      yor(4,4)=-sqtwo
      yor(5,1)= sqtwo
      yor(5,5)= sqtwo
   ELSEIF (l==3) THEN
      yoi(1,1)= sqtwo
      yoi(1,7)= sqtwo
      yoi(2,2)= sqtwo
      yoi(2,6)=-sqtwo
      yoi(3,3)= sqtwo
      yoi(3,5)= sqtwo
      yor(4,4)=1.d0
      yor(5,3)= sqtwo
      yor(5,5)=-sqtwo
      yor(6,2)= sqtwo
      yor(6,6)= sqtwo
      yor(7,1)= sqtwo
      yor(7,7)=-sqtwo
   ENDIF
END SUBROUTINE ctormt



SUBROUTINE rcl_init(l,uvalue,jvalue,rcl)
   IMPLICIT NONE

! passed variables
   INTEGER          l
   DOUBLE PRECISION rcl(4),uvalue,jvalue

! local variables
   DOUBLE PRECISION uv,jv
   INTEGER          i
   CALL dinit(rcl,4)

   uv = uvalue
   jv = jvalue
   rcl(1) = uv
   IF(l == 1) THEN
      rcl(2) = jv *5.d0
   ELSEIF(l == 2) THEN
      rcl(2) = jv * 14d0 / (1.d0 + 0.63d0)
      rcl(3) = 0.63d0 * rcl(2)
   ELSEIF(l == 3) THEN
      rcl(2) = 6435.d0 * jv / (286.d0 + 195.d0 *                    &
         451.d0 / 675.d0 + 250.d0 * 1001.d0 / 2025.d0)
      rcl(3) = 451.d0 * rcl(2) / 675.d0
      rcl(4) = 1001.d0 * rcl(2) / 2025.d0
   ENDIF
END SUBROUTINE rcl_init


SUBROUTINE dinit(array,leng)
   IMPLICIT NONE
! Passed variables:
   INTEGER leng
   DOUBLE PRECISION array(leng)
! Local variables:
   INTEGER i,m,mp1

   m = mod(leng,5)
   IF( m /= 0 ) THEN
      DO i = 1,m
         array(i) = 0.d0
      ENDDO
      IF( leng < 5 ) RETURN
   ENDIF
   mp1 = m + 1
   DO i = mp1,leng,5
      array(i) = 0.d0
      array(i + 1) = 0.d0
      array(i + 2) = 0.d0
      array(i + 3) = 0.d0
      array(i + 4) = 0.d0
   ENDDO
END SUBROUTINE dinit
