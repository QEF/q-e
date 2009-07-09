! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
       subroutine mk_u(l,mmax,uvalue,jvalue,u2)
      implicit none
! Passed variables
       double precision uvalue,jvalue
       integer          mmax,l
!Local variables
      double precision rcl(4),xu,xj,u(7,7,7,7),u2(10,10),tmp
      integer          i,j
      logical          sw2

!  calculate F^0 .... F^3
      call rcl_init(l,uvalue,jvalue,rcl)
      xu = rcl(1)
      xj = 0.d0
      if(l.eq.1) xj = rcl(2)/5.d0
      if(l.eq.2) xj = (rcl(2)+rcl(3))/14.d0
      if(l.eq.3) xj = (4.d0*rcl(2)/15.d0+2.d0*
     .                      rcl(3)/11.d0+100.d0*rcl(4)/429.d0 )/6.d0
! Produce 4index Coulomb interaction matrix
      call u4ind(u,rcl,l)
       do i = 1, mmax
          do j = 1, mmax
             u2(i,j) = u(i,j,i,j) -u(i,j,j,i)
             u2(i+mmax,j+mmax) = u(i,j,i,j) -u(i,j,j,i)
             u2(i,j+mmax) = u(i,j,i,j)
             u2(i+mmax,j) = u(j,i,j,i)
          enddo
       enddo
  
      end
      
      subroutine u4ind(u,rcl,l)
c..................................................................u4ind
c----> calculation of <m1m2|1/r12|m3m4> coulomb integrals
      implicit none
* Passed variables
      integer l
      double precision u(7,7,7,7),rcl(4)
* Local variables
      integer mmax,k,k2p1,ms1,ms2,ms3,ms4,j,
     .             ms5,ms6,ms7,ms8
      double precision cgk,cgk0,cgk1,cgk2,uc(7,7,7,7),
     .                 xk,xm1,xm2,xm3,xm,xm4,xl,
     .                 yor(7,7),yoi(7,7)
      double complex am1,am2,am3,am4,amz
      data amz/(0.d0,0.d0)/
      logical sw3
      external dinit, cgk
      
      intrinsic dfloat
      
      mmax=2*l+1
      xl=dfloat(l)

      call dinit(uc,7*7*7*7)

      do 1 k = 0, 2*l, 2
      k2p1 = k/2 + 1
      xk = dfloat(k)
      cgk0 =  cgk(xl,0.d0,xk,0.d0,xl,0.d0)
      do 1 ms1 = 1,mmax
      xm1 = dfloat(ms1-l-1)
      do 1 ms2 = 1,mmax
      xm2 = dfloat(ms2-l-1)
      do 1 ms3 = 1,mmax
      xm3 = dfloat(ms3-l-1)
      xm  = xm1 - xm3
      do 1 ms4 = 1,mmax
      if ((ms1+ms2-ms3-ms4).ne.0) go to 1
      xm4 = dfloat(ms4-l-1)
      cgk1 =  cgk(xl,xm3,xk,xm,xl,xm1)
      cgk2 =  cgk(xl,xm2,xk,xm,xl,xm4)
      uc(ms1,ms2,ms3,ms4) = uc(ms1,ms2,ms3,ms4) +
     .                       rcl(k2p1)*cgk0*cgk0*cgk1*cgk2
    1 continue
        call dinit(u,7*7*7*7)
        do 11 ms1 = 1,mmax
          do 11 ms2 = 1,mmax
            do 11 ms3 = 1,mmax
              do 11 ms4 = 1,mmax
                u(ms1,ms2,ms3,ms4) = uc(ms1,ms2,ms3,ms4)
   11   continue
      call ctormt(yor,yoi,l)
      call dinit(u,7*7*7*7)

      do 2 ms1=1,mmax
      do 2 ms2=1,mmax
      do 2 ms3=1,mmax
      do 2 ms4=1,mmax
      do 3 ms5=1,mmax
      am1 = dcmplx(yor(ms1,ms5),-yoi(ms1,ms5))
      if (am1.eq.amz) go to 3
      do 4 ms6=1,mmax
      am2 = dcmplx(yor(ms2,ms6),-yoi(ms2,ms6))
      if (am2.eq.amz) go to 4
      do 5 ms7=1,mmax
      am3 = dcmplx(yor(ms3,ms7),yoi(ms3,ms7))
      if (am3.eq.amz) go to 5
      do 6 ms8=1,mmax
      am4 = dcmplx(yor(ms4,ms8),yoi(ms4,ms8))
      if (am4.eq.amz) go to 6
      u(ms1,ms2,ms3,ms4) = u(ms1,ms2,ms3,ms4) + am1*am2*am3*am4*
     .                     uc(ms5,ms6,ms7,ms8)
    6 continue
    5 continue
    4 continue
    3 continue
    2 continue
    7 continue
 
   12   continue

333   format(30f8.4)

      end

c....................................................................cgk
      double precision function cgk(a,al,b,be,c,ga)
      implicit none

* Passed variables
      double precision a,al,b,be,c,ga

* Local variables
      integer z,zmin,zmax,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13
      double precision fa(0:20)

      data fa/1.d0, 1.d0, 2.d0, 6.d0, 24.d0, 12.d1, 72.d1, 504.d1,
     .        4032.d1, 36288.d1, 36288.d2, 399168.d2, 4790016.d2,
     .        62270208.d2, 871782912.d2, 1307674368.d3, 20922789888.d3,
     .        355687428096.d3, 6402373705728.d3, 121645100408832.d3,
     .        243290200817664.d4/

      intrinsic idint,max0,min0,dsqrt

      i1=0
      i2=idint(a+b-c)
      i3=idint(a-al)
      i4=idint(b+be)
      i5=idint(c-b+al)
      i6=idint(c-a-be)
      zmin=max0(i1,-i5,-i6)
      zmax=min0(i2, i3, i4)
      cgk=0.d0
      if (dabs(al).gt.a)   return
      if (dabs(be).gt.b)   return
      if (dabs(ga).gt.c)   return
      if ( zmin.gt.zmax )  return
      if ( (al+be).ne.ga ) return
      i7=idint(a-b+c)
      i8=idint(c+b-a)
      i9=idint(c+b+a)
      i10=idint(a+al)
      i11=idint(b-be)
      i12=idint(c+ga)
      i13=idint(c-ga)
      do 1 z=zmin,zmax
1     cgk=cgk+(-1.d0)**z/(fa(z)*fa(i2-z)*fa(i3-z)*fa(i4-z)*fa(i5+z)*
     .                    fa(i6+z))
      cgk=cgk*dsqrt(fa(i2)*fa(i7)*fa(i8)*fa(i10)*fa(i3)*
     .              fa(i4)*fa(i11)*fa(i12)*fa(i13)*(2.d0*c+1.d0)/
     .              fa(i9+1))

      end

      subroutine ctormt(yor,yoi,l)
c.................................................................ctormt
c
c---->    transformation from (ms) to real harmonics basis set
c
      implicit none
      integer l
      double precision yor(7,7),yoi(7,7),sqtwo
      external dinit
      intrinsic dsqrt
      call dinit(yor,7*7)
      call dinit(yoi,7*7)
      sqtwo=1.d0/dsqrt(2.d0)
      if (l.eq.0) then
        yor(1,1)=1.d0
      elseif (l.eq.1) then
        yoi(1,1)= sqtwo
        yoi(1,3)= sqtwo
        yor(2,2)=1.d0
        yor(3,1)= sqtwo
        yor(3,3)=-sqtwo
      elseif (l.eq.2) then
        yoi(1,1)= sqtwo
        yoi(1,5)=-sqtwo
        yoi(2,2)= sqtwo
        yoi(2,4)= sqtwo
        yor(3,3)=1.d0
        yor(4,2)= sqtwo
        yor(4,4)=-sqtwo
        yor(5,1)= sqtwo
        yor(5,5)= sqtwo
      elseif (l.eq.3) then
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
      end if
      end



      subroutine rcl_init(l,uvalue,jvalue,rcl)
      implicit none

* passed variables
      integer          l
      double precision rcl(4),uvalue,jvalue

* local variables
      double precision ev2ry,uv,jv
      integer          i
      call dinit(rcl,4)
      ev2ry = 13.6058d0

          uv = uvalue 
          jv = jvalue            
          rcl(1) = uv
          if(l .eq. 1) then
            rcl(2) = jv *5.d0
            elseif(l .eq. 2) then
              rcl(2) = jv * 14d0 / (1.d0 + 0.63d0)
              rcl(3) = 0.63d0 * rcl(2)
            elseif(l .eq. 3) then
              rcl(2) = 6435.d0 * jv / (286.d0 + 195.d0 *
     .                     451.d0 / 675.d0 + 250.d0 * 1001.d0 / 2025.d0)
              rcl(3) = 451.d0 * rcl(2) / 675.d0
              rcl(4) = 1001.d0 * rcl(2) / 2025.d0
           endif
      end

      subroutine dinit(array,leng)
      implicit none
C Passed variables:
      integer leng
      double precision array(leng)
C Local variables:
      integer i,m,mp1

      m = mod(leng,5)
      if( m .ne. 0 ) then
        do i = 1,m
          array(i) = 0.d0
        enddo
        if( leng .lt. 5 ) return
      endif
      mp1 = m + 1
      do i = mp1,leng,5
        array(i) = 0.d0
        array(i + 1) = 0.d0
        array(i + 2) = 0.d0
        array(i + 3) = 0.d0
        array(i + 4) = 0.d0
      enddo
      end
