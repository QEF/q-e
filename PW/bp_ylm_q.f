c
c-----------------------------------------------------------------------
      subroutine ylm_q(lmax,gx,g,ylm)
c-----------------------------------------------------------------------
c     REAL SPHERICAL HARMONICS,  L IS COMBINED INDEX FOR LM (L=1,2...25)
c     ORDER:  S, P_X, P_Y, P_Z, D_XY, D_XZ, D_Z^2, D_YZ, D_X^2-Y^2  ....
c     THE REAL SPHERICAL HARMONICS USED HERE FORM BASES FOR THE
c     IRRIDUCBLE REPRESENTATIONS OF THE GROUP O
c
c     SEE WEISSBLUTH 'ATOMS AND MOLECULES' PAGES 128-130
c     ERRORS IN WEISSBLUTH HAVE BEEN CORRECTED:
c        1.) ELIMINATION OF THE 7'S FROM L=20
c        2.) ADDITION OF THE FACTOR 1./sqrt(12.) TO L=25
c
      implicit none
      integer lmax
      real*8 ylm(lmax), gx(3), g
      real*8 pi, fpi, eps, c
      integer l

      PI=4.D0*DATAN(1.D0)
      fpi=4.D0*PI
      eps=1e-9

      if (lmax.ge.26) call errore
     &               (' ylm_q',' not programmed for L>',L)
      if (lmax.le.0 .or. (lmax.ne.1 .and. lmax.ne.4 .and. lmax.ne.9
     &                              .and. lmax.ne.16.and. lmax.ne.25))
     &     call errore (' ylm_q',' wrong L^2 on input',1000)

c   note :   ylm(q=0) = 1/sqrt(fpi)  WHEN L=0  AND  = 0  WHEN L>0

        ylm(1) = sqrt(1./fpi)
      if (lmax .eq. 1) return

       if(g.lt.eps) then
           do l=2,lmax
             ylm(l) = 0.0
           enddo
        return
       endif

        c=sqrt(3./fpi)

c  p_x p_y p_z

          ylm(2) = c*gx(1)/sqrt(g)   !   X
          ylm(3) = c*gx(2)/sqrt(g)   !   Y
          ylm(4) = c*gx(3)/sqrt(g)   !   Z

      if (lmax .eq. 4) return

c d_xy d_xz d_yz

        c=sqrt(15./fpi)
    
          ylm(5) = c*gx(1)*gx(2)/g   !  X*Y
          ylm(6) = c*gx(1)*gx(3)/g   !  X*Z
          ylm(8) = c*gx(2)*gx(3)/g   !  Y*Z

        c=sqrt(5.0/fpi/4.0)

          ylm(7) = c*(3.*gx(3)**2/g-1.)  ! (3.*Z*Z-1.0)

        c=sqrt(15./fpi/4.)
    
          ylm(9) = c*(gx(1)**2-gx(2)**2)/g  ! X*X-Y*Y

      if (lmax .eq. 9) return

        c=sqrt(7./fpi)*5./2.

          ylm(10) = c*gx(1)*(gx(1)**2-0.6*g)/(g*sqrt(g)) ! X(X^2-3R^2/5)
          ylm(11) = c*gx(2)*(gx(2)**2-0.6*g)/(g*sqrt(g))

        c=sqrt(7.*15./fpi)

          ylm(12) = c*gx(1)*gx(2)*gx(3)/(g*sqrt(g))

        c=sqrt(7./fpi)*5./2.

          ylm(13) = c*gx(3)*(gx(3)**2-0.6*g)/(g*sqrt(g))

        c=sqrt(7.*15./fpi)/2.

          ylm(14) = c*gx(3)*(gx(1)**2-gx(2)**2)/(g*sqrt(g))
          ylm(15) = c*gx(2)*(gx(3)**2-gx(1)**2)/(g*sqrt(g))
          ylm(16) = c*gx(1)*(gx(2)**2-gx(3)**2)/(g*sqrt(g))

      if (lmax .eq. 16) return

        c=sqrt(3.*7./fpi)*5./4.
        
          ylm(17) = c*((gx(1)**4+gx(2)**4+gx(3)**4)/(g*g)-0.6)

        c=sqrt(9.*35./fpi)/2.
        
          ylm(18) = c*gx(2)*gx(3)*(gx(2)**2-gx(3)**2)/g**2
          ylm(19) = c*gx(1)*gx(3)*(gx(3)**2-gx(1)**2)/g**2

        c=sqrt(9.*5./fpi)/4.
        
          ylm(20) = c*((gx(1)**4-gx(2)**4)-
     +      6.*gx(3)**2*(gx(1)**2-gx(2)**2))/(g*g)

        c=sqrt(9.*35./fpi)/2.
        
          ylm(21) = c*gx(1)*gx(2)*(gx(1)**2-gx(2)**2)/g**2


        c=sqrt(9.*5./fpi)*7./2.
        
          ylm(22) = c*gx(1)*gx(2)*(gx(3)**2-g/7.)/g**2
          ylm(23) = c*gx(1)*gx(3)*(gx(2)**2-g/7.)/g**2
          ylm(24) = c*gx(2)*gx(3)*(gx(1)**2-g/7.)/g**2

        c=sqrt(9.*5./fpi/3.)*7./2.
        
          ylm(25) = c*( gx(3)**4-0.5*(gx(1)**4+gx(2)**4)-
     +    6./7.*g*(gx(3)**2-0.5*(gx(1)**2+gx(2)**2) ))/ g**2


      return
      end
