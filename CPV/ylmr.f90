!-----------------------------------------------------------------------
      real(kind=8) function ylmr(l,ig)
!-----------------------------------------------------------------------
!     real spherical harmonics,  l is combined index for lm  (l=1,2...9)
!     order:  s, p_x, p_z, p_y, d_xy, d_xz, d_z^2, d_yz, d_x^2-y^2
!
      use gvec
      use constants, only: pi, fpi
!
      implicit none
      integer l, ig
      real(kind=8) x,y,z,r
!
      if (ig.gt.ng) call errore(' ylmr ',' ig.gt.ng ',ig)
      x = gx(ig,1)
      y = gx(ig,2)
      z = gx(ig,3)
!
!     yml is undefined when  g=0 and l>0
!
      r = max(sqrt(x*x+y*y+z*z),1.0d-6)
      x = x/r
      y = y/r
      z = z/r
!
!     only l=1 is ok also when  g=0
!
      if (l.eq.1) ylmr = sqrt(1.0/fpi)
      if (l.eq.2) ylmr = sqrt(3.0/fpi)*x
      if (l.eq.3) ylmr = sqrt(3.0/fpi)*z
      if (l.eq.4) ylmr = sqrt(3.0/fpi)*y
      if (l.eq.5) ylmr = sqrt(15.0/fpi)*x*y
      if (l.eq.6) ylmr = sqrt(15.0/fpi)*x*z
      if (l.eq.7) ylmr = sqrt(5.0/fpi/4.0)*(3.0*z*z-1.0)
      if (l.eq.8) ylmr = sqrt(15.0/fpi)*y*z
      if (l.eq.9) ylmr = sqrt(15.0/fpi/4.0)*(x*x-y*y)
      if (l.ge.10) call errore(' ylmr',' higher l not programmed  l=',l)
      return
      end

