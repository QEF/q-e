!
!---------------------------------------------------------------
      function exc_t(rho,rhoc,lsd)
!---------------------------------------------------------------
!
      implicit none
      integer,parameter :: dp=kind(1.d0)
      integer:: lsd
      real(kind=dp) :: exc_t, rho(2),arho,rhot, zeta,rhoc
      real(kind=dp) :: ex, ec, vx(2), vc(2)

      real(kind=dp),parameter:: e2 =2.d0

      exc_t=0.0d0

      if(lsd.eq.0) then
!
!     LDA case
!
         rhot = rho(1) + rhoc
         arho = abs(rhot)
         if (arho.gt.1.d-30) then      
            call xc(arho,ex,ec,vx,vc)
            exc_t=e2*(ex+ec)
         endif
      else
!
!     LSDA case
!
         rhot = rho(1)+rho(2)+rhoc
         arho = abs(rhot)
         if (arho.gt.1.d-30) then      
            zeta = (rho(1)-rho(2)) / arho
            call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
            exc_t=e2*(ex+ec)
         endif
      endif

      exc_t=e2*(ex+ec)

      return
      end

