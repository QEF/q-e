!     Author: Burak Himmetoglu (burakhmmtgl@gmail.com)
!     This module contains smearing methods.
!      
!     Methods:
!     w0gauss : Standard Gaussian smearing
!     w1gauss : FD smearing for n=0, Integral of Gaussian if n=1
!     sig_nk  : Variable smearing for Dirac-delta involving e_{nk} - e_{m,k+q}
!     sig0    : Variable smearing for Dirac-delta involving e_{nk} 

      module smearing_mod

      implicit none
      double precision, parameter :: sqrtpm1 = 1.0d0/1.7724538509055160d0

      contains

         double precision function w0gauss(x)
            implicit none
            double precision, intent(in) :: x
            double precision :: arg

            arg = min(200.d0,x**2.d0)
            w0gauss = sqrtpm1 * exp( -arg )

         end function w0gauss
         !
         !
         double precision function w1gauss(x,ngauss)
            implicit none
            double precision, intent(in) :: x
            integer, intent(in) :: ngauss

            ! ngauss=0 (gaussian smearing)
            ! ngauss=1 FD

            if ( ngauss .eq. 0 ) then
               w1gauss = 0.5d0 * ( 1.d0 + derf( x * sqrtpm1 ) )
            else if ( ngauss .eq. 1 ) then
               w1gauss = 1.d0 / ( exp(-x) + 1.d0)
            end if

         end function w1gauss
         !
         !
         double precision function sig_nk(nk1,nk2,nk3,vk,vkq,aa)
            implicit none
            integer, intent(in) :: nk1, nk2, nk3
            double precision, intent(in) :: vk(3), vkq(3), aa
            double precision :: dF1, dF2, dF3!, dFarr(4)

            dF1 = ( vk(1) - vkq(1) ) * 1.0/nk1
            dF2 = ( vk(2) - vkq(2) ) * 1.0/nk2
            dF3 = ( vk(3) - vkq(3) ) * 1.0/nk3
            !
            !dFarr(1) = abs(dF1+dF2+dF3)
            !dFarr(2) = abs(-dF1+dF2+dF3)
            !dFarr(3) = abs(dF1+dF2-dF3)
            !dFarr(4) = abs(dF1-dF2+dF3)
            !
            !sig_nk = aa * maxval(dFarr)
            sig_nk = aa * sqrt(abs(dF1**2+dF2**2+dF3**2))

         end function sig_nk
         !
         !
         double precision function sig0(nk1,nk2,nk3,vk,aa)
            implicit none
            integer, intent(in) :: nk1, nk2, nk3
            double precision, intent(in) :: vk(3), aa
            double precision :: dF1, dF2, dF3!, dFarr(4)

            dF1 = vk(1) * 1.0/nk1
            dF2 = vk(2) * 1.0/nk2
            dF3 = vk(3) * 1.0/nk3
            !
            !dFarr(1) = abs(dF1+dF2+dF3)
            !dFarr(2) = abs(-dF1+dF2+dF3)
            !dFarr(3) = abs(dF1+dF2-dF3)
            !dFarr(4) = abs(dF1-dF2+dF3)
            !
            !sig0 = aa * maxval(dFarr)
            sig0 = aa * sqrt(abs(dF1**2+dF2**2+dF3**2))

         end function sig0

      end module smearing_mod
