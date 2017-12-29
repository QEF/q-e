!
! Copyright (C) 2002-2020 Quantum ESPRESSO grouo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
                                                                                

!
!----------------------------------------------------------------------
SUBROUTINE gradrho(nspin,rhog,drho,d2rho,dxdyrho,dxdzrho,dydzrho)
!----------------------------------------------------------------------
!
!     calculates gradient of charge density for gradient corrections
!     in: charge density on G-space    out: gradient in R-space
!
      use cell_base
      use gvect, only: ngm, g
      USE fft_interfaces, ONLY: invfft
      USE fft_base,       ONLY: dfftp
      USE fft_helper_subroutines, ONLY: fftx_oned2threed_gamma

!
      implicit none
! input
      integer nspin
      complex(kind=8) rhog(ngm,nspin)
! output
      real(kind=8)    drho(3,dfftp%nnr), d2rho(3,dfftp%nnr),     &
     &                dxdyrho(dfftp%nnr), dxdzrho(dfftp%nnr),    &
     &                dydzrho(dfftp%nnr)
! local
      complex(kind=8), allocatable:: v(:), drhog(:,:)
      complex(kind=8) ci
      integer iss, ig, ir, j
!
!
      allocate(v(dfftp%nnr))
      allocate(drhog(dfftp%ngm,3))

      ci=(0.0d0,1.0d0)

      drho = 0.d0
      d2rho = 0.d0
      dxdyrho = 0.d0
      dxdzrho = 0.d0
      dydzrho = 0.d0

      do iss=1,nspin
         do ig=1,ngm
            drhog(ig,1) = ci*tpiba*g(1,ig)*rhog(ig,iss)
            drhog(ig,2) = ci*tpiba*g(2,ig)*rhog(ig,iss)
            drhog(ig,3) = ci*tpiba*g(3,ig)*rhog(ig,iss)
         enddo
         CALL fftx_oned2threed_gamma(dfftp, v, drhog(:,1) )
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            drho(1,ir)=drho(1,ir)+real(v(ir))
         end do
         CALL fftx_oned2threed_gamma(dfftp, v, drhog(:,2), drhog(:,3) )
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            drho(2,ir)=drho(2,ir)+real(v(ir))
            drho(3,ir)=drho(3,ir)+aimag(v(ir))
         end do

         do ig=1,ngm
            drhog(ig,1) = -1.d0*tpiba**2*g(1,ig)**2*rhog(ig,iss)
            drhog(ig,2) = -1.d0*tpiba**2*g(2,ig)**2*rhog(ig,iss)
            drhog(ig,3) = -1.d0*tpiba**2*g(3,ig)**2*rhog(ig,iss)
         enddo
         CALL fftx_oned2threed_gamma(dfftp, v, drhog(:,1) )
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            d2rho(1,ir)=d2rho(1,ir)+real(v(ir))
         end do
         CALL fftx_oned2threed_gamma(dfftp, v, drhog(:,2), drhog(:,3) )
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            d2rho(2,ir)=d2rho(2,ir)+real(v(ir))
            d2rho(3,ir)=d2rho(3,ir)+aimag(v(ir))
         end do

         do ig=1,ngm
            drhog(ig,1) = -1.d0*tpiba**2*g(1,ig)*g(2,ig)*rhog(ig,iss)
            drhog(ig,2) = -1.d0*tpiba**2*g(1,ig)*g(3,ig)*rhog(ig,iss)
            drhog(ig,3) = -1.d0*tpiba**2*g(2,ig)*g(3,ig)*rhog(ig,iss)
         enddo
         CALL fftx_oned2threed_gamma(dfftp, v, drhog(:,1) )
         CALL invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            dxdyrho(ir)=dxdyrho(ir)+real(v(ir))
         end do
         CALL fftx_oned2threed_gamma(dfftp, v, drhog(:,2), drhog(:,3) )
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            dxdzrho(ir)=dxdzrho(ir)+real(v(ir))
            dydzrho(ir)=dydzrho(ir)+aimag(v(ir))
         end do

      end do
      deallocate(v)
      deallocate(drhog)

END SUBROUTINE gradrho

