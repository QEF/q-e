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
      use gvect, only: ngm, nl, nlm, g
      USE fft_interfaces, ONLY: invfft
      USE fft_base,       ONLY: dfftp
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
      complex(kind=8), allocatable:: v(:), w(:)
      complex(kind=8) ci
      integer iss, ig, ir, j
!
!
      allocate(v(dfftp%nnr))
      allocate(w(dfftp%nnr))
      ci=(0.0d0,1.0d0)
      do ir = 1,dfftp%nnr
         do j = 1,3
            drho(j,ir) = 0.d0
            d2rho(j,ir) = 0.d0
         end do
         dxdyrho(ir) = 0.d0
         dxdzrho(ir) = 0.d0
         dydzrho(ir) = 0.d0
      end do
      do iss=1,nspin

         do ig=1,dfftp%nnr
            v(ig)=(0.0d0,0.0d0)
            w(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ngm
            v(nl(ig)) =      ci*tpiba*g(1,ig)*rhog(ig,iss)
            v(nlm(ig))=conjg(ci*tpiba*g(1,ig)*rhog(ig,iss))
            w(nl(ig)) =      -1.d0*tpiba**2*g(1,ig)**2*rhog(ig,iss)
            w(nlm(ig))=conjg(-1.d0*tpiba**2*g(1,ig)**2*rhog(ig,iss))
         end do
         call invfft('Dense',v, dfftp )
         call invfft('Dense',w, dfftp )
         do ir=1,dfftp%nnr
            drho(1,ir)=drho(1,ir)+real(v(ir))
            d2rho(1,ir)=d2rho(1,ir)+real(w(ir))
         end do
!
         do ig=1,dfftp%nnr
            v(ig)=(0.0d0,0.0d0)
            w(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ngm
            v(nl(ig))= tpiba*(      ci*g(2,ig)*rhog(ig,iss)-           &
     &                                 g(3,ig)*rhog(ig,iss) )
            v(nlm(ig))=tpiba*(conjg(ci*g(2,ig)*rhog(ig,iss))+          &
     &                              ci*conjg(ci*g(3,ig)*rhog(ig,iss)))
            w(nl(ig))= -1.d0*tpiba**2*(      g(2,ig)**2*rhog(ig,iss) + &
     &                                 ci*g(3,ig)**2*rhog(ig,iss) )
            w(nlm(ig))=-1.d0*tpiba**2*(conjg(g(2,ig)**2*rhog(ig,iss))+ &
     &                           ci*conjg(g(3,ig)**2*rhog(ig,iss)))
         end do
         call invfft('Dense',v, dfftp )
         call invfft('Dense',w, dfftp )
         do ir=1,dfftp%nnr
            drho(2,ir)=drho(2,ir)+real(v(ir))
            drho(3,ir)=drho(3,ir)+aimag(v(ir))
            d2rho(2,ir)=d2rho(2,ir)+real(w(ir))
            d2rho(3,ir)=d2rho(3,ir)+aimag(w(ir))
         end do

         do ig=1,dfftp%nnr
            v(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ngm
            v(nl(ig)) = -1.d0*tpiba**2*g(1,ig)*g(2,ig)*rhog(ig,iss)
            v(nlm(ig))=conjg(v(nl(ig)))
         end do
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            dxdyrho(ir)=dxdyrho(ir)+real(v(ir))
         end do
!
         do ig=1,dfftp%nnr
            v(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ngm
            v(nl(ig))= -1.d0*tpiba**2*(g(1,ig)*g(3,ig)*rhog(ig,iss) + &
     &                              ci*g(2,ig)*g(3,ig)*rhog(ig,iss) )
            v(nlm(ig))=-1.d0*tpiba**2*                                  &
     &                          (conjg(g(1,ig)*g(3,ig)*rhog(ig,iss))+ &
     &                        ci*conjg(g(2,ig)*g(3,ig)*rhog(ig,iss)))
         end do
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
            dxdzrho(ir)=dxdzrho(ir)+real(v(ir))
            dydzrho(ir)=dydzrho(ir)+aimag(v(ir))
         end do

      end do
      deallocate(v)
      deallocate(w)

END SUBROUTINE gradrho

