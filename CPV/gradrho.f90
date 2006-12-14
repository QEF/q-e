!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
                                                                                
#include "f_defs.h"

!
!----------------------------------------------------------------------
      subroutine gradrho(nspin,rhog,drho,d2rho,dxdyrho,dxdzrho,dydzrho)
!----------------------------------------------------------------------
!
!     calculates gradient of charge density for gradient corrections
!     in: charge density on G-space    out: gradient in R-space
!
      use cell_base
      use gvecp, only: ng => ngm
      use reciprocal_vectors
      use recvecs_indexes
      USE cp_interfaces, ONLY: fwfft, invfft
      use grid_dimensions, only : nr1, nr2, nr3, nr1x, nr2x, nr3x,      &
     &                            nnr=> nnrx
!      use grid_dimensions, only: nr1, nr2, nr3, &
!            nr1x, nr2x, nr3x, nnr => nnrx
!
      implicit none
! input
      integer nspin
      complex(kind=8) rhog(ng,nspin)
! output
      real(kind=8)    drho(3,nnr), d2rho(3,nnr),     &
     &                dxdyrho(nnr), dxdzrho(nnr),    &
     &                dydzrho(nnr)
! local
      complex(kind=8), allocatable:: v(:), w(:)
      complex(kind=8) ci
      integer iss, ig, ir, j
!
!
      allocate(v(nnr))
      allocate(w(nnr))
      ci=(0.0d0,1.0d0)
      do ir = 1,nnr
         do j = 1,3
            drho(j,ir) = 0.d0
            d2rho(j,ir) = 0.d0
         end do
         dxdyrho(ir) = 0.d0
         dxdzrho(ir) = 0.d0
         dydzrho(ir) = 0.d0
      end do
      do iss=1,nspin

         do ig=1,nnr
            v(ig)=(0.0d0,0.0d0)
            w(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ng
            v(np(ig))=      ci*tpiba*gx(1,ig)*rhog(ig,iss)
            v(nm(ig))=conjg(ci*tpiba*gx(1,ig)*rhog(ig,iss))
            w(np(ig))=      -1.d0*tpiba**2*gx(1,ig)**2*rhog(ig,iss)
            w(nm(ig))=conjg(-1.d0*tpiba**2*gx(1,ig)**2*rhog(ig,iss))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         call invfft('Dense',w,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            drho(1,ir)=drho(1,ir)+real(v(ir))
            d2rho(1,ir)=d2rho(1,ir)+real(w(ir))
         end do
!
         do ig=1,nnr
            v(ig)=(0.0d0,0.0d0)
            w(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ng
            v(np(ig))= tpiba*(      ci*gx(2,ig)*rhog(ig,iss)-           &
     &                                 gx(3,ig)*rhog(ig,iss) )
            v(nm(ig))= tpiba*(conjg(ci*gx(2,ig)*rhog(ig,iss))+          &
     &                              ci*conjg(ci*gx(3,ig)*rhog(ig,iss)))
            w(np(ig))= -1.d0*tpiba**2*(      gx(2,ig)**2*rhog(ig,iss) + &
     &                                 ci*gx(3,ig)**2*rhog(ig,iss) )
            w(nm(ig))= -1.d0*tpiba**2*(conjg(gx(2,ig)**2*rhog(ig,iss))+ &
     &                           ci*conjg(gx(3,ig)**2*rhog(ig,iss)))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         call invfft('Dense',w,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            drho(2,ir)=drho(2,ir)+real(v(ir))
            drho(3,ir)=drho(3,ir)+aimag(v(ir))
            d2rho(2,ir)=d2rho(2,ir)+real(w(ir))
            d2rho(3,ir)=d2rho(3,ir)+aimag(w(ir))
         end do

         do ig=1,nnr
            v(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ng
            v(np(ig))= -1.d0*tpiba**2*gx(1,ig)*gx(2,ig)*rhog(ig,iss)
            v(nm(ig))=conjg(v(np(ig)))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            dxdyrho(ir)=dxdyrho(ir)+real(v(ir))
         end do
!
         do ig=1,nnr
            v(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ng
            v(np(ig))= -1.d0*tpiba**2*(gx(1,ig)*gx(3,ig)*rhog(ig,iss) + &
     &                              ci*gx(2,ig)*gx(3,ig)*rhog(ig,iss) )
            v(nm(ig))= -1.d0*tpiba**2*                                  &
     &                          (conjg(gx(1,ig)*gx(3,ig)*rhog(ig,iss))+ &
     &                        ci*conjg(gx(2,ig)*gx(3,ig)*rhog(ig,iss)))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            dxdzrho(ir)=dxdzrho(ir)+real(v(ir))
            dydzrho(ir)=dydzrho(ir)+aimag(v(ir))
         end do

      end do
      deallocate(v)
      deallocate(w)

      return
      end

