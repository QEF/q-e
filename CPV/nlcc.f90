!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!=----------------------------------------------------------------------------=!
   module non_local_core_correction
!=----------------------------------------------------------------------------=!

     USE kinds

     IMPLICIT NONE
     SAVE

     PRIVATE

     PUBLIC :: add_core_charge,  core_charge_forces

!=----------------------------------------------------------------------------=!
   contains
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   subroutine add_core_charge( rhoetg, rhoetr, sfac, ps, nsp)
!=----------------------------------------------------------------------------=!

     USE fft, ONLY: pinvfft
     USE cp_types, ONLY: pseudo
     use electrons_base, only: nspin
     use gvecp, only: ngm

     implicit none

     integer :: nsp
     COMPLEX(dbl) :: rhoetg(:)
     REAL(dbl)    :: rhoetr(:,:,:)
     type (pseudo),  intent(in) :: ps
     COMPLEX(dbl), INTENT(IN) :: sfac(:,:)
          
     COMPLEX(dbl), ALLOCATABLE :: vtemp(:)
     REAL(dbl) :: fac
     integer :: is, ig

     ALLOCATE( vtemp( ngm ) )
     vtemp = CMPLX( 0.0d0 )

     fac = 1.0d0 / DBLE( nspin )
     DO is = 1, nsp
       if( ps%tnlcc( is ) ) then
         do ig = 1, ngm
           vtemp(ig) = vtemp(ig) + fac * sfac( ig, is ) * CMPLX(ps%rhoc1(ig,is),0.0d0)
         end do
       endif
     end do

     
     rhoetg( 1:ngm ) = rhoetg( 1:ngm ) + vtemp( 1:ngm )

     !  rhoetr = 1.0d0 * rhoetr + INVFFT( vtemp )
     CALL pinvfft( rhoetr, vtemp, 1.0d0 )

     DEALLOCATE( vtemp )

     RETURN
!=----------------------------------------------------------------------------=!
   end subroutine add_core_charge
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   subroutine core_charge_forces( fion, vxc, rhoc1, tnlcc, atoms, ht, ei1, ei2, ei3, kp )
!=----------------------------------------------------------------------------=!

     !   This subroutine computes the non local core correction
     !   contribution to the atomic forces

     USE constants
     USE cell_base, ONLY: tpiba
     USE cell_module, ONLY: boxdimensions
     USE brillouin, ONLY: kpoints
     USE atoms_type_module, ONLY: atoms_type
     USE grid_dimensions, ONLY: nr1, nr2, nr3
     USE reciprocal_vectors, ONLY: mill_l, gstart, gx
     USE gvecp, ONLY: ngm
     USE ions_base, ONLY: nat

     IMPLICIT NONE

     TYPE (atoms_type), INTENT(IN) :: atoms    !   atomic positions
     TYPE (boxdimensions), INTENT(IN) :: ht    !   cell parameters
     TYPE (kpoints), INTENT(IN) :: kp          !   k points
     COMPLEX(dbl) :: ei1( -nr1:nr1, nat)                  !   
     COMPLEX(dbl) :: ei2( -nr2:nr2, nat)                  !   
     COMPLEX(dbl) :: ei3( -nr3:nr3, nat)                  !   
     LOGICAL      :: tnlcc(:)                  !   NLCC flags
     REAL(dbl)    :: fion(:,:)                 !   ionic forces
     REAL(dbl)    :: rhoc1(:,:)                !   derivative of the core charge
     COMPLEX(dbl) :: vxc(:,:)                  !   XC potential

     INTEGER :: ig, ig1, ig2, ig3, isa, ia, is, ispin, nspin
     COMPLEX(dbl) :: gxc, gyc, gzc, tx, ty, tz, teigr, cxc
     COMPLEX(dbl), ALLOCATABLE :: ftmp(:,:)
     REAL(dbl) :: cost

     IF( ANY( tnlcc ) ) then

       nspin = SIZE( vxc, 2)
       ALLOCATE( ftmp( 3, atoms%nat ) )
       ftmp = CMPLX( 0.0d0, 0.0d0 )

       DO ig = gstart, ngm
         ig1  = mill_l(1,ig)
         ig2  = mill_l(2,ig)
         ig3  = mill_l(3,ig)
         GXC  = CMPLX(0.D0, gx(1,ig))
         GYC  = CMPLX(0.D0, gx(2,ig))
         GZC  = CMPLX(0.D0, gx(3,ig))
         isa = 1
         DO is = 1, atoms%nsp
           IF ( tnlcc(is) ) THEN
             CXC = 0.0_dbl
             DO ispin = 1, nspin
               CXC = CXC + rhoc1( ig, is ) * CONJG( vxc( ig, ispin ) )
             END DO
             TX = CXC * GXC / DBLE( nspin )
             TY = CXC * GYC / DBLE( nspin )
             TZ = CXC * GZC / DBLE( nspin )
             DO IA = 1, atoms%na(is)
               teigr = ei1( ig1, isa ) * ei2( ig2, isa ) * ei3( ig3, isa )
               ftmp( 1, isa ) = ftmp( 1, isa ) + TEIGR * TX
               ftmp( 2, isa ) = ftmp( 2, isa ) + TEIGR * TY
               ftmp( 3, isa ) = ftmp( 3, isa ) + TEIGR * TZ
               isa = isa + 1
             END DO
           ELSE
             isa = isa + atoms%na(is)
           END IF
         END DO
       END DO

       ! ...  each processor add its own contribution to the array FION
       IF(kp%gamma_only) THEN
         cost = 2.D0 * ht%deth * tpiba
       ELSE
         cost =        ht%deth * tpiba
       END IF

       DO isa = 1, atoms%nat
         FION(1,ISA) = FION(1,ISA) + REAL(ftmp(1,ISA)) * cost
         FION(2,ISA) = FION(2,ISA) + REAL(ftmp(2,ISA)) * cost
         FION(3,ISA) = FION(3,ISA) + REAL(ftmp(3,ISA)) * cost
       END DO

       DEALLOCATE( ftmp )

     END IF

     RETURN
!=----------------------------------------------------------------------------=!
   end subroutine core_charge_forces
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   end module non_local_core_correction
!=----------------------------------------------------------------------------=!




!-----------------------------------------------------------------------
      subroutine add_cc(rhoc,rhog,rhor)
!-----------------------------------------------------------------------
!
! add core correction to the charge density for exch-corr calculation
!
      use electrons_base,  only: nspin
      use recvecs_indexes, only: np

      ! this isn't really needed, but if I remove it, ifc 7.1
      ! gives an "internal compiler error"
      use reciprocal_vectors, only: gstart
      use gvecp,              only: ngm

      use grid_dimensions,    only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
!
      implicit none
      real(kind=8), intent(in)   :: rhoc(nnr)
      real(kind=8), intent(inout):: rhor(nnr,nspin)
      complex(kind=8), intent(inout)::  rhog(ngm,nspin)
      complex(kind=8), allocatable :: wrk1( : )
!
      integer ig, ir, iss, isup, isdw
!
! In r-space:
!
      if (nspin.eq.1) then
         iss=1
         call DAXPY(nnr,1.d0,rhoc,1,rhor(1,iss),1)
      else
         isup=1
         isdw=2
         call DAXPY(nnr,0.5d0,rhoc,1,rhor(1,isup),1)
         call DAXPY(nnr,0.5d0,rhoc,1,rhor(1,isdw),1)
      end if 
! rhoc(r) -> rhoc(g)  (wrk1 is used as work space)
      allocate( wrk1( nnr ) )
      do ir=1,nnr
         wrk1(ir)=rhoc(ir)
      end do
      call fwfft(wrk1,nr1,nr2,nr3,nr1x,nr2x,nr3x)
! In g-space:
      if (nspin.eq.1) then
         do ig=1,ngm
            rhog(ig,iss)=rhog(ig,iss)+wrk1(np(ig))
         end do
      else
         do ig=1,ngm
            rhog(ig,isup)=rhog(ig,isup)+0.5d0*wrk1(np(ig))
            rhog(ig,isdw)=rhog(ig,isdw)+0.5d0*wrk1(np(ig))
         end do
      end if

      deallocate( wrk1 )
!
      return
      end subroutine add_cc


!
!-----------------------------------------------------------------------
      subroutine force_cc(irb,eigrb,vxc,fion1)
!-----------------------------------------------------------------------
!
!     core correction force: f = \int V_xc(r) (d rhoc(r)/d R_i) dr
!     same logic as in newd - uses box grid. For parallel execution:
!     the sum over node contributions is done in the calling routine
!
      use core,            only: rhocb
      use electrons_base,  only: nspin
      use gvecb,           only: gxb, ngb, npb, nmb
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx
      use cell_base,       only: omega
      use ions_base,       only: nsp, na, nat
      use parameters,      only: natx, nsx
      use small_box,       only: tpibab
      use atom,            only: nlcc
      use reciprocal_vectors, only: gstart
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx

      implicit none

! input
      integer, intent(in)        :: irb(3,nat)
      complex(kind=8), intent(in):: eigrb(ngb,nat)
      real(kind=8), intent(in)   :: vxc(nnr,nspin)
! output
      real(kind=8), intent(inout):: fion1(3,natx)
! local
      integer iss, ix, ig, is, ia, nfft, irb3, imin3, imax3, isa
      real(kind=8) fcc(3,natx), fac, boxdotgrid
      complex(kind=8) ci, facg
      complex(kind=8), allocatable :: qv(:)
      external  boxdotgrid
!
      call start_clock( 'forcecc' )
      ci = (0.d0,1.d0)
      fac = omega/dble(nr1*nr2*nr3*nspin)
      fcc = 0.d0

      allocate( qv( nnrb ) )

      isa = 0

      do is=1,nsp
         if( .not. nlcc(is) ) go to 10
#ifdef __PARA
         do ia=1,na(is)
            nfft=1
            irb3=irb(3,ia+isa)
            call parabox(nr3b,irb3,nr3,imin3,imax3)
            if (imax3-imin3+1.le.0) go to 15
#else
         do ia=1,na(is),2
!
! two fft's on two atoms at the same time (when possible)
!
            nfft=2
            if(ia.eq.na(is)) nfft=1
#endif
            do ix=1,3
               qv(:) = (0.d0, 0.d0)
               if (nfft.eq.2) then
                  do ig=1,ngb
                     facg = tpibab*cmplx(0.d0,gxb(ix,ig))*rhocb(ig,is)
                     qv(npb(ig)) = eigrb(ig,ia+isa  )*facg                 &
     &                      + ci * eigrb(ig,ia+isa+1)*facg
                     qv(nmb(ig)) = conjg(eigrb(ig,ia+isa  )*facg)          &
     &                      + ci * conjg(eigrb(ig,ia+isa+1)*facg)
                  end do
               else
                  do ig=1,ngb
                     facg = tpibab*cmplx(0.d0,gxb(ix,ig))*rhocb(ig,is)
                     qv(npb(ig)) = eigrb(ig,ia+isa)*facg
                     qv(nmb(ig)) = conjg(eigrb(ig,ia+isa)*facg)
                  end do
               end if
!
               call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
! note that a factor 1/2 is hidden in fac if nspin=2
!
               do iss=1,nspin
                  fcc(ix,ia+isa) = fcc(ix,ia+isa) + fac *               &
     &                 boxdotgrid(irb(1,ia  +isa),1,qv,vxc(1,iss))
                  if (nfft.eq.2)                                         &
     &               fcc(ix,ia+1+isa) = fcc(ix,ia+1+isa) + fac *           &
     &                    boxdotgrid(irb(1,ia+1+isa),2,qv,vxc(1,iss))
               end do
            end do
15          continue
         end do
10       continue
         isa = isa + na(is)
      end do
!
      do ia = 1, nat
        fion1(:,ia) = fion1(:,ia) + fcc(:,ia)
      end do

      deallocate( qv )
!
      call stop_clock( 'forcecc' )
      return
      end subroutine force_cc


!
!-----------------------------------------------------------------------
      subroutine set_cc(irb,eigrb,rhoc)
!-----------------------------------------------------------------------
!
!     Calculate core charge contribution in real space, rhoc(r)
!     Same logic as for rhov: use box grid for core charges
!
      use ions_base,       only: nsp, na, nat
      use parameters,      only: natx, nsx
      use atom,            only: nlcc
      use grid_dimensions, only: nr3, nnr => nnrx
      use gvecb,           only: ngb, npb, nmb
      use control_flags,   only: iprint
      use core,            only: rhocb
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx

      implicit none
! input
      integer, intent(in)        :: irb(3,nat)
      complex(kind=8), intent(in):: eigrb(ngb,nat)
! output
      real(kind=8), intent(out)  :: rhoc(nnr)
! local
      integer nfft, ig, is, ia, irb3, imin3, imax3, isa
      complex(kind=8) ci
      complex(kind=8), allocatable :: wrk1(:)
      complex(kind=8), allocatable :: qv(:)
!
      call start_clock( 'set_cc' )
      ci=(0.,1.)
!
      allocate( qv ( nnrb ) )
      allocate( wrk1 ( nnr ) )
      wrk1 (:) = (0.d0, 0.d0)
!
      isa = 0
      do is=1,nsp
         if (.not.nlcc(is)) go to 10
#ifdef __PARA
         do ia=1,na(is)
            nfft=1
            irb3=irb(3,ia+isa)
            call parabox(nr3b,irb3,nr3,imin3,imax3)
            if (imax3-imin3+1.le.0) go to 15
#else
         do ia=1,na(is),2
            nfft=2
            if( ia.eq.na(is) ) nfft=1
!
! two ffts at the same time, on two atoms (if possible: nfft=2)
!
#endif
            qv(:) = (0.d0, 0.d0)
            if(nfft.eq.2)then
               do ig=1,ngb
                  qv(npb(ig))= eigrb(ig,ia  +isa)*rhocb(ig,is)          &
     &                    + ci*eigrb(ig,ia+1+isa)*rhocb(ig,is)
                  qv(nmb(ig))= conjg(eigrb(ig,ia  +isa)*rhocb(ig,is))   &
     &                    + ci*conjg(eigrb(ig,ia+1+isa)*rhocb(ig,is))
               end do
            else
               do ig=1,ngb
                  qv(npb(ig)) = eigrb(ig,ia+isa)*rhocb(ig,is)
                  qv(nmb(ig)) = conjg(eigrb(ig,ia+isa)*rhocb(ig,is))
               end do
            endif
!
            call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
            call box2grid(irb(1,ia+isa),1,qv,wrk1)
            if (nfft.eq.2) call box2grid(irb(1,ia+1+isa),2,qv,wrk1)
!
15          continue
         end do
10       continue
         isa = isa + na(is)
      end do
!
      call DCOPY(nnr,wrk1,2,rhoc,1)

      deallocate( qv  )
      deallocate( wrk1 )
!
      call stop_clock( 'set_cc' )
!
      return
      end subroutine set_cc

