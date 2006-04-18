!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!=----------------------------------------------------------------------------=!
  module core
!=----------------------------------------------------------------------------=!

    USE kinds
    USE splines,    ONLY: spline_data

  implicit none
  save
  !     nlcc_any = 0 no core correction on any atom
  !     rhocb  = core charge in G space (box grid)
  !     rhoc   = core charge in real space  (dense grid)
  !     rhocg  = core charge in G space  (dense grid)
  !     drhocg = derivative of core charge in G space (used for stress)
  logical :: nlcc_any
  real(8), allocatable:: rhocb(:,:)
  real(8), allocatable:: rhoc(:)
  real(8), allocatable:: rhocg(:,:)
  real(8), allocatable:: drhocg(:,:)

!=----------------------------------------------------------------------------=!
  contains
!=----------------------------------------------------------------------------=!

  subroutine allocate_core( nnrx, ngm, ngb, nsp )
     integer, intent(in) :: nnrx, ngm, ngb, nsp
     IF ( nlcc_any ) THEN
        !
        ALLOCATE( rhoc( nnrx ) )
        ALLOCATE( rhocb( ngb, nsp ) )
        ALLOCATE( rhocg( ngm, nsp ) )
        ALLOCATE( drhocg( ngm, nsp ) )
        !
     ELSE
        !
        ! ... dummy allocation required because this array appears in the
        ! ... list of arguments of some routines
        !
        ALLOCATE( rhoc( 1 ) )
        !
     END IF
  end subroutine allocate_core

  subroutine deallocate_core()
      IF( ALLOCATED( rhocb  ) ) DEALLOCATE( rhocb )
      IF( ALLOCATED( rhoc   ) ) DEALLOCATE( rhoc  )
      IF( ALLOCATED( rhocg  ) ) DEALLOCATE( rhocg  )
      IF( ALLOCATED( drhocg ) ) DEALLOCATE( drhocg )
  end subroutine deallocate_core

!=----------------------------------------------------------------------------=!
   subroutine core_charge_ftr( tpre )
!=----------------------------------------------------------------------------=!
     !
     !  Compute the fourier trasform of the core charge, from the radial
     !  mesh to the reciprocal space
     !
     use control_flags,   ONLY : program_name
     use ions_base,       ONLY : nsp
     use uspp_param,      ONLY : kkbeta
     use atom,            ONLY : nlcc, r, rab, mesh, rho_atc
     use gvecb,           ONLY : ngb, gb
     use small_box,       ONLY : omegab, tpibab
     use pseudo_base,     ONLY : compute_rhocg
     use pseudopotential, ONLY : tpstab, build_cctab, rhoc1_sp, rhocp_sp, chkpstab
     use cell_base,       ONLY : omega, tpiba2, tpiba
     USE splines,         ONLY : spline
     use reciprocal_vectors, ONLY : ngm, g, gstart
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(IN) :: tpre
     !
     INTEGER :: is, ig
     REAL(DP) :: xg, cost1
     !
     !
     IF( .NOT. nlcc_any ) RETURN
     !
     IF( tpstab ) THEN
        !
        CALL build_cctab( )
        !
     END IF
     !
     do is = 1, nsp
        !
        if( nlcc( is ) ) then
           !
           IF( program_name == 'CP90' ) THEN
              !
              CALL compute_rhocg( rhocb(:,is), rhocb(:,is), r(:,is), rab(:,is), &
                  rho_atc(:,is), gb, omegab, tpibab**2, kkbeta(is), ngb, 0 )
              !
           END IF
           !
           IF( ( program_name == 'FPMD' ) .OR. tpre ) THEN
              !
              IF( tpstab ) THEN
                 !
                 cost1 = 1.0d0/omega
                 !
                 IF( gstart == 2 ) THEN
                    rhocg (1,is) = rhoc1_sp(is)%y( 1 ) * cost1
                    drhocg(1,is) = 0.0d0
                 END IF
                 DO ig = gstart, SIZE( rhocg, 1 )
                    xg = SQRT( g(ig) ) * tpiba
                    rhocg (ig,is) = spline( rhoc1_sp(is), xg ) * cost1
                    drhocg(ig,is) = spline( rhocp_sp(is), xg ) * cost1
                 END DO
                 !
              ELSE

                 CALL compute_rhocg( rhocg(:,is), drhocg(:,is), r(:,is), rab(:,is), &
                                     rho_atc(:,is), g, omega, tpiba2, kkbeta(is), ngm, 1 )

              END IF
              !
           END IF
           !
        endif
        !
     end do

     return
   end subroutine core_charge_ftr
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   subroutine add_core_charge( rhoetg, rhoetr, sfac, rhoc, nsp)
!=----------------------------------------------------------------------------=!
     USE kinds,          ONLY: DP
     USE fft_base,       ONLY: dfftp
     use electrons_base, only: nspin
     use gvecp,          only: ngm
     use atom,           only: nlcc
     USE fft_module,     ONLY: invfft
     USE io_global,      ONLY: stdout
     USE mp_global,      ONLY: intra_image_comm
     USE cell_base,      ONLY: omega
     USE mp,             ONLY: mp_sum
     USE control_flags,  ONLY: iprsta

     implicit none

     integer :: nsp
     COMPLEX(DP) :: rhoetg(:)
     REAL(DP)    :: rhoetr(:)
     REAL(DP)    :: rhoc(:,:)
     COMPLEX(DP), INTENT(IN) :: sfac(:,:)
          
     COMPLEX(DP), ALLOCATABLE :: vtemp(:), psi(:)
     REAL(DP) :: fac
     REAL(DP) :: rsum
     integer :: is, ig

     ALLOCATE( vtemp( ngm ), psi( dfftp%nnr ) )

     vtemp = CMPLX( 0.0d0, 0.0d0 )

     fac = 1.0d0 / DBLE( nspin )
     DO is = 1, nsp
       if( nlcc( is ) ) then
         do ig = 1, ngm
           vtemp(ig) = vtemp(ig) + fac * sfac( ig, is ) * CMPLX(rhoc(ig,is),0.0d0)
         end do
       endif
     end do
     
     rhoetg( 1:ngm ) = rhoetg( 1:ngm ) + vtemp( 1:ngm )

     CALL rho2psi( 'Dense', psi, dfftp%nnr, vtemp, ngm )
     CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

     IF( SIZE( rhoetr ) /= SIZE( psi ) ) &
        CALL errore( " add_core_charge ", " inconsistent sizes ", 1 )

     IF( iprsta > 2 ) THEN
        rsum = DBLE( SUM( psi ) ) * omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
        CALL mp_sum( rsum, intra_image_comm )
        WRITE( stdout, 10 ) rsum 
10      FORMAT( 3X, 'Core Charge = ', D14.6 )
     END IF


     rhoetr(:) = rhoetr(:) + DBLE( psi )


     DEALLOCATE( vtemp, psi )

     RETURN
   end subroutine add_core_charge
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   subroutine core_charge_forces( fion, vxc, rhoc1, tnlcc, atoms, ht, ei1, ei2, ei3 )
!=----------------------------------------------------------------------------=!

     !   This subroutine computes the non local core correction
     !   contribution to the atomic forces

     USE kinds, ONLY : DP
     USE cell_base, ONLY: tpiba
     USE cell_module, ONLY: boxdimensions
     USE brillouin, ONLY: kpoints, kp
     USE atoms_type_module, ONLY: atoms_type
     USE grid_dimensions, ONLY: nr1, nr2, nr3
     USE reciprocal_vectors, ONLY: mill_l, gstart, gx, ngm
     USE ions_base, ONLY: nat

     IMPLICIT NONE

     TYPE (atoms_type), INTENT(IN) :: atoms    !   atomic positions
     TYPE (boxdimensions), INTENT(IN) :: ht    !   cell parameters
     COMPLEX(DP) :: ei1( -nr1:nr1, nat)                  !   
     COMPLEX(DP) :: ei2( -nr2:nr2, nat)                  !   
     COMPLEX(DP) :: ei3( -nr3:nr3, nat)                  !   
     LOGICAL      :: tnlcc(:)                  !   NLCC flags
     REAL(DP)    :: fion(:,:)                 !   ionic forces
     REAL(DP)    :: rhoc1(:,:)                !   derivative of the core charge
     COMPLEX(DP) :: vxc(:,:)                  !   XC potential

     INTEGER :: ig, ig1, ig2, ig3, isa, ia, is, ispin, nspin
     COMPLEX(DP) :: gxc, gyc, gzc, tx, ty, tz, teigr, cxc
     COMPLEX(DP), ALLOCATABLE :: ftmp(:,:)
     REAL(DP) :: cost

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
             CXC = 0.0_DP
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
         FION(1,ISA) = FION(1,ISA) + DBLE(ftmp(1,ISA)) * cost
         FION(2,ISA) = FION(2,ISA) + DBLE(ftmp(2,ISA)) * cost
         FION(3,ISA) = FION(3,ISA) + DBLE(ftmp(3,ISA)) * cost
       END DO

       DEALLOCATE( ftmp )

     END IF

     RETURN
!=----------------------------------------------------------------------------=!
   end subroutine core_charge_forces
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   end module core
!=----------------------------------------------------------------------------=!




!-----------------------------------------------------------------------
      subroutine add_cc( rhoc, rhog, rhor )
!-----------------------------------------------------------------------
!
! add core correction to the charge density for exch-corr calculation
!
      USE kinds,              ONLY: DP
      use electrons_base,     only: nspin
      use control_flags,      only: iprsta
      use io_global,          only: stdout
      use mp_global,          only: intra_image_comm
      use cell_base,          only: omega
      use recvecs_indexes,    only: np
      USE mp,                 ONLY: mp_sum

      ! this isn't really needed, but if I remove it, ifc 7.1
      ! gives an "internal compiler error"
      use reciprocal_vectors, only: gstart
      use gvecp,              only: ngm
      use grid_dimensions,    only: nr1, nr2, nr3, &
                                    nr1x, nr2x, nr3x, nnrx
      USE fft_module,         ONLY: fwfft
!
      implicit none
      !
      REAL(DP),    INTENT(IN)   :: rhoc( nnrx )
      REAL(DP),    INTENT(INOUT):: rhor( nnrx, nspin )
      COMPLEX(DP), INTENT(INOUT):: rhog( ngm,  nspin )
      !
      COMPLEX(DP), ALLOCATABLE :: wrk1( : )
!
      integer :: ig, ir, iss, isup, isdw
      REAL(DP) :: rsum
      !
      IF( iprsta > 2 ) THEN
         rsum = SUM( rhoc ) * omega / DBLE(nr1*nr2*nr3)
         CALL mp_sum( rsum, intra_image_comm )
         WRITE( stdout, 10 ) rsum 
10       FORMAT( 3X, 'Core Charge = ', D14.6 )
      END IF
      !
      ! In r-space:
      !
      if ( nspin .eq. 1 ) then
         iss=1
         call DAXPY(nnrx,1.d0,rhoc,1,rhor(1,iss),1)
      else
         isup=1
         isdw=2
         call DAXPY(nnrx,0.5d0,rhoc,1,rhor(1,isup),1)
         call DAXPY(nnrx,0.5d0,rhoc,1,rhor(1,isdw),1)
      end if 
      !
      ! rhoc(r) -> rhoc(g)  (wrk1 is used as work space)
      !
      allocate( wrk1( nnrx ) )

      wrk1(:) = rhoc(:)

      call fwfft('Dense',wrk1,nr1,nr2,nr3,nr1x,nr2x,nr3x)
      !
      ! In g-space:
      !
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
      USE kinds,           ONLY: DP
      use electrons_base,  only: nspin
      use gvecb,           only: gxb, ngb, npb, nmb
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx
      use cell_base,       only: omega
      use ions_base,       only: nsp, na, nat
      use parameters,      only: nsx
      use small_box,       only: tpibab
      use atom,            only: nlcc
      use core,            only: rhocb
      use fft_module,      only: invfft
      use fft_base,        only: dfftb
      use reciprocal_vectors, only: gstart
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx

      implicit none

! input
      integer, intent(in)        :: irb(3,nat)
      complex(8), intent(in):: eigrb(ngb,nat)
      real(8), intent(in)   :: vxc(nnr,nspin)
! output
      real(8), intent(inout):: fion1(3,nat)
! local
      integer iss, ix, ig, is, ia, nfft, isa
      real(8) fcc(3,nat), fac, boxdotgrid
      complex(8) ci, facg
      complex(8), allocatable :: qv(:)
      external  boxdotgrid
!
      call start_clock( 'forcecc' )
      ci = (0.d0,1.d0)
      fac = omega/DBLE(nr1*nr2*nr3*nspin)
      fcc = 0.d0

      allocate( qv( nnrb ) )

      isa = 0

      do is=1,nsp
         if( .not. nlcc(is) ) go to 10
#ifdef __PARA
         do ia=1,na(is)
            nfft=1
            if ( dfftb%np3( ia + isa ) <= 0 ) go to 15
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
                     facg = tpibab*CMPLX(0.d0,gxb(ix,ig))*rhocb(ig,is)
                     qv(npb(ig)) = eigrb(ig,ia+isa  )*facg                 &
     &                      + ci * eigrb(ig,ia+isa+1)*facg
                     qv(nmb(ig)) = CONJG(eigrb(ig,ia+isa  )*facg)          &
     &                      + ci * CONJG(eigrb(ig,ia+isa+1)*facg)
                  end do
               else
                  do ig=1,ngb
                     facg = tpibab*CMPLX(0.d0,gxb(ix,ig))*rhocb(ig,is)
                     qv(npb(ig)) = eigrb(ig,ia+isa)*facg
                     qv(nmb(ig)) = CONJG(eigrb(ig,ia+isa)*facg)
                  end do
               end if
!
               call invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,ia+isa)
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
      subroutine set_cc( irb, eigrb, rhoc )
!-----------------------------------------------------------------------
!
!     Calculate core charge contribution in real space, rhoc(r)
!     Same logic as for rhov: use box grid for core charges
!
      use ions_base,       only: nsp, na, nat
      use parameters,      only: nsx
      use atom,            only: nlcc
      use grid_dimensions, only: nr3, nnr => nnrx
      use gvecb,           only: ngb, npb, nmb
      use control_flags,   only: iprint
      use core,            only: rhocb
      use fft_module,      only: invfft
      use fft_base,        only: dfftb
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx

      implicit none
! input
      integer, intent(in)        :: irb(3,nat)
      complex(8), intent(in):: eigrb(ngb,nat)
! output
      real(8), intent(out)  :: rhoc(nnr)
! local
      integer nfft, ig, is, ia, isa
      complex(8) ci
      complex(8), allocatable :: wrk1(:)
      complex(8), allocatable :: qv(:)
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
            if ( dfftb%np3( ia + isa ) <= 0 ) go to 15
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
                  qv(nmb(ig))= CONJG(eigrb(ig,ia  +isa)*rhocb(ig,is))   &
     &                    + ci*CONJG(eigrb(ig,ia+1+isa)*rhocb(ig,is))
               end do
            else
               do ig=1,ngb
                  qv(npb(ig)) = eigrb(ig,ia+isa)*rhocb(ig,is)
                  qv(nmb(ig)) = CONJG(eigrb(ig,ia+isa)*rhocb(ig,is))
               end do
            endif
!
            call invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa+ia)
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

