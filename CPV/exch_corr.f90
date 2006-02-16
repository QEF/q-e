!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE exchange_correlation
!=----------------------------------------------------------------------------=!
#include "f_defs.h"

        USE kinds, ONLY: DP

        IMPLICIT NONE
        SAVE

        PRIVATE

! ... Gradient Correction & exchange and correlation

        REAL(DP), PARAMETER :: small_rho = 1.0d-10

        PUBLIC :: v2gc, exch_corr_energy, stress_xc

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

        SUBROUTINE v2gc( v2xc, grho, rhoer, vpot )

          USE kinds,              ONLY: DP
          USE fft_base,           ONLY: dfftp
          USE cell_base,          ONLY: tpiba
          USE reciprocal_vectors, ONLY: gstart, gx
          use grid_dimensions,    only: nnrx
          USE gvecp,              ONLY: ngm
          USE fft_module,         ONLY: fwfft, invfft
!                                                                       
          implicit none
!                                                                       
          REAL(DP) ::  vpot(:,:)
          REAL(DP), intent(in)  ::  v2xc(:,:,:)
          REAL(DP), intent(in)  ::  grho(:,:,:)
          REAL(DP), intent(in)  ::  rhoer(:,:)
!                                                                       
          integer :: ig, ipol, is, js, nspin
          COMPLEX(DP), allocatable ::  psi(:)
          COMPLEX(DP), allocatable ::  vtemp(:)
          COMPLEX(DP), allocatable ::  vtemp_pol(:)
          REAL(DP), ALLOCATABLE :: v(:)
          REAL(DP) :: fac
! ...                                                                   
          nspin = SIZE(rhoer,2)

          fac = 1.0d0

          ALLOCATE( vtemp( ngm ) )
          ALLOCATE( vtemp_pol( ngm ) )
          ALLOCATE( psi( nnrx ) )

          DO js = 1, nspin
            !
            vtemp = 0.0d0

            DO ipol = 1, 3
              DO is = 1, nspin
                !
                psi( 1:nnrx ) = fac * v2xc( 1:nnrx, js, is ) * grho( 1:nnrx, ipol, is )
                !
                CALL fwfft(   'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
                CALL psi2rho( 'Dense', psi, dfftp%nnr, vtemp_pol, ngm )
                !
                DO ig = gstart, ngm
                  vtemp(ig) = vtemp(ig) + vtemp_pol(ig) *  CMPLX( 0.d0, tpiba * gx( ipol, ig ) )
                END DO
                !
              END DO
            END DO
            !
            CALL rho2psi( 'Dense', psi, dfftp%nnr, vtemp, ngm )
            CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

            vpot( 1:nnrx, js ) = vpot( 1:nnrx, js) - DBLE( psi( 1:nnrx ) )

          END DO

          DEALLOCATE( psi )
          DEALLOCATE( vtemp_pol )
          DEALLOCATE( vtemp )

          RETURN
        END SUBROUTINE v2gc

!=----------------------------------------------------------------------------=!

    SUBROUTINE stress_gc(grho, v2xc, gcpail, omega)
!
      use grid_dimensions, only: nr1, nr2, nr3, nnrx

        IMPLICIT NONE
!
        REAL(DP) ::  v2xc(:,:,:)
        REAL(DP) ::  grho(:,:,:)
        REAL(DP) ::  gcpail(6)
        REAL(DP) ::  omega
!
        REAL(DP) :: stre, grhoi, grhoj
        INTEGER :: i, ipol, jpol, ic, is, js, nspin
        INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
        INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)
! ...
        nspin = SIZE(grho,3)

        DO ic = 1, 6
          ipol = alpha(ic)
          jpol = beta(ic)
          stre = 0.0d0
          DO is = 1, nspin
            DO js = 1, nspin
              DO i = 1, nnrx
                stre = stre + v2xc(i,is,js) * grho(i,ipol,js) * grho(i,jpol,is)
              END DO
            END DO
          END DO
          gcpail(ic) = - DBLE(nspin) * stre * omega / DBLE(nr1*nr2*nr3)
        END DO

      RETURN
    END SUBROUTINE stress_gc

!=----------------------------------------------------------------------------=!

    SUBROUTINE stress_xc( dexc, strvxc, sfac, vxc, grho, v2xc, &
        gagx_l, tnlcc, rhocp, box)

      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: nsp
      USE cell_module,        ONLY: boxdimensions
      USE cell_base,          ONLY: tpiba
      USE funct,              ONLY: dft_is_gradient
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecp,              ONLY: ngm
      USE io_global,          ONLY: stdout

      IMPLICIT NONE

      ! -- ARGUMENT

      type (boxdimensions), intent(in) :: box
      LOGICAL :: tnlcc(:)
      COMPLEX(DP) :: vxc(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      REAL(DP) :: dexc(:), strvxc
      REAL(DP) :: grho(:,:,:)
      REAL(DP) :: v2xc(:,:,:)
      REAL(DP) :: GAgx_L(:,:)
      REAL(DP) :: rhocp(:,:)

      INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)
      ! ...  dalbe(:) = delta(alpha(:),beta(:))
      REAL(DP),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)

      COMPLEX(DP) :: tex1, tex2, tex3
      REAL(DP) :: gcpail(6), omega, detmp( 3, 3 )
      REAL(DP) :: dcc( 6 )
      INTEGER :: ig, k, is, ispin, nspin
      INTEGER :: i, j

      omega = box%deth
      nspin = SIZE(vxc, 2)

      DEXC = 0.0d0
      dcc  = 0.0d0

      ! ... computes omega * \sum_{G}[ S(G)*rhopr(G)* G_{alpha} G_{beta}/|G|]
      ! ... (252) Phd thesis Dal Corso. Opposite sign.

      IF ( ANY( tnlcc ) ) THEN

        DO ig = gstart, ngm
          tex1 = (0.0_DP , 0.0_DP)
          DO is=1,nsp
            IF ( tnlcc(is) ) THEN
              tex1 = tex1 + sfac( ig, is ) * CMPLX(rhocp(ig,is), 0.d0)
            END IF
          END DO
          tex2 = 0.0_DP
          DO ispin = 1, nspin
            tex2 = tex2 + CONJG( vxc(ig, ispin) )
          END DO
          tex3 = DBLE(tex1 * tex2) / SQRT( g( ig ) ) / tpiba
          dcc = dcc + tex3 * gagx_l(:,ig)
        END DO
        dcc = dcc * 2.0_DP * omega

        ! DEBUG
        ! DO k=1,6
        !   detmp(alpha(k),beta(k)) = dcc(k)
        !   detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        ! END DO
        ! detmp = MATMUL( detmp(:,:), box%m1(:,:) )
        ! WRITE( stdout,*) "derivative of e(xc) - nlcc part"
        ! WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)


      END IF

      ! ... (E_{xc} - \int dr v_{xc}(n) n(r))/omega part of the stress
      ! ... this part of the stress is diagonal.

      dexc = strvxc * dalbe

      IF ( dft_is_gradient() ) THEN
        CALL stress_gc(grho, v2xc, gcpail, omega)
        dexc = dexc + gcpail
      END IF

      ! DEBUG
      ! DO k=1,6
      !   detmp(alpha(k),beta(k)) = dexc(k)
      !   detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
      ! END DO
      ! detmp = MATMUL( detmp(:,:), box%m1(:,:) )
      ! WRITE( stdout,*) "derivative of e(xc)"
      ! WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

      dexc = dexc + dcc

      RETURN

5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

      END SUBROUTINE stress_xc


!=----------------------------------------------------------------------------=!


     SUBROUTINE exch_corr_energy(rhoetr, rhoetg, grho, vpot, sxc, vxc, v2xc)

        USE kinds,           ONLY: DP
        use grid_dimensions, only: nnrx
        USE funct, ONLY: dft_is_gradient

        REAL (DP) :: rhoetr(:,:)
        COMPLEX(DP) :: rhoetg(:,:)
        REAL (DP) :: grho(:,:,:)
        REAL (DP) :: vpot(:,:)
        REAL (DP) :: sxc              ! E_xc   energy
        REAL (DP) :: vxc              ! SUM ( v(r) * rho(r) )
        REAL (DP) :: v2xc(:,:,:)
        !
        REAL (DP), EXTERNAL :: ddot

        INTEGER :: nspin, ispin
        logical :: is_gradient

        is_gradient =  dft_is_gradient()


        !  vpot = vxc(rhoetr); vpot(r) <-- u(r)

        nspin = SIZE( rhoetr, 2 )
        !
        IF( SIZE( vpot, 1 ) /= nnrx ) &
           CALL errore(" exch_corr_energy ", " inconsistent size for vpot ", 1 )
        !
        CALL exch_corr_wrapper( nnrx, nspin, grho(1,1,1), rhoetr(1,1), sxc, vpot(1,1), v2xc(1,1,1) )
        !
        IF( dft_is_gradient() ) THEN
          ! ... vpot additional term for gradient correction
          CALL v2gc( v2xc, grho, rhoetr, vpot )
        END If

        !
        ! vxc = SUM( vpot * rhoetr )
        !
        vxc = 0.0d0
        DO ispin = 1, nspin
           vxc = vxc + DDOT ( nnrx, vpot(1,ispin), 1, rhoetr(1,ispin), 1 )
        END DO


        RETURN
      END SUBROUTINE exch_corr_energy

!=----------------------------------------------------------------------------=!
   END MODULE exchange_correlation
!=----------------------------------------------------------------------------=!



!=----------------------------------------------------------------------------=!
!  CP subroutines
!=----------------------------------------------------------------------------=!

      subroutine exch_corr_h( nspin, rhog, rhor, rhoc, sfac, exc, dxc )
!
! calculate exch-corr potential, energy, and derivatives dxc(i,j)
! of e(xc) with respect to to cell parameter h(i,j)
!     
      use funct,           only : dft_is_gradient, dft_is_meta
      use gvecp,           only : ng => ngm
      use gvecs,           only : ngs
      use grid_dimensions, only : nr1, nr2, nr3, nnr => nnrx
      use cell_base,       only : ainv, omega
      use ions_base,       only : nsp
      use control_flags,   only : tpre
      use derho,           only : drhor
      use core,            only : drhocg, nlcc_any
      use mp,              only : mp_sum
      use metagga,         ONLY : kedtaur
      USE io_global,       ONLY : stdout
      use kinds,           ONLY : DP
!
      implicit none

      ! input     
      !
      integer nspin
      !
      ! rhog contains the charge density in G space
      ! rhor contains the charge density in R space
      !
      complex(DP) :: rhog( ng, nspin )
      complex(DP) :: sfac( ngs, nsp )
      !
      ! output
      ! rhor contains the exchange-correlation potential
      !
      real(DP) :: rhor( nnr, nspin ), rhoc( nnr )
      real(DP) :: dxc( 3, 3 ), exc
      real(DP) :: dcc( 3, 3 ), drc( 3, 3 )
      !
      ! local
      !
      integer :: i, j, ir, iss
      real(DP) :: dexc(3,3)
      real(DP), allocatable :: gradr(:,:,:)
      !
      !     filling of gradr with the gradient of rho using fft's
      !
      if ( dft_is_gradient() ) then
         !
         allocate( gradr( nnr, 3, nspin ) )
         call fillgrad( nspin, rhog, gradr )
         ! 
      end if

!
      if( dft_is_meta() ) then
            call tpssmeta( nnr, nspin, gradr, rhor, kedtaur, exc )
      else
            CALL exch_corr_cp(nnr, nspin, gradr, rhor, exc)
      end if

      call mp_sum( exc )

      exc = exc * omega / DBLE( nr1 * nr2 * nr3 )

      !
      ! exchange-correlation contribution to pressure
      !
      dxc = 0.0d0
      !
      if (tpre) then
         !
         IF( nlcc_any ) CALL  denlcc( nnr, nspin, rhor, sfac, drhocg, dcc )
         !
         ! DEBUG
         !
         ! write (stdout,*) "derivative of e(xc) - nlcc part"
         ! write (stdout,5555) ((dcc(i,j),j=1,3),i=1,3)
         !
         do iss = 1, nspin
            do j=1,3
               do i=1,3
                  do ir=1,nnr
                     dxc(i,j) = dxc(i,j) + rhor( ir, iss ) * drhor( ir, iss, i, j )
                  end do
               end do
            end do
            drc = 0.0d0
            IF( nlcc_any ) THEN
               do j=1,3
                  do i=1,3
                     do ir=1,nnr
                        drc(i,j) = drc(i,j) + rhor( ir, iss ) * rhoc( ir ) * ainv(j,i)
                     end do
                  end do
               end do
            END IF
            dxc = dxc - drc * 1.0d0 / nspin
         end do
         !
         dxc = dxc * omega / ( nr1*nr2*nr3 )
         !
         call mp_sum ( dxc )
         !
         do j=1,3
            do i=1,3
               dxc(i,j) = dxc(i,j) + exc * ainv(j,i)
            end do
         end do
         !
         ! DEBUG
         !
         ! write (stdout,*) "derivative of e(xc)"
         ! write (stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
         !
         dxc = dxc + dcc
         !
      end if
      !
      !     second part of the xc-potential
      !
      if (dft_is_gradient()) then
         !
         call gradh( nspin, gradr, rhog, rhor, dexc)
         !
         if (tpre) then
            !
            call mp_sum ( dexc )
            !
            dxc = dxc + dexc
            !
         end if
         !
         deallocate(gradr)
         !
      end if
      !
5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
      !
      return
      end subroutine exch_corr_h


!=----------------------------------------------------------------------------=!

      subroutine gradh( nspin, gradr, rhog, rhor, dexc )
!     _________________________________________________________________
!     
!     calculate the second part of gradient corrected xc potential
!     plus the gradient-correction contribution to pressure
!           
      USE kinds,              ONLY: DP
      use control_flags, only: iprint, tpre
      use reciprocal_vectors, only: gx
      use recvecs_indexes, only: np, nm
      use gvecp, only: ng => ngm
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx, nr1x, nr2x, nr3x
      use cell_base, only: ainv, tpiba, omega
      use derho, only: drhog
      USE fft_module, ONLY: fwfft, invfft
!                 
      implicit none  
! input                   
      integer nspin
      real(DP)    :: gradr( nnr, 3, nspin ), rhor( nnr, nspin ), dexc( 3, 3 )
      complex(DP) :: rhog( ng, nspin )
!
      complex(DP), allocatable:: v(:)
      complex(DP), allocatable:: x(:), vtemp(:)
      complex(DP) ::  ci, fp, fm
      integer :: iss, ig, ir, i,j
!
      allocate(v(nnr))
      allocate(x(ng))
      allocate(vtemp(ng))
      !
      ci=(0.0,1.0)
      !
      dexc = 0.0d0
      !
      do iss=1, nspin
!     _________________________________________________________________
!     second part xc-potential: 3 forward ffts
!
         do ir=1,nnr
            v(ir)=CMPLX(gradr(ir,1,iss),0.d0)
         end do
         call fwfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ig=1,ng
            x(ig)=ci*tpiba*gx(1,ig)*v(np(ig))
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ng
                     vtemp(ig) = omega*ci*CONJG(v(np(ig)))*             &
     &                    tpiba*(-rhog(ig,iss)*gx(i,ig)*ainv(j,1)+      &
     &                    gx(1,ig)*drhog(ig,iss,i,j))
                  end do
                  dexc(i,j) = dexc(i,j) + DBLE(SUM(vtemp))*2.0
               end do
            end do
         endif
!
         do ir=1,nnr
            v(ir)=CMPLX(gradr(ir,2,iss),gradr(ir,3,iss))
         end do
         call fwfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         do ig=1,ng
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*gx(2,ig)*0.5*CMPLX( DBLE(fp),AIMAG(fm))
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*gx(3,ig)*0.5*CMPLX(AIMAG(fp),-DBLE(fm))
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ng
                     fp=v(np(ig))+v(nm(ig))
                     fm=v(np(ig))-v(nm(ig))
                     vtemp(ig) = omega*ci*                              &
     &                    (0.5*CMPLX(DBLE(fp),-AIMAG(fm))*              &
     &                    tpiba*(-rhog(ig,iss)*gx(i,ig)*ainv(j,2)+      &
     &                    gx(2,ig)*drhog(ig,iss,i,j))+                  &
     &                    0.5*CMPLX(AIMAG(fp),DBLE(fm))*tpiba*          &
     &                    (-rhog(ig,iss)*gx(i,ig)*ainv(j,3)+            &
     &                    gx(3,ig)*drhog(ig,iss,i,j)))
                  end do
                  dexc(i,j) = dexc(i,j) + 2.0*DBLE(SUM(vtemp))
               end do
            end do
         endif
!     _________________________________________________________________
!     second part xc-potential: 1 inverse fft
!
         do ig=1,nnr
            v(ig)=(0.0,0.0)
         end do
         do ig=1,ng
            v(np(ig))=x(ig)
            v(nm(ig))=CONJG(x(ig))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            rhor(ir,iss)=rhor(ir,iss)-DBLE(v(ir))
         end do
      end do
!
      deallocate(vtemp)
      deallocate(x)
      deallocate(v)
!
      return
      end subroutine gradh

