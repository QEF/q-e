!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
  MODULE exchange_correlation
!=----------------------------------------------------------------------------=!

        USE kinds, ONLY: dbl

        IMPLICIT NONE
        SAVE

        PRIVATE

! ... Gradient Correction & exchange and correlation

        REAL(dbl), PARAMETER :: small_rho = 1.0d-10

        PUBLIC :: v2gc, exch_corr_energy, stress_xc

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

        SUBROUTINE v2gc(v2xc, grho, rhoer, vpot)

          USE fft
          USE fft_base, ONLY: dfftp
          USE cell_base, ONLY: tpiba
          USE cp_types
          USE mp_global
          USE reciprocal_vectors, ONLY: gstart, gx
          USE gvecp, ONLY: ngm
!                                                                       
          implicit none
!                                                                       
          REAL(dbl) ::  vpot(:,:,:,:)
          REAL(dbl), intent(in)  ::  v2xc(:,:,:,:,:)
          REAL(dbl), intent(in)  ::  grho(:,:,:,:,:)
          REAL(dbl), intent(in)  ::  rhoer(:,:,:,:)
!                                                                       
          integer :: ig, ipol, nxl, nyl, nzl, i, j, k, is, js, nspin
          integer :: ldx, ldy, ldz
          COMPLEX(dbl), allocatable ::  psi(:,:,:)
          COMPLEX(dbl), allocatable ::  vtemp(:)
          COMPLEX(dbl), allocatable ::  vtemp_pol(:)
          REAL(dbl), ALLOCATABLE :: v(:,:,:)
          REAL(dbl) :: fac
! ...                                                                   
          ldx   = dfftp%nr1x
          ldy   = dfftp%nr2x
          ldz   = dfftp%npl
          nxl   = MIN( dfftp%nr1, SIZE( grho, 1 ) )
          nyl   = MIN( dfftp%nr2, SIZE( grho, 2 ) )
          nzl   = MIN( dfftp%npl, SIZE( grho, 3 ) )
          nspin = SIZE(rhoer,4)

          !fac = REAL(nspin)
          fac = 1.0d0

          ALLOCATE( vtemp( ngm ) )
          ALLOCATE( vtemp_pol( ngm ) )

          DO js = 1, nspin
            !
            ALLOCATE( psi( ldx, ldy, ldz ) )
            !
            vtemp = 0.0d0

            DO ipol = 1, 3
              DO is = 1, nspin
                !
                DO k = 1, nzl
                  DO j = 1, nyl
                    DO i = 1, nxl
                      psi(i,j,k) = fac * v2xc(i,j,k,js,is) * grho(i,j,k,ipol,is)
                    END DO
                  END DO
                END DO
                !
                CALL pfwfft( vtemp_pol, psi )
                !
                DO ig = gstart, ngm
                  vtemp(ig) = vtemp(ig) + vtemp_pol(ig) *  CMPLX( 0.d0, tpiba * gx( ipol, ig ) )
                END DO
                !
              END DO
            END DO
            !
            DEALLOCATE( psi )

            ALLOCATE( v( ldx, ldy, ldz ) )
            !
            CALL pinvfft( v, vtemp )

            DO k = 1, nzl
              DO j = 1, nyl
                DO i = 1, nxl
                  vpot(i,j,k,js) = vpot(i,j,k,js) - v(i,j,k)
                END DO
              END DO
            END DO

            DEALLOCATE( v )

          END DO

          DEALLOCATE(vtemp_pol)
          DEALLOCATE(vtemp)

          RETURN
        END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE stress_gc(grho, v2xc, gcpail, omega)
!
      use grid_dimensions, only: nr1, nr2, nr3
      USE fft_base, ONLY: dfftp

        IMPLICIT NONE
!
        REAL(dbl) ::  v2xc(:,:,:,:,:)
        REAL(dbl) ::  grho(:,:,:,:,:)
        REAL(dbl) ::  gcpail(6)
        REAL(dbl) ::  omega
!
        REAL(dbl) :: stre, grhoi, grhoj
        INTEGER :: i, j, k, ipol, jpol, ic, nxl, nyl, nzl, is, js, nspin
        INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
        INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)
! ...
        nxl   = MIN( dfftp%nr1, SIZE( grho, 1 ) )
        nyl   = MIN( dfftp%nr2, SIZE( grho, 2 ) )
        nzl   = MIN( dfftp%npl, SIZE( grho, 3 ) )
        nspin = SIZE(grho,5)

        DO ic = 1, 6
          ipol = alpha(ic)
          jpol = beta(ic)
          stre = 0.0d0
          DO is = 1, nspin
            DO js = 1, nspin
              DO k = 1, nzl
                DO j = 1, nyl
                  DO i = 1, nxl
                    stre = stre + v2xc(i,j,k,is,js) * grho(i,j,k,ipol,js) * grho(i,j,k,jpol,is)
                  END DO
                END DO
              END DO
            END DO
          END DO
          gcpail(ic) = - REAL(nspin) * stre * omega / REAL(nr1*nr2*nr3)
        END DO

      RETURN
    END SUBROUTINE stress_gc

!=----------------------------------------------------------------------------=!

    SUBROUTINE stress_xc( dexc, strvxc, sfac, vxc, grho, v2xc, &
        gagx_l, tnlcc, rhocp, box)

      use ions_base,          only: nsp
      USE cell_module,        only: boxdimensions
      USE cell_base,          ONLY: tpiba
      USE funct,              ONLY: igcx, igcc 
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecp,              ONLY: ngm

      IMPLICIT NONE

      ! -- ARGUMENT

      type (boxdimensions), intent(in) :: box
      LOGICAL :: tnlcc(:)
      COMPLEX(dbl) :: vxc(:,:)
      COMPLEX(dbl), INTENT(IN) :: sfac(:,:)
      REAL(dbl) :: dexc(:), strvxc
      REAL(dbl) :: grho(:,:,:,:,:)
      REAL(dbl) :: v2xc(:,:,:,:,:)
      REAL(dbl) :: GAgx_L(:,:)
      REAL(dbl) :: rhocp(:,:)

      INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)
      ! ...  dalbe(:) = delta(alpha(:),beta(:))
      REAL(dbl),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_dbl, 0.0_dbl, 0.0_dbl, 1.0_dbl, 0.0_dbl, 1.0_dbl /)

      COMPLEX(dbl) :: tex1, tex2, tex3
      REAL(dbl) :: gcpail(6), omega
      INTEGER :: ig, k, is, ispin, nspin

      omega = box%deth
      nspin = SIZE(vxc, 2)

      DEXC = 0.0d0

      ! ... computes omega * \sum_{G}[ S(G)*rhopr(G)* G_{alpha} G_{beta}/|G|]
      ! ... (252) Phd thesis Dal Corso. Opposite sign.

      IF ( ANY( tnlcc ) ) THEN

        DO ig = gstart, ngm
          tex1 = (0.0_dbl , 0.0_dbl)
          DO is=1,nsp
            IF ( tnlcc(is) ) THEN
              tex1 = tex1 + sfac( ig, is ) * CMPLX(rhocp(ig,is))
            END IF
          END DO
          tex2 = 0.0_dbl
          DO ispin = 1, nspin
            tex2 = tex2 + CONJG( vxc(ig, ispin) )
          END DO
          tex3 = REAL(tex1 * tex2) / SQRT( g( ig ) ) / tpiba
          dexc = dexc + tex3 * gagx_l(:,ig)
        END DO
        dexc = dexc * 2.0_dbl * omega

      END IF

      ! ... (E_{xc} - \int dr v_{xc}(n) n(r))/omega part of the stress
      ! ... this part of the stress is diagonal.

      dexc = dexc + strvxc * dalbe

      IF ( ( igcx > 0 ) .OR. ( igcc > 0 ) ) THEN
        CALL stress_gc(grho, v2xc, gcpail, omega)
        dexc = dexc + gcpail
      END IF

      RETURN
      END SUBROUTINE


!=----------------------------------------------------------------------------=!


     SUBROUTINE exch_corr_energy(rhoetr, rhoetg, grho, vpot, sxc, vxc, v2xc)

        USE kinds, ONLY: dbl
        USE grid_dimensions, ONLY: nr1l, nr2l, nr3l
        USE funct, ONLY: igcx, igcc 

        REAL (dbl) :: rhoetr(:,:,:,:)
        COMPLEX(dbl) :: rhoetg(:,:)
        REAL (dbl) :: grho(:,:,:,:,:)
        REAL (dbl) :: vpot(:,:,:,:)
        REAL(dbl) :: sxc              ! E_xc   energy
        REAL(dbl) :: vxc              ! SUM ( v(r) * rho(r) )
        REAL (dbl) :: v2xc(:,:,:,:,:)
        REAL(dbl) :: ddot

        INTEGER :: nspin, nnr, ispin, j, k

          !  vpot = vxc(rhoetr); vpot(r) <-- u(r)

          nnr   = SIZE( rhoetr, 1 ) * SIZE( rhoetr, 2 ) * SIZE( rhoetr, 3 )
          nspin = SIZE( rhoetr, 4 )
          !
          CALL exch_corr_wrapper( nnr, nspin, grho(1,1,1,1,1), rhoetr(1,1,1,1), &
                                  sxc, vpot(1,1,1,1), v2xc(1,1,1,1,1) )

          !
          IF( ( igcx > 0 ) .OR. ( igcc > 0 ) ) THEN
            ! ... vpot additional term for gradient correction
            CALL v2gc( v2xc, grho, rhoetr, vpot )
          END If

          !
          ! vxc = SUM( vpot * rhoetr )
          !
          vxc = 0.0d0
          DO ispin = 1, nspin
            DO k = 1, nr3l
              DO j = 1, nr2l
                vxc = vxc + &
                DDOT ( nr1l, vpot(1,j,k,ispin), 1, rhoetr(1,j,k,ispin), 1 )
              END DO
            END DO
          END DO

        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!
   END MODULE exchange_correlation
!=----------------------------------------------------------------------------=!



!=----------------------------------------------------------------------------=!
!  CP subroutines
!=----------------------------------------------------------------------------=!

      subroutine exch_corr_h(nspin,rhog,rhor,exc,dxc)
!
! calculate exch-corr potential, energy, and derivatives dxc(i,j)
! of e(xc) with respect to to cell parameter h(i,j)
!     
      use funct, only: iexch, icorr, igcx, igcc
      use gvecp, only: ng => ngm
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx
      use cell_base, only: ainv, omega
      use control_flags, only: tpre
      use derho, only: drhor
      use mp, only: mp_sum
!
      implicit none
! input     
      integer nspin
! rhog contains the charge density in G space
! rhor contains the charge density in R space
      complex(kind=8) rhog(ng,nspin)
! output
! rhor contains the exchange-correlation potential
      real(kind=8) rhor(nnr,nspin), dxc(3,3), exc
! local
      integer i,j,ir
      real(kind=8) dexc(3,3)
      real(kind=8), allocatable:: gradr(:,:,:)
!
!     filling of gradr with the gradient of rho using fft's
!
      if (igcx /= 0 .or. igcc /= 0) then
         allocate(gradr(nnr,3,nspin))
         call fillgrad(nspin,rhog,gradr)
      end if
!
      CALL exch_corr_cp(nnr,nspin,gradr,rhor,exc)

      call mp_sum( exc )

      exc=exc*omega/dble(nr1*nr2*nr3)
!
! exchange-correlation contribution to pressure
!
      dxc = 0.0
      if (tpre) then
         !
         if ( nspin /= 1 ) call errore('exc-cor','spin not implemented',1)
         !
         do j=1,3
            do i=1,3
               do ir=1,nnr
                  dxc(i,j) = dxc(i,j) + rhor(ir,1)*drhor(ir,1,i,j)
               end do
               dxc(i,j)=omega/(nr1*nr2*nr3)*dxc(i,j)
            end do
         end do
         call mp_sum ( dxc )
         do j=1,3
            do i=1,3
               dxc(i,j) = dxc(i,j) + exc*ainv(j,i)
            end do
         end do
      end if
!
!     second part of the xc-potential
!
      if (igcx /= 0 .or. igcc /= 0) then
         call gradh( nspin, gradr, rhog, rhor, dexc)
         if (tpre) then
            call mp_sum ( dexc )
            dxc = dxc + dexc
         end if
         deallocate(gradr)
      end if
!
      return
      end subroutine exch_corr_h


!=----------------------------------------------------------------------------=!

      subroutine gradh(nspin,gradr,rhog,rhor,dexc)
!     _________________________________________________________________
!     
!     calculate the second part of gradient corrected xc potential
!     plus the gradient-correction contribution to pressure
!           
      use control_flags, only: iprint, tpre
      use reciprocal_vectors, only: gx
      use recvecs_indexes, only: np, nm
      use gvecp, only: ng => ngm
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx, nr1x, nr2x, nr3x
      use cell_base, only: ainv, tpiba, omega
      use derho, only: drhog
!                 
      implicit none  
! input                   
      integer nspin
      real(kind=8)  gradr(nnr,3,nspin), rhor(nnr,nspin), dexc(3,3)
      complex(kind=8) rhog(ng,nspin)
!
      complex(kind=8), allocatable:: v(:)
      complex(kind=8), allocatable:: x(:), vtemp(:)
      complex(kind=8)  ci, fp, fm
      integer iss, ig, ir, i,j
!
      allocate(v(nnr))
      allocate(x(ng))
      allocate(vtemp(ng))
      ci=(0.0,1.0)
      if (tpre .and. nspin.ne.1) &
           call errore('gradh','spin not implemented',1)
      do iss=1, nspin
!     _________________________________________________________________
!     second part xc-potential: 3 forward ffts
!
         do ir=1,nnr
            v(ir)=cmplx(gradr(ir,1,iss),0.0)
         end do
         call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ig=1,ng
            x(ig)=ci*tpiba*gx(1,ig)*v(np(ig))
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ng
                     vtemp(ig) = omega*ci*conjg(v(np(ig)))*             &
     &                    tpiba*(-rhog(ig,iss)*gx(i,ig)*ainv(j,1)+      &
     &                    gx(1,ig)*drhog(ig,iss,i,j))
                  end do
                  dexc(i,j) = real(SUM(vtemp))*2.0
               end do
            end do
         endif
!
         do ir=1,nnr
            v(ir)=cmplx(gradr(ir,2,iss),gradr(ir,3,iss))
         end do
         call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         do ig=1,ng
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*gx(2,ig)*0.5*cmplx( real(fp),aimag(fm))
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*gx(3,ig)*0.5*cmplx(aimag(fp),-real(fm))
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ng
                     fp=v(np(ig))+v(nm(ig))
                     fm=v(np(ig))-v(nm(ig))
                     vtemp(ig) = omega*ci*                              &
     &                    (0.5*cmplx(real(fp),-aimag(fm))*              &
     &                    tpiba*(-rhog(ig,iss)*gx(i,ig)*ainv(j,2)+      &
     &                    gx(2,ig)*drhog(ig,iss,i,j))+                  &
     &                    0.5*cmplx(aimag(fp),real(fm))*tpiba*          &
     &                    (-rhog(ig,iss)*gx(i,ig)*ainv(j,3)+            &
     &                    gx(3,ig)*drhog(ig,iss,i,j)))
                  end do
                  dexc(i,j) = dexc(i,j) + 2.0*real(SUM(vtemp))
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
            v(nm(ig))=conjg(x(ig))
         end do
         call invfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            rhor(ir,iss)=rhor(ir,iss)-real(v(ir))
         end do
      end do
!
      deallocate(vtemp)
      deallocate(x)
      deallocate(v)
!
      return
      end subroutine gradh

