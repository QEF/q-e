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

        LOGICAL :: tgc
        REAL(dbl), PARAMETER :: small_rho = 1.0d-10

        PUBLIC :: v2gc, exch_corr_energy, stress_xc
        PUBLIC :: exch_corr_init, exch_corr_print_info
        PUBLIC :: tgc

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

      SUBROUTINE exch_corr_init( xc_type )
        USE funct, ONLY: dft, which_dft, igcx, igcc
        CHARACTER(LEN = 20) :: xc_type
          dft = TRIM( xc_type )
          CALL which_dft( dft )
          tgc       = .FALSE.
          if( ( igcx > 0 ) .OR. ( igcc > 0 ) ) tgc = .TRUE.
        RETURN
      END SUBROUTINE exch_corr_init

!=----------------------------------------------------------------------------=!
      SUBROUTINE exch_corr_print_info(iunit)
        USE funct, ONLY: dft, iexch, icorr, igcx, igcc
        INTEGER, INTENT(IN) :: iunit
        CHARACTER(LEN = 20) :: xcgc_type
        CHARACTER(LEN = 60) :: exch_info
        CHARACTER(LEN = 60) :: corr_info
        CHARACTER(LEN = 60) :: exgc_info
        CHARACTER(LEN = 60) :: cogc_info

        WRITE(iunit,800)

          ! ...     iexch => Exchange functional form 
          ! ...     icorr => Correlation functional form
          ! ...     igcx  => Gradient Correction to the Exchange potential
          ! ...     igcc  => Gradient Correction to the Correlation potential

          SELECT CASE ( iexch )
            CASE (0)
              exch_info = 'NONE'
            CASE (1)
              exch_info = 'SLATER'
            CASE (2)
              exch_info = 'SLATER (alpha=1)'
            CASE DEFAULT
              exch_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( icorr )
            CASE (0)
              corr_info = 'NONE'
            CASE (1)
              corr_info = 'PERDEW AND ZUNGER'
            CASE (2)
              corr_info = 'VOSKO, WILK AND NUSAIR'
            CASE (3)
              corr_info = 'LEE, YANG, AND PARR'
            CASE (4)
              corr_info = 'PERDEW AND WANG'
            CASE (9)
              corr_info = 'PADE APPROXIMATION'
            CASE DEFAULT
              corr_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( igcx )
            CASE (0)
              exgc_info = 'NONE'
            CASE (1)
              exgc_info = 'BECKE'
            CASE (2)
              exgc_info = 'PERDEW'
            CASE (3)
              exgc_info = 'PERDEW BURKE ERNZERHOF'
            CASE DEFAULT
              exgc_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( igcc )
            CASE (0)
              cogc_info = 'NONE'
            CASE (1)
              cogc_info = 'PERDEW'
            CASE (2)
              cogc_info = 'LEE, YANG AND PARR'
            CASE (3)
              cogc_info = 'PERDEW AND WANG'
            CASE (4)
              cogc_info = 'PERDEW BURKE ERNZERHOF'
            CASE DEFAULT
              cogc_info = 'UNKNOWN'
          END SELECT

          WRITE(iunit,910)
          WRITE(iunit,fmt='(5X,"Exchange functional: ",A)') exch_info
          WRITE(iunit,fmt='(5X,"Correlation functional: ",A)') corr_info
          IF(tgc) THEN
            WRITE(iunit,810)
            WRITE(iunit,fmt='(5X,"Exchange functional: ",A)') exgc_info
            WRITE(iunit,fmt='(5X,"Correlation functional: ",A)') cogc_info
          END IF

        WRITE( iunit, '(5x,"Exchange-correlation      = ",a, &
       &       " (",4i1,")")') trim(dft) , iexch, icorr, igcx, igcc

800 FORMAT(//,3X,'Exchange and correlations functionals',/ &
             ,3X,'-------------------------------------')
810 FORMAT(   3X,'Using Generalized Gradient Corrections with')
910 FORMAT(   3X,'Using Local Density Approximation with')

        RETURN
      END SUBROUTINE exch_corr_print_info

!=----------------------------------------------------------------------------=!

        SUBROUTINE v2gc(v2xc, grho, rhoer, vpot, gv)

          USE fft
          USE stick, ONLY: dfftp
          USE cell_base, ONLY: tpiba
          USE cp_types
          USE mp_global
!                                                                       
          implicit none
!                                                                       
          REAL(dbl) ::  vpot(:,:,:,:)
          REAL(dbl), intent(in)  ::  v2xc(:,:,:,:,:)
          REAL(dbl), intent(in)  ::  grho(:,:,:,:,:)
          REAL(dbl), intent(in)  ::  rhoer(:,:,:,:)
          type (recvecs), intent(in)  ::  GV
!                                                                       
          integer :: ig, ipol, nxl, nyl, nzl, i, j, k, is, js, nspin
          integer :: ldx, ldy, ldz
          COMPLEX(dbl), allocatable ::  psi(:,:,:)
          COMPLEX(dbl), allocatable ::  vtemp(:)
          COMPLEX(dbl), allocatable ::  vtemp_pol(:)
          REAL(dbl), ALLOCATABLE :: v(:,:,:)
          REAL(dbl) :: fac
! ...                                                                   
          ldx   = dfftp%nr1
          ldy   = dfftp%nr2
          ldz   = dfftp%npl
          nxl   = MIN( ldx, SIZE( grho, 1 ) )
          nyl   = MIN( ldy, SIZE( grho, 2 ) )
          nzl   = MIN( ldz, SIZE( grho, 3 ) )
          nspin = SIZE(rhoer,4)

          fac = REAL(nspin) / 2.0d0

          ALLOCATE( vtemp( gv%ng_l ) )
          ALLOCATE( vtemp_pol( gv%ng_l ) )

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
                DO ig = gv%gstart, gv%ng_l
                  vtemp(ig) = vtemp(ig) + vtemp_pol(ig) *  CMPLX( 0.d0, tpiba * gv%gx_l( ipol, ig ) )
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
      USE stick, ONLY: dfftp

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
          gcpail(ic) = - REAL(nspin) / 2.0_dbl * stre * omega / REAL(nr1*nr2*nr3)
        END DO

      RETURN
    END SUBROUTINE stress_gc

!=----------------------------------------------------------------------------=!

    SUBROUTINE stress_xc( dexc, strvxc, sfac, vxc, tgc, grho, v2xc, &
        gagx_l, gv, tnlcc, rhocp, box)

      use ions_base, only: nsp
      USE cell_module, only: boxdimensions
      USE cell_base, ONLY: tpiba
      USE cp_types, ONLY: recvecs

      IMPLICIT NONE

      ! -- ARGUMENT

      type (recvecs), intent(in) :: gv
      type (boxdimensions), intent(in) :: box
      LOGICAL :: tgc, tnlcc(:)
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

        DO ig = gv%gstart, gv%ng_l
          tex1 = (0.0_dbl , 0.0_dbl)
          DO is=1,nsp
            IF ( tnlcc(is) ) THEN
              tex1 = tex1 + sfac( is, ig ) * CMPLX(rhocp(ig,is))
            END IF
          END DO
          tex2 = 0.0_dbl
          DO ispin = 1, nspin
            tex2 = tex2 + CONJG( vxc(ig, ispin) )
          END DO
          tex3 = REAL(tex1 * tex2) / SQRT(gv%hg_l(ig)) / tpiba
          dexc = dexc + tex3 * gagx_l(:,ig)
        END DO
        dexc = dexc * 2.0_dbl * omega

      END IF

      ! ... (E_{xc} - \int dr v_{xc}(n) n(r))/omega part of the stress
      ! ... this part of the stress is diagonal.

      dexc = dexc + strvxc * dalbe

      IF (tgc) THEN
        CALL stress_gc(grho, v2xc, gcpail, omega)
        dexc = dexc + gcpail
      END IF

      RETURN
      END SUBROUTINE


!=----------------------------------------------------------------------------=!


     SUBROUTINE exch_corr_energy(rhoetr, rhoetg, grho, vpot, sxc, vxc, v2xc, gv)

        USE kinds, ONLY: dbl
        USE cp_types, ONLY: recvecs

        REAL (dbl) :: rhoetr(:,:,:,:)
        COMPLEX(dbl) :: rhoetg(:,:)
        REAL (dbl) :: grho(:,:,:,:,:)
        REAL (dbl) :: vpot(:,:,:,:)
        REAL(dbl) :: sxc              ! E_xc   energy
        REAL(dbl) :: vxc              ! SUM ( v(r) * rho(r) )
        REAL (dbl) :: v2xc(:,:,:,:,:)
        TYPE (recvecs), INTENT(IN) :: gv

        INTEGER :: nspin, nnr

          !  vpot = vxc(rhoetr); vpot(r) <-- u(r)

          nnr   = SIZE( rhoetr, 1 ) * SIZE( rhoetr, 2 ) * SIZE( rhoetr, 3 )
          nspin = SIZE( rhoetr, 4 )
          !
          CALL exch_corr_wrapper( nnr, nspin, grho, rhoetr, sxc, vpot, v2xc )
          !
          vxc = SUM( vpot * rhoetr )
          !
          IF(tgc) THEN
            ! ... vpot additional term for gradient correction
            CALL v2gc( v2xc, grho, rhoetr, vpot, gv )
          END If

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
      end


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
      end

