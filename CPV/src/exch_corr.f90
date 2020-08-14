!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! ... Gradient Correction & exchange and correlation
!=----------------------------------------------------------------------------=!

   subroutine exch_corr_h( nspin, rhog, rhor, rhoc, sfac, exc, dxc, self_exc )
!
! calculate exch-corr potential, energy, and derivatives dxc(i,j)
! of e(xc) with respect to to cell parameter h(i,j)
!     
      use funct,           only : dft_is_gradient, dft_is_meta
      use fft_base,        only : dfftp, dffts
      use cell_base,       only : ainv, omega, h
      use ions_base,       only : nsp
      use control_flags,   only : tpre, iverbosity
      use core,            only : drhocg
      use gvect,           only : g
      use uspp,            only : nlcc_any
      use mp,              only : mp_sum
      use metagga_cp,      ONLY : kedtaur
      USE io_global,       ONLY : stdout
      USE mp_global,       ONLY : intra_bgrp_comm
      use kinds,           ONLY : DP
      use constants,       ONLY : au_gpa
      USE sic_module,      ONLY : self_interaction, sic_alpha
      USE cp_interfaces,   ONLY : denlcc
      use cp_main_variables,    only : drhor
!
      implicit none

      ! input     
      !
      integer nspin
      !
      ! rhog contains the charge density in G space
      ! rhor contains the charge density in R space
      !
      complex(DP) :: rhog( dfftp%ngm, nspin )
      complex(DP) :: sfac( dffts%ngm, nsp )
      !
      ! output
      ! rhor contains the exchange-correlation potential
      !
      real(DP) :: rhor( dfftp%nnr, nspin ), rhoc( dfftp%nnr )
      real(DP) :: dxc( 3, 3 ), exc
      real(DP) :: dcc( 3, 3 ), drc( 3, 3 )
      !
      ! local
      !
      integer :: i, j, ir, iss
      real(DP) :: dexc(3,3)
      real(DP), allocatable :: gradr(:,:,:)
      !
      !sic
      REAL(DP) :: self_exc
      REAL(DP), ALLOCATABLE :: self_rho( :,: ), self_gradr(:,:,:)
      complex(DP), ALLOCATABLE :: self_rhog( :,: )
      LOGICAL  :: ttsic
      real(DP) :: detmp(3,3)
      !
      !     filling of gradr with the gradient of rho using fft's
      !
      if ( dft_is_gradient() ) then
         !
         allocate( gradr( 3, dfftp%nnr, nspin ) )
         do iss = 1, nspin
            CALL fft_gradient_g2r ( dfftp, rhog(1,iss), g, gradr(1,1,iss) )
         end do
         ! 
      else
         ! 
         allocate( gradr( 3, 1, 2 ) )
         !
      end if


      ttsic = (self_interaction /= 0 )
      !
      IF ( ttsic ) THEN
         !
         IF ( dft_is_meta() ) CALL errore ('exch_corr_h', &
                               'SIC and meta-GGA not together', 1)
         IF ( tpre ) CALL errore( 'exch_corr_h', 'SIC and stress not implemented', 1)

         !  allocate the sic_arrays
         !
         ALLOCATE( self_rho( dfftp%nnr, nspin ) )
         ALLOCATE( self_rhog(dfftp%ngm, nspin ) )

         self_rho(:, :) = rhor( :, :)
         IF( dft_is_gradient() ) THEN
            ALLOCATE( self_gradr( 3, dfftp%nnr, nspin ) )
            self_gradr(:, :, :) = gradr(:, :, :)
         ENDIF
         self_rhog(:, :) = rhog( :, :)
         !
      END IF
!
      self_exc = 0.d0
!
      if( dft_is_meta() ) then
         !
         call tpssmeta( dfftp%nnr, nspin, gradr, rhor, kedtaur, exc )
         !
      else
         !
         CALL exch_corr_cp(dfftp%nnr, nspin, gradr, rhor, exc)
         !
         IF ( ttsic ) THEN
            CALL exch_corr_cp(dfftp%nnr, nspin, self_gradr, self_rho, self_exc)
            self_exc = sic_alpha * (exc - self_exc)
            exc = exc - self_exc
         END IF
!
      end if

      call mp_sum( exc, intra_bgrp_comm )
      IF ( ttsic ) call mp_sum( self_exc, intra_bgrp_comm )

      exc = exc * omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
      IF ( ttsic ) self_exc = self_exc * omega/DBLE(dfftp%nr1 * dfftp%nr2 *dfftp%nr3 )

      !     WRITE(*,*) 'Debug: calcolo exc', exc, 'eself', self_exc
      !
      ! exchange-correlation contribution to pressure
      !
      dxc = 0.0d0
      !
      if ( tpre ) then
         !
         !  Add term: Vxc( r ) * Drhovan( r )_ij - Vxc( r ) * rho( r ) * ((H^-1)^t)_ij
         !
         do iss = 1, nspin
            do j=1,3
               do i=1,3
                  do ir=1,dfftp%nnr
                     dxc(i,j) = dxc(i,j) + rhor( ir, iss ) * drhor( ir, iss, i, j )
                  end do
               end do
            end do
         end do
         !
         dxc = dxc * omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
         !
         call mp_sum ( dxc, intra_bgrp_comm )
         !
         do j = 1, 3
            do i = 1, 3
               dxc( i, j ) = dxc( i, j ) + exc * ainv( j, i )
            end do
         end do
         !
         ! DEBUG
         !
         ! write (stdout,*) "derivative of e(xc)"
         ! write (stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
         !
         IF( iverbosity > 1 ) THEN
            DO i=1,3
               DO j=1,3
                  detmp(i,j)=exc*ainv(j,i)
               END DO
            END DO
            WRITE( stdout,*) "derivative of e(xc) - diag - kbar"
            detmp = -1.0d0 * MATMUL( detmp, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
         END IF
         !
      end if
      !
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !
         call gradh( nspin, gradr, rhog, rhor, dexc)
         !
         if (tpre) then
            !
            call mp_sum ( dexc, intra_bgrp_comm )
            !
            dxc = dxc + dexc
            !
         end if
         !
      end if
      !

      IF( ttsic ) THEN
!
         IF (dft_is_gradient()) then
      
            call gradh( nspin, self_gradr, self_rhog, self_rho, dexc)
          
            gradr(:,:, 1) = (1.d0 - sic_alpha ) * gradr(:,:, 1)
            gradr(:,:, 2) = (1.d0 - sic_alpha ) * gradr(:,:, 2) + &
                          &  sic_alpha * ( self_gradr(:,:,1) + self_gradr(:,:,2) )
         ENDIF

         rhor(:, 1) = (1.d0 - sic_alpha ) * rhor(:, 1)
         rhor(:, 2) = (1.d0 - sic_alpha ) * rhor(:, 2) + &
                    &  sic_alpha * ( self_rho(:,1) + self_rho(:,2) )

         IF(ALLOCATED(self_gradr)) DEALLOCATE(self_gradr)
         IF(ALLOCATED(self_rhog)) DEALLOCATE(self_rhog)
         IF(ALLOCATED(self_rho)) DEALLOCATE(self_rho)
!
      ENDIF

      IF( tpre ) THEN
         !
         dcc = 0.0d0
         !
         IF( nlcc_any ) CALL  denlcc( dfftp%nnr, nspin, rhor, sfac, drhocg, dcc )
         !
         ! DEBUG
         !
         !  write (stdout,*) "derivative of e(xc) - nlcc part"
         !  write (stdout,5555) ((dcc(i,j),j=1,3),i=1,3)
         !
         dxc = dxc + dcc
         !
         do iss = 1, nspin
            drc = 0.0d0
            IF( nlcc_any ) THEN
               do j=1,3
                  do i=1,3
                     do ir=1,dfftp%nnr
                        drc(i,j) = drc(i,j) + rhor( ir, iss ) * rhoc( ir ) * ainv(j,i)
                     end do
                  end do
               end do
               call mp_sum ( drc, intra_bgrp_comm )
            END IF
            dxc = dxc - drc * ( 1.0d0 / nspin ) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
         end do
         !
      END IF
      !

      IF( ALLOCATED( gradr ) ) DEALLOCATE( gradr )

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
      use gvect, only: g
      use cell_base, only: ainv, tpiba, omega
      use cp_main_variables, only: drhog
      USE fft_interfaces, ONLY: fwfft, invfft
      USE fft_base,       ONLY: dfftp
      USE fft_helper_subroutines, ONLY: fftx_threed2oned, fftx_oned2threed
!                 
      implicit none  
! input                   
      integer nspin
      real(DP)    :: gradr( 3, dfftp%nnr, nspin ), rhor( dfftp%nnr, nspin ), dexc( 3, 3 )
      complex(DP) :: rhog( dfftp%ngm, nspin )
!
      complex(DP), allocatable:: v(:), vp(:), vm(:)
      complex(DP), allocatable:: x(:), vtemp(:)
      complex(DP) ::  ci, fp, fm
      integer :: iss, ig, ir, i,j
!
      allocate(v(dfftp%nnr))
      allocate(x(dfftp%ngm))
      allocate(vp(dfftp%ngm))
      allocate(vm(dfftp%ngm))
      allocate(vtemp(dfftp%ngm))
      !
      ci=(0.0d0,1.0d0)
      !
      dexc = 0.0d0
      !
      do iss=1, nspin
!     _________________________________________________________________
!     second part xc-potential: 3 forward ffts
!
         do ir=1,dfftp%nnr
            v(ir)=CMPLX(gradr(1,ir,iss),0.d0,kind=DP)
         end do
         call fwfft('Rho',v, dfftp )
         CALL fftx_threed2oned( dfftp, v, vp )
         do ig=1,dfftp%ngm
            x(ig)=ci*tpiba*g(1,ig)*vp(ig)
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,dfftp%ngm
                     vtemp(ig) = omega*ci*CONJG(vp(ig))*             &
     &                    tpiba*(-rhog(ig,iss)*g(i,ig)*ainv(j,1)+      &
     &                    g(1,ig)*drhog(ig,iss,i,j))
                  end do
                  dexc(i,j) = dexc(i,j) + DBLE(SUM(vtemp))*2.0d0
               end do
            end do
         endif
!
         do ir=1,dfftp%nnr
            v(ir)=CMPLX(gradr(2,ir,iss),gradr(3,ir,iss),kind=DP)
         end do
         call fwfft('Rho',v, dfftp )
         CALL fftx_threed2oned( dfftp, v, vp, vm )
!
         do ig=1,dfftp%ngm
            x(ig) = x(ig) + ci*tpiba*g(2,ig)*vp(ig)
            x(ig) = x(ig) + ci*tpiba*g(3,ig)*vm(ig)
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,dfftp%ngm
                     vtemp(ig) = omega*ci*( &
     &                    CONJG(vp(ig))*tpiba*(-rhog(ig,iss)*g(i,ig)*ainv(j,2)+g(2,ig)*drhog(ig,iss,i,j))+ &
     &                    CONJG(vm(ig))*tpiba*(-rhog(ig,iss)*g(i,ig)*ainv(j,3)+g(3,ig)*drhog(ig,iss,i,j))  )
                  end do
                  dexc(i,j) = dexc(i,j) + 2.0d0*DBLE(SUM(vtemp))
               end do
            end do
         endif
!     _________________________________________________________________
!     second part xc-potential: 1 inverse fft
!
         CALL fftx_oned2threed( dfftp, v, x )
         call invfft('Rho',v, dfftp )
         do ir=1,dfftp%nnr
            rhor(ir,iss)=rhor(ir,iss)-DBLE(v(ir))
         end do
      end do
!
      deallocate(vtemp)
      deallocate(x)
      deallocate(v)
      deallocate(vp)
      deallocate(vm)
!
      return
   end subroutine gradh

!=----------------------------------------------------------------------------=!
!
!  For CP we need a further small interface subroutine
!
!=----------------------------------------------------------------------------=!

subroutine exch_corr_cp(nnr,nspin,grhor,rhor,etxc)
  use kinds,       only: DP
  use funct,       only: dft_is_gradient, get_igcc
  use xc_lda_lsda, only: xc
  use xc_gga,      only: xc_gcx, change_threshold_gga
  implicit none
  integer, intent(in) :: nnr
  integer, intent(in) :: nspin
  real(DP) :: grhor( 3, nnr, nspin )
  real(DP) :: rhor( nnr, nspin )
  real(DP) :: etxc

  real(DP), parameter :: epsr = 1.0d-10
  real(DP), parameter :: e2=1.0_dp
  integer :: ir, is, k, ipol, neg(3)
  real(DP) ::  grup, grdw
  real(DP), allocatable :: v(:,:)
  real(DP), allocatable :: h(:,:,:)
  real(DP), allocatable :: rhox (:,:)!^
  real(DP), allocatable :: ex(:), ec(:)
  real(DP), allocatable :: vx(:,:), vc(:,:)
  REAL(DP), allocatable :: sx(:), sc(:)
  REAL(DP), allocatable :: v1x(:,:), v2x(:,:), v1c(:,:), v2c(:,:)
  real(dp), allocatable :: v2c_ud(:)
  real(dp) :: zetas
  !
  logical :: debug_xc = .false.
  logical :: igcc_is_lyp
  !
  allocate( v( nnr, nspin ) )
  if( dft_is_gradient() ) then
    allocate( h( nnr, nspin, nspin ) )
  else
    allocate( h( 1, 1, 1 ) )
  endif
  !
  igcc_is_lyp = (get_igcc() == 3)
  !
  etxc = 0.0d0
  !
  allocate ( ex(nnr), ec(nnr), vx(nnr,nspin), vc(nnr,nspin) )
  IF ( nspin == 1 ) THEN
     !
     ! spin-unpolarized case
     !
     CALL xc( nnr, nspin, nspin, rhor, ex, ec, vx, vc )
     !
     v(:,nspin) = e2 * (vx(:,1) + vc(:,1) )
     etxc = e2 * SUM( (ex + ec)*rhor(:,nspin) )
     !
  ELSE
     !
     ! spin-polarized case
     !
     neg(1) = 0
     neg(2) = 0
     neg(3) = 0
     !
     allocate ( rhox(nnr,2) ) !^
     rhox(:,1) = rhor(:,1) + rhor(:,2)
     rhox(:,2) = rhor(:,1) - rhor(:,2)
     !
     CALL xc( nnr, 2, 2, rhox, ex, ec, vx, vc )
     !
     DO ir = 1, nnr
        !
        DO is = 1, nspin
           v(ir,is) = e2 * (vx(ir,is) + vc(ir,is))
        ENDDO
        etxc = etxc + e2 * (ex(ir) + ec(ir)) * rhox(ir,1)
        !
        zetas =  rhox(ir,2) / rhox(ir,1)
        IF (rhor(ir,1) < 0.d0) neg(1) = neg(1) + 1
        IF (rhor(ir,2) < 0.d0) neg(2) = neg(2) + 1
        IF (ABS(zetas) > 1.d0) neg(3) = neg(3) + 1
        !
     ENDDO
     !
     deallocate ( rhox ) !^
     !
  ENDIF
  deallocate ( vc, vx, ec, ex )
  !
  if( debug_xc ) then
    open(unit=17,form='unformatted')
    write(17) nnr, nspin
    write(17) rhor
    write(17) grhor
    close(17)
    debug_xc = .false.
  end if
  !
  ! gradient corrections
  !
  if ( dft_is_gradient() ) then
    !
    call change_threshold_gga( epsr )
    !
    allocate ( sx(nnr), sc(nnr), v1x(nnr,nspin), v1c(nnr,nspin), &
               v2x(nnr,nspin), v2c(nnr,nspin) )
    if (nspin == 1) then
       !
       ! ... This is the spin-unpolarised case
       !
       call xc_gcx( nnr, nspin, rhor, grhor, sx, sc, v1x, v2x, v1c, v2c )
       !
       do k = 1, nnr
          ! first term of the gradient correction: D(rho*Exc)/D(rho)
          v(k,1) = v(k,1) + e2 * (v1x(k,1) + v1c(k,1))
          ! HERE h contains D(rho*Exc)/D(|grad rho|) / |grad rho|
          h(k, 1, 1) = e2 * (v2x(k,1) + v2c(k,1))
          etxc = etxc + e2 * (sx(k) + sc(k))
       enddo
       !
    else
       !
       ! ... Spin-polarised case
       !
       allocate (v2c_ud(nnr))
       call xc_gcx( nnr, 2, rhor, grhor, sx, sc, v1x, v2x, v1c, v2c, v2c_ud )
       !
       ! first term of the gradient correction : D(rho*Exc)/D(rho)
       !
       v = v + e2*( v1x + v1c )
       !
       ! HERE h contains D(rho*Exc)/D(|grad rho|) / |grad rho|
       !
       h(:,1,1) = e2 * (v2x(:,1) + v2c(:,1))  ! Spin UP-UP
       h(:,1,2) = e2 * v2c_ud(:)              ! Spin UP-DW
       h(:,2,1) = e2 * v2c_ud(:)              ! Spin DW-UP
       h(:,2,2) = e2 * (v2x(:,2) + v2c(:,2))  ! Spin DW-DW
       !
       etxc = etxc + e2 * SUM( sx(:)+sc(:) )
       !
       deallocate (v2c_ud)
    endif
    !
    deallocate ( v2c, v2x, v1c, v1x, sc, sx )
  end if

  if( dft_is_gradient() ) then
     !
     if( nspin == 1 ) then
        !
        ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
!$omp parallel default(none), shared(nnr,grhor,h), private(ipol,k)
        do ipol = 1, 3
!$omp do
           do k = 1, nnr
              grhor (ipol, k, 1) = h (k, 1, 1) * grhor (ipol, k, 1)
           enddo
!$omp end do
        end do
!$omp end parallel
        !
        !
     else
        !
!$omp parallel default(none), shared(nnr,grhor,h), private(ipol,k,grup,grdw)
        do ipol = 1, 3
!$omp do
           do k = 1, nnr
             grup = grhor (ipol, k, 1)
             grdw = grhor (ipol, k, 2)
             grhor (ipol, k, 1) = h (k, 1, 1) * grup + h (k, 1, 2) * grdw
             grhor (ipol, k, 2) = h (k, 2, 2) * grdw + h (k, 2, 1) * grup
           enddo
!$omp end do
        enddo
!$omp end parallel
        !
     end if
     !
  end if

  rhor = v 

  deallocate( v )
  deallocate( h )

  return
end subroutine exch_corr_cp
