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
      use gvect,           only : ngm
      use gvecs,           only : ngms
      use fft_base,        only : dfftp
      use cell_base,       only : ainv, omega, h
      use ions_base,       only : nsp
      use control_flags,   only : tpre, iverbosity
      use core,            only : drhocg
      use uspp,            only : nlcc_any
      use mp,              only : mp_sum
      use metagga,         ONLY : kedtaur
      USE io_global,       ONLY : stdout
      USE mp_global,       ONLY : intra_bgrp_comm
      use kinds,           ONLY : DP
      use constants,       ONLY : au_gpa
      USE sic_module,      ONLY : self_interaction, sic_alpha
      USE cp_interfaces,   ONLY : fillgrad, denlcc
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
      complex(DP) :: rhog( ngm, nspin )
      complex(DP) :: sfac( ngms, nsp )
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
         allocate( gradr( dfftp%nnr, 3, nspin ) )
         call fillgrad( nspin, rhog, gradr )
         ! 
      else
         ! 
         allocate( gradr( 1, 3, 2 ) )
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
         ALLOCATE( self_rhog(ngm, nspin ) )
         IF( dft_is_gradient() ) ALLOCATE( self_gradr( dfftp%nnr, 3, nspin ) )

         self_rho(:, 1) = rhor( :, 2)
         self_rho(:, 2) = rhor( :, 2)

         IF( dft_is_gradient() ) THEN
            self_gradr(:, :, 1) = gradr(:, :, 2)
            self_gradr(:, :, 2) = gradr(:, :, 2)
         ENDIF

         self_rhog(:, 1) = rhog( :, 2)
         self_rhog(:, 2) = rhog( :, 2)
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
      use gvect, only: ngm, nl, nlm
      use cell_base, only: ainv, tpiba, omega
      use cp_main_variables, only: drhog
      USE fft_interfaces, ONLY: fwfft, invfft
      USE fft_base,       ONLY: dfftp
!                 
      implicit none  
! input                   
      integer nspin
      real(DP)    :: gradr( dfftp%nnr, 3, nspin ), rhor( dfftp%nnr, nspin ), dexc( 3, 3 )
      complex(DP) :: rhog( ngm, nspin )
!
      complex(DP), allocatable:: v(:)
      complex(DP), allocatable:: x(:), vtemp(:)
      complex(DP) ::  ci, fp, fm
      integer :: iss, ig, ir, i,j
!
      allocate(v(dfftp%nnr))
      allocate(x(ngm))
      allocate(vtemp(ngm))
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
            v(ir)=CMPLX(gradr(ir,1,iss),0.d0,kind=DP)
         end do
         call fwfft('Dense',v, dfftp )
         do ig=1,ngm
            x(ig)=ci*tpiba*g(1,ig)*v(nl(ig))
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ngm
                     vtemp(ig) = omega*ci*CONJG(v(nl(ig)))*             &
     &                    tpiba*(-rhog(ig,iss)*g(i,ig)*ainv(j,1)+      &
     &                    g(1,ig)*drhog(ig,iss,i,j))
                  end do
                  dexc(i,j) = dexc(i,j) + DBLE(SUM(vtemp))*2.0d0
               end do
            end do
         endif
!
         do ir=1,dfftp%nnr
            v(ir)=CMPLX(gradr(ir,2,iss),gradr(ir,3,iss),kind=DP)
         end do
         call fwfft('Dense',v, dfftp )
!
         do ig=1,ngm
            fp=v(nl(ig))+v(nlm(ig))
            fm=v(nl(ig))-v(nlm(ig))
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*g(2,ig)*0.5d0*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*g(3,ig)*0.5d0*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ngm
                     fp=v(nl(ig))+v(nlm(ig))
                     fm=v(nl(ig))-v(nlm(ig))
                     vtemp(ig) = omega*ci*                              &
     &                    (0.5d0*CMPLX(DBLE(fp),-AIMAG(fm),kind=DP)*              &
     &                    tpiba*(-rhog(ig,iss)*g(i,ig)*ainv(j,2)+      &
     &                    g(2,ig)*drhog(ig,iss,i,j))+                  &
     &                    0.5d0*CMPLX(AIMAG(fp),DBLE(fm),kind=DP)*tpiba*          &
     &                    (-rhog(ig,iss)*g(i,ig)*ainv(j,3)+            &
     &                    g(3,ig)*drhog(ig,iss,i,j)))
                  end do
                  dexc(i,j) = dexc(i,j) + 2.0d0*DBLE(SUM(vtemp))
               end do
            end do
         endif
!     _________________________________________________________________
!     second part xc-potential: 1 inverse fft
!
         do ig=1,dfftp%nnr
            v(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ngm
            v(nl(ig))=x(ig)
            v(nlm(ig))=CONJG(x(ig))
         end do
         call invfft('Dense',v, dfftp )
         do ir=1,dfftp%nnr
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

!=----------------------------------------------------------------------------=!
!
!  This wrapper interface CP/FPMD to the PW xc and gga functionals
!
!  tested with PP/xctest.f90 code
!
!=----------------------------------------------------------------------------=!

subroutine exch_corr_wrapper(nnr, nspin, grhor, rhor, etxc, v, h)
  use kinds, only: DP
  use funct, only: dft_is_gradient, get_igcc, &
                   xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
  implicit none
  integer, intent(in) :: nnr
  integer, intent(in) :: nspin
  real(DP), intent(in) :: grhor( nnr, 3, nspin )
  real(DP) :: h( nnr, nspin, nspin )
  real(DP), intent(in) :: rhor( nnr, nspin )
  real(DP) :: v( nnr, nspin )
  real(DP) :: etxc
  integer :: ir, is, k
  real(DP) :: rup, rdw, ex, ec, vx(2), vc(2)
  real(DP) :: rh, grh2, zeta 
  real(DP) :: sx, sc, v1x, v2x, v1c, v2c
  real(DP) :: rhox, arhox, e2
  real(DP) :: grho2(2), arho, segno
  real(DP) :: v1xup, v1xdw, v2xup, v2xdw
  real(DP) :: v1cup, v1cdw
  real(DP) :: grhoup, grhodw, grhoud
  real(DP) :: v2cup, v2cdw, v2cud
  integer :: neg(3)
  real(DP), parameter :: epsr = 1.0d-10, epsg = 1.0d-10
  logical :: debug_xc = .false.
  logical :: igcc_is_lyp

  igcc_is_lyp = (get_igcc() == 3)
  !
  e2  = 1.0d0
  etxc = 0.0d0
  if( nspin == 1 ) then
     !
     ! spin-unpolarized case
     !
!$omp parallel do private( rhox, arhox, ex, ec, vx, vc ), reduction(+:etxc)
     do ir = 1, nnr
        rhox = rhor (ir, nspin)
        arhox = abs (rhox)
        if (arhox.gt.1.d-30) then
           CALL xc( arhox, ex, ec, vx(1), vc(1) )
           v(ir,nspin) = e2 * (vx(1) + vc(1) )
           etxc = etxc + e2 * (ex + ec) * rhox
        else
           v(ir,nspin) = 0.0D0
        endif
     enddo
!$omp end parallel do
     !
  else
     !
     ! spin-polarized case
     !
     neg (1) = 0
     neg (2) = 0
     neg (3) = 0
     do ir = 1, nnr
        rhox = rhor(ir,1) + rhor(ir,2)
        arhox = abs(rhox)
        if (arhox.gt.1.d-30) then
           zeta = ( rhor(ir,1) - rhor(ir,2) ) / arhox
           if (abs(zeta) .gt.1.d0) then
              neg(3) = neg(3) + 1
              zeta = sign(1.d0,zeta)
           endif
           ! WRITE(6,*) rhox, zeta
           if (rhor(ir,1) < 0.d0) neg(1) = neg(1) + 1
           if (rhor(ir,2) < 0.d0) neg(2) = neg(2) + 1
           call xc_spin (arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           do is = 1, nspin
              v(ir,is) = e2 * (vx(is) + vc(is) )
           enddo
           etxc = etxc + e2 * (ex + ec) * rhox
        else
           do is = 1, nspin
              v(ir,is) = 0.0D0
           end do
        endif
     enddo
  endif

  if( debug_xc ) then
    open(unit=17,form='unformatted')
    write(17) nnr, nspin
    write(17) rhor
    write(17) grhor
    close(17)
    debug_xc = .false.
  end if

  ! now come the corrections

  if( dft_is_gradient() ) then

    if (nspin == 1) then
       !
       !    This is the spin-unpolarised case
       !
!$omp parallel do &
!$omp private( is, grho2, arho, segno, sx, sc, v1x, v2x, v1c, v2c  ), reduction(+:etxc)
       do k = 1, nnr
          !
          grho2 (1) = grhor(k, 1, 1)**2 + grhor(k, 2, 1)**2 + grhor(k, 3, 1)**2
          arho = abs (rhor (k, 1) )
          segno = sign (1.d0, rhor (k, 1) )
          if (arho > epsr .and. grho2 (1) > epsg) then

             call gcxc (arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c)
             !
             ! first term of the gradient correction : D(rho*Exc)/D(rho)

             v (k, 1) = v (k, 1) + e2 * (v1x + v1c)

             ! HERE h contains D(rho*Exc)/D(|grad rho|) / |grad rho|
             !
             h (k, 1, 1) = e2 * (v2x + v2c) 
             etxc = etxc + e2 * (sx + sc) * segno

          else
             h (k, 1, 1) = 0.d0
          endif
          !
       end do
!$omp end parallel do
       !
    else
       !
       !    spin-polarised case
       !
       do k = 1, nnr
          do is = 1, nspin
             grho2 (is) = grhor(k, 1, is)**2 + grhor(k, 2, is)**2 + grhor(k, 3, is)**2
          enddo
          rup = rhor (k, 1)
          rdw = rhor (k, 2)
          call gcx_spin ( rup, rdw, grho2 (1), grho2 (2), sx, v1xup, v1xdw, v2xup, v2xdw)
          !
          rh = rhor (k, 1) + rhor (k, 2)
          !
          if (rh.gt.epsr) then
             if( igcc_is_lyp ) then
                grhoup = grhor(k,1,1)**2 + grhor(k,2,1)**2 + grhor(k,3,1)**2
                grhodw = grhor(k,1,2)**2 + grhor(k,2,2)**2 + grhor(k,3,2)**2
                grhoud =          grhor(k,1,1)* grhor(k,1,2)
                grhoud = grhoud + grhor(k,2,1)* grhor(k,2,2)
                grhoud = grhoud + grhor(k,3,1)* grhor(k,3,2)
                call gcc_spin_more(rup, rdw, grhoup, grhodw, grhoud, sc, &
                     v1cup, v1cdw, v2cup, v2cdw, v2cud)
             else
                zeta = (rhor (k, 1) - rhor (k, 2) ) / rh
                !
                grh2 = (grhor (k, 1, 1) + grhor (k, 1, 2) ) **2 + &
                       (grhor (k, 2, 1) + grhor (k, 2, 2) ) **2 + &
                       (grhor (k, 3, 1) + grhor (k, 3, 2) ) **2
                call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
                v2cup = v2c
                v2cdw = v2c
                v2cud = v2c
             end if
          else
             sc = 0.d0
             v1cup = 0.d0
             v1cdw = 0.d0
             v2c = 0.d0
             v2cup = 0.0d0
             v2cdw = 0.0d0
             v2cud = 0.0d0
          endif
          !
          ! first term of the gradient correction : D(rho*Exc)/D(rho)
          !
          v (k, 1) = v (k, 1) + e2 * (v1xup + v1cup)
          v (k, 2) = v (k, 2) + e2 * (v1xdw + v1cdw)
          !
          ! HERE h contains D(rho*Exc)/D(|grad rho|) / |grad rho|
          !
          h (k, 1, 1) = e2 * (v2xup + v2cup)  ! Spin UP-UP
          h (k, 1, 2) = e2 *          v2cud   ! Spin UP-DW
          h (k, 2, 1) = e2 *          v2cud   ! Spin DW-UP
          h (k, 2, 2) = e2 * (v2xdw + v2cdw)  ! Spin DW-DW
          !
          etxc = etxc + e2 * (sx + sc)
          !
          !
       enddo
       !
    endif
    !
  end if

  return
end subroutine exch_corr_wrapper


!=----------------------------------------------------------------------------=!
!
!  For CP we need a further small interface subroutine
!
!=----------------------------------------------------------------------------=!

subroutine exch_corr_cp(nnr,nspin,grhor,rhor,etxc)
  use kinds, only: DP
  use funct, only: dft_is_gradient
  implicit none
  integer, intent(in) :: nnr
  integer, intent(in) :: nspin
  real(DP) :: grhor( nnr, 3, nspin )
  real(DP) :: rhor( nnr, nspin )
  real(DP) :: etxc
  integer :: k, ipol
  real(DP) ::  grup, grdw
  real(DP), allocatable :: v(:,:)
  real(DP), allocatable :: h(:,:,:)
  !
  allocate( v( nnr, nspin ) )
  if( dft_is_gradient() ) then
    allocate( h( nnr, nspin, nspin ) )
  else
    allocate( h( 1, 1, 1 ) )
  endif
  !
  call exch_corr_wrapper(nnr,nspin,grhor,rhor,etxc,v,h)

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
              grhor (k, ipol, 1) = h (k, 1, 1) * grhor (k, ipol, 1)
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
             grup = grhor (k, ipol, 1)
             grdw = grhor (k, ipol, 2)
             grhor (k, ipol, 1) = h (k, 1, 1) * grup + h (k, 1, 2) * grdw
             grhor (k, ipol, 2) = h (k, 2, 2) * grdw + h (k, 2, 1) * grup
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
