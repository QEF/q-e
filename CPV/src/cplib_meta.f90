!
! Copyright (C) 2005-2010 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      subroutine dforce_meta (c,ca,df,da, psi,iss1,iss2,fi,fip)
!-----------------------------------------------------------------------
      !! This subroutine computes the generalized force \(\text{df}=
      !! \text{cmplx}(\text{dfr},\text{dfi})\) acting on the i-th electron
      !! state at the gamma point of the Brillouin zone represented by the
      !! vector c=cmplx(cr,ci)\(c=\text{cmplx}(\text{cr},\text{ci})\).
      !! Contribution from metaGGA.
      !
!
!          contribution from metaGGA
      use kinds, only: dp
      use gvect,     only : g
      use gvecw,                  only : ngw
      use cell_base,              only : tpiba2
      USE metagga_cp,             ONLY : kedtaus
      USE fft_interfaces,         ONLY : fwfft, invfft
      USE fft_base,               ONLY : dffts
      USE fft_helper_subroutines, ONLY : fftx_c2psi_gamma, fftx_psi2c_gamma
!
      implicit none
!
      complex(dp) c(ngw), ca(ngw), df(ngw), da(ngw),psi(dffts%nnr)
      complex(dp), allocatable ::  dc(:,:), dca(:)
      integer iss1, iss2
      real(dp) fi, fip
! local variables
      integer ir,ig, ipol !metagga
      complex(dp) fp,fm,ci
!
!
      ci=(0.0d0,1.0d0)
      allocate( dc( ngw,1 ) )
      allocate( dca( ngw ) )
!
         do ipol = 1, 3
            
            dc(:,1)  = ci*g(ipol,1:ngw)*c(:)
            dca(:) = ci*g(ipol,1:ngw)*ca(:)
            CALL fftx_c2psi_gamma( dffts, psi, dc, dca )
            CALL invfft( 'Wave', psi, dffts )

!           on smooth grids--> grids for charge density

            do ir=1, dffts%nnr
               psi(ir) = CMPLX (kedtaus(ir,iss1)*DBLE(psi(ir)), &
                                kedtaus(ir,iss2)*AIMAG(psi(ir)),kind=DP)
            end do
            call fwfft('Wave',psi, dffts )
            CALL fftx_psi2c_gamma( dffts, psi, dc(:,1:1), vout2=dca )
            do ig=1,ngw
               df(ig)= df(ig) - ci*fi *tpiba2*g(ipol,ig) * dc(ig,1)
               da(ig)= da(ig) - ci*fip*tpiba2*g(ipol,ig) * dca(ig) 
            end do
         end do

      deallocate( dc )
      deallocate( dca )
!
      return
    end subroutine dforce_meta
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      subroutine kedtauofr_meta (c)
!-----------------------------------------------------------------------
      !! Calculates the kinetic energy density (metaGGA case).
      !
      use kinds, only: dp
      use control_flags, only: tpre
      use gvecw, only: ngw
      use gvect, only: g
      use cell_base, only : omega, tpiba, ainv
      use electrons_base, only: nx => nbspx, n => nbsp, f, ispin, nspin
      use constants, only: pi, fpi
!
      use dener
      use metagga_cp, ONLY : kedtaur, kedtaus, kedtaug, crosstaus, gradwfc, &
                             dkedtaus
      USE fft_interfaces, ONLY: fwfft, invfft
      USE fft_base,       ONLY: dffts, dfftp
      USE fft_helper_subroutines, ONLY : fftx_c2psi_gamma
      USE fft_rho
      
      implicit none

      complex(dp) :: c(ngw,nx)
      complex(dp), allocatable :: psis( : )
      complex(dp), allocatable :: dc( :, : ), dca( : )

! local variables
      integer iss, isup, isdw, iss1, iss2, ios, i, ir, ig
      integer ipol, ix,iy, ipol2xy(3,3)
      real(dp) sa1, sa2
      complex(dp) ci,fp,fm
!
      ALLOCATE( psis( dffts%nnr ) )
      ALLOCATE( dc( ngw, 1 ) )
      ALLOCATE( dca( ngw ) )
!
      ci=(0.0d0,1.0d0)
      kedtaur(:,:)=0.d0
      kedtaus(:,:)=0.d0
      kedtaug(:,:)=(0.d0,0.d0)
      if(tpre) crosstaus(:,:,:)=0.d0

!
!    
!    warning! trhor and thdyn are not compatible yet!   
!
!     important: if n is odd then nx must be .ge.n+1 and c(*,n+1)=0.
! 
      if (mod(n,2).ne.0) then
         c(1:ngw,n+1)=(0.d0,0.d0)
      endif
         !
      do i=1,n,2
         iss1=ispin(i)
         sa1=f(i)/omega
         if (i.ne.n) then
            iss2=ispin(i+1)
            sa2=f(i+1)/omega
         else
            iss2=iss1
            sa2=0.0d0
         end if

         do ipol = 1, 3
            psis( : ) = (0.d0,0.d0)
            ! gradient of wfc in real space
            do ig=1,ngw
               dc( ig, 1 )  = ci * tpiba * g(ipol,ig) * c(ig,i)
               dca( ig ) = ci * tpiba * g(ipol,ig) * c(ig,i+1)
            end do
            CALL fftx_c2psi_gamma( dffts, psis, dc, dca )
            call invfft('Wave',psis, dffts )
            ! on smooth grids--> grids for charge density
            do ir=1, dffts%nnr
               kedtaus(ir,iss1)=kedtaus(ir,iss1)+0.5d0*sa1*DBLE(psis(ir))**2
               kedtaus(ir,iss2)=kedtaus(ir,iss2)+0.5d0*sa2*AIMAG(psis(ir))**2
            end do
            if(tpre) then
               do ir=1, dffts%nnr
                  gradwfc(ir,ipol)=psis(ir)
               end do
            end if
         end do
         if(tpre) then
            ipol=1
            do ix=1,3
               do iy=1,ix
                  ipol2xy(ix,iy)=ipol
                  ipol2xy(iy,ix)=ipol
                  do ir=1,dffts%nnr
                     crosstaus(ir,ipol,iss1) = crosstaus(ir,ipol,iss1) +&
                          sa1*DBLE(gradwfc(ir,ix))*DBLE(gradwfc(ir,iy))
                     crosstaus(ir,ipol,iss2) = crosstaus(ir,ipol,iss2) +&
                          sa2*AIMAG(gradwfc(ir,ix))*AIMAG(gradwfc(ir,iy))
                  end do
                  ipol=ipol+1
               end do
            end do
         end if

            !        d kedtaug / d h
         if(tpre) then
            do iss=1,nspin
               do ix=1,3
                  do iy=1,3
                     do ir=1,dffts%nnr
                        dkedtaus(ir,ix,iy,iss)=-kedtaus(ir,iss)*ainv(iy,ix)&
                             -crosstaus(ir,ipol2xy(1,ix),iss)*ainv(iy,1)&
                             -crosstaus(ir,ipol2xy(2,ix),iss)*ainv(iy,2)&
                             -crosstaus(ir,ipol2xy(3,ix),iss)*ainv(iy,3)
                     end do
                  end do
               end do
            end do
         end if  !end metagga
         !
      end do

      DEALLOCATE( dca )
      DEALLOCATE( dc )
      DEALLOCATE( psis )

!     kinetic energy density (kedtau) in g-space (kedtaug)

      CALL rho_r2g( dffts, kedtaus, kedtaug )
!
      kedtaug(dffts%ngm+1:,:) = 0.0d0

      CALL rho_g2r( dfftp, kedtaug, kedtaur )
!
      return
    end subroutine kedtauofr_meta
!
!
!-----------------------------------------------------------------------
      subroutine vofrho_meta ( )
!-----------------------------------------------------------------------
      !! This subroutine computes the one-particle potential v in real space,
      !! the total energy etot, the forces fion acting on the ions, the
      !! derivative of total energy to cell parameters h.
      !
      !! rhor input : electronic charge on dense real space grid
      !! (plus core charge if present);  
      !! rhog input : electronic charge in g space (up to density cutoff);  
      !! rhos input : electronic charge on smooth real space grid;  
      !! rhor output: total potential on dense real space grid;  
      !! rhos output: total potential on smooth real space grid.
      !
      use kinds, only: dp
      use control_flags, only: thdyn, tpre, tfor, tprnfor
      use io_global, only: stdout
      use ions_base, only: nsp, na, nat
      use cell_base, only: omega
      use electrons_base, only: nspin
      use constants, only: pi, fpi
      use energies, only: etot, eself, enl, ekin, epseu, esr, eht, exc
      use local_pseudo, only: vps, rhops
      use core
      use smallbox_gvec
      use dener
      use mp,      ONLY : mp_sum
      use mp_global, ONLY : intra_bgrp_comm
      use metagga_cp, ONLY : kedtaur, kedtaug, kedtaus, dkedtaus
      USE fft_interfaces, ONLY: fwfft, invfft
      USE fft_base,       ONLY: dffts, dfftp
      USE fft_rho
!
      implicit none
!
      integer iss, isup, isdw, ig, ir,i,j,k,is, ia
      real(dp) dkedxc(3,3) !metagga
      complex(dp)  fp, fm, ci
!
      ci=(0.d0,1.d0)
!
!     ===================================================================
!      calculation exchange and correlation energy and potential
!     -------------------------------------------------------------------
!      if (nlcc.gt.0) call add_cc(rhoc,rhog,rhor)
!
#ifdef VARIABLECELL
!      call exch_corr_h(nspin,rhog,rhor,exc,dxc)
#else
!      call exch_corr(nspin,rhog,rhor,exc)
#endif
!
!     rhor contains the xc potential in r-space
!
!     ===================================================================
!     fourier transform of xc potential to g-space (dense grid)
!     -------------------------------------------------------------------
!
      CALL rho_r2g( dfftp, kedtaur, kedtaug )
!
      CALL rho_g2r ( dffts, kedtaug, kedtaus )

      !calculate dkedxc in real space on smooth grids  !metagga
      if(tpre) then
         do iss=1,nspin
            do j=1,3
               do i=1,3
                  dkedxc(i,j)=0.d0
                  do ir=1,dffts%nnr
                     !2.d0 : because kedtau = 0.5d0 d_Exc/d_kedtau
                      dkedxc(i,j)= dkedxc(i,j)+kedtaus(ir,iss)*2.d0*&
                           dkedtaus(ir,i,j,iss)
                   end do
                end do
             end do
          end do
          !
          call mp_sum( dkedxc, intra_bgrp_comm )
          !
          do j=1,3
             do i=1,3
                dxc(i,j) = dxc(i,j) + omega/(dffts%nr1*dffts%nr2*dffts%nr3)*dkedxc(i,j)
             end do
          end do
       end if        
       return
     end subroutine vofrho_meta
!-----------------------------------------------------------------------
