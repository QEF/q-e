!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      subroutine vofrho2(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,           &
     &     ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion,v0s,vhxcs)
!-----------------------------------------------------------------------
!     computes: the one-particle potential v in real space,
!               the total energy etot, the total free energy atot
!               the forces fion acting on the ions,
!               the derivative of total energy to cell parameters h
!     rhor input : electronic charge on dense real space grid
!                  (plus core charge if present)
!     rhog input : electronic charge in g space (up to density cutoff)
!     rhos input : electronic charge on smooth real space grid
!     rhor output: total potential on dense real space grid
!     rhos output: total potential on smooth real space grid
!     v0s output : total local pseudopotential on smooth real space grid
!     vhxcs out  : hartree-xc potential on smooth real space grid
!
      use control_flags, only: iprint, tvlocw, iprsta, thdyn, tpre, tfor
      use io_global, only: stdout
      use parameters, only: natx, nsx
      use ions_base, only: nas => nax, nsp, na, nat
      use gvecs
      use reciprocal_vectors, only: g
      use recvecs_indexes, only: np, nm
      use gvecp, only: ngm
      use cell_base, only: omega, tpiba2
      use cell_base, only: a1, a2, a3
      use reciprocal_vectors, only: ng0 => gstart
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base, only: nspin
      use constants, only: pi, fpi
      use energies, only: etot, eself, enl, ekin, epseu, esr, eht, exc,   &
     &        atot, egrand, entropy 

      use local_pseudo, only: vps, rhops
      use core
      use gvecb
!
      use dener
      use derho
      use mp, only: mp_sum
!
      implicit none
!
      logical tlast,tfirst
      integer nfi
      real(kind=8) rhor(nnr,nspin), rhos(nnrsx,nspin), fion(3,natx),&
     &             v0s(nnrsx), vhxcs(nnrsx,nspin)
      real(kind=8)  rhoc(nnr), tau0(3,natx)
      complex(kind=8) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),     &
     &                ei3(-nr3:nr3,nat), eigrb(ngb,nat),        &
     &                rhog(ngm,nspin), sfac(ngs,nsp)
!
      integer irb(3,nat), iss, isup, isdw, ig, ir,i,j,k,is, ia
      real(kind=8) fion1(3,natx), vave, ebac, wz, eh
      complex(kind=8)  fp, fm, ci
      complex(kind=8), allocatable :: v(:), vs(:)
      complex(kind=8), allocatable :: rhotmp(:), vtemp(:), drhotmp(:,:,:)
      
!     complex(kind=8), allocatable:: vtemp1(:,:)
!
      ci=(0.,1.)
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz = 2.0
      allocate( v ( nnr   ) )
      allocate( vs( nnrsx ) )

      allocate(vtemp(ngm))
      allocate(rhotmp(ngm))
      if (tpre) allocate(drhotmp(ngm,3,3))
!
!     first routine in which fion is calculated: annihilation
!
      fion =0.d0
      fion1=0.d0
!
!     ===================================================================
!     forces on ions, ionic term in real space
!     -------------------------------------------------------------------
      if(tfor.or.tfirst.or.thdyn) call force_ion(tau0,esr,fion,dsr)
!
      if(nspin.eq.1) then
         iss=1
         do ig=1,ngm
            rhotmp(ig)=rhog(ig,iss)
         end do
         if(tpre)then
            do j=1,3
               do i=1,3
                  do ig=1,ngm
                     drhotmp(ig,i,j)=drhog(ig,iss,i,j)
                  enddo
               enddo
            enddo
         endif
      else
         isup=1
         isdw=2
         do ig=1,ngm
            rhotmp(ig)=rhog(ig,isup)+rhog(ig,isdw)
         end do
         if(tpre)then
            do i=1,3
               do j=1,3
                  do ig=1,ngm
                     drhotmp(ig,i,j) = drhog(ig,isup,i,j) +           &
     &                                 drhog(ig,isdw,i,j)
                  enddo
               enddo
            enddo
         endif
      end if
!     ===================================================================
!     calculation local potential energy
!     -------------------------------------------------------------------
      vtemp=(0.,0.)
      do is=1,nsp
         do ig=1,ngs
            vtemp(ig)=vtemp(ig)+conjg(rhotmp(ig))*sfac(ig,is)*vps(ig,is)
         end do
      end do
!
      epseu = wz * real( SUM( vtemp( 1:ngs ) ) )
      if (ng0.eq.2) epseu=epseu-vtemp(1)

      call mp_sum( epseu )

      epseu=epseu*omega
!
      if(tpre) call denps(rhotmp,drhotmp,sfac,vtemp,dps)
!
!     ===================================================================
!     calculation hartree energy
!     -------------------------------------------------------------------
      do is=1,nsp
         do ig=1,ngs
            rhotmp(ig)=rhotmp(ig)+sfac(ig,is)*rhops(ig,is)
         end do
      end do
      if (ng0.eq.2) vtemp(1)=0.0
      do ig=ng0,ngm
         vtemp(ig)=conjg(rhotmp(ig))*rhotmp(ig)/g(ig)
      end do
!
      eh=real( SUM( vtemp( 1:ngm ) ) ) *wz*0.5*fpi/tpiba2

      call mp_sum( eh )

      if(tpre) call denh(rhotmp,drhotmp,sfac,vtemp,eh,dh)
      if(tpre) deallocate(drhotmp)
!     ===================================================================
!     forces on ions, ionic term in reciprocal space
!     -------------------------------------------------------------------
      if(tfor.or.thdyn)                                                  &
     &    call force_ps(rhotmp,rhog,vtemp,ei1,ei2,ei3,fion1)
!     ===================================================================
!     calculation hartree + local pseudo potential
!     -------------------------------------------------------------------
!
      if (ng0.eq.2) vtemp(1)=(0.,0.)

      do ig=ng0,ngm
         vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      end do

      do is=1,nsp
         do ig=1,ngs
           vtemp(ig)=vtemp(ig)+sfac(ig,is)*vps(ig,is)
         end do
      end do
     
      vs = 0.0d0
      do is=1,nsp
         do ig=1,ngs
           vs(nms(ig))=vs(nms(ig))+conjg(sfac(ig,is)*vps(ig,is))
           vs(nps(ig))=vs(nps(ig))+sfac(ig,is)*vps(ig,is)
         end do
      end do
!
      call ivffts(vs,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!
      do ir=1,nnrsx
         v0s(ir)=real(vs(ir))
      end do


!
!     vtemp = v_loc(g) + v_h(g)
!
!     ===================================================================
!      calculation exchange and correlation energy and potential
!     -------------------------------------------------------------------
      if ( nlcc_any ) call add_cc(rhoc,rhog,rhor)
!
      call exch_corr_h(nspin,rhog,rhor,exc,dxc)
!
!     rhor contains the xc potential in r-space
!
!     ===================================================================
!     fourier transform of xc potential to g-space (dense grid)
!     -------------------------------------------------------------------
!
      if(nspin.eq.1) then
         iss=1
         do ir=1,nnr
            v(ir)=cmplx(rhor(ir,iss),0.0)
         end do
!
!     v_xc(r) --> v_xc(g)
!
         call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         do ig=1,ngm
            rhog(ig,iss)=vtemp(ig)+v(np(ig))
         end do
!
!     v_tot(g) = (v_tot(g) - v_xc(g)) +v_xc(g)
!     rhog contains the total potential in g-space
!
      else
         isup=1
         isdw=2
         do ir=1,nnr
            v(ir)=cmplx(rhor(ir,isup),rhor(ir,isdw))
         end do
         call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ig=1,ngm
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            rhog(ig,isup)=vtemp(ig)+0.5*cmplx( real(fp),aimag(fm))
            rhog(ig,isdw)=vtemp(ig)+0.5*cmplx(aimag(fp),-real(fm))
         end do
      endif
!
!     rhog contains now the total (local+Hartree+xc) potential in g-space
!
      if(tfor) then
         if ( nlcc_any ) call force_cc(irb,eigrb,rhor,fion1)
   
         call mp_sum( fion1 )
!
!    add g-space ionic and core correction contributions to fion
!
         fion = fion + fion1

      end if
!     ===================================================================
!     fourier transform of total potential to r-space (dense grid)
!     -------------------------------------------------------------------
      v = 0.0d0
      if(nspin.eq.1) then
         iss=1
         do ig=1,ngm
            v(np(ig))=rhog(ig,iss)
            v(nm(ig))=conjg(rhog(ig,iss))
         end do
!
!     v(g) --> v(r)
!
         call invfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         do ir=1,nnr
            rhor(ir,iss)=real(v(ir))
         end do
!
!     calculation of average potential
!
         vave= SUM( rhor(1:nnr,iss) ) /dfloat(nr1*nr2*nr3)
      else
         isup=1
         isdw=2
         do ig=1,ngm
            v(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            v(nm(ig))=conjg(rhog(ig,isup)) +ci*conjg(rhog(ig,isdw))
         end do
!
         call invfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            rhor(ir,isup)= real(v(ir))
            rhor(ir,isdw)=aimag(v(ir))
         end do
!
!     calculation of average potential
!
         vave=( SUM( rhor(1:nnr,isup) ) + SUM( rhor(1:nnr,isdw) ) )      &
     &        /2.0/dfloat(nr1*nr2*nr3)
      endif

      call mp_sum( vave )

!     ===================================================================
!     fourier transform of total potential to r-space (smooth grid)
!     -------------------------------------------------------------------
    
      vs = 0.0d0

      if(nspin.eq.1)then
         iss=1
         do ig=1,ngs
            vs(nms(ig))=conjg(rhog(ig,iss))
            vs(nps(ig))=rhog(ig,iss)
         end do
!
         call ivffts(vs,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!
         do ir=1,nnrsx
            vhxcs(ir,iss)=real(vs(ir))-v0s(ir)
            rhos(ir,iss)=real(vs(ir))
         end do
      else
         isup=1
         isdw=2
         do ig=1,ngs
            vs(nps(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            vs(nms(ig))=conjg(rhog(ig,isup)) +ci*conjg(rhog(ig,isdw))
         end do 
         call ivffts(vs,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
         do ir=1,nnrsx
            
            vhxcs(ir,isup)= real(vs(ir))-v0s(ir)
            vhxcs(ir,isdw)=aimag(vs(ir))-v0s(ir)

            rhos(ir,isup)= real(vs(ir))
            rhos(ir,isdw)=aimag(vs(ir))
         end do
      endif

      ebac=0.0
!
      eht=eh*omega+esr-eself
!
!     etot is the total energy ; ekin, enl were calculated in rhoofr
!
      etot=ekin+eht+epseu+enl+exc+ebac
      atot=etot+entropy+egrand

      if(tpre) detot=dekin+dh+dps+denl+dxc+dsr
!
      if(tvlocw.and.tlast)then
#ifdef __PARA
         call write_rho(46,nspin,rhor)
#else
         write(46) ((rhor(ir,iss),ir=1,nnr),iss=1,nspin)
#endif
      endif
!
      deallocate(rhotmp)
      deallocate(vtemp)
      deallocate( vs )
      deallocate( v  )
!
!
      if((nfi.eq.0).or.tfirst.or.tlast) goto 999
      if(mod(nfi-1,iprint).ne.0 ) return
!
999  write(stdout,1) etot,atot,ekin,eht,esr,eself,epseu,enl,exc,entropy,&
               egrand,vave
    1 format(//' VOFRHO2:                                    '/         &
     &         '                total energy = ',f18.9,' a.u.'/         &
     &         '           total free energy = ',f18.9,' a.u.'/         &
     &         '              kinetic energy = ',f18.9,' a.u.'/         &
     &         '        electrostatic energy = ',f18.9,' a.u.'/         &
     &         '                         esr = ',f18.9,' a.u.'/         &
     &         '                       eself = ',f18.9,' a.u.'/         &
     &         '      pseudopotential energy = ',f18.9,' a.u.'/         &
     &         '  n-l pseudopotential energy = ',f18.9,' a.u.'/         &
     &         ' exchange-correlation energy = ',f18.9,' a.u.'/         &
     &         '               entropy (-TS) = ',f18.9,' a.u.'/         &
     &         '                      egrand = ',f18.9,' a.u.'/         &
     &         '           average potential = ',f18.9,' a.u.'//)
!
      if(tpre)then
         write (stdout,*) "cell parameters h"
         write (stdout,5555) (a1(i),a2(i),a3(i),i=1,3)
         write (stdout,*)
         write (stdout,*) "derivative of e(tot)"
         write (stdout,5555) ((detot(i,j),j=1,3),i=1,3)
         write (stdout,*)
         if(tpre.and.iprsta.ge.2) then
            write (stdout,*) "derivative of e(kin)"
            write (stdout,5555) ((dekin(i,j),j=1,3),i=1,3)
            write (stdout,*) "derivative of e(electrostatic)"
            write (stdout,5555) (((dh(i,j)+dsr(i,j)),j=1,3),i=1,3)
            write (stdout,*) "derivative of e(h)"
            write (stdout,5555) ((dh(i,j),j=1,3),i=1,3)
            write (stdout,*) "derivative of e(sr)"
            write (stdout,5555) ((dsr(i,j),j=1,3),i=1,3)
            write (stdout,*) "derivative of e(ps)"
            write (stdout,5555) ((dps(i,j),j=1,3),i=1,3)
            write (stdout,*) "derivative of e(nl)"
            write (stdout,5555) ((denl(i,j),j=1,3),i=1,3)
            write (stdout,*) "derivative of e(xc)"
            write (stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
         endif
      endif
5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
!
      return
      end subroutine vofrho2

