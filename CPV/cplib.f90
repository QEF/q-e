!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!
!-----------------------------------------------------------------------
      subroutine atomic_wfc(eigr,n_atomic_wfc,wfc)
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space
!
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart, g, gx
      use ions_base, only: nsp, na, nat
      use cell_base, only: tpiba
      use atom, only: nchi, lchi, mesh, r, chi, rab
!
      implicit none
      integer, intent(in) :: n_atomic_wfc
      complex(kind=8), intent(in) ::  eigr(ngw,nat)
      complex(kind=8), intent(out):: wfc(ngw,n_atomic_wfc)
!
      integer :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc, isa
      real(kind=8), allocatable::  ylm(:,:), q(:), jl(:), vchi(:),      &
     &     chiq(:)
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      do is = 1,nsp
         do nb = 1, nchi(is)
            lmax_wfc = max (lmax_wfc, lchi (nb, is) )
         enddo
      enddo
      allocate(ylm(ngw,(lmax_wfc+1)**2))
      call ylmr2 ((lmax_wfc+1)**2, ngw, gx, g, ylm)
      ndm = MAXVAL(mesh(1:nsp))
      allocate(jl(ndm), vchi(ndm))
      allocate(q(ngw), chiq(ngw))
!
      do i=1,ngw
         q(i) = sqrt(g(i))*tpiba
      end do
!
      natwfc=0
      isa   = 0
      do is=1,nsp
!
!   radial fourier transform of the chi functions
!   NOTA BENE: chi is r times the radial part of the atomic wavefunction
!
         do nb = 1,nchi(is)
            l = lchi(nb,is)
            do i=1,ngw
               call sph_bes (mesh(is), r(1,is), q(i), l, jl)
               do ir=1,mesh(is)
                  vchi(ir) = chi(ir,nb,is)*r(ir,is)*jl(ir)
               enddo
               call simpson_cp90(mesh(is),vchi,rab(1,is),chiq(i))
            enddo
!
!   multiply by angular part and structure factor
!   NOTA BENE: the factor i^l MUST be present!!!
!
            do m = 1,2*l+1
               lm = l**2 + m
               do ia = 1 + isa, na(is) + isa
                  natwfc = natwfc + 1
                  wfc(:,natwfc) = (0.d0,1.d0)**l * eigr(:,ia)* ylm(:,lm)*chiq(:)
               enddo
            enddo
         enddo
         isa = isa + na(is)
      enddo
!
      if (natwfc.ne.n_atomic_wfc)                                       &
     &     call errore('atomic_wfc','unexpected error',natwfc)
!
      deallocate(q, chiq, vchi, jl, ylm)
!
      return
      end subroutine atomic_wfc
!
!-----------------------------------------------------------------------
      subroutine box2grid(irb,nfft,qv,vr)
!-----------------------------------------------------------------------
!
! add array qv(r) on box grid to array vr(r) on dense grid
! irb   : position of the box in the dense grid
! nfft=1  add      real part of qv(r) to real part of array vr(r) 
! nfft=2  add imaginary part of qv(r) to real part of array vr(r) 
!
      use parameters, only: natx, nsx
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nnrx
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      use para_mod
      implicit none
      integer, intent(in):: nfft, irb(3)
      real(kind=8), intent(in):: qv(2,nnrb)
      complex(kind=8), intent(inout):: vr(nnr)
!
      integer ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig

      if(nfft.le.0.or.nfft.gt.2) call errore('box2grid','wrong data',nfft)

      do ir3=1,nr3b
         ibig3=irb(3)+ir3-1
         ibig3=1+mod(ibig3-1,nr3)
         if(ibig3.lt.1.or.ibig3.gt.nr3)                                 &
     &        call errore('box2grid','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%ipp(me)
         if ( ibig3 .gt. 0 .and. ibig3 .le. ( dfftp%npp(me) ) ) then
            do ir2=1,nr2b
               ibig2=irb(2)+ir2-1
               ibig2=1+mod(ibig2-1,nr2)
               if(ibig2.lt.1.or.ibig2.gt.nr2)                           &
     &              call errore('box2grid','ibig2 wrong',ibig2)
               do ir1=1,nr1b
                  ibig1=irb(1)+ir1-1
                  ibig1=1+mod(ibig1-1,nr1)
                  if(ibig1.lt.1.or.ibig1.gt.nr1)                        &
     &                 call errore('box2grid','ibig1 wrong',ibig1)
                  ibig=ibig1+(ibig2-1)*nr1x+(ibig3-1)*nr1x*nr2x
                  ir=ir1+(ir2-1)*nr1bx+(ir3-1)*nr1bx*nr2bx
                  vr(ibig) = vr(ibig)+qv(nfft,ir)
               end do
            end do
         end if
      end do
!
      return
      end subroutine box2grid
!
!-----------------------------------------------------------------------
      subroutine box2grid2(irb,qv,v)
!-----------------------------------------------------------------------
!
! add array qv(r) on box grid to array v(r) on dense grid
! irb   : position of the box in the dense grid
!
      use parameters, only: nsx, natx
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nnrx
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      use para_mod
      implicit none
      integer, intent(in):: irb(3)
      complex(kind=8), intent(in):: qv(nnrb)
      complex(kind=8), intent(inout):: v(nnr)
!
      integer ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig

      do ir3=1,nr3b
         ibig3=irb(3)+ir3-1
         ibig3=1+mod(ibig3-1,nr3)
         if(ibig3.lt.1.or.ibig3.gt.nr3)                                 &
     &        call errore('box2grid2','ibig3 wrong',ibig3)
         ibig3=ibig3-dfftp%ipp(me)
         if (ibig3.gt.0.and.ibig3.le. dfftp%npp(me) ) then
            do ir2=1,nr2b
               ibig2=irb(2)+ir2-1
               ibig2=1+mod(ibig2-1,nr2)
               if(ibig2.lt.1.or.ibig2.gt.nr2)                           &
     &              call errore('box2grid2','ibig2 wrong',ibig2)
               do ir1=1,nr1b
                  ibig1=irb(1)+ir1-1
                  ibig1=1+mod(ibig1-1,nr1)
                  if(ibig1.lt.1.or.ibig1.gt.nr1)                        &
     &                 call errore('box2grid2','ibig1 wrong',ibig1)
                  ibig=ibig1+(ibig2-1)*nr1x+(ibig3-1)*nr1x*nr2x
                  ir=ir1+(ir2-1)*nr1bx+(ir3-1)*nr1bx*nr2bx
                  v(ibig) = v(ibig)+qv(ir)
               end do
            end do
         end if
      end do

      return
      end subroutine box2grid2
!
!-----------------------------------------------------------------------
      real(kind=8) function boxdotgrid(irb,nfft,qv,vr)
!-----------------------------------------------------------------------
!
! Calculate \sum_i qv(r_i)*vr(r_i)  with r_i on box grid
! array qv(r) is defined on box grid, array vr(r)on dense grid
! irb   : position of the box in the dense grid
! nfft=1 (2): use real (imaginary) part of qv(r)
! Parallel execution: remember to sum the contributions from other nodes
!
      use parameters, only: nsx, natx
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nnr => nnrx
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nnrb => nnrbx
      use para_mod
      implicit none
      integer, intent(in):: nfft, irb(3)
      real(kind=8), intent(in):: qv(2,nnrb), vr(nnr)
!
      integer ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig
!
!
      if(nfft.le.0.or.nfft.gt.2) call errore('box2grid','wrong data',nfft)

      boxdotgrid=0.d0

      do ir3=1,nr3b
         ibig3=irb(3)+ir3-1
         ibig3=1+mod(ibig3-1,nr3)
         ibig3=ibig3-dfftp%ipp(me)
         if (ibig3.gt.0.and.ibig3.le. dfftp%npp(me) ) then
            do ir2=1,nr2b
               ibig2=irb(2)+ir2-1
               ibig2=1+mod(ibig2-1,nr2)
               do ir1=1,nr1b
                  ibig1=irb(1)+ir1-1
                  ibig1=1+mod(ibig1-1,nr1)
                  ibig=ibig1 + (ibig2-1)*nr1x + (ibig3-1)*nr1x*nr2x
                  ir  =ir1 + (ir2-1)*nr1bx + (ir3-1)*nr1bx*nr2bx
                  boxdotgrid = boxdotgrid + qv(nfft,ir)*vr(ibig)
               end do
            end do
         endif
      end do

      return
      end function boxdotgrid
!
!-----------------------------------------------------------------------
      subroutine calbec (nspmn,nspmx,eigr,c,bec)
!-----------------------------------------------------------------------
!     this routine calculates array bec
!
!        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
!
!     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
!
      use ions_base, only: na, nas => nax, nat
      use io_global, only: stdout
      use cvan, only: ish
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use control_flags, only: iprint, iprsta
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb
!
      implicit none
      integer nspmn, nspmx
      real(kind=8)  bec(nhsa,n)
      complex(kind=8) c(ngw,n), eigr(ngw,nat)
! local variables
      integer is, ia, i , iv
!
!
      call start_clock( 'calbec' )
      call nlsm1(n,nspmn,nspmx,eigr,c,bec)
!
      if (iprsta.gt.2) then
         WRITE( stdout,*)
         do is=1,nspmx
            if(nspmx.gt.1) then
               WRITE( stdout,'(33x,a,i4)') ' calbec: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &              ((bec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
            else
               do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calbec: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &             ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
               end do
            end if
         end do
      endif
      call stop_clock( 'calbec' )
!
      return
      end subroutine calbec
!-------------------------------------------------------------------------
      subroutine calphi(c0,ema0bg,bec,betae,phi)
!-----------------------------------------------------------------------
!     input: c0 (orthonormal with s(r(t)), bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
!     where s'=s(r(t))  
!
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use cvan, only: ish, nvb
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use gvecw, only: ngw
      use electrons_base, only: n => nbsp
      use constants, only: pi, fpi
      use control_flags, only: iprint, iprsta
      use mp, only: mp_sum
!
      implicit none
      complex(kind=8) c0(ngw,n), phi(ngw,n), betae(ngw,nhsa)
      real(kind=8)    ema0bg(ngw), bec(nhsa,n), emtot
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j
      real(kind=8) qtemp(nhsavb,n) ! automatic array
!
      call start_clock( 'calphi' )
      phi(:,:) = (0.d0, 0.d0)
!
      if (nvb.gt.0) then
         qtemp (:,:) = 0.d0
         do is=1,nvb
            do iv=1,nh(is)
               do jv=1,nh(is)
                  if(abs(qq(iv,jv,is)) > 1.e-5) then
                     do ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        do i=1,n
                           qtemp(inl,i) = qtemp(inl,i) +                &
     &                                    qq(iv,jv,is)*bec(jnl,i)
                        end do
                     end do
                  endif
               end do
            end do
         end do
!
         call MXMA                                                     &
     &       (betae,1,2*ngw,qtemp,1,nhsavb,phi,1,2*ngw,2*ngw,nhsavb,n)
      end if
!
      do j=1,n
         do i=1,ngw
            phi(i,j)=(phi(i,j)+c0(i,j))*ema0bg(i)
         end do
      end do
!     =================================================================
      if(iprsta > 2) then
         emtot=0.
         do j=1,n
            do i=1,ngw
               emtot=emtot                                              &
     &        +2.*real(phi(i,j)*conjg(c0(i,j)))*ema0bg(i)**(-2.)
            end do
         end do
         emtot=emtot/n

         call mp_sum( emtot )

         WRITE( stdout,*) 'in calphi sqrt(emtot)=',sqrt(emtot)
         WRITE( stdout,*)
         do is=1,nsp
            if(nsp > 1) then
               WRITE( stdout,'(33x,a,i4)') ' calphi: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &            ((bec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
            else
               do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calphi: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &               ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
               end do
            end if
         end do
      endif
      call stop_clock( 'calphi' )
!
      return
      end subroutine calphi
!-----------------------------------------------------------------------
      real(kind=8) function cscnorm(bec,cp,i)
!-----------------------------------------------------------------------
!     requires in input the updated bec(i)
!
      use ions_base, only: na
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use electrons_base, only: n => nbsp
      use cvan, only: ish, nvb
      use uspp_param, only: nh
      use uspp, only: nhsa=>nkb, nhsavb=>nkbus, qq
      use mp, only: mp_sum
!
      implicit none
      integer i
      real(kind=8) bec(nhsa,n)
      complex(kind=8) cp(ngw,n)
!
      integer ig, is, iv, jv, ia, inl, jnl
      real(kind=8) rsum
      real(kind=8), allocatable:: temp(:)
!
!
      allocate(temp(ngw))
      do ig=1,ngw
         temp(ig)=real(conjg(cp(ig,i))*cp(ig,i))
      end do
      rsum=2.*SUM(temp)
      if (gstart == 2) rsum=rsum-temp(1)

      call mp_sum( rsum )

      deallocate(temp)
!
      do is=1,nvb
         do iv=1,nh(is)
            do jv=1,nh(is)
               if(abs(qq(iv,jv,is)).gt.1.e-5) then 
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     rsum = rsum +                                        &
     &                    qq(iv,jv,is)*bec(inl,i)*bec(jnl,i)
                  end do
               endif
            end do
         end do
      end do
!
      cscnorm=sqrt(rsum)
!
      return
      end function cscnorm
!
!-----------------------------------------------------------------------
      subroutine denkin(c,dekin)
!-----------------------------------------------------------------------
!
      use constants, only: pi, fpi
      use electrons_base, only: n => nbsp, nx => nbspx, f
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart, g, gx
      use cell_base, only: ainv, tpiba2
      use gvecw, only: ggp, ecutz, ecsig, ecfix
      use mp, only: mp_sum
!
      implicit none
! input
      complex(kind=8) c(ngw,nx)
! output
      real(kind=8) dekin(3,3)
! local
      integer j, k, ig, i
      real(kind=8), allocatable:: gtmp(:)
      real(kind=8) sk(n)  ! automatic array
      real(kind=8) :: ga, dggp, efac
!
      allocate (gtmp(ngw))
      dekin=0.d0
      do j=1,3
         do k=1,3
            do ig=1,ngw
               efac     = 2.d0 * ecutz / ecsig / sqrt(pi)
               dggp     = 1.d0 + efac * exp( - ( tpiba2 * g(ig) - ecfix ) * ( tpiba2 * g(ig) - ecfix ) / ecsig / ecsig )
               ga       = gx(1,ig) * ainv(k,1) + gx(2,ig) * ainv(k,2) + gx(3,ig) * ainv(k,3)
               gtmp(ig) = gx(j,ig) * ga * dggp
            end do
            do i=1,n
               sk(i)=0.d0
               do ig=gstart,ngw
                  sk(i)=sk(i)+real(conjg(c(ig,i))*c(ig,i))*gtmp(ig)
               end do
            end do
            do i=1,n
               dekin(j,k)=dekin(j,k)-2.d0*tpiba2*(f(i)*sk(i))
            end do
         end do
      end do
      deallocate (gtmp)

      call mp_sum( dekin( 1:3, 1:3 ) )
!
      return
      end subroutine denkin
!
!-----------------------------------------------------------------------
      subroutine denh(rhotmp,drhotmp,sfac,vtemp,eh,dh)
!-----------------------------------------------------------------------
!
! derivative of hartree energy wrt cell parameters h
! Output in dh
!
! rhotmp input : total electronic + ionic broadened charge (G)
! drhotmp input and work space
! sfac   input : structure factors
! wtemp work space
! eh input: hartree energy
!
      use constants, only: pi, fpi
      use ions_base, only: nsp
      use gvecs
      use gvecp, only: ng => ngm
      use reciprocal_vectors, only: gstart, gx, g
      use cell_base, only: omega
      use cell_base, only: ainv, tpiba2
      use local_pseudo, only: rhops, drhops
      use mp, only: mp_sum

      implicit none
! input
      complex(kind=8) rhotmp(ng), drhotmp(ng,3,3), vtemp(ng), sfac(ngs,nsp)
      real(kind=8) eh
! output
      real(kind=8) dh(3,3)
! local
      integer i, j, ig, is
      real(kind=8) wz
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz=2.d0
      do j=1,3
         do i=1,3
            do is=1,nsp
               do ig=1,ngs
                  drhotmp(ig,i,j) = drhotmp(ig,i,j) -                   &
     &                    sfac(ig,is)*drhops(ig,is)*                    &
     &                    2.d0*tpiba2*gx(i,ig)*(gx(1,ig)*ainv(j,1)+     &
     &                     gx(2,ig)*ainv(j,2)+gx(3,ig)*ainv(j,3))-      &
     &                    sfac(ig,is)*rhops(ig,is)*ainv(j,i)
               enddo
            enddo
            if (gstart == 2) vtemp(1)=(0.d0,0.d0)
            do ig=gstart,ng
               vtemp(ig)=conjg(rhotmp(ig))*rhotmp(ig)/(tpiba2*g(ig))**2 &
     &                 * tpiba2*gx(i,ig)*(gx(1,ig)*ainv(j,1)+           &
     &                   gx(2,ig)*ainv(j,2)+gx(3,ig)*ainv(j,3)) +       &
     &                 conjg(rhotmp(ig))/(tpiba2*g(ig))*drhotmp(ig,i,j)
            enddo
            dh(i,j)=fpi*omega*real(SUM(vtemp))*wz
         enddo
      enddo

      call mp_sum( dh( 1:3, 1:3 ) )

      do i=1,3
         do j=1,3
            dh(i,j)=dh(i,j)+omega*eh*ainv(j,i)
         end do
      end do

      return
      end subroutine denh
!
!-----------------------------------------------------------------------
      subroutine dennl(bec,denl)
!-----------------------------------------------------------------------
!
      use cvan, only: ish
      use uspp_param, only: nh
      use uspp, only: nhsa=>nkb, dvan
      use cdvan
      use electrons_base, only: n => nbsp, ispin => fspin, f, nspin
      use reciprocal_vectors, only: gstart
      use ions_base, only: nsp, na
      implicit none
! input
      real(kind=8) bec(nhsa,n)
! output
      real(kind=8) denl(3,3)
! local
      real(kind=8) dsum(3,3),dsums(2,3,3)
      integer is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss, i,j,k
!
      denl=0.d0
      do is=1,nsp
         ijv=0
         do iv=1,nh(is)
            do jv=iv,nh(is)
               ijv=ijv+1
               isa=0
               do ism=1,is-1
                  isa=isa+na(ism)
               end do
               do ia=1,na(is)
                  inl=ish(is)+(iv-1)*na(is)+ia
                  jnl=ish(is)+(jv-1)*na(is)+ia
                  isa=isa+1
                  dsums=0.d0
                  do i=1,n
                     iss=ispin(i) 
                     do k=1,3
                        do j=1,3
                           dsums(iss,k,j)=dsums(iss,k,j)+f(i)*       &
     &                          (dbec(inl,i,k,j)*bec(jnl,i)          &
     &                          + bec(inl,i)*dbec(jnl,i,k,j))
                        enddo
                     enddo
                  end do
                  dsum=0.d0
                  do iss=1,nspin
                     do k=1,3
                        do j=1,3
                           drhovan(ijv,isa,iss,j,k)=dsums(iss,j,k)
                           dsum(j,k)=dsum(j,k)+dsums(iss,j,k)
                        enddo
                     enddo
                  end do
                  if(iv.ne.jv) dsum=2.d0*dsum
                  denl = denl + dsum*dvan(jv,iv,is)
               end do
            end do
         end do
      end do
!
      return
      end subroutine dennl
!
!-----------------------------------------------------------------------
      subroutine denps(rhotmp,drhotmp,sfac,vtemp,dps)
!-----------------------------------------------------------------------
!
! derivative of local potential energy wrt cell parameters h
! Output in dps
!
! rhotmp input : rho(G) (up and down spin components summed)
! drhotmp input
! sfac   input : structure factors
! wtemp work space
!
      use ions_base, only: nsp
      use gvecs, only: ngs
      use gvecp, only: ng => ngm
      use reciprocal_vectors, only: gstart, gx
      use cell_base, only: omega
      use cell_base, only: ainv, tpiba2
      use local_pseudo, only: vps, dvps
      use mp, only: mp_sum

      implicit none
! input
      complex(kind=8) rhotmp(ng), drhotmp(ng,3,3), vtemp(ng), sfac(ngs,nsp)
! output
      real(kind=8) dps(3,3)
! local
      integer i, j, ig, is
      real(kind=8) wz
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz=2.d0
      do i=1,3
         do j=1,3
            do ig=1,ngs
               vtemp(ig)=(0.,0.)
            enddo
            do is=1,nsp
               do ig=1,ngs
                  vtemp(ig)=vtemp(ig)-conjg(rhotmp(ig))*sfac(ig,is)*    &
     &                    dvps(ig,is)*2.d0*tpiba2*gx(i,ig)*             &
     &                    (gx(1,ig)*ainv(j,1) +                         &
     &                     gx(2,ig)*ainv(j,2) +                         &
     &                     gx(3,ig)*ainv(j,3) ) +                       &
     &                    conjg(drhotmp(ig,i,j))*sfac(ig,is)*vps(ig,is)
               enddo
            enddo
            dps(i,j)=omega*real(wz*SUM(vtemp))
            if (gstart == 2) dps(i,j)=dps(i,j)-omega*real(vtemp(1))
         enddo
      enddo

      call mp_sum( dps( 1:3, 1:3 ) )

      return
      end subroutine denps
!
!-------------------------------------------------------------------------
      subroutine dforce (bec,betae,i,c,ca,df,da,v)
!-----------------------------------------------------------------------
!computes: the generalized force df=cmplx(dfr,dfi) acting on the i-th
!          electron state at the gamma point of the brillouin zone
!          represented by the vector c=cmplx(cr,ci)
!
!     d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
!              sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
!                                 e^-ig.r_i < beta_i,j | c_n >}
      use control_flags, only: iprint, tbuff
      use gvecs
      use gvecw, only: ngw
      use cvan, only: ish
      use uspp, only: nhsa=>nkb, dvan, deeq
      use uspp_param, only: nhm, nh
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base, only: n => nbsp, ispin => fspin, f, nspin
      use constants, only: pi, fpi
      use ions_base, only: nsp, na, nat
      use gvecw, only: ggp
      use cell_base, only: tpiba2
      use ensemble_dft, only: tens
!
      implicit none
!
      complex(kind=8) betae(ngw,nhsa), c(ngw), ca(ngw), df(ngw), da(ngw)
      real(kind=8) bec(nhsa,n), v(nnrsx,nspin)
      integer i
! local variables
      integer iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
      real(kind=8) fi, fip, dd
      complex(kind=8) fp,fm,ci
      real(kind=8) af(nhsa), aa(nhsa) ! automatic arrays
      complex(kind=8)  dtemp(ngw)    !
      complex(kind=8), allocatable :: psi(:)
!
!
      call start_clock( 'dforce' ) 
      !
      allocate( psi( nnrsx ) )
!
!     important: if n is odd => c(*,n+1)=0.
! 
      if (mod(n,2).ne.0.and.i.eq.n) then
         do ig=1,ngw
            ca(ig)=(0.,0.)
         end do
      endif
!
      ci=(0.0,1.0)
!
      if (.not.tbuff) then
!
         psi (:) = (0.d0, 0.d0)
         do ig=1,ngw
            psi(nms(ig))=conjg(c(ig)-ci*ca(ig))
            psi(nps(ig))=c(ig)+ci*ca(ig)
         end do
!
         call ivfftw(psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!     
      else
!
!     read psi from buffer 21
!
#if defined(__CRAYY)
         buffer in(21,0) (psi(1),psi(nnrsx))
         ios = unit(21)
#else
         read(21,iostat=ios) psi
#endif
         if(ios.ne.0) call errore                                        &
     &       (' dforce',' error in reading unit 21',ios)
!
      endif
! 
      iss1=ispin(i)
!
! the following avoids a potential out-of-bounds error
!
      if (i.ne.n) then
         iss2=ispin(i+1)
      else
         iss2=iss1
      end if
!
      do ir=1,nnrsx
         psi(ir)=cmplx(v(ir,iss1)* real(psi(ir)),                       &
     &                 v(ir,iss2)*aimag(psi(ir)) )
      end do
!
      call fwfftw(psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!
!     note : the factor 0.5 appears 
!       in the kinetic energy because it is defined as 0.5*g**2
!       in the potential part because of the logics
!
   
     if (tens) then
        fi =-0.5
        fip=-0.5
      else
        fi =-  f(i)*0.5
        fip=-f(i+1)*0.5
      end if

      do ig=1,ngw
         fp= psi(nps(ig)) + psi(nms(ig))
         fm= psi(nps(ig)) - psi(nms(ig))
         df(ig)= fi*(tpiba2*ggp(ig)* c(ig)+cmplx(real(fp), aimag(fm)))
         da(ig)=fip*(tpiba2*ggp(ig)*ca(ig)+cmplx(aimag(fp),-real(fm)))
      end do
!
!     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
! 
      if(nhsa.gt.0)then
         do inl=1,nhsa
            af(inl)=0.
            aa(inl)=0.
         end do
!
         do is=1,nsp
            do iv=1,nh(is)
               do jv=1,nh(is)
                  isa=0
                  do ism=1,is-1
                     isa=isa+na(ism)
                  end do
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     isa=isa+1
                     dd = deeq(iv,jv,isa,iss1)+dvan(iv,jv,is)
                     if(tens) then
                      af(inl)=af(inl)-dd*bec(jnl,  i)
                     else
                      af(inl)=af(inl)- f(i)*dd*bec(jnl,  i)
                     end if
                     dd = deeq(iv,jv,isa,iss2)+dvan(iv,jv,is)
                     if(tens) then
                      if (i.ne.n) aa(inl)=aa(inl)-dd*bec(jnl,i+1)
                     else
                      if (i.ne.n) aa(inl)=aa(inl)-f(i+1)*dd*bec(jnl,i+1)
                     end if
                  end do
               end do
            end do
         end do
!
         do ig=1,ngw
            dtemp(ig)=(0.,0.)
         end do
         call MXMA                                                      &
     &        (betae,1,2*ngw,af,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         do ig=1,ngw
            df(ig)=df(ig)+dtemp(ig)
         end do
!
         do ig=1,ngw
            dtemp(ig)=(0.,0.)
         end do
         call MXMA                                                      &
     &        (betae,1,2*ngw,aa,1,nhsa,dtemp,1,2*ngw,2*ngw,nhsa,1)
         do ig=1,ngw
            da(ig)=da(ig)+dtemp(ig)
         end do
      endif

      deallocate( psi )
!
      call stop_clock( 'dforce' ) 
!
      return
      end subroutine dforce
!
!-----------------------------------------------------------------------
      subroutine dotcsc(eigr,cp)
!-----------------------------------------------------------------------
!
      use ions_base, only: nas => nax, na, nsp, nat
      use io_global, only: stdout
      use gvecw, only: ngw
      use electrons_base, only: n => nbsp
      use reciprocal_vectors, only: gstart
      use cvan, only: ish, nvb
      use uspp, only: nhsa=>nkb, qq
      use uspp_param, only: nh
      use mp, only: mp_sum
!
      implicit none
!
      complex(kind=8)  eigr(ngw,nat), cp(ngw,n)
! local variables
      real(kind=8) rsum, csc(n) ! automatic array
      complex(kind=8) temp(ngw) ! automatic array
 
      real(kind=8), allocatable::  becp(:,:)
      integer i,kmax,nnn,k,ig,is,ia,iv,jv,inl,jnl
!
      allocate(becp(nhsa,n))
!
!     < beta | phi > is real. only the i lowest:
!
      nnn=min(12,n)
      do i=nnn,1,-1
         kmax=i
         call nlsm1(i,1,nvb,eigr,cp,becp)
!
         do k=1,kmax
            do ig=1,ngw
               temp(ig)=conjg(cp(ig,k))*cp(ig,i)
            end do
            csc(k)=2.*real(SUM(temp))
            if (gstart == 2) csc(k)=csc(k)-real(temp(1))
         end do

         call mp_sum( csc( 1:kmax ) )

         do k=1,kmax
            rsum=0.
            do is=1,nvb
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     do ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        rsum = rsum +                                    &
     &                   qq(iv,jv,is)*becp(inl,i)*becp(jnl,k)
                     end do
                  end do
               end do
            end do
            csc(k)=csc(k)+rsum
         end do
!
         WRITE( stdout,'(a,12f18.15)')' dotcsc = ',(csc(k),k=1,i)
!
      end do
      WRITE( stdout,*)
!
      deallocate(becp)
!
      return
      end subroutine dotcsc
!-----------------------------------------------------------------------
      subroutine drhov(irb,eigrb,rhovan,rhog,rhor)
!-----------------------------------------------------------------------
!     this routine calculates arrays drhog drhor, derivatives wrt h of:
!
!        n_v(g) = sum_i,ij rho_i,ij q_i,ji(g) e^-ig.r_i
!
!     Same logic as in routine rhov.
!     On input rhor and rhog must contain the smooth part only !!!
!     Output in module derho (drhor, drhog)
!
      use control_flags, only: iprint
      use parameters, only: natx, nsx
      use ions_base, only: na, nsp, nat, nas => nax
      use cvan
      use uspp_param, only: nhm, nh
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use electrons_base, only: nspin
      use gvecb
      use gvecp, only: ng => ngm
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      use cell_base, only: ainv
      use qgb_mod
      use para_mod
      use cdvan
      use derho
      use dqgb_mod
      use recvecs_indexes, only: nm, np

      implicit none
! input
      integer, intent(in) ::  irb(3,nat)
      real(kind=8), intent(in)::  rhor(nnr,nspin)
      real(kind=8) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      complex(kind=8), intent(in)::  eigrb(ngb,nat), rhog(ng,nspin)
! local
      integer i, j, isup, isdw, nfft, ifft, iv, jv, ig, ijv, is, iss,   &
     &     isa, ia, ir, irb3, imin3, imax3
      real(kind=8) sum, dsum
      complex(kind=8) fp, fm, ci
      complex(kind=8), allocatable :: v(:)
      complex(kind=8), allocatable:: dqgbt(:,:)
      complex(kind=8), allocatable :: qv(:)
!
!
      do j=1,3
         do i=1,3
            do iss=1,nspin
               do ir=1,nnr
                  drhor(ir,iss,i,j)=-rhor(ir,iss)*ainv(j,i)
               end do
               do ig=1,ng
                  drhog(ig,iss,i,j)=-rhog(ig,iss)*ainv(j,i)
               end do
            end do
         end do
      end do
!
      if (nvb.eq.0) return
!
      allocate( v( nnr ) )
      allocate( qv( nnrb ) )
      allocate( dqgbt( ngb, 2 ) )

      ci=(0.,1.)
!
      if(nspin.eq.1) then
!     ------------------------------------------------------------------
!     nspin=1 : two fft at a time, one per atom, if possible
!     ------------------------------------------------------------------
         do i=1,3
            do j=1,3
!
               v(:) = (0.d0, 0.d0)
!
               iss=1
               isa=1
               do is=1,nvb
#ifdef __PARA
                  do ia=1,na(is)
                     nfft=1
                     irb3=irb(3,isa)
                     call parabox(nr3b,irb3,nr3,imin3,imax3)
                     if (imax3-imin3+1.le.0) go to 15
#else
                  do ia=1,na(is),2
                     nfft=2
#endif
                     dqgbt(:,:) = (0.d0, 0.d0) 
                     if (ia.eq.na(is)) nfft=1
!
!  nfft=2 if two ffts at the same time are performed
!
                     do ifft=1,nfft
                        ijv=0
                        do iv=1,nh(is)
                           do jv=iv,nh(is)
                              ijv=ijv+1
                              sum = rhovan(ijv,isa+ifft-1,iss)
                              dsum=drhovan(ijv,isa+ifft-1,iss,i,j)
                              if(iv.ne.jv) then
                                 sum =2.*sum
                                 dsum=2.*dsum
                              endif
                              do ig=1,ngb
                                 dqgbt(ig,ifft)=dqgbt(ig,ifft) +        &
     &                                (sum*dqgb(ig,ijv,is,i,j) +        &
     &                                dsum*qgb(ig,ijv,is) )
                              end do
                           end do
                        end do
                     end do
!     
! add structure factor
!
                     qv(:) = (0.d0, 0.d0)
                     if(nfft.eq.2) then
                        do ig=1,ngb
                           qv(npb(ig)) = eigrb(ig,isa   )*dqgbt(ig,1)  &
     &                        + ci*      eigrb(ig,isa+1 )*dqgbt(ig,2)
                           qv(nmb(ig))=                                 &
     &                             conjg(eigrb(ig,isa  )*dqgbt(ig,1)) &
     &                        + ci*conjg(eigrb(ig,isa+1)*dqgbt(ig,2))
                        end do
                     else
                        do ig=1,ngb
                           qv(npb(ig)) = eigrb(ig,isa)*dqgbt(ig,1)
                           qv(nmb(ig)) =                                &
     &                             conjg(eigrb(ig,isa)*dqgbt(ig,1))
                        end do
                     endif
!
                     call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
!  qv = US contribution in real space on box grid
!       for atomic species is, real(qv)=atom ia, imag(qv)=atom ia+1
!
!  add qv(r) to v(r), in real space on the dense grid
!
                     call box2grid(irb(1,isa),1,qv,v)
                     if (nfft.eq.2) call box2grid(irb(1,isa+1),2,qv,v)
  15                 isa=isa+nfft
!
                  end do
               end do
!
               do ir=1,nnr
                  drhor(ir,iss,i,j)=drhor(ir,iss,i,j)+real(v(ir))
               end do
!
               call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
               do ig=1,ng
                  drhog(ig,iss,i,j)=drhog(ig,iss,i,j)+v(np(ig))
               end do
!
            enddo
         enddo
!
      else
!     ------------------------------------------------------------------
!     nspin=2: two fft at a time, one for spin up and one for spin down
!     ------------------------------------------------------------------
         isup=1
         isdw=2
         do i=1,3
            do j=1,3
               v(:) = (0.d0, 0.d0)
               isa=1
               do is=1,nvb
                  do ia=1,na(is)
#ifdef __PARA
                     irb3=irb(3,isa)
                     call parabox(nr3b,irb3,nr3,imin3,imax3)
                     if (imax3-imin3+1.le.0) go to 25
#endif
                     do iss=1,2
                        dqgbt(:,iss) = (0.d0, 0.d0)
                        ijv=0
                        do iv= 1,nh(is)
                           do jv=iv,nh(is)
                              ijv=ijv+1
                              sum=rhovan(ijv,isa,iss)
                              dsum =drhovan(ijv,isa,iss,i,j)
                              if(iv.ne.jv) then
                                 sum =2.*sum
                                 dsum=2.*dsum
                              endif
                              do ig=1,ngb
                                 dqgbt(ig,iss)=dqgbt(ig,iss)  +         &
     &                               (sum*dqgb(ig,ijv,is,i,j) +         &
     &                               dsum*qgb(ig,ijv,is))
                              end do
                           end do
                        end do
                     end do
!     
! add structure factor
!
                     qv(:) = (0.d0, 0.d0)
                     do ig=1,ngb
                        qv(npb(ig))= eigrb(ig,isa)*dqgbt(ig,1)        &
     &                    + ci*      eigrb(ig,isa)*dqgbt(ig,2)
                        qv(nmb(ig))= conjg(eigrb(ig,isa)*dqgbt(ig,1)) &
     &                    +       ci*conjg(eigrb(ig,isa)*dqgbt(ig,2))
                     end do
!
                     call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
!  qv is the now the US augmentation charge for atomic species is
!  and atom ia: real(qv)=spin up, imag(qv)=spin down
!
!  add qv(r) to v(r), in real space on the dense grid
!
                     call box2grid2(irb(1,isa),qv,v)
  25                 isa=isa+1
                  end do
               end do
!
               do ir=1,nnr
                  drhor(ir,isup,i,j) = drhor(ir,isup,i,j) + real(v(ir))
                  drhor(ir,isdw,i,j) = drhor(ir,isdw,i,j) +aimag(v(ir))
               enddo
!
               call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
               do ig=1,ng
                  fp=v(np(ig))+v(nm(ig))
                  fm=v(np(ig))-v(nm(ig))
                  drhog(ig,isup,i,j) = drhog(ig,isup,i,j) +             &
     &                 0.5*cmplx( real(fp),aimag(fm))
                  drhog(ig,isdw,i,j) = drhog(ig,isdw,i,j) +             &
     &                 0.5*cmplx(aimag(fp),-real(fm))
               end do
!
            end do
         end do
      endif
      deallocate(dqgbt)
      deallocate( v )
      deallocate( qv )
!
      return
      end subroutine drhov

!
!-----------------------------------------------------------------------
      real(kind=8) function enkin(c)
!-----------------------------------------------------------------------
!
! calculation of kinetic energy term
!
      use constants, only: pi, fpi
      use electrons_base, only: nx => nbspx, n => nbsp, f
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use gvecw, only: ggp
      use mp, only: mp_sum
      use cell_base, only: tpiba2

      implicit none
! input
      complex(kind=8) c(ngw,nx)
! local
      integer ig, i
      real(kind=8) sk(n)  ! automatic array
!
!
      do i=1,n
         sk(i)=0.0
         do ig=gstart,ngw
            sk(i)=sk(i)+real(conjg(c(ig,i))*c(ig,i))*ggp(ig)
         end do
      end do

      call mp_sum( sk(1:n) )

      enkin=0.0
      do i=1,n
         enkin=enkin+f(i)*sk(i)
      end do
      enkin=enkin*tpiba2
!
      return
      end function enkin
!
!-----------------------------------------------------------------------
      real(kind=8) function ennl(rhovan, bec)
!-----------------------------------------------------------------------
!
! calculation of nonlocal potential energy term
!
      use cvan, only: ish
      use uspp_param, only: nhm, nh
      use uspp, only :nhsa=>nkb, dvan
      use electrons_base, only: n => nbsp, nspin, ispin => fspin, f
      use ions_base, only: nsp, nat, na
      implicit none
! input
      real(kind=8) bec(nhsa,n)
      real(kind=8) rhovan(nhm*(nhm+1)/2,nat,nspin)
! local
      real(kind=8) sum, sums(2)
      integer is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss, i
!
!
      ennl=0.d0
      do is=1,nsp
         ijv=0
         do iv= 1,nh(is)
            do jv=iv,nh(is)
               ijv=ijv+1
               isa=0
               do ism=1,is-1
                  isa=isa+na(ism)
               end do
               do ia=1,na(is)
                  inl=ish(is)+(iv-1)*na(is)+ia
                  jnl=ish(is)+(jv-1)*na(is)+ia
                  isa=isa+1
                  sums=0.d0
                  do i=1,n
                     iss=ispin(i) 
                     sums(iss) = sums(iss) +f(i)*bec(inl,i)*bec(jnl,i)
                  end do
                  sum=0.d0
                  do iss=1,nspin
                     rhovan(ijv,isa,iss) = sums(iss)
                     sum=sum+sums(iss)
                  end do
                  if(iv.ne.jv) sum=2.d0*sum
                  ennl=ennl+sum*dvan(jv,iv,is)
               end do
            end do
         end do
      end do
!
      return
      end function ennl

!
!
!-----------------------------------------------------------------------
      subroutine force_ion(tau0,esr,fion,dsr)
!-----------------------------------------------------------------------
!
!     forces on ions, ionic term in real space (also stress if requested)
!
      use parameters, only: nsx, natx
      use control_flags, only: iprint, tpre
      use constants, only: pi, fpi
      use cell_base, only: ainv, a1, a2, a3
      use ions_base, only: nsp, na, rcmax, zv
      implicit none
! input
      real(kind=8) tau0(3,natx)
! output
      real(kind=8) fion(3,natx), dsr(3,3), esr
! local variables
      integer i,j,k,l,m, ii, lax, inf, isak, isaj
      real(kind=8) rlm(3), rckj, rlmn, arg, addesr, addpre, repand, fxx
      real(kind=8), external :: erfc
!
!
      esr=0.d0
      if(tpre) dsr=0.d0
!
      isak = 0
      do k=1,nsp
         isaj = 0
         do j = 1, k-1
           isaj = isaj + na(j)
         end do
         do j=k,nsp
            rckj=sqrt(rcmax(k)**2+rcmax(j)**2)
            lax=na(k)
            if(k.eq.j) lax=lax-1
!
            do l=1,lax
               inf=1
               if(k.eq.j) inf=l+1
!
               do m=inf,na(j)
                  rlm(1) = tau0(1,l + isak) - tau0(1,m + isaj)
                  rlm(2) = tau0(2,l + isak) - tau0(2,m + isaj)
                  rlm(3) = tau0(3,l + isak) - tau0(3,m + isaj)
                  call pbc(rlm,a1,a2,a3,ainv,rlm)
!
                  rlmn=sqrt(rlm(1)**2+rlm(2)**2+rlm(3)**2)
!
                  arg=rlmn/rckj
                  addesr=zv(k)*zv(j)*erfc(arg)/rlmn
                  esr=esr+addesr
                  addpre=2.d0*zv(k)*zv(j)*exp(-arg*arg)/rckj/sqrt(pi)
                  repand=(addesr+addpre)/rlmn/rlmn
!
                  do i=1,3
                     fxx=repand*rlm(i)
                     fion(i,l+isak)=fion(i,l+isak)+fxx
                     fion(i,m+isaj)=fion(i,m+isaj)-fxx
                     if(tpre)then
                        do ii=1,3
                           dsr(i,ii)=dsr(i,ii)-                         &
     &                             repand*rlm(i)*rlm(1)*ainv(ii,1)-     &
     &                             repand*rlm(i)*rlm(2)*ainv(ii,2)-     &
     &                             repand*rlm(i)*rlm(3)*ainv(ii,3)
                        end do
                     endif
                  end do
               end do
            end do
            isaj = isaj + na(j)
         end do
         isak = isak + na(k)
      end do

      return
      end subroutine force_ion
!
!-----------------------------------------------------------------------
      subroutine force_ps(rhotemp,rhog,vtemp,ei1,ei2,ei3,fion1)
!-----------------------------------------------------------------------
!
! Contribution to ionic forces from local pseudopotential
!
      use constants, only: pi, fpi
      use electrons_base, only: nspin
      use gvecs
      use gvecp, only: ng => ngm
      use reciprocal_vectors, only: gstart, gx, mill_l, g
      use cell_base, only: omega, tpiba, tpiba2
      use ions_base, only: nsp, na, nas => nax, nat
      use grid_dimensions, only: nr1, nr2, nr3
      use parameters, only: nsx, natx
      use local_pseudo, only: vps, rhops
!
      implicit none
! input
      complex(kind=8) rhotemp(ng), rhog(ng,nspin), vtemp(ng),           &
     &           ei1(-nr1:nr1,nat),                                 &
     &           ei2(-nr2:nr2,nat),                                 &
     &           ei3(-nr3:nr3,nat)
! output
      real(kind=8) fion1(3,natx)
! local
      integer ig, is, isa, ism, ia, ix, iss, isup, isdw
      integer i, j, k
      real(kind=8)  wz
      complex(kind=8) eigrx, vcgs, cnvg, cvn
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz=2.0
      do is=1,nsp
         isa=0
         do ism=1,is-1
            isa=isa+na(ism)
         end do
         do ia=1,na(is)
            isa=isa+1
            do ix=1,3
               if(nspin.eq.1)then
                  iss=1
                  if (gstart == 2) vtemp(1)=0.0
                  do ig=gstart,ngs
                     vcgs=conjg(rhotemp(ig))*fpi/(tpiba2*g(ig))
                     cnvg=rhops(ig,is)*vcgs
                     cvn=vps(ig,is)*conjg(rhog(ig,iss))
                     i = mill_l(1,ig)
                     j = mill_l(2,ig)
                     k = mill_l(3,ig)
                     eigrx=ei1(i,isa)*ei2(j,isa)*ei3(k,isa)
                     vtemp(ig)=eigrx*(cnvg+cvn)*cmplx(0.0,gx(ix,ig)) 
                  end do
               else
                  isup=1
                  isdw=2
                  if (gstart == 2) vtemp(1)=0.0
                  do ig=gstart,ngs
                     vcgs=conjg(rhotemp(ig))*fpi/(tpiba2*g(ig))
                     cnvg=rhops(ig,is)*vcgs
                     cvn=vps(ig,is)*conjg(rhog(ig,isup)                 &
     &                                   +rhog(ig,isdw))
                     i = mill_l(1,ig)
                     j = mill_l(2,ig)
                     k = mill_l(3,ig)
                     eigrx=ei1(i,isa)*ei2(j,isa)*ei3(k,isa)
                     vtemp(ig)=eigrx*(cnvg+cvn)*cmplx(0.0,gx(ix,ig)) 
                  end do
               endif
               fion1(ix,isa) = fion1(ix,isa) + tpiba*omega* wz*real(SUM(vtemp))
            end do
         end do
      end do
!
      return
      end subroutine force_ps
!
!-----------------------------------------------------------------------
      subroutine gausin(eigr,cm)
!-----------------------------------------------------------------------
!
! initialize wavefunctions with gaussians - edit to fit your system
!
      use ions_base, only: nas => nax, na, nsp, nat
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use reciprocal_vectors, only: gx, g
!
      implicit none
!
      complex(kind=8) eigr(ngw,nat), cm(ngw,n)
      real(kind=8)    sigma, auxf
      integer nband, is, ia, ig, isa
!
      sigma=12.0
      nband=0
!!!      do is=1,nsp
      isa = 0
      is=1
         do ia=1,na(is)
! s-like gaussians
            nband=nband+1
            do ig=1,ngw
               auxf=exp(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)
            end do
! px-like gaussians
            nband=nband+1
            do ig=1,ngw
               auxf=exp(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)
            end do
! py-like gaussians
            nband=nband+1
            do ig=1,ngw
               auxf=exp(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)
            end do
! pz-like gaussians
            nband=nband+1
            do ig=1,ngw
               auxf=exp(-g(ig)/sigma**2)
               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(3,ig)
            end do
         end do
      isa = isa + na(is)
      is=2
         do ia=1,na(is)
! s-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)
!            end do
! px-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)
!            end do
! py-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)
!            end do
! pz-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(3,ig)
!            end do
! dxy-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)*gx(2,ig)
!            end do
! dxz-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)*gx(3,ig)
!            end do
! dxy-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)*gx(3,ig)
!            end do
! dx2-y2-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*                        &
!     &              (gx(1,ig)**2-gx(2,ig)**2)
!            end do
         end do
!!!      end do
      return
      end subroutine gausin
!            

!-------------------------------------------------------------------------
      subroutine gracsc(bec,betae,cp,i,csc)
!-----------------------------------------------------------------------
!     requires in input the updated bec(k) for k<i
!     on output: bec(i) is recalculated
!
      use ions_base, only: na
      use cvan, only :nvb, ish
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use uspp_param, only:  nh
      use electrons_base, only: n => nbsp, ispin => fspin, nx => nbspx
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use mp, only: mp_sum
!
      implicit none
!
      integer i
      complex(kind=8) betae(ngw,nhsa)
      real(kind=8)  bec(nhsa,n), cp(2,ngw,n)
      real(kind=8)  csc(nx)
      integer k, kmax,ig, is, iv, jv, ia, inl, jnl
      real(kind=8) rsum, temp(ngw) ! automatic array
!
!     calculate csc(k)=<cp(i)|cp(k)>,  k<i
!
      kmax=i-1
      do k=1,kmax
         csc(k)=0.
         if (ispin(i).eq.ispin(k)) then
            do ig=1,ngw
               temp(ig)=cp(1,ig,k)*cp(1,ig,i)+cp(2,ig,k)*cp(2,ig,i)
            end do
            csc(k)=2.*SUM(temp)
            if (gstart == 2) csc(k)=csc(k)-temp(1)
         endif
      end do

      call mp_sum( csc( 1:kmax ) )

!
!     calculate bec(i)=<cp(i)|beta>
!
      do inl=1,nhsavb
         do ig=1,ngw
            temp(ig)=cp(1,ig,i)* real(betae(ig,inl))+             &
     &               cp(2,ig,i)*aimag(betae(ig,inl))
         end do
         bec(inl,i)=2.*SUM(temp)
         if (gstart == 2) bec(inl,i)= bec(inl,i)-temp(1)
      end do

      call mp_sum( bec( 1:nhsavb, i ) )
!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i
!
      do k=1,kmax
         if (ispin(i).eq.ispin(k)) then
            rsum=0.
            do is=1,nvb
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     if(abs(qq(iv,jv,is)).gt.1.e-5) then 
                        do ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           rsum = rsum + qq(iv,jv,is)*bec(inl,i)*bec(jnl,k)
                        end do
                     endif
                  end do
               end do
            end do
            csc(k)=csc(k)+rsum
         endif
      end do
!
!     orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
!
!     corresponing bec:  bec(i)=<cp(i)|beta>-csc(k)<cp(k)|beta>
!
      do k=1,kmax
         do inl=1,nhsavb
            bec(inl,i)=bec(inl,i)-csc(k)*bec(inl,k)
         end do
      end do
!
      return
      end subroutine gracsc
!-------------------------------------------------------------------------
      subroutine gram(betae,bec,cp)
!-----------------------------------------------------------------------
!     gram-schmidt orthogonalization of the set of wavefunctions cp
!
      use uspp, only :nhsa=>nkb, nhsavb=> nkbus
      use electrons_base, only: nx => nbspx, n => nbsp
      use gvecw, only: ngw
!
      implicit none
!
      real(kind=8)  bec(nhsa,n)
      complex(kind=8)   cp(ngw,n), betae(ngw,nhsa)
!
      real(kind=8) :: anorm, cscnorm
      real(kind=8), allocatable :: csc( : )
      integer :: i,k
      external cscnorm
!
      call start_clock( 'gram' )

      allocate( csc( nx ) )
!
      do i=1,n
         call gracsc(bec,betae,cp,i,csc)
!
! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
!
         do k=1,i-1
            call DAXPY(2*ngw,-csc(k),cp(1,k),1,cp(1,i),1)
         end do
         anorm =cscnorm(bec,cp,i)
         call DSCAL(2*ngw,1.0/anorm,cp(1,i),1)
!
!         these are the final bec's
!
         call DSCAL(nhsavb,1.0/anorm,bec(1,i),1)
      end do
!
      deallocate( csc )

      call stop_clock( 'gram' )
!
      return
      end subroutine gram
!
!-----------------------------------------------------------------------
      subroutine herman_skillman_grid(mesh,z,cmesh,r)
!-----------------------------------------------------------------------
!
      implicit none
!
      integer mesh
      real(kind=8) z, cmesh, r(mesh)
!
      real(kind=8) deltax
      integer nblock,i,j,k
!
      nblock = mesh/40
      i=1
      r(i)=0.0
      cmesh=0.88534138/z**(1.0/3.0)
      deltax=0.0025*cmesh
      do j=1,nblock
         do k=1,40
            i=i+1
            r(i)=r(i-1)+deltax
         end do
         deltax=deltax+deltax
      end do
!
      return
      end subroutine herman_skillman_grid
!
!-----------------------------------------------------------------------
      subroutine herman_skillman_int(mesh,cmesh,func,asum)
!-----------------------------------------------------------------------
!     simpsons rule integration for herman skillman mesh
!     mesh - # of mesh points
!     c    - 0.8853418/z**(1/3.)
!
      implicit none
      integer mesh
      real(kind=8) cmesh, func(mesh), asum
!
      integer i, j, k, i1, nblock
      real(kind=8) a1, a2e, a2o, a2es, h
!
      a1=0.0
      a2e=0.0
      asum=0.0
      h=0.0025*cmesh
      nblock=mesh/40
      i=1
      func(1)=0.0
      do j=1,nblock
         do k=1,20
            i=i+2
            i1=i-1
            a2es=a2e
            a2o=func(i1)/12.0
            a2e=func(i)/12.0
            a1=a1+5.0*a2es+8.0*a2o-a2e
            func(i1)=asum+a1*h
            a1=a1-a2es+8.0*a2o+5.0*a2e
            func(i)=asum+a1*h
         end do
         asum=func(i)
         a1=0.0
         h=h+h
      end do
!
      return
      end subroutine herman_skillman_int
!
!-----------------------------------------------------------------------
      subroutine initbox ( tau0, taub, irb )
!-----------------------------------------------------------------------
!
!     sets the indexes irb and positions taub for the small boxes 
!     around atoms
!
      use parameters, only: natx, nsx
      use ions_base, only: nsp, na, nat
      use grid_dimensions, only: nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
      use control_flags, only: iprsta
      use cvan, only: nvb
      use io_global, only: stdout

      implicit none
! input
      real(kind=8), intent(in):: tau0(3,natx)
! output
      integer, intent(out):: irb(3,nat)
      real(kind=8), intent(out):: taub(3,natx)
! local
      real(kind=8) x(3), xmod
      integer nr(3), nrb(3), xint, is, ia, i, isa
!
      nr (1)=nr1
      nr (2)=nr2
      nr (3)=nr3
      nrb(1)=nr1b
      nrb(2)=nr2b
      nrb(3)=nr3b
!
      isa = 0
      do is=1,nsp
         do ia=1,na(is)
           isa = isa + 1
!
            do i=1,3
!
! bring atomic positions to crystal axis
!
               x(i) = ainv(i,1)*tau0(1,isa) +                         &
     &                ainv(i,2)*tau0(2,isa) +                         &
     &                ainv(i,3)*tau0(3,isa)
!
! bring x in the range between 0 and 1
!
               x(i) = mod(x(i),1.d0)
               if (x(i).lt.0.d0) x(i)=x(i)+1.d0
!
! case of nrb(i) even
!
               if (mod(nrb(i),2).eq.0) then
!
! find irb = index of the grid point at the corner of the small box
!           (the indices of the small box run from irb to irb+nrb-1)
!
                  xint=int(x(i)*nr(i))
                  irb (i,isa)=xint+1-nrb(i)/2+1
                  if(irb(i,isa).lt.1) irb(i,isa)=irb(i,isa)+nr(i)
!
! x(i) are the atomic positions in crystal coordinates, where the
! "crystal lattice" is the small box lattice and the origin is at
! the corner of the small box. Used to calculate phases exp(iG*taub)
!
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+nrb(i)/2-1)/nr(i)
               else
!
! case of nrb(i) odd - see above for comments
!
                  xint=nint(x(i)*nr(i))
                  irb (i,isa)=xint+1-(nrb(i)-1)/2
                  if(irb(i,isa).lt.1) irb(i,isa)=irb(i,isa)+nr(i)
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+(nrb(i)-1)/2)/nr(i)
               end if
            end do
!
! bring back taub in cartesian coordinates
!
            do i=1,3
               taub(i,isa)= x(1)*a1(i) + x(2)*a2(i) + x(3)*a3(i)
            end do
         end do
      end do

      if( iprsta > 2 ) then
           isa = 1
           do is=1,nvb
              WRITE( stdout,'(/,2x,''species= '',i2)') is
              do ia=1,na(is)
                 WRITE( stdout,2000) ia, (irb(i,isa),i=1,3)
2000             format(2x,'atom= ',i3,' irb1= ',i3,' irb2= ',i3,' irb3= ',i3)
                 isa = isa + 1
               end do
            end do
      endif

!
      return
      end subroutine initbox
!
!-------------------------------------------------------------------------
      subroutine newd(vr,irb,eigrb,rhovan,fion)
!-----------------------------------------------------------------------
!     this routine calculates array deeq:
!         deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!     and the corresponding term in forces
!         fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!     where
!         rho_lm = \sum_j f_j <psi_j|beta_l><beta_m|psi_j>
!
      use uspp_param, only: nh, nhm
      use uspp, only: deeq
      use cvan, only: nvb
      use ions_base, only: nas => nax, nat, nsp, na
      use parameters, only: natx, nsx
      use constants, only: pi, fpi
      use grid_dimensions, only: nr3, nnr => nnrx
      use gvecb
      use small_box, only: omegab, tpibab
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      use qgb_mod
      use electrons_base, only: nspin
      use control_flags, only: iprint, thdyn, tfor, tprnfor
      use para_mod
      use mp, only: mp_sum
!
      implicit none
! input
      integer irb(3,nat)
      real(kind=8) rhovan(nhm*(nhm+1)/2,nat,nspin)
      complex(kind=8) eigrb(ngb,nat)
      real(kind=8)  vr(nnr,nspin)
! output
      real(kind=8)  fion(3,natx)
! local
      integer isup,isdw,iss, iv,ijv,jv, ik, nfft, isa, ia, is, ig
      integer irb3, imin3, imax3
      real(kind=8)  fvan(3,natx,nsx), fac, fac1, fac2, boxdotgrid
      complex(kind=8) ci, facg1, facg2
      complex(kind=8), allocatable :: qv(:)
      external boxdotgrid
!
      call start_clock( 'newd' )
      ci=(0.d0,1.d0)
      fac=omegab/float(nr1b*nr2b*nr3b)
      deeq (:,:,:,:) = 0.d0
      fvan (:,:,:) = 0.d0

      allocate( qv( nnrb ) )
!
! calculation of deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!
      isa=1
      do is=1,nvb
#ifdef __PARA
         do ia=1,na(is)
            nfft=1
            irb3=irb(3,isa)
            call parabox(nr3b,irb3,nr3,imin3,imax3)
            if (imax3-imin3+1.le.0) go to 15
#else
         do ia=1,na(is),2
            nfft=2
#endif
            if(ia.eq.na(is)) nfft=1
!
! two ffts at the same time, on two atoms (if possible: nfft=2)
!
            ijv=0
            do iv=1,nh(is)
               do jv=iv,nh(is)
                  ijv=ijv+1
                  qv(:) = (0.d0, 0.d0)
                  if (nfft.eq.2) then
                     do ig=1,ngb
                        qv(npb(ig))= eigrb(ig,isa  )*qgb(ig,ijv,is)   &
     &                          + ci*eigrb(ig,isa+1)*qgb(ig,ijv,is)
                        qv(nmb(ig))= conjg(                             &
     &                               eigrb(ig,isa  )*qgb(ig,ijv,is))  &
     &                          + ci*conjg(                             &
     &                               eigrb(ig,isa+1)*qgb(ig,ijv,is))
                     end do
                  else
                     do ig=1,ngb
                        qv(npb(ig)) = eigrb(ig,isa)*qgb(ig,ijv,is)
                        qv(nmb(ig)) = conjg(                            &
     &                                eigrb(ig,isa)*qgb(ig,ijv,is))
                     end do
                  end if
!
                  call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
                  do iss=1,nspin
                     deeq(iv,jv,isa,iss) = fac *                        &
     &                    boxdotgrid(irb(1,isa),1,qv,vr(1,iss))
                     if (iv.ne.jv)                                      &
     &                    deeq(jv,iv,isa,iss)=deeq(iv,jv,isa,iss)
!
                     if (nfft.eq.2) then
                        deeq(iv,jv,isa+1,iss) = fac*                    &
     &                       boxdotgrid(irb(1,isa+1),2,qv,vr(1,iss))
                        if (iv.ne.jv)                                   &
     &                       deeq(jv,iv,isa+1,iss)=deeq(iv,jv,isa+1,iss)
                     end if
                  end do
               end do
            end do
  15        isa=isa+nfft
         end do
      end do

      call reduce(nat*nhm*nhm*nspin,deeq)

      if (.not.( tfor .or. thdyn .or. tprnfor ) ) go to 10
!
! calculation of fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!
      isa=1
      if(nspin.eq.1) then
!     =================================================================
!     case nspin=1: two ffts at the same time, on two atoms (if possible)
!     -----------------------------------------------------------------
         iss=1
         isa=1
         do is=1,nvb
#ifdef __PARA
            do ia=1,na(is)
               nfft=1
               irb3=irb(3,isa)
               call parabox(nr3b,irb3,nr3,imin3,imax3)
               if (imax3-imin3+1.le.0) go to 20
#else
            do ia=1,na(is),2
               nfft=2
#endif
               if( ia.eq.na(is)) nfft=1
               do ik=1,3
                  qv(:) = (0.d0, 0.d0)
                  ijv=0
                  do iv=1,nh(is)
                     do jv=iv,nh(is)
                        ijv=ijv+1
                        if(iv.ne.jv) then
                           fac1=2.d0*fac*tpibab*rhovan(ijv,isa,iss)
                           if (nfft.eq.2) fac2=2.d0*fac*tpibab*         &
     &                                           rhovan(ijv,isa+1,iss)
                        else
                           fac1=     fac*tpibab*rhovan(ijv,isa,iss)
                           if (nfft.eq.2) fac2=     fac*tpibab*        &
     &                                           rhovan(ijv,isa+1,iss)
                        endif
                        if (nfft.eq.2) then
                           do ig=1,ngb
                              facg1 = cmplx(0.d0,-gxb(ik,ig)) *         &
     &                                   qgb(ig,ijv,is) * fac1
                              facg2 = cmplx(0.d0,-gxb(ik,ig)) *         &
     &                                   qgb(ig,ijv,is) * fac2
                              qv(npb(ig)) = qv(npb(ig))                 &
     &                                    +    eigrb(ig,isa  )*facg1  &
     &                                    + ci*eigrb(ig,isa+1)*facg2
                              qv(nmb(ig)) = qv(nmb(ig))                 &
     &                                +   conjg(eigrb(ig,isa  )*facg1)&
     &                                +ci*conjg(eigrb(ig,isa+1)*facg2)
                           end do
                        else
                           do ig=1,ngb
                              facg1 = cmplx(0.d0,-gxb(ik,ig)) *         &
     &                                   qgb(ig,ijv,is)*fac1
                              qv(npb(ig)) = qv(npb(ig))                 &
     &                                    +    eigrb(ig,isa)*facg1
                              qv(nmb(ig)) = qv(nmb(ig))                 &
     &                               +  conjg( eigrb(ig,isa)*facg1)
                           end do
                        end if
                     end do
                  end do
!
                  call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
                  fvan(ik,ia,is) =                                      &
     &                    boxdotgrid(irb(1,isa),1,qv,vr(1,iss))
!
                  if (nfft.eq.2) fvan(ik,ia+1,is) =                     &
     &                    boxdotgrid(irb(1,isa+1),2,qv,vr(1,iss))
               end do
 20            isa = isa+nfft
            end do
         end do
      else
!     =================================================================
!     case nspin=2: up and down spin fft's combined into a single fft
!     -----------------------------------------------------------------
         isup=1
         isdw=2
         isa=1
         do is=1,nvb
            do ia=1,na(is)
#ifdef __PARA
               irb3=irb(3,isa)
               call parabox(nr3b,irb3,nr3,imin3,imax3)
               if (imax3-imin3+1.le.0) go to 25
#endif
               do ik=1,3
                  qv(:) = (0.d0, 0.d0)
                  ijv=0
!
                  do iv=1,nh(is)
                     do jv=iv,nh(is)
                        ijv=ijv+1
                        if(iv.ne.jv) then
                           fac1=2.d0*fac*tpibab*rhovan(ijv,isa,isup)
                           fac2=2.d0*fac*tpibab*rhovan(ijv,isa,isdw)
                        else
                           fac1=     fac*tpibab*rhovan(ijv,isa,isup)
                           fac2=     fac*tpibab*rhovan(ijv,isa,isdw)
                        end if
                        do ig=1,ngb
                           facg1 = fac1 * cmplx(0.d0,-gxb(ik,ig)) *     &
     &                                qgb(ig,ijv,is) * eigrb(ig,isa)
                           facg2 = fac2 * cmplx(0.d0,-gxb(ik,ig)) *     &
     &                                qgb(ig,ijv,is) * eigrb(ig,isa)
                           qv(npb(ig)) = qv(npb(ig))                    &
     &                                    + facg1 + ci*facg2
                           qv(nmb(ig)) = qv(nmb(ig))                    &
     &                                    +conjg(facg1)+ci*conjg(facg2)
                        end do
                     end do
                  end do
!
                  call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
                  fvan(ik,ia,is) =                                      &
     &                    boxdotgrid(irb(1,isa),isup,qv,vr(1,isup)) + &
     &                    boxdotgrid(irb(1,isa),isdw,qv,vr(1,isdw))
               end do
25             isa = isa+1
            end do
         end do
      end if

      call reduce(3*natx*nvb,fvan)

      isa = 0
      DO is = 1, nvb
        DO ia = 1, na(is)
          isa = isa + 1
          fion(:,isa) = fion(:,isa) - fvan(:,ia,is)
        END DO
      END DO

      deallocate( qv )
!
  10  call stop_clock( 'newd' )
!
      return
      end subroutine newd
!-------------------------------------------------------------------------
      subroutine nlfl(bec,becdr,lambda,fion)
!-----------------------------------------------------------------------
!     contribution to fion due to the orthonormality constraint
! 
!
      use io_global, only: stdout
      use ions_base, only: na, nsp
      use parameters, only: natx
      use uspp, only :nhsa=>nkb, qq
      use uspp_param, only: nhm, nh
      use cvan, only: ish, nvb
      use electrons_base, only: nx => nbspx, n => nbsp
      use constants, only: pi, fpi
!
      implicit none
      real(kind=8) bec(nhsa,n), becdr(nhsa,n,3), lambda(nx,nx)
      real(kind=8) fion(3,natx)
!
      integer k, is, ia, iv, jv, i, j, inl, isa
      real(kind=8) temp(nx,nx), tmpbec(nhm,nx),tmpdr(nx,nhm) ! automatic arrays
!
      call start_clock( 'nlfl' )
      do k=1,3
         isa = 0
         do is=1,nvb
            do ia=1,na(is)
               isa = isa + 1
!
               tmpbec = 0.d0
               tmpdr  = 0.d0
!
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     inl=ish(is)+(jv-1)*na(is)+ia
                     if(abs(qq(iv,jv,is)).gt.1.e-5) then
                        do i=1,n
                           tmpbec(iv,i)=tmpbec(iv,i)                    &
     &                          + qq(iv,jv,is)*bec(inl,i)
                        end do
                     endif
                  end do
               end do
!
               do iv=1,nh(is)
                  inl=ish(is)+(iv-1)*na(is)+ia
                  do i=1,n
                     tmpdr(i,iv)=becdr(inl,i,k)
                  end do
               end do
!
               if(nh(is).gt.0)then
                  temp = 0.d0
!
                  call MXMA                                             &
     &                 (tmpdr,1,nx,tmpbec,1,nhm,temp,1,nx,n,nh(is),n)
!
                  do j=1,n
                     do i=1,n
                        temp(i,j)=temp(i,j)*lambda(i,j)
                     end do
                  end do
!
                  fion(k,isa)=fion(k,isa)+2.*SUM(temp)
               endif
!
            end do
         end do
      end do
!
!     end of x/y/z loop
!
      call stop_clock( 'nlfl' )
      return
      end subroutine nlfl
!-----------------------------------------------------------------------
      subroutine nlfq(c,eigr,bec,becdr,fion)
!-----------------------------------------------------------------------
!     contribution to fion due to nonlocal part
!
      use uspp, only :nhsa=>nkb, dvan, deeq
      use uspp_param, only: nhm, nh
      use cvan, only: ish, nvb
      use ions_base, only: nas => nax, nat, nsp, na
      use parameters, only: natx, nsx
      use electrons_base, only: n => nbsp, ispin => fspin, f
      use gvecw, only: ngw
      use constants, only: pi, fpi
      !use parm
! 
      implicit none
      real(kind=8) bec(nhsa,n), becdr(nhsa,n,3), c(2,ngw,n)
      complex(kind=8) eigr(ngw,nat)
      real(kind=8) fion(3,natx)
!
      integer k, is, ia, isa, iss, inl, iv, jv, i
      real(kind=8) tmpbec(nhm,n), tmpdr(nhm,n) ! automatic arrays
      real(kind=8) temp
!
!     nlsm2 fills becdr
!
      call start_clock( 'nlfq' )
      call nlsm2(eigr,c,becdr)
!
      do k=1,3
!
         isa=0
         do is=1,nsp
            do ia=1,na(is)
               isa=isa+1
!
               tmpbec = 0.d0
               tmpdr  = 0.d0
!
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     inl=ish(is)+(jv-1)*na(is)+ia
                     do i=1,n
                        iss=ispin(i)
                        temp=dvan(iv,jv,is)+deeq(jv,iv,isa,iss)
                        tmpbec(iv,i)=tmpbec(iv,i)+temp*bec(inl,i)
                     end do
                  end do
               end do
!  
               do iv=1,nh(is)
                  inl=ish(is)+(iv-1)*na(is)+ia
                  do i=1,n
                     tmpdr(iv,i)=f(i)*becdr(inl,i,k)
                  end do
               end do
!
               do i=1,n
                  do iv=1,nh(is)
                     tmpdr(iv,i)=tmpdr(iv,i)*tmpbec(iv,i)
                  end do
               end do
!
               fion(k,isa)=fion(k,isa)-2.*SUM(tmpdr)
!
            end do
         end do
      end do
!
!     end of x/y/z loop
!
      call stop_clock( 'nlfq' )
!
      return
      end subroutine nlfq
!-----------------------------------------------------------------------
      subroutine nlsm1 (n,nspmn,nspmx,eigr,c,becp)
!-----------------------------------------------------------------------
!     computes: the array becp
!     becp(ia,n,iv,is)=
!         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
!         = delta_l0 beta(g=0,iv,is) c(g=0,n)
!          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
!
!     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
!     input : beta(ig,l,is), eigr, c
!     output: becp as parameter
!
      use ions_base, only: na, nas => nax, nat
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use constants, only: pi, fpi
      use uspp, only :nhsa=>nkb, nhtol, beta
      use cvan, only: ish
      use uspp_param, only: nh
!
      implicit none
      integer n, nspmn, nspmx
      real(kind=8)  eigr(2,ngw,nat), c(2,ngw,n)
      real(kind=8)  becp(nhsa,n)
      complex(kind=8), allocatable :: wrk2(:,:)
!
      integer isa, ig, is, iv, ia, l, ixr, ixi, inl, i
      real(kind=8) signre, signim, arg
!
      call start_clock( 'nlsm1' )

      allocate( wrk2( ngw, nas ) )

      isa = 0
      do is = 1, nspmn - 1
        isa = isa + na(is)
      end do

      do is=nspmn,nspmx
         do iv=1,nh(is)
            l=nhtol(iv,is)
            if (l == 0) then
               ixr = 1
               ixi = 2
               signre =  1.0
               signim =  1.0
            else if (l == 1) then
               ixr = 2
               ixi = 1
               signre =  1.0
               signim = -1.0
            else if (l == 2) then
               ixr = 1
               ixi = 2
               signre = -1.0
               signim = -1.0
            else if (l == 3) then
               ixr = 2
               ixi = 1
               signre = -1.0
               signim =  1.0
            endif
!
            do ia=1,na(is)
               if (gstart == 2) then
!                   q = 0   component (with weight 1.0)
                  wrk2(1,ia)= cmplx(                                   &
     &               signre*beta(1,iv,is)*eigr(ixr,1,ia+isa),           &
     &               signim*beta(1,iv,is)*eigr(ixi,1,ia+isa) )
!                   q > 0   components (with weight 2.0)
               end if
               do ig=gstart,ngw
                  arg = 2.0*beta(ig,iv,is)
                  wrk2(ig,ia) = cmplx(                                 &
     &                  signre*arg*eigr(ixr,ig,ia+isa),                 &
     &                  signim*arg*eigr(ixi,ig,ia+isa) )
               end do
            end do
            inl=ish(is)+(iv-1)*na(is)+1
            call MXMA(wrk2,2*ngw,1,c,1,2*ngw,becp(inl,1),1,nhsa,       &
     &           na(is),2*ngw,n)
         end do

#ifdef __PARA
         inl=ish(is)+1
         do i=1,n
            call reduce(na(is)*nh(is),becp(inl,i))
         end do
#endif

        isa = isa + na(is)

      end do

      deallocate( wrk2 )
      call stop_clock( 'nlsm1' )

      return
      end subroutine nlsm1
!-------------------------------------------------------------------------
      subroutine nlsm2(eigr,c,becdr)
!-----------------------------------------------------------------------
!     computes: the array becdr 
!     becdr(ia,n,iv,is,k)
!      =2.0 sum_g> g_k beta(g,iv,is) re[ (i)**(l+1) e^(ig.r_ia) c(g,n)]
! 
!     routine makes use of  c*(g)=c(-g)  (g> see routine ggen)
!     input : eigr, c
!     output: becdr
!
      use ions_base, only: nas => nax, nsp, na, nat
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use constants, only: pi, fpi
      use uspp, only :nhsa=>nkb, nhtol, beta
      use cvan, only: ish
      use uspp_param, only: nh
      use reciprocal_vectors, only: gx
      use cell_base, only: tpiba
!
      implicit none
      real(kind=8)  eigr(2,ngw,nat),c(2,ngw,n), becdr(nhsa,n,3)
      integer ig, is, iv, ia, k, l, ixr, ixi, inl, isa
      real(kind=8) signre, signim, arg
      real(kind=8), allocatable:: gk(:)
      complex(kind=8), allocatable :: wrk2(:,:)
!
      call start_clock( 'nlsm2' )

      allocate( gk( ngw ) )
      allocate( wrk2( ngw, nas ) )

      becdr = 0.d0
!
      do k=1,3
         do ig=1,ngw
            gk(ig)=gx(k,ig)*tpiba
         end do
!
         isa = 0
         do is=1,nsp
            do iv=1,nh(is)
!
!     order of states:  s_1  p_x1  p_z1  p_y1  s_2  p_x2  p_z2  p_y2
!
               l=nhtol(iv,is)
               if (l.eq.0) then
                  ixr = 2
                  ixi = 1
                  signre =  1.0
                  signim = -1.0
               else if (l.eq.1) then
                  ixr = 1
                  ixi = 2
                  signre = -1.0
                  signim = -1.0
               else if (l.eq.2) then
                  ixr = 2
                  ixi = 1
                  signre = -1.0
                  signim =  1.0
               else if (l == 3) then
                  ixr = 1
                  ixi = 2
                  signre =  1.0
                  signim =  1.0
               endif
!     
               do ia=1,na(is)
                  if (gstart == 2) then
!                             q = 0   component (with weight 1.0)
                     wrk2(1,ia) = cmplx (                               &
     &                  signre*gk(1)*beta(1,iv,is)*eigr(ixr,1,ia+isa),   &
     &                  signim*gk(1)*beta(1,iv,is)*eigr(ixi,1,ia+isa) )
!                            q > 0   components (with weight 2.0)
                  end if
                  do ig=gstart,ngw
                     arg = 2.0*gk(ig)*beta(ig,iv,is)
                     wrk2(ig,ia) = cmplx (                              &
    &                     signre*arg*eigr(ixr,ig,ia+isa),                &
    &                     signim*arg*eigr(ixi,ig,ia+isa) )
                  end do
               end do
               inl=ish(is)+(iv-1)*na(is)+1
               call MXMA(wrk2,2*ngw,1,c,1,2*ngw,becdr(inl,1,k),1,       &
     &                   nhsa,na(is),2*ngw,n)
            end do
 
            isa = isa + na(is)

         end do
      end do

      call reduce(3*nhsa*n,becdr)

      deallocate( gk )
      deallocate( wrk2 )

      call stop_clock( 'nlsm2' )
!
      return
      end subroutine nlsm2
!-----------------------------------------------------------------------
      subroutine ortho                                                  &
     &      (eigr,cp,phi,x0,diff,iter,ccc,eps,max,delt,bephi,becp)
!-----------------------------------------------------------------------
!     input = cp (non-orthonormal), beta
!     input = phi |phi>=s'|c0>
!     output= cp (orthonormal with s( r(t+dt) ) )
!     output= bephi, becp
!     the method used is similar to the version in les houches 1988
!     'simple molecular systems at..'  p. 462-463  (18-22)
!      xcx + b x + b^t x^t + a = 1
!     where c = <s'c0|s|s'c0>   b = <s'c0|s cp>   a = <cp|s|cp>
!     where s=s(r(t+dt)) and s'=s(r(t))  
!     for vanderbilt pseudo pot - kl & ap
!
      use ions_base, only: na, nsp, nas => nax, nat
      use cvan, only: ish, nvb
      use uspp, only :nhsa=>nkb, qq
      use uspp_param, only: nh
      use electrons_base, only: n => nbsp, nx => nbspx, nspin, nupdwn, iupdwn, f
      use gvecw, only: ngw
      use control_flags, only: iprint, iprsta
      use io_global, only: stdout
!
      implicit none
!
      complex(kind=8)   cp(ngw,n), phi(ngw,n), eigr(ngw,nat)
      real(kind=8) x0(nx,nx), diff, ccc, eps, delt
      integer iter, max
      real(kind=8) bephi(nhsa,n), becp(nhsa,n)
!
      real(kind=8), allocatable :: diag(:), work1(:), work2(:), xloc(:,:), &
                                   tmp1(:,:), tmp2(:,:), dd(:,:), x1(:,:), &
                                   rhos(:,:), rhor(:,:), con(:,:), u(:,:), &
                                   sig(:,:), rho(:,:), tau(:,:)

! the above are all automatic arrays
      integer istart, nss, ifail, i, j, iss, iv, jv, ia, is, inl, jnl
      real(kind=8), allocatable:: qbephi(:,:), qbecp(:,:)

      allocate( diag(nx), work1(nx), work2(nx), xloc(nx,nx), tmp1(nx,nx),    &
                tmp2(nx,nx), dd(nx,nx), x1(nx,nx), rhos(nx,nx), rhor(nx,nx), &
                con(nx,nx),  u(nx,nx), sig(nx,nx), rho(nx,nx), tau(nx,nx) )

!
!     calculation of becp and bephi
!
      call start_clock( 'ortho' )
      call nlsm1(n,1,nvb,eigr, cp, becp)
      call nlsm1(n,1,nvb,eigr,phi,bephi)
!
!     calculation of qbephi and qbecp
!
      allocate(qbephi(nhsa,n))
      allocate(qbecp (nhsa,n))
      qbephi = 0.d0
      qbecp  = 0.d0
!
      do is=1,nvb
         do iv=1,nh(is)
            do jv=1,nh(is)
               if(abs(qq(iv,jv,is)).gt.1.e-5) then
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     do i=1,n
                        qbephi(inl,i)= qbephi(inl,i)                    &
     &                       +qq(iv,jv,is)*bephi(jnl,i)
                        qbecp (inl,i)=qbecp (inl,i)                     &
     &                       +qq(iv,jv,is)*becp (jnl,i)
                     end do
                  end do
               endif
            end do
         end do
      end do
!
      do iss=1,nspin
         nss=nupdwn(iss)
         istart=iupdwn(iss)
!
!     rho = <s'c0|s|cp>
!     sig = 1-<cp|s|cp>
!     tau = <s'c0|s|s'c0>
!
         call rhoset(cp,phi,bephi,qbecp,nss,istart,rho)
         call sigset(cp,becp,qbecp,nss,istart,sig)
         call tauset(phi,bephi,qbephi,nss,istart,tau)
!
         if(iprsta.gt.4) then
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    rho '
            do i=1,nss
               WRITE( stdout,'(7f11.6)') (rho(i,j),j=1,nss)
            end do
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    sig '
            do i=1,nss
               WRITE( stdout,'(7f11.6)') (sig(i,j),j=1,nss)
            end do
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    tau '
            do i=1,nss
               WRITE( stdout,'(7f11.6)') (tau(i,j),j=1,nss)
            end do
         endif
!
!
!----------------------------------------------------------------by ap--
! 
         do j=1,nss
            do i=1,nss
               xloc(i,j) = x0(istart-1+i,istart-1+j)*ccc
               dd(i,j) = 0.d0
               x1(i,j) = 0.d0
               tmp1(i,j)=0.d0
               rhos(i,j)=0.5d0*( rho(i,j)+rho(j,i) )
!
! on some machines (IBM RS/6000 for instance) the following test allows
! to distinguish between Numbers and Sodium Nitride (NaN, Not a Number).
! If a matrix of Not-Numbers is passed to rs, the most likely outcome is
! that the program goes on forever doing nothing and writing nothing.
!
               if (rhos(i,j).ne.rhos(i,j))                                &
     &             call errore('ortho','ortho went bananas',1)
               rhor(i,j)=rho(i,j)-rhos(i,j)
            end do
         end do
!     
         do i=1,nss
            tmp1(i,i)=1.d0
         end do
         ifail=0
         call start_clock( 'rsg' )
         call rs(nx,nss,rhos,diag,1,u,work1,work2,ifail) 
         call stop_clock( 'rsg' )
!
!                calculation of lagranges multipliers
!
         do iter=1,max
!
!       the following 4 MXMA-calls do the following matrix 
!       multiplications:
!                       tmp1 = x0*rhor    (1st call)
!                       dd   = x0*tau*x0  (2nd and 3rd call)
!                       tmp2 = x0*rhos    (4th call)
!
            call MXMA( xloc,1,nx,rhor,1,nx,tmp1,1,nx,nss,nss,nss)
            call MXMA( tau ,1,nx,xloc,1,nx,tmp2,1,nx,nss,nss,nss)
            call MXMA( xloc,1,nx,tmp2,1,nx,  dd,1,nx,nss,nss,nss)
            call MXMA( xloc,1,nx,rhos,1,nx,tmp2,1,nx,nss,nss,nss)
            do i=1,nss
               do j=1,nss
                  x1(i,j) = sig(i,j)-tmp1(i,j)-tmp1(j,i)-dd(i,j)       
                  con(i,j)= x1(i,j)-tmp2(i,j)-tmp2(j,i)
               end do
            end do
!
!         x1      = sig      -x0*rho    -x0*rho^t  -x0*tau*x0
!
            diff=0.d0
            do i=1,nss
               do j=1,nss
                  if(abs(con(i,j)).gt.diff) diff=abs(con(i,j))
               end do
            end do
!
            if( diff.le.eps ) go to 20
!     
!     the following two MXMA-calls do:   
!                       tmp1 = x1*u
!                       tmp2 = ut*x1*u
!
            call MXMA(x1,1,nx,   u,1,nx,tmp1,1,nx,nss,nss,nss)
            call MXMA(u ,nx,1,tmp1,1,nx,tmp2,1,nx,nss,nss,nss)   
!
!       g=ut*x1*u/d  (g is stored in tmp1)
! 
            do i=1,nss
               do j=1,nss
                  tmp1(i,j)=tmp2(i,j)/(diag(i)+diag(j))
               end do
            end do
!     
!       the following two MXMA-calls do:   
!                       tmp2 = g*ut
!                       x0 = u*g*ut
!
            call MXMA(tmp1,1,nx,  u,nx,1,tmp2,1,nx,nss,nss,nss)
            call MXMA(   u,1,nx,tmp2,1,nx,xloc,1,nx,nss,nss,nss)
         end do
         WRITE( stdout,*) ' diff= ',diff,' iter= ',iter
         call errore('ortho','max number of iterations exceeded',iter)
!
 20      continue
!
!-----------------------------------------------------------------------
!
         if(iprsta.gt.4) then
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    lambda '
            do i=1,nss
               WRITE( stdout,'(7f11.6)') (xloc(i,j)/f(i+istart-1),j=1,nss)
            end do
         endif
!     
         if(iprsta.gt.2) then
            WRITE( stdout,*) ' diff= ',diff,' iter= ',iter
         endif
!     
!     lagrange multipliers
!
         do i=1,nss
            do j=1,nss
               x0(istart-1+i,istart-1+j)=xloc(i,j)/ccc
               if (xloc(i,j).ne.xloc(i,j))                                &
     &             call errore('ortho','ortho went bananas',2)
            end do
         end do
!
      end do
!
      deallocate(qbecp )
      deallocate(qbephi)
      deallocate( diag, work1, work2, xloc, tmp1, tmp2, dd, x1, rhos, rhor, &
                  con, u, sig, rho, tau )
!
      call stop_clock( 'ortho' )
      return
      end subroutine ortho
!
!-----------------------------------------------------------------------
      subroutine pbc(rin,a1,a2,a3,ainv,rout)
!-----------------------------------------------------------------------
!
!     brings atoms inside the unit cell
!
      implicit none
! input
      real(kind=8) rin(3), a1(3),a2(3),a3(3), ainv(3,3)
! output
      real(kind=8) rout(3)
! local
      real(kind=8) x,y,z
!
! bring atomic positions to crystal axis
!
      x = ainv(1,1)*rin(1)+ainv(1,2)*rin(2)+ainv(1,3)*rin(3)
      y = ainv(2,1)*rin(1)+ainv(2,2)*rin(2)+ainv(2,3)*rin(3)
      z = ainv(3,1)*rin(1)+ainv(3,2)*rin(2)+ainv(3,3)*rin(3)
!
! bring x,y,z in the range between -0.5 and 0.5
!
      x = x - nint(x)
      y = y - nint(y)
      z = z - nint(z)
!
! bring atomic positions back in cartesian axis
!
      rout(1) = x*a1(1)+y*a2(1)+z*a3(1)
      rout(2) = x*a1(2)+y*a2(2)+z*a3(2)
      rout(3) = x*a1(3)+y*a2(3)+z*a3(3)
!
      return
      end subroutine pbc

!
!-------------------------------------------------------------------------
      subroutine prefor(eigr,betae)
!-----------------------------------------------------------------------
!
!     input :        eigr =  e^-ig.r_i
!     output:        betae_i,i(g) = (-i)**l beta_i,i(g) e^-ig.r_i 
!
      use ions_base, only: nas => nax, nsp, na, nat
      use gvecw, only: ngw
      use cvan, only: ish
      use uspp, only :nhsa=>nkb, beta, nhtol
      use uspp_param, only: nh
!
      implicit none
      complex(kind=8) eigr(ngw,nat)
      complex(kind=8) betae(ngw,nhsa)
!
      integer is, iv, ia, inl, ig, isa
      complex(kind=8) ci
!
      call start_clock( 'prefor' )
      isa = 0
      do is=1,nsp
         do iv=1,nh(is)
            ci=(0.,-1.)**nhtol(iv,is)
            do ia=1,na(is)
               inl=ish(is)+(iv-1)*na(is)+ia
               do ig=1,ngw
                  betae(ig,inl)=ci*beta(ig,iv,is)*eigr(ig,ia+isa)
               end do
            end do
         end do
         isa = isa + na(is)
      end do
      call stop_clock( 'prefor' )
!
      return
      end subroutine prefor
!
!-----------------------------------------------------------------------
      subroutine projwfc(c,eigr,betae)
!-----------------------------------------------------------------------
!
! Projection on atomic wavefunctions
!
      use io_global, only: stdout
      use electrons_base, only: nx => nbspx, n => nbsp
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use ions_base, only: nsp, na, nas => nax, nat
      use uspp, only: nhsa => nkb
      use atom
!
      implicit none
      complex(kind=8), intent(in) :: c(ngw,nx), eigr(ngw,nat),      &
     &                               betae(ngw,nhsa)
!
      complex(kind=8), allocatable:: wfc(:,:), swfc(:,:), becwfc(:,:)
      real(kind=8), allocatable   :: overlap(:,:), e(:), z(:,:),        &
     &                               proj(:,:), temp(:)
      real(kind=8)                :: somma
      integer n_atomic_wfc
      integer is, ia, nb, l, m, k, i
!
! calculate number of atomic states
!
      n_atomic_wfc=0
      do is=1,nsp
         do nb = 1,nchi(is)
            l = lchi(nb,is)
            n_atomic_wfc = n_atomic_wfc + (2*l+1)*na(is)
         end do
      end do
      if (n_atomic_wfc.eq.0) return
!
      allocate(wfc(ngw,n_atomic_wfc))
!
! calculate wfc = atomic states
!
      call atomic_wfc(eigr,n_atomic_wfc,wfc)
!
! calculate bec = <beta|wfc>
!
      allocate(becwfc(nhsa,n_atomic_wfc))
      call nlsm1 (n_atomic_wfc,1,nsp,eigr,wfc,becwfc)

! calculate swfc = S|wfc>
!
      allocate(swfc(ngw,n_atomic_wfc))
      call s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!
! calculate overlap(i,j) = <wfc_i|S|wfc_j> 
!
      allocate(overlap(n_atomic_wfc,n_atomic_wfc))
!
      call MXMA(wfc,2*ngw,1,swfc,1,2*ngw,overlap,1,                     &
     &          n_atomic_wfc,n_atomic_wfc,2*ngw,n_atomic_wfc)

      call reduce(n_atomic_wfc**2,overlap)

      overlap=overlap*2.d0
      if (gstart == 2) then
         do l=1,n_atomic_wfc
            do m=1,n_atomic_wfc
               overlap(m,l)=overlap(m,l)-real(wfc(1,m))*real(swfc(1,l))
            end do
         end do
      end if
!
! calculate (overlap)^(-1/2)(i,j). An orthonormal set of vectors |wfc_i>
! is obtained by introducing |wfc_j>=(overlap)^(-1/2)(i,j)*S|wfc_i>
!
      allocate(z(n_atomic_wfc,n_atomic_wfc))
      allocate(e(n_atomic_wfc))
      call rdiag(n_atomic_wfc,overlap,n_atomic_wfc,e,z)
      overlap=0.d0
      do l=1,n_atomic_wfc
         do m=1,n_atomic_wfc
            do k=1,n_atomic_wfc
               overlap(l,m)=overlap(l,m)+z(m,k)*z(l,k)/sqrt(e(k))
            end do
         end do
      end do
      deallocate(e)
      deallocate(z)
!
! calculate |wfc_j>=(overlap)^(-1/2)(i,j)*S|wfc_i>   (note the S matrix!)
!
      wfc=0.d0
      do m=1,n_atomic_wfc
         do l=1,n_atomic_wfc
            wfc(:,m)=wfc(:,m)+overlap(l,m)*swfc(:,l)
         end do
      end do
      deallocate(overlap)
      deallocate(swfc)
      deallocate(becwfc)
!
! calculate proj = <c|S|wfc> 
!
      allocate(proj(n,n_atomic_wfc))
      allocate(temp(ngw))
      do m=1,n
         do l=1,n_atomic_wfc
            temp(:)=real(conjg(c(:,m))*wfc(:,l))
            proj(m,l)=2.d0*SUM(temp)
            if (gstart == 2) proj(m,l)=proj(m,l)-temp(1)
         end do
      end do
      deallocate(temp)

      call reduce(n*n_atomic_wfc,proj)

      i=0
      WRITE( stdout,'(/''Projection on atomic states:'')')
      do is=1,nsp
         do nb = 1,nchi(is)
            l=lchi(nb,is)
            do m = -l,l
               do ia=1,na(is)
                  i=i+1
                  WRITE( stdout,'(''atomic state # '',i3,'': atom # '',i3,    &
     &                      ''  species # '',i2,''  wfc # '',i2,        &
     &                      '' (l='',i1,'' m='',i2,'')'')')             &
     &                 i, ia, is, nb, l, m
               end do
            end do
         end do
      end do

      WRITE( stdout,*)
      do m=1,n
         somma=0.d0
         do l=1,n_atomic_wfc
            somma=somma+proj(m,l)**2
         end do
         WRITE( stdout,'(''state # '',i4,''    sum c^2 ='',f7.4)') m,somma
         WRITE( stdout,'(10f7.4)') (abs(proj(m,l)),l=1,n_atomic_wfc)
      end do
!
      deallocate(proj)
      deallocate(wfc)
      return
      end subroutine projwfc
!-----------------------------------------------------------------------
      subroutine raddrizza(nspin,nx,nupdwn,iupdwn,f,lambda,ngw,c)
!-----------------------------------------------------------------------
!
!     transform wavefunctions into eigenvectors of the hamiltonian
!     via diagonalization of the constraint matrix lambda
!
      implicit none
      integer, intent(in)           :: nspin, nx, ngw, nupdwn(nspin),   &
     &                                 iupdwn(nspin)
      real   (kind=8), intent(in)   :: lambda(nx,nx), f(nx)
      complex(kind=8), intent(inout):: c(ngw,nx)

      real(kind=8)                :: lambdar(nx,nx), wr(nx), zr(nx,nx)
      complex(kind=8), allocatable:: csave(:,:)
      integer                     :: iss, n, j, i, i0
!
      do iss=1,nspin
         n=nupdwn(iss)
         i0=iupdwn(iss)-1
         allocate(csave(ngw,n))
         do i=1,n
            do j=1,n
               lambdar(j,i)=lambda(i0+j,i0+i)
            end do
         end do

         call rdiag(n,lambdar,nx,wr,zr)

         csave=0.d0
         do i=1,n
            do j=1,n
               csave(:,i) = csave(:,i) + zr(j,i)*c(:,i0+j)
            end do
         end do
         do i=1,n
            c(:,i0+i)=csave(:,i)
         end do
         deallocate(csave)

!     uncomment to print out eigenvalues
!         do i=1,n
!            if (f(i0+i).gt.1.e-6) then
!               wr(i)=27.212*wr(i)/f(i0+i)
!            else
!               wr(i)=0.0
!            end if
!         end do
!         WRITE( stdout,'(/10f8.2/)') (wr(i),i=1,nupdwn(iss))
      end do
      return
      end subroutine raddrizza
!
!---------------------------------------------------------------------
      subroutine randin(nmin,nmax,gstart,ngw,ampre,c)
!---------------------------------------------------------------------
!
      use wave_functions, only: wave_rand_init
      implicit none

! input
      integer nmin, nmax, gstart, ngw
      real(kind=8) ampre
! output
      complex(kind=8) c(ngw,nmax)
! local
      integer i,j
      real(kind=8) ranf1, randy, ranf2, ampexp
!
      CALL wave_rand_init( c )
!      do i=nmin,nmax
!         do j=gstart,ngw
!            ranf1=.5-randy()
!            ranf2=.5-randy()
!            ampexp=ampre*exp(-(4.*j)/ngw)
!            c(j,i)=c(j,i)+ampexp*cmplx(ranf1,ranf2)
!         end do
!      end do
!
      return
      end subroutine randin
!
!-----------------------------------------------------------------------
      subroutine rdiag (n,h,ldh,e,v)
!-----------------------------------------------------------------------
!
!   calculates all the eigenvalues and eigenvectors of a complex
!   hermitean matrix H . On output, the matrix H is destroyed
!
      implicit none
      integer, intent(in)           :: n, ldh
      complex(kind=8), intent(inout):: h(ldh,n)
      real   (kind=8), intent(out)  :: e(n)
      complex(kind=8), intent(out)  :: v(ldh,n)
!
      real(kind=8) fv1(n), fv2(n)
      integer ierr
!
      call rs(ldh,n,h,e,1,v,fv1,fv2,ierr)
!
      return
      end subroutine rdiag
!-----------------------------------------------------------------------
   subroutine rhoofr (nfi,c,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
!-----------------------------------------------------------------------
!     the normalized electron density rhor in real space
!     the kinetic energy ekin
!     subroutine uses complex fft so it computes two ft's
!     simultaneously
!
!     rho_i,ij = sum_n < beta_i,i | psi_n >< psi_n | beta_i,j >
!     < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                   2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
!
!     e_v = sum_i,ij rho_i,ij d^ion_is,ji
!
      use control_flags, only: iprint, tbuff, iprsta, thdyn, tpre, trhor
      use ions_base, only: nat, nas => nax, nsp
      use parameters, only: natx, nsx
      use gvecp, only: ng => ngm
      use gvecs
      use gvecb, only: ngb
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use recvecs_indexes, only: np, nm
      use uspp, only: nhsa => nkb
      use uspp_param, only: nh, nhm
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use cell_base, only: omega
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base, only: nx => nbspx, n => nbsp, f, ispin => fspin, nspin
      use constants, only: pi, fpi
      use mp, ONLY: mp_sum
      ! use local_pseudo
!
      use cdvan
      use dener
      use io_global, only: stdout
!
      implicit none
      real(kind=8) bec(nhsa,n), rhovan(nhm*(nhm+1)/2,nat,nspin)
      real(kind=8) rhor(nnr,nspin), rhos(nnrsx,nspin)
      real(kind=8) enl, ekin
      complex(kind=8) eigrb(ngb,nat), c(ngw,nx), rhog(ng,nspin)
      integer irb(3,nat), nfi
! local variables
      integer iss, isup, isdw, iss1, iss2, ios, i, ir, ig
      real(kind=8) rsumr(2), rsumg(2), sa1, sa2
      real(kind=8) rnegsum, rmin, rmax, rsum
      real(kind=8), external :: enkin, ennl
      complex(kind=8) ci,fp,fm
      complex(kind=8), allocatable :: psi(:), psis(:)
!
!
      call start_clock( 'rhoofr' )
      allocate( psi( nnr ) ) 
      allocate( psis( nnrsx ) ) 
      ci=(0.0,1.0)
      do iss=1,nspin
         rhor(:,iss) = 0.d0
         rhos(:,iss) = 0.d0
         rhog(:,iss) = (0.d0, 0.d0)
      end do
!
!     ==================================================================
!     calculation of kinetic energy ekin
!     ==================================================================
      ekin=enkin(c)
      if(tpre) call denkin(c,dekin)
!
!     ==================================================================
!     calculation of non-local energy
!     ==================================================================
      enl=ennl(rhovan, bec)
      if(tpre) call dennl(bec,denl)
!    
!    warning! trhor and thdyn are not compatible yet!   
!
      if(trhor.and.(.not.thdyn))then
!     ==================================================================
!     charge density is read from unit 47
!     ==================================================================
#ifdef __PARA
         call read_rho(47,nspin,rhor)
#else
         read(47) ((rhor(ir,iss),ir=1,nnr),iss=1,nspin)
#endif
         rewind 47
!
         if(nspin.eq.1)then
            iss=1
            do ir=1,nnr
               psi(ir)=cmplx(rhor(ir,iss),0.)
            end do
            call fwfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ig=1,ng
               rhog(ig,iss)=psi(np(ig))
            end do
         else
            isup=1
            isdw=2
            do ir=1,nnr
               psi(ir)=cmplx(rhor(ir,isup),rhor(ir,isdw))
            end do
            call fwfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ig=1,ng
               fp=psi(np(ig))+psi(nm(ig))
               fm=psi(np(ig))-psi(nm(ig))
               rhog(ig,isup)=0.5*cmplx( real(fp),aimag(fm))
               rhog(ig,isdw)=0.5*cmplx(aimag(fp),-real(fm))
            end do
         endif
!
      else

         !     ==================================================================
         !     self-consistent charge
         !     ==================================================================
         !
         !     important: if n is odd then nx must be .ge.n+1 and c(*,n+1)=0.
         ! 
         if (mod(n,2).ne.0) then
            do ig=1,ngw
               c(ig,n+1)=(0.,0.)
            end do
         endif
         !
         do i=1,n,2
            psis (:) = (0.d0, 0.d0)
            do ig=1,ngw
               psis(nms(ig))=conjg(c(ig,i))+ci*conjg(c(ig,i+1))
               psis(nps(ig))=c(ig,i)+ci*c(ig,i+1)
               ! write(6,'(I6,4F15.10)') ig, psis(nms(ig)), psis(nps(ig))
            end do

            call ivfftw(psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)

            !     wavefunctions in unit 21
            !
#if defined(__CRAYY)
            if(tbuff) buffer out(21,0) (psis(1),psis(nnrsx))
#else
            if(tbuff) write(21,iostat=ios) psis
#endif
            iss1=ispin(i)
            sa1=f(i)/omega
            if (i.ne.n) then
               iss2=ispin(i+1)
               sa2=f(i+1)/omega
            else
               iss2=iss1
               sa2=0.0
            end if
            do ir=1,nnrsx
               rhos(ir,iss1)=rhos(ir,iss1) + sa1*( real(psis(ir)))**2
               rhos(ir,iss2)=rhos(ir,iss2) + sa2*(aimag(psis(ir)))**2
            end do

            !
            !       buffer 21
            !     
            if(tbuff) then
#if defined(__CRAYY)
               ios=unit(21)
#endif
               if(ios.ne.0) call errore(' rhoofr',' error in writing unit 21',ios)
            endif
            !
         end do
         !
         if(tbuff) rewind 21
         !
         !     smooth charge in g-space is put into rhog(ig)
         !
         if(nspin.eq.1)then
            iss=1
            do ir=1,nnrsx
               psis(ir)=cmplx(rhos(ir,iss),0.)
            end do
            call fwffts(psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
            do ig=1,ngs
               rhog(ig,iss)=psis(nps(ig))
            end do
         else
            isup=1
            isdw=2
             do ir=1,nnrsx
               psis(ir)=cmplx(rhos(ir,isup),rhos(ir,isdw))
            end do
            call fwffts(psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
            do ig=1,ngs
               fp= psis(nps(ig)) + psis(nms(ig))
               fm= psis(nps(ig)) - psis(nms(ig))
               rhog(ig,isup)=0.5*cmplx( real(fp),aimag(fm))
               rhog(ig,isdw)=0.5*cmplx(aimag(fp),-real(fm))
            end do
         endif
!
         if(nspin.eq.1) then
            ! 
            !     case nspin=1
            ! 
            iss=1
            psi (:) = (0.d0, 0.d0)
            do ig=1,ngs
               psi(nm(ig))=conjg(rhog(ig,iss))
               psi(np(ig))=      rhog(ig,iss)
            end do
            call invfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ir=1,nnr
               rhor(ir,iss)=real(psi(ir))
            end do
         else 
            !
            !     case nspin=2
            !
            isup=1
            isdw=2
            psi (:) = (0.d0, 0.d0)
            do ig=1,ngs
               psi(nm(ig))=conjg(rhog(ig,isup))+ci*conjg(rhog(ig,isdw))
               psi(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            end do
            call invfft(psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            do ir=1,nnr
               rhor(ir,isup)= real(psi(ir))
               rhor(ir,isdw)=aimag(psi(ir))
            end do
         endif
!
         if(iprsta.ge.3)then
            do iss=1,nspin
               rsumg(iss)=omega*real(rhog(1,iss))
               rsumr(iss)=SUM(rhor(:,iss))*omega/dble(nr1*nr2*nr3)
            end do

            if ( gstart /= 2 ) then
               !
               !    in the parallel case, only one processor has G=0 ! 
               !
               do iss=1,nspin
                  rsumg(iss)=0.0
               end do
            end if
            call mp_sum( rsumg( 1:nspin ) )
            call mp_sum( rsumr( 1:nspin ) )

            if ( nspin == 1 ) then
              WRITE( stdout, 10) rsumg(1), rsumr(1)
            else
              WRITE( stdout, 20) rsumg(1), rsumr(1), rsumg(2), rsumr(2)
            endif

         endif
         !
         !     add vanderbilt contribution to the charge density
         !     drhov called before rhov because input rho must be the smooth part
         !
         if (tpre) call drhov(irb,eigrb,rhovan,rhog,rhor)
         !
         call rhov(irb,eigrb,rhovan,rhog,rhor)

      endif

!     ======================================endif for trhor=============
!
!     here to check the integral of the charge density
!
!
      if(iprsta.ge.2) then
         call checkrho(nnr,nspin,rhor,rmin,rmax,rsum,rnegsum)
         rnegsum=rnegsum*omega/dble(nr1*nr2*nr3)
         rsum=rsum*omega/dble(nr1*nr2*nr3)
         WRITE( stdout,'(a,4(1x,f12.6))')                                     &
     &     ' rhoofr: rmin rmax rnegsum rsum  ',rmin,rmax,rnegsum,rsum
      end if
!
      if( nfi == 0 .or. mod(nfi, iprint) == 0 ) then

         do iss=1,nspin
            rsumg(iss)=omega*real(rhog(1,iss))
            rsumr(iss)=SUM(rhor(:,iss),1)*omega/dble(nr1*nr2*nr3)
         end do

         if (gstart.ne.2) then
            ! in the parallel case, only one processor has G=0 ! 
            do iss=1,nspin
               rsumg(iss)=0.0
            end do
         end if

         call mp_sum( rsumg( 1:nspin ) )
         call mp_sum( rsumr( 1:nspin ) )

         if ( nspin == 1 ) then
           WRITE( stdout, 10) rsumg(1), rsumr(1)
         else
           WRITE( stdout, 20) rsumg(1), rsumr(1), rsumg(2), rsumr(2)
         endif

      endif

      deallocate( psi ) 
      deallocate( psis ) 

10    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 )
20    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'spin up', &
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 , &
            & /, 3X, 'spin down', &
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 )
!
      call stop_clock( 'rhoofr' )

!
      return
      end subroutine rhoofr
!
!-----------------------------------------------------------------------
      subroutine rhoset(cp,phi,bephi,qbecp,nss,ist,rho)
!-----------------------------------------------------------------------
!     input: cp (non-orthonormal), phi, bephi, qbecp
!     computes the matrix
!       rho = <s'c0|s cp> = <phi|s cp>
!     where  |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
!     where s=s(r(t+dt)) and s'=s(r(t))  
!     routine makes use of  c(-q)=c*(q)
!
      use parameters, only: nsx, natx
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use uspp, only: nhsa => nkb, nhsavb=>nkbus
      use cvan, only: nvb
      use electrons_base, only: nx => nbspx, n => nbsp ! , f, ispin => fspin, nspin
!
      implicit none
!
      integer nss, ist
      complex(kind=8)   cp(ngw,n), phi(ngw,n)
      real(kind=8)       bephi(nhsa,n), qbecp(nhsa,n), rho(nx,nx)
      integer i, j
      real(kind=8)    tmp1(nx,nx) ! automatic array
!
      rho (:,:) = 0.d0
!
!     <phi|cp>
!
      call MXMA(phi(1,ist),2*ngw,1,cp(1,ist),1,2*ngw,                   &
     &                              rho,1,nx,nss,2*ngw,nss)
!
!     q >= 0  components with weight 2.0
!
      do j=1,nss
         do i=1,nss
            rho(i,j)=2.*rho(i,j)
         end do
      end do
!
      if (gstart == 2) then
!
!     q = 0  components has weight 1.0
!
         do j=1,nss
            do i=1,nss
               rho(i,j) = rho(i,j) -                                    &
     &              real(phi(1,i+ist-1))*real(cp(1,j+ist-1))
            end do
         end do
      end if

      call reduce(nx*nss,rho)
!
      if(nvb.gt.0)then
         tmp1 (:,:) = 0.d0
!
         call MXMA(bephi(1,ist),nhsa,1,qbecp(1,ist),1,nhsa,               &
     &                                tmp1,1,nx,nss,nhsavb,nss)
!
         do j=1,nss
            do i=1,nss
               rho(i,j)=rho(i,j)+tmp1(i,j)
            end do
         end do
      endif
!
      return
      end subroutine rhoset
!
!-----------------------------------------------------------------------
      subroutine rhov(irb,eigrb,rhovan,rhog,rhor)
!-----------------------------------------------------------------------
!     Add Vanderbilt contribution to rho(r) and rho(g)
!
!        n_v(g) = sum_i,ij rho_i,ij q_i,ji(g) e^-ig.r_i
!
!     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
!
      use ions_base, only: nas => nax, nat, na, nsp
      use io_global, only: stdout
      use parameters, only: natx, nsx
      use cvan, only: nvb
      use uspp_param, only: nh, nhm
      use uspp, only: deeq
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use electrons_base, only: nspin
      use gvecb
      use gvecp, only: ng => ngm
      use cell_base, only: omega
      use small_box, only: omegab
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx
      use control_flags, only: iprint, iprsta
      use qgb_mod
      use para_mod
      use recvecs_indexes, only: np, nm
!
      implicit none
!
      real(kind=8) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      integer, intent(in) :: irb(3,nat)
      complex(kind=8), intent(in):: eigrb(ngb,nat)
      real(kind=8), intent(inout):: rhor(nnr,nspin)
      complex(kind=8),  intent(inout):: rhog(ng,nspin)
!
      integer isup, isdw, nfft, ifft, iv, jv, ig, ijv, is, iss,           &
     &     isa, ia, ir, irb3, imin3, imax3
      real(kind=8) sumrho
      complex(kind=8) ci, fp, fm, ca
      complex(kind=8), allocatable::  qgbt(:,:)
      complex(kind=8), allocatable:: v(:)
      complex(kind=8), allocatable:: qv(:)
!
      if (nvb.eq.0) return
      call start_clock( 'rhov' )
      ci=(0.,1.)
!
!
      allocate( v( nnr ) )
      allocate( qv( nnrb ) )
      v (:) = (0.d0, 0.d0)
      allocate( qgbt( ngb, 2 ) )

!
      if(nspin.eq.1) then
         ! 
         !     nspin=1 : two fft at a time, one per atom, if possible
         !
         iss=1
         isa=1

         do is = 1, nvb

#ifdef __PARA

            do ia=1,na(is)
               nfft=1
               irb3=irb(3,isa)
               call parabox(nr3b,irb3,nr3,imin3,imax3)
               if (imax3-imin3+1.le.0) go to 15
#else

            do ia = 1, na(is), 2
               nfft = 2
               if( ia .eq. na(is) ) nfft = 1

#endif

               !
               !  nfft=2 if two ffts at the same time are performed
               !
               do ifft=1,nfft
                  qgbt(:,ifft) = (0.d0, 0.d0)
                  ijv=0
                  do iv= 1,nh(is)
                     do jv=iv,nh(is)
                        ijv=ijv+1
                        sumrho=rhovan(ijv,isa+ifft-1,iss)
                        if(iv.ne.jv) sumrho=2.*sumrho
                        do ig=1,ngb
                           qgbt(ig,ifft)=qgbt(ig,ifft) + sumrho*qgb(ig,ijv,is)
                        end do
                     end do
                  end do
               end do
               !
               ! add structure factor
               !
               qv(:) = (0.d0, 0.d0)
               if(nfft.eq.2)then
                  do ig=1,ngb
                     qv(npb(ig))=  &
                                   eigrb(ig,isa  )*qgbt(ig,1)  &
                        + ci*      eigrb(ig,isa+1)*qgbt(ig,2)
                     qv(nmb(ig))=                                       &
                             conjg(eigrb(ig,isa  )*qgbt(ig,1))        &
                        + ci*conjg(eigrb(ig,isa+1)*qgbt(ig,2))
                  end do
               else
                  do ig=1,ngb
                     qv(npb(ig)) = eigrb(ig,isa)*qgbt(ig,1)
                     qv(nmb(ig)) = conjg(eigrb(ig,isa)*qgbt(ig,1))
                  end do
               endif

               call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)

               !
               !  qv = US augmentation charge in real space on box grid
               !       for atomic species is, real(qv)=atom ia, imag(qv)=atom ia+1
 
               if(iprsta.gt.2) then
                  ca = SUM(qv)
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom g-sp = ',         &
     &                 omegab*real(qgbt(1,1))
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom r-sp = ',         &
     &                 omegab*real(ca)/(nr1b*nr2b*nr3b)
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom g-sp = ',         &
     &                 omegab*real(qgbt(1,2))
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom r-sp = ',         &
     &                 omegab*aimag(ca)/(nr1b*nr2b*nr3b)
               endif
               !
               !  add qv(r) to v(r), in real space on the dense grid
               !
               call  box2grid(irb(1,isa),1,qv,v)
               if (nfft.eq.2) call  box2grid(irb(1,isa+1),2,qv,v)
  15           isa=isa+nfft
!
            end do
         end do
         !
         !  rhor(r) = total (smooth + US) charge density in real space
         !
         do ir=1,nnr
            rhor(ir,iss)=rhor(ir,iss)+real(v(ir))        
         end do
!
         if(iprsta.gt.2) then
            ca = SUM(v)

            call reduce(2,ca)

            WRITE( stdout,'(a,2f12.8)')                                  &
     &           ' rhov: int  n_v(r)  dr = ',omega*ca/(nr1*nr2*nr3)
         endif
!
         call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         if(iprsta.gt.2) then
            WRITE( stdout,*) ' rhov: smooth ',omega*rhog(1,iss)
            WRITE( stdout,*) ' rhov: vander ',omega*v(1)
            WRITE( stdout,*) ' rhov: all    ',omega*(rhog(1,iss)+v(1))
         endif
         !
         !  rhog(g) = total (smooth + US) charge density in G-space
         !
         do ig=1,ng
            rhog(ig,iss)=rhog(ig,iss)+v(np(ig))
         end do
!
         if(iprsta.gt.1) WRITE( stdout,'(a,2f12.8)')                          &
     &        ' rhov: n_v(g=0) = ',omega*real(rhog(1,iss))
!
      else
         !
         !     nspin=2: two fft at a time, one for spin up and one for spin down
         !
         isup=1
         isdw=2
         isa=1
         do is=1,nvb
            do ia=1,na(is)
#ifdef __PARA
               irb3=irb(3,isa)
               call parabox(nr3b,irb3,nr3,imin3,imax3)
               if (imax3-imin3+1.le.0) go to 25
#endif
               do iss=1,2
                  qgbt(:,iss) = (0.d0, 0.d0)
                  ijv=0
                  do iv=1,nh(is)
                     do jv=iv,nh(is)
                        ijv=ijv+1
                        sumrho=rhovan(ijv,isa,iss)
                        if(iv.ne.jv) sumrho=2.*sumrho
                        do ig=1,ngb
                           qgbt(ig,iss)=qgbt(ig,iss)+sumrho*qgb(ig,ijv,is)
                        end do
                     end do
                  end do
               end do
!     
! add structure factor
!
               qv(:) = (0.d0, 0.d0)
               do ig=1,ngb
                  qv(npb(ig)) =    eigrb(ig,isa)*qgbt(ig,1)           &
     &                  + ci*      eigrb(ig,isa)*qgbt(ig,2)
                  qv(nmb(ig)) = conjg(eigrb(ig,isa)*qgbt(ig,1))       &
     &                  + ci*   conjg(eigrb(ig,isa)*qgbt(ig,2))
               end do
!
               call ivfftb(qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!
!  qv is the now the US augmentation charge for atomic species is
!  and atom ia: real(qv)=spin up, imag(qv)=spin down
!
               if(iprsta.gt.2) then
                  ca = SUM(qv)
                  WRITE( stdout,'(a,f12.8)') ' rhov: up   g-space = ',        &
     &                 omegab*real(qgbt(1,1))
                  WRITE( stdout,'(a,f12.8)') ' rhov: up r-sp = ',             &
     &                 omegab*real(ca)/(nr1b*nr2b*nr3b)
                  WRITE( stdout,'(a,f12.8)') ' rhov: dw g-space = ',          &
     &                 omegab*real(qgbt(1,2))
                  WRITE( stdout,'(a,f12.8)') ' rhov: dw r-sp = ',             &
     &                 omegab*aimag(ca)/(nr1b*nr2b*nr3b)
               endif
!
!  add qv(r) to v(r), in real space on the dense grid
!
               call box2grid2(irb(1,isa),qv,v)
  25           isa=isa+1
!
            end do
         end do
!
         do ir=1,nnr
            rhor(ir,isup)=rhor(ir,isup)+real(v(ir)) 
            rhor(ir,isdw)=rhor(ir,isdw)+aimag(v(ir)) 
         end do
!
         if(iprsta.gt.2) then
            ca = SUM(v)
            call reduce(2,ca)
            WRITE( stdout,'(a,2f12.8)') 'rhov:in n_v  ',omega*ca/(nr1*nr2*nr3)
         endif
!
         call fwfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         if(iprsta.gt.2) then
            WRITE( stdout,*) 'rhov: smooth up',omega*rhog(1,isup)
            WRITE( stdout,*) 'rhov: smooth dw',omega*rhog(1,isdw)
            WRITE( stdout,*) 'rhov: vander up',omega*real(v(1))
            WRITE( stdout,*) 'rhov: vander dw',omega*aimag(v(1))
            WRITE( stdout,*) 'rhov: all up',                                  &
     &           omega*(rhog(1,isup)+real(v(1)))
            WRITE( stdout,*) 'rhov: all dw',                                  &
     &           omega*(rhog(1,isdw)+aimag(v(1)))
         endif
!
         do ig=1,ng
            fp=  v(np(ig)) + v(nm(ig))
            fm=  v(np(ig)) - v(nm(ig))
            rhog(ig,isup)=rhog(ig,isup) + 0.5*cmplx(real(fp),aimag(fm))
            rhog(ig,isdw)=rhog(ig,isdw) + 0.5*cmplx(aimag(fp),-real(fm))
         end do
!
         if(iprsta.gt.2) WRITE( stdout,'(a,2f12.8)')                          &
     &        ' rhov: n_v(g=0) up   = ',omega*real (rhog(1,isup))
         if(iprsta.gt.2) WRITE( stdout,'(a,2f12.8)')                          &
     &        ' rhov: n_v(g=0) down = ',omega*real(rhog(1,isdw))
!
      endif

      deallocate(qgbt)
      deallocate( v )
      deallocate( qv )

      call stop_clock( 'rhov' )
!
      return
      end subroutine rhov
!
!
!-------------------------------------------------------------------------
      subroutine s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!-----------------------------------------------------------------------
!
!     input: wfc, becwfc=<wfc|beta>, betae=|beta>
!     output: swfc=S|wfc>
!
      use ions_base, only: na
      use cvan, only: nvb, ish
      use uspp, only: nhsa => nkb, nhsavb=>nkbus, qq
      use uspp_param, only: nh
      use gvecw, only: ngw
      !use parm
      use constants, only: pi, fpi
      implicit none
! input
      integer, intent(in)         :: n_atomic_wfc
      complex(kind=8), intent(in) :: betae(ngw,nhsa),                   &
     &                               wfc(ngw,n_atomic_wfc)
      real(kind=8), intent(in)    :: becwfc(nhsa,n_atomic_wfc)
! output
      complex(kind=8), intent(out):: swfc(ngw,n_atomic_wfc)
! local
      integer is, iv, jv, ia, inl, jnl, i
      real(kind=8) qtemp(nhsavb,n_atomic_wfc)
!
      swfc=0.d0
!
      if (nvb.gt.0) then
         qtemp=0.d0
         do is=1,nvb
            do iv=1,nh(is)
               do jv=1,nh(is)
                  if(abs(qq(iv,jv,is)).gt.1.e-5) then
                     do ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        do i=1,n_atomic_wfc
                           qtemp(inl,i) = qtemp(inl,i) +                &
     &                                    qq(iv,jv,is)*becwfc(jnl,i)
                        end do
                     end do
                  endif
               end do
            end do
         end do
!
         call MXMA (betae,1,2*ngw,qtemp,1,nhsavb,swfc,1,                &
     &              2*ngw,2*ngw,nhsavb,n_atomic_wfc)
      end if
!
      swfc=swfc+wfc
!
      return
      end subroutine s_wfc

!
!-------------------------------------------------------------------------
      subroutine sigset(cp,becp,qbecp,nss,ist,sig)
!-----------------------------------------------------------------------
!     input: cp (non-orthonormal), becp, qbecp
!     computes the matrix
!       sig = 1 - a ,  a = <cp|s|cp> = <cp|cp> + sum q_ij <cp|i><j|cp>
!     where s=s(r(t+dt)) 
!     routine makes use of c(-q)=c*(q)
!
      use parameters, only: natx, nsx
      use uspp, only: nhsa => nkb, nhsavb=>nkbus
      use cvan, only : nvb
      use electrons_base, only: nx => nbspx, n => nbsp
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
!
      implicit none
!
      integer nss, ist
      complex(kind=8)  cp(ngw,n)
      real(kind=8) becp(nhsa,n), qbecp(nhsa,n), sig(nx,nx)
!
      integer i, j
      real(kind=8)    tmp1(nx,nx) ! automatic array
!
      sig = 0.d0
      call MXMA(cp(1,ist),2*ngw,1,cp(1,ist),1,2*ngw,                    &
     &                                  sig,1,nx,nss,2*ngw,nss)
!
!     q >= 0  components with weight 2.0
!
      do j=1,nss
         do i=1,nss
            sig(i,j)=-2.*sig(i,j)
         end do
      end do
      if (gstart == 2) then
!
!     q = 0  components has weight 1.0
!
         do j=1,nss
            do i=1,nss
               sig(i,j) = sig(i,j) +                                    &
     &              real(cp(1,i+ist-1))*real(cp(1,j+ist-1))
            end do
         end do
      end if
      call reduce(nx*nss,sig)
      do i=1,nss
         sig(i,i) = sig(i,i)+1.
      end do
!
      if(nvb.gt.0)then
         tmp1 = 0.d0
!
         call MXMA(becp(1,ist),nhsa,1,qbecp(1,ist),1,nhsa,                &
     &                              tmp1,1,nx,nss,nhsavb,nss)
!
         do j=1,nss
            do i=1,nss
               sig(i,j)=sig(i,j)-tmp1(i,j)
            end do
         end do
      endif
!
      return
      end subroutine sigset
!
!-----------------------------------------------------------------------
      subroutine spinsq (c,bec,rhor)
!-----------------------------------------------------------------------
!
!     estimate of <S^2>=s(s+1) in two different ways.
!     1) using as many-body wavefunction a single Slater determinant
!        constructed with Kohn-Sham orbitals:
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw - 
!                \sum_up\sum_dw < psi_up | psi_dw >
!
!        where Nup, Ndw = number of up and down states, the sum is over 
!        occupied states. Not suitable for fractionary occupancy.
!        In the ultrasoft scheme (c is the smooth part of \psi): 
!
!        < psi_up | psi_dw > = \sum_G c*_up(G) c_dw(G) +
!                              \int Q_ij <c_up|beta_i><beta_j|c_dw>
!
!        This is the usual formula, unsuitable for fractionary occupancy.
!     2) using the "LSD model" of Wang, Becke, Smith, JCP 102, 3477 (1995):
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw -
!                \int max(rhoup(r),rhodw(r)) dr
!
!     Requires on input: c=psi, bec=<c|beta>, rhoup(r), rhodw(r)
!     Assumes real psi, with only half G vectors.
!
      use electrons_base, only: nx => nbspx, n => nbsp, iupdwn, nupdwn, f, nel, nspin
      use io_global, only: stdout
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use grid_dimensions, only: nr1, nr2, nr3, &
            nnr => nnrx
      use cell_base, only: omega
      use cvan, only: nvb, ish
      use uspp, only: nhsa => nkb, nhsavb=>nkbus, qq
      use uspp_param, only: nh
      use ions_base, only: na
!
      implicit none
! input
      real(kind=8) bec(nhsa,n), rhor(nnr,nspin)
      complex(kind=8) c(ngw,nx)
! local variables
      integer nup, ndw, ir, i, j, jj, ig, ia, is, iv, jv, inl, jnl
      real(kind=8) spin0, spin1, spin2, fup, fdw
      real(kind=8), allocatable:: overlap(:,:), temp(:)
      logical frac
!
!
      if (nspin.eq.1) return
!
! find spin-up and spin-down states
!
      fup = 0.0
      do i=iupdwn(1),nupdwn(1)
         fup = fup + f(i)
      end do
      nup = nint(fup)
      ndw = nel(1)+nel(2) - nup
!
! paranoid checks
!
      frac= abs(fup-nup).gt.1.0e-6
      fup = 0.0
      do i=1,nup
         fup = fup + f(i)
      end do
      frac=frac.or.abs(fup-nup).gt.1.0e-6
      fdw = 0.0
      do j=iupdwn(2),iupdwn(2)-1+ndw
         fdw = fdw + f(j)
      end do
      frac=frac.or.abs(fdw-ndw).gt.1.0e-6
!
      spin0 = abs(fup-fdw)/2.d0 * ( abs(fup-fdw)/2.d0 + 1.d0 ) + fdw
!
!     Becke's formula for spin polarization
!
      spin1 = 0.0
      do ir=1,nnr
         spin1 = spin1 - min(rhor(ir,1),rhor(ir,2))
      end do
      call reduce(1,spin1)
      spin1 = spin0 + omega/(nr1*nr2*nr3)*spin1
      if (frac) then
         WRITE( stdout,'(/'' Spin contamination: s(s+1)='',f5.2,'' (Becke) '',&
     &                             f5.2,'' (expected)'')')              &
     &          spin1, abs(fup-fdw)/2.d0*(abs(fup-fdw)/2.d0+1.d0)
         return
      end if
!
!     Slater formula, smooth contribution to  < psi_up | psi_dw >
!
      allocate (overlap(nup,ndw))
      allocate (temp(ngw))
      do j=1,ndw
         jj=j+iupdwn(2)-1
         do i=1,nup
            overlap(i,j)=0.d0
            do ig=1,ngw
               temp(ig)=real(conjg(c(ig,i))*c(ig,jj))
            end do
            overlap(i,j) = 2.d0*SUM(temp)
            if (gstart == 2) overlap(i,j) = overlap(i,j) - temp(1)
         end do
      end do
      deallocate (temp)
      call reduce(nup*ndw,overlap)
      do j=1,ndw
         jj=j+iupdwn(2)-1
         do i=1,nup
!
!     vanderbilt contribution to  < psi_up | psi_dw >
!
            do is=1,nvb
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     if(abs(qq(iv,jv,is)).gt.1.e-5) then 
                        do ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           overlap(i,j) = overlap(i,j) +                &
     &                          qq(iv,jv,is)*bec(inl,i)*bec(jnl,jj)
                        end do
                     endif
                  end do
               end do
            end do
         end do
      end do
!
      spin2 = spin0
      do j=1,ndw
         do i=1,nup
            spin2 = spin2 - overlap(i,j)**2
         end do
      end do
!
      deallocate (overlap)
!
      WRITE( stdout,'(/" Spin contamination: s(s+1)=",f5.2," (Slater) ",  &
     &          f5.2," (Becke) ",f5.2," (expected)")')              &
     &     spin2,spin1, abs(fup-fdw)/2.d0*(abs(fup-fdw)/2.d0+1.d0)
!
      return
      end subroutine spinsq
!-------------------------------------------------------------------------
      subroutine tauset(phi,bephi,qbephi,nss,ist,tau)
!-----------------------------------------------------------------------
!     input: phi
!     computes the matrix
!        tau = <s'c0|s|s'c0> = <phi|s|phi>,  where  |phi> = s'|c0>
!     where s=s(r(t+dt)) and s'=s(r(t))  
!     routine makes use of c(-q)=c*(q)
!
      use parameters, only: nsx, natx
      use cvan, only: nvb
      use uspp, only: nhsa => nkb, nhsavb=>nkbus
      use electrons_base, only: nx => nbspx, n => nbsp
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
!
      implicit none
      integer nss, ist
      complex(kind=8) phi(ngw,n)
      real(kind=8)  bephi(nhsa,n), qbephi(nhsa,n), tau(nx,nx)
      integer i, j
      real(kind=8)    tmp1(nx,nx) ! automatic array
!
      tau = 0.d0
      call MXMA(phi(1,ist),2*ngw,1,phi(1,ist),1,2*ngw,                  &
     &                                   tau,1,nx,nss,2*ngw,nss)
!
!     q >= 0  components with weight 2.0
!
      do j=1,nss
         do i=1,nss
            tau(i,j)=2.*tau(i,j)
         end do
      end do
      if (gstart == 2) then
!
!     q = 0  components has weight 1.0
!
         do j=1,nss
            do i=1,nss
               tau(i,j) = tau(i,j) -                                    &
     &              real(phi(1,i+ist-1))*real(phi(1,j+ist-1))
            end do
         end do
      end if
      call reduce(nx*nss,tau)
!
      if(nvb.gt.0)then
         tmp1 = 0.d0
!
         call MXMA(bephi(1,ist),nhsa,1,qbephi(1,ist),1,nhsa,              &
     &                                     tmp1,1,nx,nss,nhsavb,nss)
!
         do j=1,nss
            do i=1,nss
               tau(i,j)=tau(i,j)+tmp1(i,j)
            end do
         end do
      endif
!
      return
      end subroutine tauset
!
!-------------------------------------------------------------------------
      subroutine updatc(ccc,x0,phi,bephi,becp,bec,cp)
!-----------------------------------------------------------------------
!     input ccc : dt**2/emass (unchanged in output)
!     input x0  : converged lambdas from ortho-loop (unchanged in output)
!     input cp  : non-orthonormal cp=c0+dh/dc*ccc
!     input bec : <cp|beta_i>
!     input phi 
!     output cp : orthonormal cp=cp+lambda*phi
!     output bec: bec=becp+lambda*bephi
!
      use ions_base, only: nsp, na
      use io_global, only: stdout
      use cvan, only: nvb, ish
      use uspp, only: nhsa => nkb, nhsavb=>nkbus
      use uspp_param, only: nh
      use gvecw, only: ngw
      use electrons_base, only: nx => nbspx, n => nbsp
      use control_flags, only: iprint, iprsta
!
      implicit none
!
      complex(kind=8) cp(ngw,n), phi(ngw,n)
      real(kind=8)   bec(nhsa,n), x0(nx,nx), ccc
      real(kind=8)   bephi(nhsa,n), becp(nhsa,n)
! local variables
      integer i, j, ig, is, iv, ia, inl
      real(kind=8) wtemp(n,nhsa) ! automatic array
      complex(kind=8), allocatable :: wrk2(:,:)
!
!     lagrange multipliers
!
      call start_clock( 'updatc' )
      
      allocate( wrk2( ngw, n ) )

      wrk2 = (0.d0, 0.d0)
      do j=1,n
         call DSCAL(n,ccc,x0(1,j),1)
      end do
!
!     wrk2 = sum_m lambda_nm s(r(t+dt))|m>
!
      call MXMA(phi,1,2*ngw,x0,nx,1,wrk2,1,2*ngw,2*ngw,n,n)
!
      do i=1,n
         do ig=1,ngw
            cp(ig,i)=cp(ig,i)+wrk2(ig,i)
         end do
      end do
!    
!     updating of the <beta|c(n,g)>
!
!     bec of vanderbilt species are updated 
!
      if(nvb.gt.0)then
         call MXMA(x0,1,nx,bephi,nhsa,1,wtemp,1,n,n,n,nhsavb)
!
         do i=1,n
            do inl=1,nhsavb
               bec(inl,i)=wtemp(i,inl)+becp(inl,i)
            end do
         end do
      endif
!
      if (iprsta.gt.2) then
         WRITE( stdout,*)
         do is=1,nsp
            if(nsp.gt.1) then
               WRITE( stdout,'(33x,a,i4)') ' updatc: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &            ((bec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
            else
               do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' updatc: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &            ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
               end do
            end if
            WRITE( stdout,*)
         end do
      endif
!
      do j=1,n
         call DSCAL(n,1.0/ccc,x0(1,j),1)
      end do

      deallocate( wrk2 )
!
      call stop_clock( 'updatc' )
!
      return
      end subroutine updatc
!
!-----------------------------------------------------------------------
      subroutine vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,           &
     &     ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!-----------------------------------------------------------------------
!     computes: the one-particle potential v in real space,
!               the total energy etot,
!               the forces fion acting on the ions,
!               the derivative of total energy to cell parameters h
!     rhor input : electronic charge on dense real space grid
!                  (plus core charge if present)
!     rhog input : electronic charge in g space (up to density cutoff)
!     rhos input : electronic charge on smooth real space grid
!     rhor output: total potential on dense real space grid
!     rhos output: total potential on smooth real space grid
!
      use control_flags, only: iprint, tvlocw, iprsta, thdyn, tpre, tfor, tprnfor
      use io_global, only: stdout
      use parameters, only: natx, nsx
      use ions_base, only: nas => nax, nsp, na, nat
      use gvecs
      use gvecp, only: ng => ngm
      use cell_base, only: omega
      use cell_base, only: a1, a2, a3, tpiba2
      use reciprocal_vectors, only: gstart, g
      use recvecs_indexes, only: np, nm
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base, only: nspin
      use constants, only: pi, fpi
      use energies, only: etot, eself, enl, ekin, epseu, esr, eht, exc 
      use local_pseudo, only: vps, rhops
      use core, only: nlcc_any
      use gvecb
      use dener
      use derho
      use mp, only: mp_sum
!
      implicit none
!
      logical tlast,tfirst
      integer nfi
      real(kind=8)  rhor(nnr,nspin), rhos(nnrsx,nspin), fion(3,natx)
      real(kind=8)  rhoc(nnr), tau0(3,natx)
      complex(kind=8) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),     &
     &                ei3(-nr3:nr3,nat), eigrb(ngb,nat),        &
     &                rhog(ng,nspin), sfac(ngs,nsp)
!
      integer irb(3,nat), iss, isup, isdw, ig, ir,i,j,k,is, ia
      real(kind=8) fion1(3,natx), vave, ebac, wz, eh
      complex(kind=8)  fp, fm, ci
      complex(kind=8), allocatable :: v(:), vs(:)
      complex(kind=8), allocatable :: rhotmp(:), vtemp(:), drhotmp(:,:,:)
!
      call start_clock( 'vofrho' )
      ci=(0.,1.)
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz = 2.0
      allocate( v( nnr ) )
      allocate( vs( nnrsx ) )
      allocate(vtemp(ng))
      allocate(rhotmp(ng))
      if (tpre) allocate(drhotmp(ng,3,3))
!
!     first routine in which fion is calculated: annihilation
!
      fion =0.d0
      fion1=0.d0
!
!     ===================================================================
!     forces on ions, ionic term in real space
!     -------------------------------------------------------------------
      if( tprnfor .or. tfor .or. tfirst .or. thdyn ) then
        call force_ion(tau0,esr,fion,dsr)
      end if
!
      if(nspin.eq.1) then
         iss=1
         do ig=1,ng
            rhotmp(ig)=rhog(ig,iss)
         end do
         if(tpre)then
            do j=1,3
               do i=1,3
                  do ig=1,ng
                     drhotmp(ig,i,j)=drhog(ig,iss,i,j)
                  enddo
               enddo
            enddo
         endif
      else
         isup=1
         isdw=2
         do ig=1,ng
            rhotmp(ig)=rhog(ig,isup)+rhog(ig,isdw)
         end do
         if(tpre)then
            do i=1,3
               do j=1,3
                  do ig=1,ng
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
      epseu=wz*real(SUM(vtemp))
      if (gstart == 2) epseu=epseu-vtemp(1)
      call reduce(1,epseu)
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
      if (gstart == 2) vtemp(1)=0.0
      do ig=gstart,ng
         vtemp(ig)=conjg(rhotmp(ig))*rhotmp(ig)/g(ig)
      end do
!
      eh=real(SUM(vtemp))*wz*0.5*fpi/tpiba2
      call reduce(1,eh)
      if(tpre) call denh(rhotmp,drhotmp,sfac,vtemp,eh,dh)
      if(tpre) deallocate(drhotmp)
!     ===================================================================
!     forces on ions, ionic term in reciprocal space
!     -------------------------------------------------------------------
      if( tprnfor .or. tfor .or. thdyn)                                                  &
     &    call force_ps(rhotmp,rhog,vtemp,ei1,ei2,ei3,fion1)
!     ===================================================================
!     calculation hartree + local pseudo potential
!     -------------------------------------------------------------------
!
      if (gstart == 2) vtemp(1)=(0.,0.)
      do ig=gstart,ng
         vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      end do
!
      do is=1,nsp
         do ig=1,ngs
            vtemp(ig)=vtemp(ig)+sfac(ig,is)*vps(ig,is)
         end do
      end do
!
!     vtemp = v_loc(g) + v_h(g)
!
!     ===================================================================
!      calculation exchange and correlation energy and potential
!     -------------------------------------------------------------------
      if (nlcc_any) call add_cc(rhoc,rhog,rhor)
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
         do ig=1,ng
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
         do ig=1,ng
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            rhog(ig,isup)=vtemp(ig)+0.5*cmplx( real(fp),aimag(fm))
            rhog(ig,isdw)=vtemp(ig)+0.5*cmplx(aimag(fp),-real(fm))
         end do
      endif
!
!     rhog contains now the total (local+Hartree+xc) potential in g-space
!
      if( tprnfor .or. tfor ) then

         if (nlcc_any) call force_cc(irb,eigrb,rhor,fion1)

         call mp_sum( fion1 )
!
!    add g-space ionic and core correction contributions to fion
!
         fion = fion + fion1

      end if
!     ===================================================================
!     fourier transform of total potential to r-space (dense grid)
!     -------------------------------------------------------------------
      v(:) = (0.d0, 0.d0)
      if(nspin.eq.1) then
         iss=1
         do ig=1,ng
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
         vave=SUM(rhor(:,iss))/dble(nr1*nr2*nr3)
      else
         isup=1
         isdw=2
         do ig=1,ng
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
         vave=(SUM(rhor(:,isup))+SUM(rhor(:,isdw)))       &
     &        /2.0/dble(nr1*nr2*nr3)
      endif
      call reduce(1,vave)
!     ===================================================================
!     fourier transform of total potential to r-space (smooth grid)
!     -------------------------------------------------------------------
      vs (:) = (0.d0, 0.d0)
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
      deallocate( v )
      deallocate( vs )
!
!
      call stop_clock( 'vofrho' )
      if((nfi.eq.0).or.tfirst.or.tlast) goto 999
      if(mod(nfi-1,iprint).ne.0 ) return
!
 999  WRITE( stdout,1) etot,ekin,eht,esr,eself,epseu,enl,exc,vave
    1 format(//'                total energy = ',f14.5,' a.u.'/         &
     &         '              kinetic energy = ',f14.5,' a.u.'/         &
     &         '        electrostatic energy = ',f14.5,' a.u.'/         &
     &         '                         esr = ',f14.5,' a.u.'/         &
     &         '                       eself = ',f14.5,' a.u.'/         &
     &         '      pseudopotential energy = ',f14.5,' a.u.'/         &
     &         '  n-l pseudopotential energy = ',f14.5,' a.u.'/         &
     &         ' exchange-correlation energy = ',f14.5,' a.u.'/         &
     &         '           average potential = ',f14.5,' a.u.'//)
!
      if(tpre)then
         WRITE( stdout,*) "cell parameters h"
         WRITE( stdout,5555) (a1(i),a2(i),a3(i),i=1,3)
         WRITE( stdout,*)
         WRITE( stdout,*) "derivative of e(tot)"
         WRITE( stdout,5555) ((detot(i,j),j=1,3),i=1,3)
         WRITE( stdout,*)
         if(tpre.and.iprsta.ge.2) then
            WRITE( stdout,*) "derivative of e(kin)"
            WRITE( stdout,5555) ((dekin(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(electrostatic)"
            WRITE( stdout,5555) (((dh(i,j)+dsr(i,j)),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(h)"
            WRITE( stdout,5555) ((dh(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(sr)"
            WRITE( stdout,5555) ((dsr(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(ps)"
            WRITE( stdout,5555) ((dps(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(nl)"
            WRITE( stdout,5555) ((denl(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(xc)"
            WRITE( stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
         endif
      endif
5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
!
      return
      end subroutine vofrho

!
!----------------------------------------------------------------------
      subroutine checkrho(nnr,nspin,rhor,rmin,rmax,rsum,rnegsum)
!----------------------------------------------------------------------
!
!     check \int rho(r)dr and the negative part of rho
!
      implicit none
      integer nnr, nspin
      real(kind=8) rhor(nnr,nspin), rmin, rmax, rsum, rnegsum
!
      real(kind=8) roe
      integer ir, iss
!
      rsum   =0.0
      rnegsum=0.0
      rmin   =100.
      rmax   =0.0
      do iss=1,nspin
         do ir=1,nnr
            roe=rhor(ir,iss)
            rsum=rsum+roe
            if (roe.lt.0.0) rnegsum=rnegsum+roe
            rmax=max(rmax,roe)
            rmin=min(rmin,roe)
         end do
      end do
      call reduce(1,rsum)
      call reduce(1,rnegsum)
      return
      end subroutine checkrho
!______________________________________________________________________


!-----------------------------------------------------------------------
      subroutine vofrho_wf(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,           &
     &     ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!-----------------------------------------------------------------------
!     computes: the one-particle potential v in real space,
!               the total energy etot,
!               the forces fion acting on the ions,
!               the derivative of total energy to cell parameters h
!     rhor input : electronic charge on dense real space grid
!                  (plus core charge if present)
!     rhog input : electronic charge in g space (up to density cutoff)
!     rhos input : electronic charge on smooth real space grid
!     rhor output: total potential on dense real space grid
!     rhos output: total potential on smooth real space grid
!
      use control_flags, only: iprint, tvlocw, iprsta, thdyn, tpre, tfor, tprnfor
      use io_global, only: stdout
      use parameters, only: natx, nsx
      use ions_base, only: nas => nax, nsp, na, nat
      use gvecs
      use cell_base, only: omega, tpiba2
      use cell_base, only: a1, a2, a3, alat
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
            nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base, only: nspin, qbac
      use constants, only: pi, fpi
      use energies, only: etot, eself, enl, ekin, epseu, esr, eht, exc 
      use local_pseudo, only: rhops, vps
      use core, only: nlcc_any
      use gvecb
      use atom, only: nlcc
      use reciprocal_vectors, only: g
      use reciprocal_vectors, only: ng0 => gstart
      use recvecs_indexes, only: np, nm
      use gvecp, only: ng => ngm
!
      use dener
      use derho
!
      implicit none
!
      logical tlast,tfirst
      integer nfi
      real(kind=8)  rhor(nnr,nspin), rhos(nnrsx,nspin), fion(3,natx)
      real(kind=8)  rhoc(nnr), tau0(3,natx)
      complex(kind=8) ei1(-nr1:nr1,nat), ei2(-nr2:nr2,nat),     &
     &                ei3(-nr3:nr3,nat), eigrb(ngb,nat),        &
     &                rhog(ng,nspin), sfac(ngs,nsp)
!
      integer irb(3,nat), iss, isup, isdw, ig, ir,i,j,k,is, ia
      real(kind=8) fion1(3,natx), vave, ebac, wz, eh
      complex(kind=8)  fp, fm, ci
      complex(kind=8), allocatable :: v(:), vs(:)
      complex(kind=8), allocatable:: rhotmp(:), vtemp(:), drhotmp(:,:,:)

! Makov Payne Variables
!
      real(kind=8) dipole,quadrupole
      real(kind=8) E_dip,E_quad,en1,en2
      real(kind=8), allocatable:: rhortot(:)
      real(kind=8) alpha

!
      call start_clock( 'vofrho_wf' )
      
      ci=(0.,1.)
!
!     wz = factor for g.neq.0 because of c*(g)=c(-g)
!
      wz = 2.0
      allocate( v( nnr ) )
      allocate( vs( nnrsx ) )
      allocate(vtemp(ng))
!      write(6,*) 'Allocated vtemp'
      allocate(rhotmp(ng))
!      write(6,*) 'Allocated rhotmp'
      allocate(rhortot(nnr))                ! for Makov Payne
!      write(6,*) 'Allocated rhortot'
      if (tpre) allocate(drhotmp(ng,3,3))
!      write(6,*) 'Allocated all'
!
!     first routine in which fion is calculated: annihilation
!
      fion =0.d0
      fion1=0.d0

!      write(6,*) 'Annihilation'
!
!     ===================================================================
!     forces on ions, ionic term in real space
!     -------------------------------------------------------------------
      if( tprnfor .or. tfor .or. tfirst .or. thdyn ) then
        call force_ion(tau0,esr,fion,dsr)
      end if
!
      if(nspin.eq.1) then
         iss=1
         do ig=1,ng
            rhotmp(ig)=rhog(ig,iss)
         end do
         if(tpre)then
            do j=1,3
               do i=1,3
                  do ig=1,ng
                     drhotmp(ig,i,j)=drhog(ig,iss,i,j)
                  enddo
               enddo
            enddo
         endif
      else
         isup=1
         isdw=2
         do ig=1,ng
            rhotmp(ig)=rhog(ig,isup)+rhog(ig,isdw)
         end do
         if(tpre)then
            do i=1,3
               do j=1,3
                  do ig=1,ng
                     drhotmp(ig,i,j) = drhog(ig,isup,i,j) +           &
     &                                 drhog(ig,isdw,i,j)
                  enddo
               enddo
            enddo
         endif
      end if
!      write(6,*) 'fion'
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
      epseu=wz*real(SUM(vtemp(1:ngs)))
      if (ng0.eq.2) epseu=epseu-vtemp(1)
      call reduce(1,epseu)
      epseu=epseu*omega
!
      if(tpre) call denps(rhotmp,drhotmp,sfac,vtemp,dps)

!      write(6,*) 'Local Energy'
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
      do ig=ng0,ng
         vtemp(ig)=conjg(rhotmp(ig))*rhotmp(ig)/g(ig)
      end do
!
      eh=real(SUM(vtemp(1:ng)))*wz*0.5*fpi/tpiba2
      call reduce(1,eh)
      if(tpre) call denh(rhotmp,drhotmp,sfac,vtemp,eh,dh)
      if(tpre) deallocate(drhotmp)
!      write(6,*) 'Hartree Energy'
!     ===================================================================
!     forces on ions, ionic term in reciprocal space
!     -------------------------------------------------------------------
      if( tprnfor .or. tfor .or. thdyn)                                                  &
     &    call force_ps(rhotmp,rhog,vtemp,ei1,ei2,ei3,fion1)
!     ===================================================================
!     calculation hartree + local pseudo potential
!     -------------------------------------------------------------------
!
      if (ng0.eq.2) vtemp(1)=(0.,0.)
      do ig=ng0,ng
         vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      end do
!
      do is=1,nsp
         do ig=1,ngs
            vtemp(ig)=vtemp(ig)+sfac(ig,is)*vps(ig,is)
         end do
      end do
!
!     vtemp = v_loc(g) + v_h(g)
!  
!       write(6,*) 'Hartree + Local'
! Makov-Payne corrections, by Filippo
!
      if(tlast) then
!     ===================================================================
!     fourier transform of total density to r-space (dense grid)
!     -------------------------------------------------------------------
      v(:) = (0.d0, 0.d0)
         do ig=1,ng
            v(nm(ig))=conjg(rhotmp(ig))
            v(np(ig))=rhotmp(ig)
         end do
!
!     v(g) --> v(r)
!
         call invfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!
         do ir=1,nnr
          rhortot(ir)=real(v(ir))
         end do
!
       call poles(rhortot,dipole,quadrupole)
!
!      Madelung constant for cubic lattice (NaCl)

!
       alpha=1.7476
!
       en1=qbac**2.*alpha/(2.*alat)
       en2=2.*pi*qbac*quadrupole/(3.*alat**3)
!
       write (6,*) "en1: ", en1
       write (6,*) "en2: ", en2
!
       E_quad= en1 + en2
!
!      The interaction energy of the background charge (minus the
!      molecular charge) with itself on a lattice (Madelung energy).
!      +
!      The interaction energy of the background charge with the nuclear
!      quadupole moment on a lattice, with reversed sign due to the fact
!      that the electron density is assumed to be positive.
!
      end if
! END of Makov-Payne corrections, written by Filippo
!
!
!     ===================================================================
!      calculation exchange and correlation energy and potential
!     -------------------------------------------------------------------
      if ( ANY( nlcc ) ) call add_cc(rhoc,rhog,rhor)
!
!      write(6,*) 'add_cc'

      call exch_corr_h(nspin,rhog,rhor,exc,dxc)
!
!     rhor contains the xc potential in r-space

!      write(6,*) 'XC R Space'
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
         do ig=1,ng
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
         do ig=1,ng
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            rhog(ig,isup)=vtemp(ig)+0.5*cmplx( real(fp),aimag(fm))
            rhog(ig,isdw)=vtemp(ig)+0.5*cmplx(aimag(fp),-real(fm))
         end do
      endif
!
!     rhog contains now the total (local+Hartree+xc) potential in g-space
!
!     write(6,*) 'XC G-Space'

      if( tprnfor .or. tfor ) then
         if ( ANY( nlcc ) ) call force_cc(irb,eigrb,rhor,fion1)
         call reduce(3*natx,fion1)
!
!    add g-space ionic and core correction contributions to fion
!
          fion = fion + fion1
      end if
!     ===================================================================
!     fourier transform of total potential to r-space (dense grid)
!     -------------------------------------------------------------------
      v(:) = (0.d0, 0.d0)
      if(nspin.eq.1) then
         iss=1
         do ig=1,ng
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
         vave=SUM(rhor(1:nnr,iss))/dfloat(nr1*nr2*nr3)
      else
         isup=1
         isdw=2
         do ig=1,ng
            v(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            v(nm(ig))=conjg(rhog(ig,isup)) +ci*conjg(rhog(ig,isdw))
         end do
!
         call invfft(v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            rhor(ir,isup)= real(v(ir))
            rhor(ir,isdw)=aimag(v(ir))
         end do

!       write(6,*) 'Average Potential'
!
!     calculation of average potential
!
         vave=(SUM(rhor(1:nnr,isup))+SUM(rhor(1:nnr,isdw)))       &
     &        /2.0/dfloat(nr1*nr2*nr3)
      endif
      call reduce(1,vave)
!     ===================================================================
!     fourier transform of total potential to r-space (smooth grid)
!     -------------------------------------------------------------------
      vs (:) = (0.d0, 0.d0)
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
            rhos(ir,isup)= real(vs(ir))
            rhos(ir,isdw)=aimag(vs(ir))
         end do
      endif


!      write(6,*) 'Total Potential r-space'

      ebac=0.0
!
      eht=eh*omega+esr-eself
!
!     etot is the total energy ; ekin, enl were calculated in rhoofr
!
      etot=ekin+eht+epseu+enl+exc+ebac
      if(tpre) detot=dekin+dh+dps+denl+dxc+dsr

     if(tlast) then
         write (6,*)'MAKOV-PAYNE CORRECTED TOTAL ENERGY',etot+E_quad
         write (6,*)'THIS CORRECTION IS VALID ONLY FOR CUBIC LATTICES'
      end if


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
      deallocate(rhortot)                ! Makov Payne Variable - M.S
      deallocate( v )
      deallocate( vs )

!      write(6,*) 'Deallocations'
!
!
      call stop_clock( 'vofrho_wf' )
      if((nfi.eq.0).or.tfirst.or.tlast) goto 999
      if(mod(nfi-1,iprint).ne.0 ) return
!
 999  WRITE( stdout,1) etot,ekin,eht,esr,eself,epseu,enl,exc,vave
    1 format(//'                total energy = ',f14.5,' a.u.'/         &
     &         '              kinetic energy = ',f14.5,' a.u.'/         &
     &         '        electrostatic energy = ',f14.5,' a.u.'/         &
     &         '                         esr = ',f14.5,' a.u.'/         &
     &         '                       eself = ',f14.5,' a.u.'/         &
     &         '      pseudopotential energy = ',f14.5,' a.u.'/         &
     &         '  n-l pseudopotential energy = ',f14.5,' a.u.'/         &
     &         ' exchange-correlation energy = ',f14.5,' a.u.'/         &
     &         '           average potential = ',f14.5,' a.u.'//)
!
      if(tpre)then
         WRITE( stdout,*) "cell parameters h"
         WRITE( stdout,5555) (a1(i),a2(i),a3(i),i=1,3)
         WRITE( stdout,*)
         WRITE( stdout,*) "derivative of e(tot)"
         WRITE( stdout,5555) ((detot(i,j),j=1,3),i=1,3)
         WRITE( stdout,*)
         if(tpre.and.iprsta.ge.2) then
            WRITE( stdout,*) "derivative of e(kin)"
            WRITE( stdout,5555) ((dekin(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(electrostatic)"
            WRITE( stdout,5555) (((dh(i,j)+dsr(i,j)),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(h)"
            WRITE( stdout,5555) ((dh(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(sr)"
            WRITE( stdout,5555) ((dsr(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(ps)"
            WRITE( stdout,5555) ((dps(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(nl)"
            WRITE( stdout,5555) ((denl(i,j),j=1,3),i=1,3)
            WRITE( stdout,*) "derivative of e(xc)"
            WRITE( stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
         endif
      endif
5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
!
      return
      end subroutine vofrho_wf

!------------------------------------------------------------------------
      subroutine poles(rhortot,dipole,quadrupole)
!------------------------------------------------------------------------
!
      use para_mod
!      use parm
      use grid_dimensions, only : nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr=> nnrx
      use cell_base, only : a1, a2, a3, omega
      use electrons_base, only: qbac
!
      implicit none
      real(kind=8), parameter :: debye=1./0.39344, angs=1./0.52917726
!
      real(kind=8)  dipole,quadrupole,mu(3),quad(6)
      real(kind=8)  ax,ay,az,XG0,YG0,ZG0,X,Y,Z,D,s,rzero,x0,y0,z0
      real(kind=8)  en1,en2, pass1, pass2, pass3
      real(kind=8)  rhortot(nnr)
!     real(kind=8), allocatable:: x(:),y(:),z(:)
      real(kind=8), allocatable:: dip(:)
      integer (kind=4) ix,ir, i, j, k
!
      allocate(dip(nnr))

!     compute the dipole moment
!
        ax=a1(1)
        ay=a2(2)
        az=a3(3)
!
        XG0 = -ax/2.
        YG0 = -ay/2.
        ZG0 = -az/2.
        pass1=ax/nr1
        pass2=ax/nr2
        pass3=ax/nr3
!        pass1 = ax / (nr1-1)
!        pass2 = ay / (nr2-1)
!        pass3 = az / (nr3-1)
!
        do ix=1,3
        ir=1
!
        do k = dfftp%ipp(me)+1, dfftp%ipp(me)+ dfftp%npp(me)
         do j=1,nr2x
          do i=1,nr1x
            X=XG0+(i-1)*pass1
            Y=YG0+(j-1)*pass2
            Z=ZG0+(k-1)*pass3
            if (ix.eq.1) D=X
            if (ix.eq.2) D=Y
            if (ix.eq.3) D=Z
            dip(ir)=D*rhortot(ir)
            ir=ir+1
           end do
          end do
         end do
!
         mu(ix)=sum(dip(1:nnr))
!
         end do !!!!!!! ix
!
         call reduce(3,mu)
!
        do ix=1,3
         mu(ix)=mu(ix)*omega/dfloat(nr1*nr2*nr3)
        end do
!
        dipole=sqrt(mu(1)**2+mu(2)**2+mu(3)**2)
!
!
!       compute the coordinates which put the dipole moment to zero
!
        if (abs(qbac).gt.1.d-05) then
         x0=mu(1)/abs(qbac)
         y0=mu(2)/abs(qbac)
         z0=mu(3)/abs(qbac)
         rzero=x0**2+y0**2+z0**2
        else
         rzero=0.
        end if
!
!       compute the quadrupole moment
!
        do ix=1,6
!
         ir=1
         do k=dfftp%ipp(me)+1, dfftp%ipp(me) + dfftp%npp(me)
          do j=1,nr2x
           do i=1,nr1x
!
            X=XG0+(i-1)*pass1
            Y=YG0+(j-1)*pass2
            Z=ZG0+(k-1)*pass3
!
            if (ix.eq.1) D=X*X
            if (ix.eq.2) D=Y*Y
            if (ix.eq.3) D=Z*Z
            if (ix.eq.4) D=X*Y
            if (ix.eq.5) D=X*Z
            if (ix.eq.6) D=Y*Z
!
            dip(ir)=D*rhortot(ir)
!
            ir=ir+1
           end do
          end do
         end do
!
        quad(ix)=SUM(dip(1:nnr))
        end do
!
         call reduce(6,quad)

        do ix=1,6
         quad(ix)=quad(ix)*omega/dfloat(nr1*nr2*nr3)
        end do
!
        quadrupole=quad(1)+quad(2)+quad(3)-rzero*qbac
!
!  only the diagonal elements contribute to the inetaction energy
!  the term rzero*qbac is subtracted to zero the dipole moment
!
        write (*,1001)(mu(ix),ix=1,3)
        write (*,1002) dipole
        write (*,*) ' '
        write (*,1003)(quad(ix),ix=1,3)
        write (*,1004)(quad(ix),ix=4,6)
        write (*,1005) quadrupole,rzero*qbac
!
1001  format('DIPOLE XYZ-COMPONENTS (A.U.)',f10.4,2x,f10.4,2x,f10.4)
1002  format('DIPOLE MOMENT         (A.U.)',f10.4)
1003  format('QUADRUPOLE XX-YY-ZZ COMPONENTS (A.U.)',             &
     &f9.4,2x,f9.4,2x,f9.4)
1004  format('QUADRUPOLE XY-XZ-YZ COMPONENTS (A.U.)',             &
     &f9.4,2x,f9.4,2x,f9.4)
1005  format('QUADRUPOLE MOMENT              (A.U.)',2f9.4)
!
      deallocate(dip)
!
      return
      end subroutine poles

