!
! Copyright (C) 2002-2004 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
     subroutine caldbec(nspmn,nspmx,eigr,c)
!-----------------------------------------------------------------------
!     this routine calculates array dbec, derivative of bec:
!
!        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i) 
!
!     with respect to cell parameters h
!
!     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
!
      use ions_base, only: na, nas => nax
      use elct, only: n
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use constants, only: pi, fpi
      use cvan, only: ish
      use uspp, only: nhtol, nhsa => nkb
      use uspp_param, only: nh, nhm
      use cdvan
      use work, only: wrk2
!
      implicit none
      integer nspmn, nspmx
      complex(kind=8) c(ngw,n)
      real(kind=8) eigr(2,ngw,nas,nspmx)
!
      integer ig, is, iv, ia, l, ixr, ixi, inl, i, j, ii
      real(kind=8) signre, signim, arg
!
!
      do j=1,3
         do i=1,3

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
                        wrk2(1,ia)= cmplx(                             &
     &                  signre*dbeta(1,iv,is,i,j)*eigr(ixr,1,ia,is),   &
     &                  signim*dbeta(1,iv,is,i,j)*eigr(ixi,1,ia,is) )
!                   q > 0   components (with weight 2.0)
                     end if
                     do ig=gstart,ngw
                        arg = 2.0*dbeta(ig,iv,is,i,j)
                        wrk2(ig,ia) = cmplx(                           &
     &                         signre*arg*eigr(ixr,ig,ia,is),          &
     &                         signim*arg*eigr(ixi,ig,ia,is) )
                     end do
                  end do
                  inl=ish(is)+(iv-1)*na(is)+1
                  call MXMA(wrk2,2*ngw,1,c,1,2*ngw,dbec(inl,1,i,j),1,  &
     &                      nhsa,na(is),2*ngw,n)
               end do
#ifdef __PARA
               inl=ish(is)+1
               do ii=1,n
                  call reduce(na(is)*nh(is),dbec(inl,ii,i,j))
               end do
#endif
            end do
         end do
      end do
!
      return
      end
!
!---------------------------------------------------------------------
      subroutine exch_corr_h(nspin,rhog,rhor,exc,dxc)
!---------------------------------------------------------------------
!
! calculate exch-corr potential, energy, and derivatives dxc(i,j)
! of e(xc) with respect to to cell parameter h(i,j)
!
      use funct, only: iexch, icorr, igcx, igcc
      use gvecp, only: ng => ngm
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx
      use cell_base, only: ainv
      !use cell_module
      use cell_base, only: omega
      use control_flags, only: tpre
      use derho
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
      exc=0.0
      if (iexch==1.and.icorr==1.and.igcx==0.and.igcc==0) then
         ! LDA (Perdew-Zunger)
         call expxc(nnr,nspin,rhor,exc)
      else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
         ! PW91
         call ggapwold(nspin,rhog,gradr,rhor,exc)
      else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
         ! BLYP
         call ggablyp4(nspin,rhog,gradr,rhor,exc)
      else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
         ! PBE
         call ggapbe(nspin,rhog,gradr,rhor,exc)
      else
         call errore('exc-cor','no such exch-corr',1)
      end if
      exc=exc*omega/dble(nr1*nr2*nr3)
!
! exchange-correlation contribution to pressure
!
      dxc(:,:) = 0.0
      if (tpre) then
         if (nspin.ne.1) call errore('exc-cor','spin not implemented',1)
!
         do j=1,3
            do i=1,3
               do ir=1,nnr
                  dxc(i,j) = dxc(i,j) + rhor(ir,1)*drhor(ir,1,i,j)
               end do
               dxc(i,j)=omega/(nr1*nr2*nr3)*dxc(i,j)
            end do
         end do
#ifdef __PARA
         call reduce (9,dxc)
#endif
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
         call gradh(nspin,gradr,rhog,rhor,dexc)
         if (tpre) then
#ifdef __PARA
            call reduce (9,dexc)
#endif
            do j=1,3
               do i=1,3
                  dxc(i,j) = dxc(i,j) + dexc(i,j)
               end do
            end do
         end if
         deallocate(gradr)
      end if
!
      return
      end
!-----------------------------------------------------------------------
      subroutine formf(tfirst, eself)
!-----------------------------------------------------------------------
!computes (a) the self-energy eself of the ionic pseudocharges;
!         (b) the form factors of: (i) pseudopotential (vps),
!             (ii) ionic pseudocharge (rhops)
!         all quantities are returned in common /pseu/
!         also calculated the derivative of vps with respect to
!         g^2 (dvps)
! 
      use control_flags, only: iprint, tpre, iprsta
      use io_global, only: stdout
      use bhs
      use gvec, only: tpiba, tpiba2, g
      use gvecp, only: ng => ngm ! , ngl => ngml, ng_g => ngmt
      use gvecs
      use cell_base, only: omega
      use constants, only: pi, fpi
      use ions_base, only: rcmax,  zv, nsp, na
      use cvan, only: ipp
      use pseu
      use reciprocal_vectors, only: gstart
      use atom, only: r, rab, mesh
      use uspp_param, only: vloc_at
      use qrl_mod, only: cmesh
!
      use dpseu
      use dener
!
      implicit none
      logical tfirst
      real(kind=8) :: eself
!
      real(kind=8), allocatable:: f(:),vscr(:), figl(:)
      real(kind=8) el, ql, par, sp, e1, e2, emax, vpsum, rhopsum, fint, &
     &             fpibg, gps, sfp, xg, dsfp, dgps, r2new, r2max, r21,  &
     &             r22, r2l
      real(kind=8), external :: erf
      integer is, irmax, ir, ig, ib, ndm
      real(kind=8), allocatable:: df(:), dfigl(:)
!
!     ==================================================================
!     calculation of gaussian selfinteraction
!     ==================================================================
      call start_clock( 'formf' )
      eself=0.
      do is=1,nsp
         eself=eself+float(na(is))*zv(is)*zv(is)/rcmax(is)
      end do
      eself=eself/sqrt(2.*pi)
      if(tfirst.or.iprsta.ge.4)then
         WRITE( stdout,1200) eself
      endif
 1200 format(2x,'formf: eself=',f10.5)
!
      ndm = MAXVAL (mesh(1:nsp))
      allocate(figl(ngs))
      allocate(f(ndm))
      allocate(vscr(ndm))
      if (tpre) then
         allocate(dfigl(ngs))
         allocate(df(ndm))
      end if
!
!     ==================================================================
!     fourier transform of local pp and gaussian nuclear charge
!     ==================================================================
      do is=1,nsp
         if(ipp(is).ne.3) then
!!!        if (numeric(is)) then
!     ==================================================================
!     local potential given numerically on logarithmic mesh 
!     ==================================================================
!
!     ------------------------------------------------------------------
!     g=0
!     ------------------------------------------------------------------
!
!     definition of irmax: gridpoint beyond which potential is zero
!
            irmax=0
            do ir=1,mesh(is)
               if(r(ir,is).le.10.0)then
                  irmax=ir
               endif
            end do
!
            do ir=1,irmax
               vscr(ir)=0.5*r(ir,is)*vloc_at(ir,is) +                   &
     &                  zv(is)*erf(r(ir,is)/rcmax(is))
               f(ir)=vscr(ir)*r(ir,is)
            end do
            do ir=irmax+1,mesh(is)
               vscr(ir)=0.0
               f(ir)=0.0
            end do
!!!         if (oldpseudo(is)) then
            if (ipp(is).eq.0) then
               call herman_skillman_int(mesh(is),cmesh(is),f,fint)
            else
               call simpson_cp90(mesh(is),f,rab(1,is),fint)
            end if
!
            if (gstart == 2) then
               vps(1,is)=fpi*fint/omega
               rhops(1,is)=-zv(is)/omega
               vpsum=vps(1,is)
               rhopsum=rhops(1,is)
            else
               vpsum=0.0
               rhopsum=0.0
            end if
            r2new=0.25*tpiba2*rcmax(is)**2
!
!     ------------------------------------------------------------------
!     g>0
!     ------------------------------------------------------------------
            do ig=gstart,ngs
               xg=sqrt(g(ig))*tpiba
               do ir=1,mesh(is)
                  f(ir)=vscr(ir)*sin(r(ir,is)*xg)
                  if(tpre) then
                     df(ir)=vscr(ir)*cos(r(ir,is)*xg)*.5*r(ir,is)/xg
                  endif
               end do
!!!         if (oldpseudo(is)) then
               if (ipp(is).eq.0) then
                  call herman_skillman_int                              &
     &                 (mesh(is),cmesh(is),f,figl(ig))
                  if(tpre) call herman_skillman_int                     &
     &                          (mesh(is),cmesh(is),df,dfigl(ig))
               else
                  call simpson_cp90(mesh(is),f,rab(1,is),figl(ig))
                  if(tpre) call simpson_cp90(mesh(is),df,rab(1,is),dfigl(ig))
               end if
            end do
!
            do ig=gstart,ngs
               xg=sqrt(g(ig))*tpiba
               rhops(ig,is)=-zv(is)*exp(-r2new*g(ig))/omega
               vps(ig,is)=fpi*figl(ig)/(omega*xg)
               if(tpre)then
                  drhops(ig,is)=-rhops(ig,is)*r2new/tpiba2
                  dvps(ig,is)=fpi*dfigl(ig)/(omega*xg)-                 &
     &                 0.5*vps(ig,is)/(xg*xg)
               endif
               rhopsum=rhopsum+rhops(ig,is)
               vpsum=vpsum+vps(ig,is)
            end do
!
         else
!     ==================================================================
!     bhs pseudopotentials can be fourier transformed analytically
!     ==================================================================
            r2new=0.25*tpiba2*rcmax(is)**2
            r2max=rcmax(is)**2
            r21=rc1(is)**2
            r22=rc2(is)**2
!
!     ------------------------------------------------------------------
!     g=0
!     ------------------------------------------------------------------
            if (gstart == 2) then
               rhops(1,is)=-zv(is)/omega
               gps=-zv(is)*pi*(-wrc2(is)*r22-wrc1(is)*r21+r2max)/omega
               sfp=0.
               do ib=1,3
                  r2l=rcl(ib,is,lloc(is))**2
                  ql=0.25*r2l*g(1)*tpiba2
                  el=exp(-ql)
                  par=al(ib,is,lloc(is))+bl(ib,is,lloc(is))*r2l*(1.5-ql)
                  sp=(pi*r2l)**1.5*el/omega
                  sfp=sp*par+sfp
               end do
               vps(1,is)=gps+sfp
               vpsum=vps(1,is)
               rhopsum=rhops(1,is)
            else
               vpsum=0.0
               rhopsum=0.0
            end if
!
!     ------------------------------------------------------------------
!     g>0
!     ------------------------------------------------------------------
            do ig=gstart,ngs
               rhops(ig,is)=-zv(is)*exp(-r2new*g(ig))/omega
               if(tpre) drhops(ig,is)=-rhops(ig,is)*r2new/tpiba2
               emax=exp(-0.25*r2max*g(ig)*tpiba2)
               e1=exp(-0.25*r21*g(ig)*tpiba2)
               e2=exp(-0.25*r22*g(ig)*tpiba2)
               fpibg=fpi/(g(ig)*tpiba2)
               gps=-zv(is)*(wrc1(is)*e1-emax+wrc2(is)*e2)/omega
               gps=gps*fpibg
               if(tpre) dgps=-gps/(tpiba2*g(ig)) +                      &
     &                       fpibg*zv(is)*(wrc1(is)*r21*e1-             &
     &                       r2max*emax+wrc2(is)*r22*e2)*0.25/omega
               sfp=0.
               dsfp=0.
               do ib=1,3
                  r2l=rcl(ib,is,lloc(is))**2
                  ql=0.25*r2l*g(ig)*tpiba2
                  par=al(ib,is,lloc(is))+bl(ib,is,lloc(is))*r2l*(1.5-ql)
                  sp=(pi*r2l)**1.5*exp(-ql)/omega
                  sfp=sp*par+sfp
                  if(tpre) dsfp = dsfp -                                &
     &                 sp*(par+bl(ib,is,lloc(is))*r2l)*ql/(tpiba2*g(ig))
               end do
               vps(ig,is)=sfp+gps
               if(tpre) dvps(ig,is)=dsfp+dgps
               rhopsum=rhopsum+rhops(ig,is)
               vpsum=vpsum+vps(ig,is)
            end do
! 
         endif
!
         if(tfirst.or.(iprsta.ge.4))then
#ifdef __PARA
            call reduce(1,vpsum)
            call reduce(1,rhopsum)
#endif
            WRITE( stdout,1250) vps(1,is),rhops(1,is)
            WRITE( stdout,1300) vpsum,rhopsum
         endif
!
      end do
!
      if (tpre) then
         deallocate(df)
         deallocate(dfigl)
      end if
      deallocate(vscr)
      deallocate(f)
      deallocate(figl)
      call stop_clock( 'formf' )
!
 1250 format(2x,'formf:     vps(g=0)=',f12.7,'     rhops(g=0)=',f12.7)
 1300 format(2x,'formf: sum_g vps(g)=',f12.7,' sum_g rhops(g)=',f12.7)
!
      return
      end
!______________________________________________________________________
      subroutine gradh(nspin,gradr,rhog,rhor,dexc)
!     _________________________________________________________________
!
!     calculate the second part of gradient corrected xc potential
!     plus the gradient-correction contribution to pressure
!
      use control_flags, only: iprint, tpre
      use gvec
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx, &
            nr1x, nr2x, nr3x
      use cell_base, only: ainv
      use work, only: wrk1
      !use cell_module
      use cell_base, only: omega
      use derho
!
      implicit none
! input
      integer nspin
      real(kind=8)  gradr(nnr,3,nspin), rhor(nnr,nspin), dexc(3,3)
      complex(kind=8) rhog(ng,nspin)
!
      complex(kind=8), pointer:: v(:)
      complex(kind=8), allocatable:: x(:), vtemp(:)
      complex(kind=8)  ci, fp, fm
      integer iss, ig, ir, i,j
!
      v => wrk1
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
!
      return
      end
!-----------------------------------------------------------------------
      subroutine init (ibrav,celldm, ecut, ecutw,ndr,nbeg,  &
                       tfirst,tau0,taus,delt,tps,iforce)
!-----------------------------------------------------------------------
!
!     initialize G-vectors and related quantities
!     use ibrav=0 for generic cell vectors given by the matrix h(3,3)
!
      use control_flags, only: iprint, thdyn
      use io_global, only: stdout
      use gvecw, only: ngw
      use ions_base, only: na, pmass, nsp
      use cell_base, only: ainv, a1, a2, a3, r_to_s, s_to_r
      use constants, only: pi, fpi
      use cell_base, only: hold, h
      use gvecw, only: agg => ecutz, sgg => ecsig, e0gg => ecfix
      use betax, only: mmx, refg
      use restart, only: readfile
      use parameters, only: nacx, nsx, natx

      implicit none
! input/output
      integer ibrav, ndr, nbeg
      logical tfirst
      real(kind=8) tau0(3,natx), taus(3,natx)
      integer iforce(3,natx)
      real(kind=8) celldm(6), ecut, ecutw
      real(kind=8) delt
! local
      real(kind=8) randy
      integer i, j, ia, is, nfi, isa, isat
! present in the call to read(p)file, not actually used
      complex(kind=8) c0(1,1,1,1),cm(1,1,1,1)
      real(kind=8) taum(1,1),vel(1,1),velm(1,1),acc(nacx)
      real(kind=8) lambda(1,1),lambdam(1,1)
      real(kind=8) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp, ekincm
      real(kind=8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      real(kind=8) fion(1,1), tps
!
!
! taus = scaled, tau0 = alat units
!
      CALL r_to_s( tau0, taus, na, nsp, ainv )
!
      refg = 1.0d0 * ecut / ( mmx - 1 )
      WRITE( stdout,*) '   NOTA BENE: refg, mmx = ',refg,mmx
!
      if( nbeg >= 0 ) then
!
! read only h and hold from file ndr
!
         call readfile                                              &
     &     (-1,ndr,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),tau0,taum,vel,velm,acc,             &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion, tps)
!
         WRITE( stdout,344) ibrav
         do i=1,3
            WRITE( stdout,345) (h(i,j),j=1,3)
         enddo
         WRITE( stdout,*)

      else
!
! with variable-cell we use h to describe the cell
!
         do i = 1, 3
            h(i,1) = a1(i)
            h(i,2) = a2(i)
            h(i,3) = a3(i)
         enddo

         hold = h

      end if
!
!     ==============================================================
!     ==== generate true g-space                                ==== 
!     ==============================================================
!
      call newinit( ibrav )
!
      !
 344  format(' ibrav = ',i4,'       cell parameters ',/)
 345  format(3(4x,f10.5))
      return
      end
!
!-----------------------------------------------------------------------
      subroutine newinit(ibrav)
!-----------------------------------------------------------------------
!     re-initialization of lattice parameters and g-space vectors.
!     Note that direct and reciprocal lattice primitive vectors
!     a1,a2,a3, ainv, and corresponding quantities for small boxes
!     are recalculated according to the value of cell parameter h
!
      use control_flags, only: iprint, iprsta
      use io_global, only: stdout
      use gvec
      use grid_dimensions, only: nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3
      use cell_base, only: omega, alat
      use constants, only: pi, fpi
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b
      use small_box, only: a1b, a2b, a3b, ainvb, omegab, tpibab
      use cell_base, only: h, deth
      use gvecw, only: agg => ecutz, sgg => ecsig, e0gg => ecfix
!
      implicit none
      integer ibrav
!
! local
      integer i, j
      real(kind=8) alatb, gmax, b1(3),b2(3),b3(3), b1b(3),b2b(3),b3b(3)
      real(kind=8) ddum
!
!
      alat = sqrt( h(1,1)*h(1,1) + h(2,1)*h(2,1) + h(3,1)*h(3,1) )

!     ==============================================================
      tpiba  = 2.d0 * pi / alat
      tpiba2 = tpiba * tpiba

!     ==============================================================
!     ==== generate g-space                                     ==== 
!     ==============================================================
      call invmat (3, h, ainv, deth)
      omega = deth
!
      do i = 1, 3
         a1(i) = h(i,1)
         a2(i) = h(i,2)
         a3(i) = h(i,3)
      enddo
!
      call recips( a1, a2, a3, b1, b2, b3 )
      b1 = b1 * alat
      b2 = b2 * alat
      b3 = b3 * alat
      call gcal( b1, b2, b3, gmax )
!
!     ==============================================================
!     generation of little box g-vectors
!     ==============================================================
!
      call newgb( a1, a2, a3, omega, alat )

!     ==============================================================
      if(iprsta.ge.4)then
         WRITE( stdout,34) ibrav,alat,omega
         if(ibrav.eq.0) then
            WRITE( stdout,344)
            do i=1,3
               WRITE( stdout,345) (h(i,j),j=1,3)
            enddo
            WRITE( stdout,*)
         endif
      endif
! 
 34   format(' initialization ',//,                                     &
     &       ' ibrav=',i3,' alat=',f7.3,' omega=',f10.4,//)
 344  format(' cell parameters ',/)
 345  format(3(4x,f10.5))
!
      return
      end
!-----------------------------------------------------------------------
      subroutine newnlinit
!-----------------------------------------------------------------------
!     this routine calculates arrays beta, qradb, qq, qgb, rhocb
!     and derivatives w.r.t. cell parameters dbeta, dqrad 
!     See also comments in nlinit
!
      use control_flags, only: iprint, tpre, iprsta
      use io_global, only: stdout
      !use gvec
      use gvecw, only: ngw
      use reciprocal_vectors, only: g, gx
      use reciprocal_vectors, only: gstart
      use cell_base, only: omega, tpiba2, tpiba
      use cell_base, only: ainv
      use cvan, only: nvb
      use uspp, only: qq, nhtolm, beta
      use uspp_param, only: nh
      use core
      use constants, only: pi, fpi
      use ions_base, only: nsp
      use elct
      use uspp_param, only: lmaxq, nqlc, lmaxkb, kkbeta, nbrx, nbeta
      use atom, only: nlcc, r, rab, mesh, rho_atc
      use qradb_mod
      use qgb_mod
      use gvecb
      use small_box,  only: omegab, tpibab
      use cdvan
      use dqrad_mod
      use dqgb_mod
      use betax
!
      implicit none
      integer  is, l, lp, ig, ir, iv, jv, ijv, i,j, jj
      real(kind=8), allocatable:: fint(:), jl(:), dqradb(:,:,:,:,:)
      real(kind=8), allocatable:: ylmb(:,:), ylm(:,:), &
                                  dylmb(:,:,:,:), dylm(:,:,:,:)
      complex(kind=8), allocatable:: dqgbs(:,:,:)
      real(kind=8) xg, c, betagl, dbetagl, gg
!
! 
      allocate(ylmb(ngb,lmaxq*lmaxq))
!
      qradb(:,:,:,:,:) = 0.d0
      call ylmr2 (lmaxq*lmaxq, ngb, gxb, gb, ylmb)

!     ===============================================================
!     initialization for vanderbilt species
!     ===============================================================
      do is=1,nvb
!     ---------------------------------------------------------------
!     calculation of array qradb(igb,iv,jv,is)
!     ---------------------------------------------------------------
         if(iprsta.ge.4) WRITE( stdout,*)  '  qradb  '
         c=fpi/omegab
!
         do l=1,nqlc(is)
            do iv= 1,nbeta(is)
               do jv=iv,nbeta(is)
                  do ig=1,ngb
                     gg=gb(ig)*tpibab*tpibab/refg
                     jj=int(gg)+1
                     if(jj.ge.mmx) then
                        qradb(ig,iv,jv,l,is)=0.
                        qradb(ig,jv,iv,l,is)=qradb(ig,iv,jv,l,is)
                     else
                        qradb(ig,iv,jv,l,is)=                           &
     &                       c*qradx(jj+1,iv,jv,l,is)*(gg-real(jj-1))+  &
     &                       c*qradx(jj,iv,jv,l,is)*(real(jj)-gg)
                        qradb(ig,jv,iv,l,is)=qradb(ig,iv,jv,l,is)
                     endif
                  enddo
               enddo
            enddo
         enddo
!
!     ---------------------------------------------------------------
!     stocking of qgb(igb,ijv,is) and of qq(iv,jv,is)
!     ---------------------------------------------------------------
         ijv=0
         do iv= 1,nh(is)
            do jv=iv,nh(is)
!
!       compact indices because qgb is symmetric:
!       ivjv:  11 12 13 ... 22 23...
!       ijv :   1  2  3 ...  
!
               ijv=ijv+1
               call qvan2b(ngb,iv,jv,is,ylmb,qgb(1,ijv,is) )
!
               qq(iv,jv,is)=omegab*real(qgb(1,ijv,is))
               qq(jv,iv,is)=qq(iv,jv,is)
!
            end do
         end do
      end do
!
      if (tpre) then
!     ---------------------------------------------------------------
!     arrays required for stress calculation, variable-cell dynamics
!     ---------------------------------------------------------------
         allocate(dqradb(ngb,nbrx,nbrx,lmaxq,nsp))
         allocate(dylmb(ngb,lmaxq*lmaxq,3,3))
         allocate(dqgbs(ngb,3,3))
         dqrad(:,:,:,:,:,:,:) = 0.d0
         !
         call dylmr2_(lmaxq*lmaxq, ngb, gxb, gb, ainv, dylmb)
         !
         do is=1,nvb
            !
            do l=1,nqlc(is)
               do iv= 1,nbeta(is)
                  do jv=iv,nbeta(is) 
                    do ig=1,ngb
                        gg=gb(ig)*tpibab*tpibab/refg
                        jj=int(gg)+1
                        if(jj.ge.mmx) then
                           dqradb(ig,iv,jv,l,is) = 0.
                        else
                           dqradb(ig,iv,jv,l,is) =                      &
     &                       dqradx(jj+1,iv,jv,l,is)*(gg-real(jj-1))+   &
     &                       dqradx(jj,iv,jv,l,is)*(real(jj)-gg)
                        endif
                     enddo
                     do i=1,3
                        do j=1,3
                           dqrad(1,iv,jv,l,is,i,j) = &
                                -qradb(1,iv,jv,l,is) * ainv(j,i)
                           dqrad(1,jv,iv,l,is,i,j) = &
                                dqrad(1,iv,jv,l,is,i,j)
                           do ig=2,ngb
                              dqrad(ig,iv,jv,l,is,i,j) =                &
     &                          -qradb(ig,iv,jv,l,is)*ainv(j,i)         &
     &                          -c*dqradb(ig,iv,jv,l,is)*               &
     &                          gxb(i,ig)/gb(ig)*                       &
     &                          (gxb(1,ig)*ainv(j,1)+                   &
     &                           gxb(2,ig)*ainv(j,2)+                   &
     &                           gxb(3,ig)*ainv(j,3)) 
                              dqrad(ig,jv,iv,l,is,i,j) =                &
     &                          dqrad(ig,iv,jv,l,is,i,j)
                           enddo
                        enddo
                     enddo
                  end do
               enddo
            enddo
            !
            ijv=0
            !
            do iv= 1,nh(is)
               do jv=iv,nh(is)
                  !
                  !       compact indices because qgb is symmetric:
                  !       ivjv:  11 12 13 ... 22 23...
                  !       ijv :   1  2  3 ...  
                  !
                  ijv=ijv+1
                  call dqvan2b(ngb,iv,jv,is,ylmb,dylmb,dqgbs )
                  do i=1,3
                     do j=1,3
                        do ig=1,ngb
                           dqgb(ig,ijv,is,i,j)=dqgbs(ig,i,j)
                        enddo
                     enddo
                  enddo
               end do
            end do
         end do
         deallocate(dqgbs)
         deallocate(dylmb)
         deallocate(dqradb)
      end if
      deallocate(ylmb)
!
!     ===============================================================
!     initialization that is common to all species
!     ===============================================================
!
      allocate(ylm(ngw,(lmaxkb+1)**2))
      call ylmr2 ((lmaxkb+1)**2, ngw, gx, g, ylm)
!
      do is=1,nsp
!     ---------------------------------------------------------------
!     calculation of array  beta(ig,iv,is)
!     ---------------------------------------------------------------
         if(iprsta.ge.4) WRITE( stdout,*)  '  beta  '
         c=fpi/sqrt(omega)
         do iv=1,nh(is)
            lp=nhtolm(iv,is)
            do ig=1,ngw
               gg=g(ig)*tpiba*tpiba/refg
               jj=int(gg)+1
               betagl=betagx(jj+1,iv,is)*(gg-real(jj-1))+               &
     &              betagx(jj,iv,is)*(real(jj)-gg)
               beta(ig,iv,is)=c*ylm(ig,lp)*betagl
            end do
         end do
      end do
!
      if (tpre) then
!     ---------------------------------------------------------------
!     calculation of array dbeta required for stress, variable-cell
!     ---------------------------------------------------------------
         allocate(dylm(ngw,(lmaxkb+1)**2,3,3))
         !
         call dylmr2_((lmaxkb+1)**2, ngw, gx, g, ainv, dylm)
         !
         do is=1,nsp
            if(iprsta.ge.4) WRITE( stdout,*)  '  dbeta  '
            c=fpi/sqrt(omega)
            do iv=1,nh(is)
               lp=nhtolm(iv,is)
               betagl=betagx(1,iv,is)
               do i=1,3
                  do j=1,3
                     dbeta(1,iv,is,i,j)=-0.5*beta(1,iv,is)*ainv(j,i)    &
     &                                 +c*dylm(1,lp,i,j)*betagl
                  enddo
               enddo
               do ig=gstart,ngw
                  gg=g(ig)*tpiba*tpiba/refg
                  jj=int(gg)+1
                  betagl = betagx(jj+1,iv,is)*(gg-real(jj-1)) +         &
     &                     betagx(jj,iv,is)*(real(jj)-gg)
                  dbetagl= dbetagx(jj+1,iv,is)*(gg-real(jj-1)) +        &
     &                     dbetagx(jj,iv,is)*(real(jj)-gg)
                  do i=1,3
                     do j=1,3
                        dbeta(ig,iv,is,i,j)=                            &
     &                    -0.5*beta(ig,iv,is)*ainv(j,i)                 &
     &                    +c*dylm(ig,lp,i,j)*betagl                     &
     &                    -c*ylm (ig,lp)*dbetagl*gx(i,ig)/g(ig)         &
     &                    *(gx(1,ig)*ainv(j,1)+                         &
     &                      gx(2,ig)*ainv(j,2)+                         &
     &                      gx(3,ig)*ainv(j,3))
                     end do
                  end do
               end do
            end do
         end do
         !
         deallocate(dylm)
      end if
      !
      deallocate(ylm)
!     ---------------------------------------------------------------
!     non-linear core-correction   ( rhocb(ig,is) )
!     ---------------------------------------------------------------
      do is=1,nsp
         if(nlcc(is)) then
            allocate(fint(kkbeta(is)))
            allocate(jl(kkbeta(is)))
            c=fpi/omegab
            l=1
            do ig=1,ngb
               xg=sqrt(gb(ig))*tpibab
               call sph_bes (kkbeta(is), r(1,is), xg, l-1, jl)
               do ir=1,kkbeta(is)
                  fint(ir)=r(ir,is)**2*rho_atc(ir,is)*jl(ir)
               end do
               call simpson_cp90(kkbeta(is),fint,rab(1,is),rhocb(ig,is))
            end do
            do ig=1,ngb
               rhocb(ig,is)=c*rhocb(ig,is)
            end do
            if(iprsta.ge.4) WRITE( stdout,'(a,f12.8)')                     &
      &              ' integrated core charge= ',omegab*rhocb(1,is)
            deallocate(jl)
            deallocate(fint)
         endif
      end do
!
!
      return
      end
!-----------------------------------------------------------------------
      subroutine nlfh(bec,dbec,lambda)
!-----------------------------------------------------------------------
!     contribution to the internal stress tensor due to the constraints
!
      use gvec
      use cvan, only: nvb, ish
      use uspp, only: nhsa => nkb, qq
      use uspp_param, only: nh, nhm
      use ions_base, only: na
      use elct
      use cell_base, only: omega, h
      use constants, only: pi, fpi
      use stre
!
      implicit none
      real(kind=8) bec(nhsa,n), dbec(nhsa,n,3,3), lambda(nx,nx)
!
      integer i, j, ii, jj, inl, iv, jv, ia, is
      real(kind=8) fpre(3,3), tmpbec(nhm,nx), tmpdh(nx,nhm), temp(nx,nx)
!
      fpre(:,:) = 0.d0
      do ii=1,3
         do jj=1,3
            do is=1,nvb
               do ia=1,na(is)
!
                  tmpbec(:, 1:n) = 0.d0
                  tmpdh (1:n, :) = 0.d0
!
                  do iv=1,nh(is)
                     do jv=1,nh(is)
                        inl=ish(is)+(jv-1)*na(is)+ia
                        if(abs(qq(iv,jv,is)).gt.1.e-5) then
                           do i=1,n
                              tmpbec(iv,i) = tmpbec(iv,i) +             &
     &                             qq(iv,jv,is)*bec(inl,i)
                           end do
                        endif
                     end do
                  end do
!
                  do iv=1,nh(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        tmpdh(i,iv)=dbec(inl,i,ii,jj)
                     end do
                  end do
!
                  if(nh(is).gt.0)then
                     temp(:, 1:n) = 0.d0
!
                     call MXMA                                          &
     &                    (tmpdh,1,nx,tmpbec,1,nhm,temp,1,nx,n,nh(is),n)
!
                     do j=1,n
                        do i=1,n
                           temp(i,j)=temp(i,j)*lambda(i,j)
                        end do
                     end do
!
                     fpre(ii,jj)=fpre(ii,jj)+2.*SUM(temp(1:n,1:n))
                  endif
!
               end do
            end do
         end do
      end do
      do i=1,3
         do j=1,3
            stress(i,j)=stress(i,j)+(fpre(i,1)*h(j,1)+                  &
     &           fpre(i,2)*h(j,2)+fpre(i,3)*h(j,3))/omega
         enddo
      enddo
!
      return
      end
!-----------------------------------------------------------------------
      subroutine nlinit
!-----------------------------------------------------------------------
!
!     this routine allocates and initalizes arrays beta, qradb, qq, qgb,
!     rhocb, and derivatives w.r.t. cell parameters dbeta, dqrad 
!
!       beta(ig,l,is) = 4pi/sqrt(omega) y^r(l,q^)
!                               int_0^inf dr r^2 j_l(qr) betar(l,is,r)
!
!       Note that beta(g)_lm,is = (-i)^l*beta(ig,l,is) (?)
!
!       qradb(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
!
!       qq_ij=int_0^r q_ij(r)=omega*qg(g=0)
!
!     beta and qradb are first calculated on a fixed linear grid in |G|
!     (betax, qradx) then calculated on the box grid by interpolation
!     (this is done in routine newnlinit)
!     
      use parameters, only: lmaxx
      use control_flags, only: iprint, tpre
      use io_global, only: stdout
      use gvec
      use gvecw, only: ngw
      use cvan, only: ish, nvb, ipp
      use core
      use constants, only: pi, fpi
      use ions_base, only: na, nsp
      use elct
      use uspp, only: aainit, beta, qq, dvan, nhtol, nhtolm, indv, &
           nhsa => nkb, nhsavb=>nkbus
      use uspp_param, only: kkbeta, qqq, nqlc, betar, nbrx, lmaxq, dion, &
           nbeta, lmaxkb, lll, nhm, nh
      use qrl_mod, only: qrl, cmesh
      use atom, only: mesh, r, rab, nlcc
      use qradb_mod
      use qgb_mod
      use gvecb
      use cdvan
      use dqrad_mod
      use dqgb_mod
      use betax
!
      implicit none
!
      integer  is, il, l, ir, iv, jv, lm, ind, ltmp, i0
      real(kind=8), allocatable:: fint(:), jl(:),  jltmp(:), djl(:),    &
     &              dfint(:)
      real(kind=8) xg, xrg, fac
!     ------------------------------------------------------------------
!     find  number of beta functions per species, max dimensions,
!     total number of beta functions (all and Vanderbilt only)
!     ------------------------------------------------------------------
      lmaxkb=-1
      nhm=0
      nhsa=0
      nhsavb=0
      nlcc_any=.false.
      do is=1,nsp
         ind=0
         do iv=1,nbeta(is)
            lmaxkb = max(lmaxkb,lll(iv,is))
            ind=ind+2*lll(iv,is)+1
         end do
         nh(is)=ind
         nhm=max(nhm,nh(is))
         ish(is)=nhsa
         nhsa=nhsa+na(is)*nh(is)
         if(ipp(is).le.1) nhsavb=nhsavb+na(is)*nh(is)
         !!! if(tvanp(is)) nhsavb=nhsavb+na(is)*nh(is)
         nlcc_any = nlcc_any .OR. nlcc(is)
      end do
      if (lmaxkb > lmaxx) call errore('nlinit ',' l > lmax ',lmaxkb)
      lmaxq = 2*lmaxkb + 1
      if (nhsa.le.0) call errore('nlinit ','not implemented ?',nhsa)
!
!     initialize array ap
!
      call aainit(lmaxkb+1)
!
      allocate(beta(ngw,nhm,nsp))
      allocate(qradb(ngb,nbrx,nbrx,lmaxq,nsp))
      allocate(qgb(ngb,nhm*(nhm+1)/2,nsp))
      allocate(qq(nhm,nhm,nsp))
      allocate(dvan(nhm,nhm,nsp))
      if (nlcc_any) allocate(rhocb(ngb,nsp))
      allocate(nhtol(nhm,nsp))
      allocate(indv (nhm,nsp))
      allocate(nhtolm(nhm,nsp))
!
      allocate(dqrad(ngb,nbrx,nbrx,lmaxq,nsp,3,3))
      allocate(dqgb(ngb,nhm*(nhm+1)/2,nsp,3,3))
      allocate(dbeta(ngw,nhm,nsp,3,3))
      allocate(betagx(mmx,nhm,nsp))
      allocate(dbetagx(mmx,nhm,nsp))
      allocate(qradx(mmx,nbrx,nbrx,lmaxq,nsp))
      allocate(dqradx(mmx,nbrx,nbrx,lmaxq,nsp))
!
      qradb(:,:,:,:,:) = 0.d0
      qq  (:,:,:) =0.d0
      dvan(:,:,:) =0.d0
      if(tpre) dqrad(:,:,:,:,:,:,:) = 0.d0
!
!     ------------------------------------------------------------------
!     definition of indices nhtol, indv, nhtolm
!     ------------------------------------------------------------------
      do is=1,nsp
         ind=0
         do iv=1,nbeta(is)
            lm = lll(iv,is)**2
            do il=1,2*lll(iv,is)+1
               lm=lm+1
               ind=ind+1
               nhtolm(ind,is)=lm
               nhtol(ind,is)=lll(iv,is)
               indv(ind,is)=iv
            end do
         end do
      end do
!
!     ===============================================================
!     initialization for vanderbilt species
!     ===============================================================
      do is=1,nvb
         if (tpre) then
            allocate(dfint(kkbeta(is)))
            allocate(djl(kkbeta(is)))
            allocate(jltmp(kkbeta(is)))
         end if
         allocate(fint(kkbeta(is)))
         allocate(jl(kkbeta(is)))
!
!     qqq and beta are now indexed and taken in the same order
!     as vanderbilts ppot-code prints them out
!
!     ---------------------------------------------------------------
!     calculation of array qradx(igb,iv,jv,is)
!     ---------------------------------------------------------------
         WRITE( stdout,*) ' nlinit  nh(is),ngb,is,kkbeta,lmaxq = ',             &
     &        nh(is),ngb,is,kkbeta(is),nqlc(is)
         do l=1,nqlc(is)
            do il=1,mmx
               xg=sqrt(refg*(il-1))
               call sph_bes (kkbeta(is), r(1,is), xg, l-1, jl)
!
               if(tpre) then
                  ltmp=l-1
                  !
                  ! r(i0) is the first point such that r(i0) >0
                  !
                  i0 = 1
                  if ( r(1,is) < 1.0d-8 ) i0 = 2  
                  ! special case q=0
                  if (xg < 1.0d-8) then
                     if (l == 1) then
                        ! Note that dj_1/dx (x=0) = 1/3
                        jltmp(:) = 1.0d0/3.d0
                     else
                        jltmp(:) = 0.0d0
                     end if
                  else
                     call sph_bes &
                          (kkbeta(is)+1-i0, r(i0,is), xg, ltmp-1, jltmp )
                  end if
                  do ir=i0, kkbeta(is)
                     xrg=r(ir,is)*xg
                     djl(ir)=jltmp(ir)*xrg-l*jl(ir)
                  end do
                  if (i0.eq.2) djl(1) = djl(2)
               endif
!
               do iv= 1,nbeta(is)
                  do jv=iv,nbeta(is)
!
!      note qrl(r)=r^2*q(r)
!
                     do ir=1,kkbeta(is)
                        fint(ir)=qrl(ir,iv,jv,l,is)*jl(ir)
                     end do
!!!         if (oldpseudo(is)) then
                     if (ipp(is).eq.0) then
                        call herman_skillman_int                        &
     &                    (kkbeta(is),cmesh(is),fint,qradx(il,iv,jv,l,is))
                     else
                        call simpson_cp90                               &
     &                    (kkbeta(is),fint,rab(1,is),qradx(il,iv,jv,l,is))
                     end if
                     qradx(il,jv,iv,l,is)=qradx(il,iv,jv,l,is)
!
                     if(tpre) then
                        do ir=1,kkbeta(is)
                           dfint(ir)=qrl(ir,iv,jv,l,is)*djl(ir)
                        end do
!!!         if (oldpseudo(is)) then
                        if (ipp(is).eq.0) then
                           call herman_skillman_int                     &
     &                          (kkbeta(is),cmesh(is),dfint,            &
     &                          dqradx(il,iv,jv,l,is))
                        else
                           call simpson_cp90                            &
     &                          (kkbeta(is),dfint,rab(1,is),            &
     &                          dqradx(il,iv,jv,l,is))
                        end if
                     end if
!
                  end do
               end do
            end do
         end do
!
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    qqq '
         do iv=1,nbeta(is)
            WRITE( stdout,'(8f9.4)') (qqq(iv,jv,is),jv=1,nbeta(is))
         end do
         WRITE( stdout,*)
!
         deallocate(jl)
         deallocate(fint)
         if (tpre) then
            deallocate(jltmp)
            deallocate(djl)
            deallocate(dfint)
         end if
!
      end do
!
!     ===============================================================
!     initialization that is common to all species
!     ===============================================================
      do is=1,nsp
!!!         if (.not.numeric(is)) then
         if (ipp(is).eq.3) then
            fac=1.0
         else
!     fac converts ry to hartree
            fac=0.5
         end if
         if (tpre) then
            allocate(dfint(kkbeta(is)))
            allocate(djl(kkbeta(is)))
         end if
         allocate(fint(kkbeta(is)))
         allocate(jl(kkbeta(is)))
         allocate(jltmp(kkbeta(is)))
!     ---------------------------------------------------------------
!     calculation of array  betagx(ig,iv,is)
!     ---------------------------------------------------------------
         WRITE( stdout,*)  '  betagx  '
         do iv=1,nh(is)
            l=nhtol(iv,is)+1
            do il=1,mmx
               xg=sqrt(refg*(il-1))
               call sph_bes (kkbeta(is), r(1,is), xg, l-1, jl )
!
               if(tpre)then
                  ltmp=l-1
                  !
                  ! r(i0) is the first point such that r(i0) >0
                  !
                  i0 = 1
                  if ( r(1,is) < 1.0d-8 ) i0 = 2  
                  ! special case q=0
                  if (xg < 1.0d-8) then
                     if (l == 1) then
                        ! Note that dj_1/dx (x=0) = 1/3
                        jltmp(:) = 1.0d0/3.d0
                     else
                        jltmp(:) = 0.0d0
                     end if
                  else
                     call sph_bes &
                          (kkbeta(is)+1-i0, r(i0,is), xg, ltmp-1, jltmp )
                  end if
                  do ir=i0, kkbeta(is)
                     xrg=r(ir,is)*xg
                     djl(ir)=jltmp(ir)*xrg-l*jl(ir)
                  end do
                  if (i0.eq.2) djl(1) = djl(2)
!
               endif
!
!     beta(ir)=r*beta(r)
!
               do ir=1,kkbeta(is)
                  fint(ir)=r(ir,is)*betar(ir,indv(iv,is),is)*jl(ir)
               end do
!!!         if (oldpseudo(is)) then
               if (ipp(is).eq.0) then
                  call herman_skillman_int                              &
     &                 (kkbeta(is),cmesh(is),fint,betagx(il,iv,is))
               else
                  call simpson_cp90                                     &
     &                 (kkbeta(is),fint,rab(1,is),betagx(il,iv,is))
               endif
!
               if(tpre) then
                  do ir=1,kkbeta(is)
                     dfint(ir)=r(ir,is)*betar(ir,indv(iv,is),is)*djl(ir)
                  end do
!!!         if (oldpseudo(is)) then
                  if (ipp(is).eq.0) then
                     call herman_skillman_int                           &
     &                 (kkbeta(is),cmesh(is),dfint,dbetagx(il,iv,is))
                  else
                     call simpson_cp90                                  &
     &                 (kkbeta(is),dfint,rab(1,is),dbetagx(il,iv,is))
                  end if
               endif
!
            end do
         end do
! 
!     ---------------------------------------------------------------
!     calculate array  dvan(iv,jv,is)
!     ---------------------------------------------------------------
         do iv=1,nh(is)
            do jv=1,nh(is)
               if ( nhtolm(iv,is) == nhtolm(jv,is) ) then
                  dvan(iv,jv,is)=fac*dion(indv(iv,is),indv(jv,is),is)
               endif 
            end do
         end do
!
         do iv=1,nh(is)
            WRITE( stdout,901) iv,indv(iv,is),nhtol(iv,is)
         end do
 901     format(2x,i2,'  indv= ',i2,'   ang. mom= ',i2)
!
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    dion '
         do iv=1,nbeta(is)
            WRITE( stdout,'(8f9.4)') (fac*dion(iv,jv,is),jv=1,nbeta(is))
         end do
!
         deallocate(jltmp)
         deallocate(jl)
         deallocate(fint)
         if (tpre) then
            deallocate(djl)
            deallocate(dfint)
         end if
      end do
!
! newnlinit stores qgb and qq, calculates arrays  beta  qradb  rhocb
! and derivatives wrt cell    dbeta dqrad
!
      call newnlinit
!
      return
      end

!-------------------------------------------------------------------------
      subroutine qvan2b(ngy,iv,jv,is,ylm,qg)
!--------------------------------------------------------------------------
!     q(g,l,k) = sum_lm (-i)^l ap(lm,l,k) yr_lm(g^) qrad(g,l,l,k)
!
      use control_flags, only: iprint, tpre
      use qradb_mod
      use uspp, only: nlx, lpx, lpl, ap, indv, nhtolm
      use gvecb
      use cdvan
      use uspp_param, only: lmaxq
! 
      implicit none
      integer ngy, iv, jv, is
      complex(kind=8)   qg(ngb)
!
      integer ivs, jvs, ivl, jvl, i, ii, ij, l, lp, ig
      complex(kind=8) sig
      real(kind=8) :: ylm(ngb,lmaxq*lmaxq), dylm(ngb,lmaxq*lmaxq,3,3)
! 
!       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
!       ivs = 1..4     s_1 s_2 p_1 p_2
!       ivl = 1..4     s p_x p_z p_y
! 
      ivs=indv(iv,is)
      jvs=indv(jv,is)
      ivl=nhtolm(iv,is)
      jvl=nhtolm(jv,is)
      if(ivl > nlx)  call errore(' qvan2b ',' ivl out of bounds ',ivl)
      if(jvl > nlx)  call errore(' qvan2b ',' jvl out of bounds ',jvl)
!
      qg(:) = (0.d0, 0.d0)
!
!     lpx = max number of allowed y_lm
!     lp  = composite lm to indentify them
!
      do i=1,lpx(ivl,jvl)
         lp=lpl(ivl,jvl,i)
         if (lp > lmaxq*lmaxq) call errore(' qvan2b ',' lp out of bounds ',lp)
!
!     extraction of angular momentum l from lp:  
!     l = int ( sqrt( float(l-1) + epsilon) ) + 1
!
         if (lp == 1) then
            l=1         
         else if ((lp >= 2) .and. (lp <= 4)) then
            l=2
         else if ((lp >= 5) .and. (lp <= 9)) then
            l=3
         else if ((lp >= 10).and.(lp <= 16)) then
            l=4
         else if ((lp >= 17).and.(lp <= 25)) then
            l=5
         else if ((lp >= 26).and.(lp <= 36)) then 
            l=6
         else if ((lp >= 37).and.(lp <= 49)) then 
            l=7
         else
            call errore(' qvan2b ',' not implemented ',lp)
         endif
!     
!       sig= (-i)^l
!
         sig=(0.,-1.)**(l-1)
         sig=sig*ap(lp,ivl,jvl)
         do ig=1,ngy
            qg(ig)=qg(ig)+sig*ylm(ig,lp)*qradb(ig,ivs,jvs,l,is)
         end do
      end do
!
      return
      end

!-------------------------------------------------------------------------
      subroutine dqvan2b(ngy,iv,jv,is,ylm,dylm,dqg)
!--------------------------------------------------------------------------
!
!     dq(i,j) derivatives wrt to h(i,j) of q(g,l,k) calculated in qvan2b
!
      use control_flags, only: iprint, tpre
      use qradb_mod
      use uspp, only: nlx, lpx, lpl, ap, indv, nhtolm
      use gvecb
      use dqrad_mod
      use cdvan
      use uspp_param, only: lmaxq
! 
      implicit none
      integer ngy, iv, jv, is
      complex(kind=8) dqg(ngb,3,3)
!
      integer ivs, jvs, ivl, jvl, i, ii, ij, l, lp, ig
      complex(kind=8) sig
      real(kind=8) :: ylm(ngb,lmaxq*lmaxq), dylm(ngb,lmaxq*lmaxq,3,3)
! 
!       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
!       ivs = 1..4     s_1 s_2 p_1 p_2
!       ivl = 1..4     s p_x p_z p_y
! 
      ivs=indv(iv,is)
      jvs=indv(jv,is)
      ivl=nhtolm(iv,is)
      jvl=nhtolm(jv,is)
      if(ivl > nlx)  call errore(' dqvan2b ',' ivl out of bounds ',ivl)
      if(jvl > nlx)  call errore(' dqvan2b ',' jvl out of bounds ',jvl)
!
      dqg(:,:,:) = (0.d0, 0.d0)
!
!     lpx = max number of allowed y_lm
!     lp  = composite lm to indentify them
!
      do i=1,lpx(ivl,jvl)
         lp=lpl(ivl,jvl,i)
         if (lp > lmaxq*lmaxq) call errore(' dqvan2b ',' lp out of bounds ',lp)
!
!     extraction of angular momentum l from lp:  
!     l = int ( sqrt( float(l-1) + epsilon) ) + 1
!
         if (lp == 1) then
            l=1         
         else if ((lp >= 2) .and. (lp <= 4)) then
            l=2
         else if ((lp >= 5) .and. (lp <= 9)) then
            l=3
         else if ((lp >= 10).and.(lp <= 16)) then
            l=4
         else if ((lp >= 17).and.(lp <= 25)) then
            l=5
         else if ((lp >= 26).and.(lp <= 36)) then 
            l=6
         else if ((lp >= 37).and.(lp <= 49)) then 
            l=7
         else
            call errore(' qvan2b ',' not implemented ',lp)
         endif
!     
!       sig= (-i)^l
!
         sig=(0.,-1.)**(l-1)
         sig=sig*ap(lp,ivl,jvl)
         do ij=1,3
            do ii=1,3
               do ig=1,ngy
                  dqg(ig,ii,ij) = dqg(ig,ii,ij) +  sig *                &
     &                    ( ylm(ig,lp) * dqrad(ig,ivs,jvs,l,is,ii,ij) + &
     &                     dylm(ig,lp,ii,ij)*qradb(ig,ivs,jvs,l,is)   )
               end do
            end do
         end do
      end do
!
      return
      end
!-----------------------------------------------------------------------
subroutine dylmr2_(nylm, ngy, g, gg, ainv, dylm)
  !-----------------------------------------------------------------------
  !
  ! temporary CP interface for PW routine dylmr2
  ! dylmr2  calculates d Y_{lm} /d G_ipol
  ! dylmr2_ calculates G_ipol \sum_k h^(-1)(jpol,k) (dY_{lm} /dG_k)
  !
  USE kinds
  implicit none
  !
  integer, intent(IN) :: nylm, ngy
  real(kind=DP), intent(IN) :: g (3, ngy), gg (ngy), ainv(3,3)
  real(kind=DP), intent(OUT) :: dylm (ngy, nylm, 3, 3)
  !
  integer :: ipol, jpol, lm
  real(kind=DP), allocatable :: dylmaux (:,:,:)
  !
  allocate ( dylmaux(ngy,nylm,3) )
  dylmaux(:,:,:) = 0.d0
  do ipol =1,3
     call dylmr2 (nylm, ngy, g, gg, dylmaux(1,1,ipol), ipol)
  enddo
  !
  do ipol =1,3
     do jpol =1,3
        do lm=1,nylm
           dylm (:,lm,ipol,jpol) = (dylmaux(:,lm,1) * ainv(jpol,1) + &
                                    dylmaux(:,lm,2) * ainv(jpol,2) + &
                                    dylmaux(:,lm,3) * ainv(jpol,3) ) &
                                    * g(ipol,:)
        end do
     end do
  end do
  deallocate ( dylmaux )
end subroutine dylmr2_
