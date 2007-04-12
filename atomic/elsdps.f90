!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
      subroutine elsdps 
!---------------------------------------------------------------
!
!   atomic total energy in the local-spin-density scheme
!   atomic pseudopotentials with nonlinear core correction are allowed
!   gradient correction  allowed (A. Dal Corso fecit AD 1993)
!
use kinds, only: DP
use constants, only: fpi
use ld1inc
use funct, only : dft_is_gradient
implicit none
      real(DP) :: &
             int_0_inf_dr,  &   ! the integral function
             rh(2),rh0(2),  &   ! the charge in a given point
             rhc,           &   ! core charge in a given point
             vxc,           &   ! the exchange correlation potential
             exc,           &   ! the exchange correlation energy
             exc_t,         &   ! the exchange correlation energy
             vxcp(2),       &   
             rho_tot,       &   !
             gi(ndm),       &  ! 
             work(nwfsx)        ! 

      real(DP),allocatable :: &
             f1(:,:),     &   ! auxiliary
             f2(:),      &   ! auxiliary
             f3(:),      &   ! auxiliary
             f4(:),      &   ! auxiliary
             f8(:),      &   ! auxiliary
             vgc(:,:),    &   ! the gga potential
             egc(:),      &   ! the gga energy
             egcc(:)          ! the gga core energy

      integer :: &
             n,i,ns,is,nst,lam,n1,n2,ikl,ierr, ind

      logical :: &
             gga                ! if true it is a gga calculation


      gga=dft_is_gradient()
      allocate(vgc(ndm,2),stat=ierr)
      allocate(egc(ndm),stat=ierr)
      if (nlcc) allocate(egcc(ndm), stat=ierr)

      allocate(f1(ndm,2), stat=ierr)
      allocate(f2(mesh), stat=ierr)
      allocate(f3(mesh), stat=ierr)
      allocate(f4(mesh), stat=ierr)
      allocate(f8(mesh), stat=ierr)
      vgc=0.0_DP
      egc=0.0_DP
      if (gga.and.nlcc) then
         f1=0.0_DP
         call vxcgc(ndm,mesh,nspin,r,r2,f1,rhoc,vgc,egcc)
      endif

      if (gga) call vxcgc(ndm,mesh,nspin,r,r2,rhos,rhoc,vgc,egc)

      rh0(1)=0.0_DP
      rh0(2)=0.0_DP
      rhc=0.0_DP
      do i=1,mesh
         rho_tot=rhos(i,1)
         if (lsd.eq.1) rho_tot=rho_tot+rhos(i,2)
         f1(i,1)= vpsloc(i)*rho_tot
         f4(i)= vxt(i)*rho_tot
         vh(i)= vh(i)*rho_tot
         !
         do is=1,nspin
            rh(is) = rhos(i,is)/r2(i)/fpi
         enddo
         if (nlcc) then
            rhc= rhoc(i)/r2(i)/fpi
            call vxc_t(rh,rhc,lsd,vxcp)
            if (gga) then
               f3(i) = exc_t(rh,rhc,lsd)*(rho_tot+rhoc(i)) &
                     + egc(i)*r2(i)*fpi 
               f8(i) = exc_t(rh0,rhc,lsd)*rhoc(i) + &
                       egcc(i)*r2(i)*fpi
               f2(i) =-(vgc(i,1)+vxcp(1))*rhos(i,1) &
                      -f1(i,1)-vh(i)-f4(i)
               if (lsd.eq.1) f2(i)=f2(i)-  &
                              (vgc(i,2)+vxcp(2))*rhos(i,2)
            else
               f3(i) = exc_t(rh,rhc,lsd) * (rho_tot+rhoc(i)) 
               f8(i) = exc_t(rh0,rhc,lsd)*rhoc(i)
               f2(i) =-vxcp(1)*rhos(i,1)-f1(i,1)-vh(i)-f4(i)
               if (lsd.eq.1) f2(i)=f2(i)-vxcp(2)*rhos(i,2)
            endif
         else
            call vxc_t(rh,rhc,lsd,vxcp)
            if (gga) then
               f3(i) = exc_t(rh,rhc,lsd)*rho_tot + &
                       egc(i)*r2(i)*fpi
               f2(i) =-(vgc(i,1)+vxcp(1))*rhos(i,1) &
                                 -f1(i,1)-vh(i)-f4(i)
               if (lsd.eq.1) f2(i)=f2(i)  &
                           -(vgc(i,2)+vxcp(2))*rhos(i,2)
            else
               f3(i) = exc_t(rh,rhc,lsd)* rho_tot
               f2(i) =-vxcp(1)*rhos(i,1)-f1(i,1)-vh(i)-f4(i)
               if (lsd.eq.1) f2(i)=f2(i)-vxcp(2)*rhos(i,2)
            endif
         endif
      enddo

      encl=    int_0_inf_dr(f1,r,r2,dx,mesh,1)
      ehrt=0.5_DP*int_0_inf_dr(vh,r,r2,dx,mesh,2)
      ecxc=    int_0_inf_dr(f3,r,r2,dx,mesh,2)
      evxt=    int_0_inf_dr(f4,r,r2,dx,mesh,2)
      if (nlcc) then
         ecc=    int_0_inf_dr(f8,r,r2,dx,mesh,2)
      else
         ecc= 0.d0
      endif
!
      epseu=0.0_DP
      if (pseudotype == 1) then
         do ns=1,nwfts
            if ( rel < 2 .or. llts(ns) == 0 .or. &
                 abs(jjts(ns)-llts(ns)+0.5_dp) < 0.001_dp) then
               ind=1
            else if ( rel == 2 .and. llts(ns) > 0 .and. &
                 abs(jjts(ns)-llts(ns)-0.5_dp) < 0.001_dp) then
               ind=2
            endif
            f1=0.0_DP
            lam=llts(ns)
            if (octs(ns).gt.0.0_DP) then
               do n=1, mesh
                  f1(n,1) = f1(n,1) + phis(n,ns)**2 * octs(ns)
               enddo
            endif
            do n=1,mesh
               f1(n,1) = f1(n,1) * vnl(n,lam,ind)
            end do
            if (ikk(ns) > 0) &
                epseu = epseu + int_0_inf_dr(f1,r,r2,dx,ikk(ns),2*(lam+1))
         enddo
      else
         do ns=1,nwfts
            if (octs(ns).gt.0.0_DP) then
               do n1=1,nbeta
                  if ( llts(ns).eq.lls(n1).and.   &
                       abs(jjts(ns)-jjs(n1)).lt.1.e-7_DP) then
                     nst=(llts(ns)+1)*2
                     ikl=ikk(n1)
                     do n=1,ikl
                        gi(n)=betas(n,n1)*phis(n,ns)
                     enddo
                     work(n1)=int_0_inf_dr(gi,r,r2,dx,ikl,nst)
                  else
                     work(n1)=0.0_DP
                  endif
               enddo
               do n1=1,nbeta
                  do n2=1,nbeta
                     epseu=epseu  &
                       + bmat(n1,n2)*work(n1)*work(n2)*octs(ns)
                  enddo
               enddo
            endif
         enddo
      endif 
      ekin = int_0_inf_dr(f2,r,r2,dx,mesh,2) - epseu
      do ns=1,nwfts
         if (octs(ns).gt.0.0_DP) then
            ekin=ekin+octs(ns)*enls(ns)
         endif
      end do

      etots= ekin + encl + epseu + ehrt + ecxc + evxt

      deallocate(f8)
      deallocate(f4)
      deallocate(f3)
      deallocate(f2)
      deallocate(f1)
      if (nlcc) deallocate(egcc)
      deallocate(egc)
      deallocate(vgc)

      return
      end subroutine elsdps
