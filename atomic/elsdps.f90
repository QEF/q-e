!
!---------------------------------------------------------------
      subroutine elsdps 
!---------------------------------------------------------------
!
!   atomic total energy in the local-spin-density scheme
!   atomic pseudopotentials with nonlinear core correction are allowed
!   gradient correction  allowed (A. Dal Corso fecit AD 1993)
!
use ld1inc
use funct
implicit none
      real(kind=dp) :: &
             int_0_inf_dr,  &   ! the integral function
             rh(2),rh0(2),  &   ! the charge in a given point
             rhc,           &   ! core charge in a given point
             vxc,           &   ! the exchange correlation potential
             exc,           &   ! the exchange correlation energy
             exc_t,         &   ! the exchange correlation energy
             ecc,           &   ! the core correction energy
             vxcp(2),       &   
             rho_tot,       &   !
             gi(ndm),       &  ! 
             work(nwfsx)        ! 

      real(kind=dp),allocatable :: &
             f1(:,:),     &   ! auxiliary
             f2(:),      &   ! auxiliary
             f3(:),      &   ! auxiliary
             f4(:),      &   ! auxiliary
             f8(:),      &   ! auxiliary
             vgc(:,:),    &   ! the gga potential
             egc(:),      &   ! the gga energy
             egcc(:)          ! the gga core energy

      integer :: &
             n,i,ns,is,nst,lam,n1,n2,ikl,ierr

      logical :: &
             gga                ! if true it is a gga calculation

      real(kind=dp), parameter :: fourpi = 4.d0 * 3.141592653589793d0  

      gga=igcx.ne.0.or.igcc.ne.0
      allocate(vgc(ndm,2),stat=ierr)
      allocate(egc(ndm),stat=ierr)
      if (nlcc) allocate(egcc(ndm), stat=ierr)

      allocate(f1(ndm,2), stat=ierr)
      allocate(f2(mesh), stat=ierr)
      allocate(f3(mesh), stat=ierr)
      allocate(f4(mesh), stat=ierr)
      allocate(f8(mesh), stat=ierr)
      vgc=0.d0
      egc=0.d0
      if (gga.and.nlcc) then
         f1=0.d0
         call vxcgc(ndm,mesh,nspin,r,r2,f1,rhoc,vgc,egcc)
      endif

      if (gga) call vxcgc(ndm,mesh,nspin,r,r2,rhos,rhoc,vgc,egc)

      rh0(1)=0.d0
      rh0(2)=0.d0
      rhc=0.d0
      do i=1,mesh
         rho_tot=rhos(i,1)
         if (lsd.eq.1) rho_tot=rho_tot+rhos(i,2)
         f1(i,1)= vpsloc(i)*rho_tot
         f4(i)= vxt(i)*rho_tot
         vh(i)= vh(i)*rho_tot
         do is=1,nspin
            rh(is) = rhos(i,is)/r2(i)/fourpi
         enddo
         if (nlcc) then
            rhc= rhoc(i)/r2(i)/fourpi
            call vxc_t(rh,rhc,lsd,vxcp)
            if (gga) then
               f2(i) =-(vgc(i,1)+vxcp(1))*rhos(i,1) &
                      -f1(i,1)-vh(i)-f4(i)
               f3(i) = exc_t(rh,rhc,lsd)*(rho_tot+rhoc(i)) &
                     + egc(i)*r2(i)*fourpi  &
                     - exc_t(rh0,rhc,lsd)*rhoc(i) &
                     - egcc(i)*r2(i)*fourpi
               f8(i) = exc_t(rh0,rhc,lsd)*rhoc(i) + &
                       egcc(i)*r2(i)*fourpi
               if (lsd.eq.1) f2(i)=f2(i)-  &
                              (vgc(i,2)+vxcp(2))*rhos(i,2)
            else
               f2(i) =-vxcp(1)*rhos(i,1)-f1(i,1)-vh(i)-f4(i)
               f3(i) = exc_t(rh,rhc,lsd) * (rho_tot+rhoc(i)) &
                     - exc_t(rh0,rhc,lsd)*rhoc(i)
               f8(i) = exc_t(rh0,rhc,lsd)*rhoc(i)
               if (lsd.eq.1) f2(i)=f2(i)-vxcp(2)*rhos(i,2)
            endif
         else
            call vxc_t(rh,rhc,lsd,vxcp)
            if (gga) then
               f2(i) =-(vgc(i,1)+vxcp(1))*rhos(i,1) &
                                 -f1(i,1)-vh(i)-f4(i)
               f3(i) = exc_t(rh,rhc,lsd)*rho_tot + &
                       egc(i)*r2(i)*fourpi
               if (lsd.eq.1) f2(i)=f2(i)  &
                           -(vgc(i,2)+vxcp(2))*rhos(i,2)
            else
               f2(i) =-vxcp(1)*rhos(i,1)-f1(i,1)-vh(i)-f4(i)
               f3(i) = exc_t(rh,rhc,lsd)* rho_tot
               if (lsd.eq.1) f2(i)=f2(i)-vxcp(2)*rhos(i,2)
            endif
         endif
      enddo

      encl=    int_0_inf_dr(f1,r,r2,dx,mesh,1)
      ehrt=0.5d0*int_0_inf_dr(vh,r,r2,dx,mesh,2)
      ecxc=    int_0_inf_dr(f3,r,r2,dx,mesh,2)
      evxt=    int_0_inf_dr(f4,r,r2,dx,mesh,2)
      if (nlcc) then
         ecc=    int_0_inf_dr(f8,r,r2,dx,mesh,2)
         write(6,'(5x,'' Core only energy '',f15.8 )') ecc
      endif
!
      epseu=0.d0
      if (pseudotype.eq.1) then
         do ns=1,nwfts
            f1=0.d0
            lam=llts(ns)
            if (octs(ns).gt.0.d0) then
               do n=1, mesh
                  f1(n,1) = f1(n,1) + phis(n,ns)**2 * octs(ns)
               enddo
            endif
            do n=1,mesh
               f1(n,1) = f1(n,1) * vnl(n,lam)
            end do
            if (ikk(ns) > 0) &
                epseu = epseu + int_0_inf_dr(f1,r,r2,dx,ikk(ns),2*(lam+1))
         enddo
      elseif ((pseudotype.eq.2).or.(pseudotype.eq.3)) then
         do ns=1,nwfts
            if (octs(ns).gt.0.d0) then
               do n1=1,nbeta
                  if (llts(ns).eq.lls(n1).and.   &
                                   abs(jjts(ns)-jjs(n1)).lt.1.d-7) then
                     nst=(llts(ns)+1)*2
                     ikl=ikk(n1)
                     do n=1,ikl
                        gi(n)=betas(n,n1)*phis(n,ns)
                     enddo
                     work(n1)=int_0_inf_dr(gi,r,r2,dx,ikl,nst)
                  else
                     work(n1)=0.d0
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
         if (octs(ns).gt.0.d0) then
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
      end
