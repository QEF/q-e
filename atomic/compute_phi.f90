!
! Copyright (C) 2001-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine compute_phi(lam,ik,nwf0,ns,xc,iflag,iok,occ)
  !--------------------------------------------------------------------------
  !
  !     This routine computes the phi functions by pseudizing the
  !     all_electron chi functions. In input it receives, the point
  !     ik where the cut is done, the angular momentum lam,
  !     and the correspondence with the all electron wavefunction
  !
  !
  !      
  use ld1inc
  implicit none

  integer ::    &
       lam,  &   ! the angular momentum
       ik,   &   ! the point corresponding to rc
       nwf0, &   ! 
       ns,   &   ! the function to pseudize
       iflag,&   ! if 1 print
       iok       ! output: number of nodes of the pseudized wavefunction

  real(kind=dp) :: &
       xc(8), occ
  !
  real(kind=dp), parameter :: pi=3.14159265358979d0
  real(kind=dp) :: &
       fae,    & ! the value of the all-electron function
       f1ae,   & ! its first derivative
       f2ae,   & ! the second derivative
       wmax,   & !
       den,    & ! denominator
       faenor    ! the norm of the function

  integer ::    &
       isign,   &! sign of the max of the ae-wfc
       m, n, nst, nnode, nc, nc1, ij, imax, iq, i

  real(kind=dp) :: &
       chir(ndm,nwfx), rab(ndm), chi_dir(ndm,2), gi(ndm), j1(ndm,4), &
       f1aep1, f1aem1, jnor, psnor, fact(4), &
       cm(10), bm(4), ze2, cn(6), c2, &
       delta, a, b, c, deter, gamma, &
       lamda0, lamda3, lamda4, mu0, mu4, s0, s4, t0, t4

  real(kind=dp), external :: deriv_7pts, deriv2_7pts, int_0_inf_dr, pr
  integer ::  &
       nbes   ! number of Bessel functions to be used

  !
  !   decide whether to use 4 Bessel functions or 3 (default)
  !
  if ( (rho0 /= 0.d0) .and. (lam == 0) ) then
     nbes = 4
  else
     nbes = 3
  end if
  !
  nst=(lam+1)*2
  !
  !   compute the reference wavefunction
  !
  if (new(ns)) then
     if (rel.eq.1) then
        call lschps(3,zed,exp(dx),dx,mesh,mesh,mesh,  &
             1,lam,enls(ns),chir(1,ns),r,vpot)
     elseif (rel.eq.2) then
        do i=1,mesh
           rab(i)=r(i)*dx
        enddo
        call dir_outward(ndm,mesh,lam,jjs(ns),enls(ns),dx,chi_dir,r,rab,vpot)
        chir(:,ns)=chi_dir(:,2)
     else
        ze2=-zed*2.d0
        call intref(lam,enls(ns),mesh,dx,r,r2,sqr,vpot,ze2,chir(1,ns))
     endif
     !
     !    fix arbitrarily the norm at the cut-off radius equal to 0.5
     !
     jnor=chir(ik,ns)
     do n=1,mesh
        chir(n,ns)=chir(n,ns)*0.5d0/jnor
    enddo
  else 
     do n=1,mesh
        chir(n,ns)=psi(n,nwf0)
     enddo
  endif
  !
  !   compute the first and second derivative of all-electron function
  !
  fae=chir(ik,ns)
  f1ae=deriv_7pts(chir(1,ns),ik,r(ik),dx)
  f2ae=deriv2_7pts(chir(1,ns),ik,r(ik),dx)
  !
  !   compute the norm of the all-electron function
  !
  do n=1,ik+1
     gi(n)=chir(n,ns)**2  
  enddo
  faenor=int_0_inf_dr(gi,r,r2,dx,ik,nst)
  !
  !
  if (tm) then
     !
     ! TM: the pseudo-wavefunction is written as polynomial times exponential
     !
     call find_coefficients &
          (ik, chir(1,ns), enls(ns), r, dx, faenor, vpot, cn, c2, lam, mesh)
      do i=1,ik
         phis(i,ns) =sign( r(i)**(lam+1)*exp(pr(cn,c2,r(i))), &
                          chir(ik,ns) )
      end do
      xc(1:6) = cn(:)
      xc(7) = c2
  else
     !   RRKJ: the pseudo-wavefunction is written as an expansion into (3 or 4) 
     !         spherical bessel functions for r < r_c
     !   find q_i with the correct log derivatives
     
     call find_qi(f1ae/fae,xc(nbes+1),ik,lam,nbes,1,iok)
     if (iok.ne.0) return
     !
     !    compute the bessel functions
     !
     do nc=1,nbes
        call sph_besr(ik+5,r,xc(nbes+nc),lam,j1(1,nc))
        jnor=j1(ik,nc)
        fact(nc)=chir(ik,ns)/jnor
        do n=1,ik+5
           j1(n,nc)=j1(n,nc)*chir(ik,ns)/jnor
        enddo
     enddo
  !
  !    compute the bm functions (second derivative of the j1)
  !    and the integrals for cm
  !
     ij=0
     do nc=1, nbes
        bm(nc)=deriv2_7pts(j1(1,nc),ik,r(ik),dx)
        do nc1=1,nc
           ij=ij+1
           do n=1,ik
              gi(n)=j1(n,nc)*j1(n,nc1)
           enddo
           cm(ij)=int_0_inf_dr(gi,r,r2,dx,ik,nst)
        enddo
     enddo
     !
     !    solve the second order equation
     !
     if (nbes == 4) then
        wmax=0.d0
        do n=1,mesh
           if (abs(chir(n,ns)) > wmax .and. r(n) < 4.d0) then
              wmax=abs(chir(n,ns))
              if(chir(n,ns).lt.0.d0)then
                 isign=-1
              else
                 isign=+1
              endif
           endif
        enddo
        lamda0=(f2ae-bm(1))/(bm(2)-bm(1))
        lamda3=(bm(3)-bm(1))/(bm(2)-bm(1))
        lamda4=(bm(4)-bm(1))/(bm(2)-bm(1))
        if (new(ns)) then
           den=1.d0
        else
           den=abs(occ)
        endif
        mu0=(isign*sqrt(rho0*4.d0*pi/den) &
             -fact(1)+fact(1)*lamda0-fact(2)*lamda0)  &
             /(fact(1)*lamda3-fact(1)-fact(2)*lamda3+fact(3))
        mu4=(fact(1)*lamda4-fact(1)-fact(2)*lamda4+fact(4)) &
             /(fact(1)*lamda3-fact(1)-fact(2)*lamda3+fact(3))
        s0=lamda0-mu0*lamda3
        s4=lamda4-mu4*lamda3
        t0=1.d0-lamda0+mu0*lamda3-mu0
        t4=1.d0+mu4*lamda3-mu4-lamda4

        a=cm(1)*t4**2+cm(3)*s4**2+cm(6)*mu4**2+cm(10)  &
             +2.d0*cm(4)*t4*mu4+2.d0*cm(2)*t4*s4+2.d0*cm(5)*s4*mu4 &
             -2.d0*cm(7)*t4-2.d0*cm(8)*s4-2.d0*cm(9)*mu4
        
        b=-2.d0*cm(1)*t0*t4-2.d0*cm(3)*s0*s4-2.d0*cm(6)*mu0*mu4  &
             -2.d0*cm(4)*t0*mu4-2.d0*cm(4)*t4*mu0-2.d0*cm(2)*t0*s4  &
             -2.d0*cm(2)*t4*s0-2.d0*cm(5)*s0*mu4-2.d0*cm(5)*s4*mu0  &
             +2.d0*cm(7)*t0+2.d0*cm(8)*s0+2.d0*cm(9)*mu0
        
        c=cm(1)*t0**2+cm(3)*s0**2+cm(6)*mu0**2+2.d0*cm(4)*t0*mu0 &
             +2.d0*cm(2)*t0*s0+2.d0*cm(5)*s0*mu0-faenor
        
        deter=b**2-4.d0*a*c
        if (deter < 0.d0) then
           call errore('compute phi','negative determinant',-1) 
           write(6,100) ns,f1ae/fae, f2ae, faenor
100        format(/5x,i5,' ld= ',f10.6,' f2ae',f10.6,' faenor',f10.6)
           iok=1
           return
        endif
        
        xc(4)=(-b+sqrt(deter))/(2.d0*a)
        xc(3)=mu0-mu4*xc(4)
        xc(2)=s0-s4*xc(4)
        xc(1)=t0-t4*xc(4)
        !
        
        do n=1,ik
           phis(n,ns)= xc(1)*j1(n,1) + xc(2)*j1(n,2)  &
                     + xc(3)*j1(n,3) + xc(4)*j1(n,4)
        enddo
     else
        gamma=(bm(3)-bm(1))/(bm(2)-bm(1))
        delta=(f2ae-bm(1))/(bm(2)-bm(1))
       
        a=(gamma-1.d0)**2*cm(1)+gamma**2*cm(3)  &
             -2.d0*gamma*(gamma-1.d0)*cm(2) &
             +2.d0*(gamma-1.d0)*cm(4)-2.d0*gamma*cm(5)+cm(6)
        
        b=2.d0*(1.d0-delta)*(gamma-1.d0)*cm(1) &
             -2.d0*gamma*delta*cm(3) &
             -2.d0*gamma*(1.d0-delta)*cm(2) &
             +2.d0*delta*(gamma-1.d0)*cm(2) & 
             +2.d0*(1.d0-delta)*cm(4)+2.d0*delta*cm(5)
        
        c=-faenor+(1.d0-delta)**2*cm(1)+delta**2*cm(3) &
             +2.d0*delta*(1.d0-delta)*cm(2)
        
        deter=b**2-4.d0*a*c
        if (deter < 0.d0) then
           call errore('compute phi','negative determinant',-1) 
           write(6,110) ns,f1ae/fae, f2ae, faenor
110        format (/5x,i5,' ld= ',f10.6,' f2ae',f10.6, ' faenor',f10.6)
           iok=1
           return
        endif
        
        xc(3)=(-b+sqrt(deter))/(2.d0*a)
        xc(2)=-xc(3)*gamma+delta
        xc(1)=1.d0-xc(2)-xc(3)
        
        do n=1,ik
           phis(n,ns)= xc(1)*j1(n,1) + xc(2)*j1(n,2) &
                     + xc(3)*j1(n,3)
        enddo

     endif
     do nc=1,nbes
        xc(nc)=xc(nc)*fact(nc)
     enddo
  end if
  !
  !      for r > r(ik) the pseudo and all-electron psi(r) coincide
  !
  do n=ik+1,mesh
     phis(n,ns)= chir(n,ns)
  enddo
  !
  !      check for the norm of the pseudo wavefunction
  !
  do n=1,ik
     gi(n)=phis(n,ns)**2
  enddo
  psnor=int_0_inf_dr(gi,r,r2,dx,ik,nst)

  if (iflag == 1) then
     if (tm) then
        write(6,120) els(ns), rcut(ns)
120     format (/ /5x, ' Wfc  ',a3,'  rcut=',f6.3, &
          '  Using Troullier-Martins method ')
     else
        write(6,130) els(ns),rcut(ns),2.d0*xc(6)**2 
130     format (/ /5x, ' Wfc  ',a3,'  rcut=',f6.3, &
          '  Estimated cut-off energy= ', f11.2,' Ry')
        if (nbes == 4) write(6,140) rho0
140        format (5x,' Using 4 Bessel functions for this wfc,', &
                ' rho(0) =',f6.3)
     end if
     ! write(6,'(5x," AE norm = ",f6.3,"  PS norm = ",f6.3)') faenor, psnor
  end if
  !
  !      check for absence of nodes in the pseudo wavefunction
  !
  nnode=0  
  do n=1,ik+1
     if ( phis(n,ns) .ne. sign(phis(n,ns),phis(n+1,ns)) ) then
        write(6,150) lam,ns,r(n)
150     format (5x,'l=',i4,' ns=',i4,' Node at ',f10.8)
        nnode=nnode+1
     endif
  enddo
  iok=nnode
  if (iflag == 1) write(6,160) nnode,r(ik)
160 format (5x,' This function has ',i4,' nodes', ' for 0 < r < ',f8.3)

  return
end subroutine compute_phi
