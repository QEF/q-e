!
!--------------------------------------------------------------------------
      subroutine compute_phi(lam,ik,nwf0,ns,xc,iflag,iok,occ)
!--------------------------------------------------------------------------
!
!     This routine computes the phi functions by pseusizing the
!     all_electron chi functions. In input it receives, the point
!     ik where the cut is done, the angualar momentum lam,
!     and the correspondence with the all electron wavefunction
!
!
!      
use ld1inc
implicit none

real(kind=dp) :: &
        fae,    & ! the value of the all-electron function
        f1ae,   & ! its first derivative
        f2ae,   & ! the second derivative
        ff,     & ! compute the second derivative
        wmax,   & !
        den,    & ! denominator
        faenor    ! the norm of the function
      
real(kind=dp), parameter :: pi=3.14159265358979d0

integer ::    &
        ik,   &   ! the point corresponding to rc
        isign,&   ! sign of the max of the ae-wfc
        ns,   &   ! the function to pseudize
        iflag,&   ! if 1 print
        iok,  &   ! gives the number of nodes
        lam       ! the angualar momentum

real(kind=dp) :: &
        f1aep1,f1aem1,jnor,psnor,fact(4), occ,  &
        gi(ndm),j1(ndm,4),cm(10),bm(4),ze2, &
        xc(8),int_0_inf_dr,delta,chir(ndm,nwfx), &
        a,b,c,deter,gamma, &
        lamda0,lamda3,lamda4,mu0,mu4,s0,s4,t0,t4
     
real(kind=dp) :: &
        deriv_7pts, deriv2_7pts, rab(ndm), chi_dir(ndm,2)

integer :: &
        m, n, nwf0, nst, nnode, nc, nc1, ij, imax, iq, i

logical  ::  &
        lbes4  !decide if use 3 or 4 Bessel function

!
!   decide wether use 4 Bessel functions or 3 (default)
!
lbes4=rho0.ne.0.d0.and.lam.eq.0
!
!   compute the pseudowavefunctions by expansion in spherical
!   bessel function before r_c
!
ff=1.d0-dx**2/48.d0
ze2=-zed*2.d0
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
      call intref(lam,enls(ns),mesh,dx,r,r2,sqr,vpot,ze2,chir(1,ns))
   endif
!
!    fix arbitrarely the norm at the cut-off radius equal to 0.5
!
   jnor=chir(ik,ns)
   do n=1,mesh
      chir(n,ns)=chir(n,ns)*0.5d0/jnor
   enddo
else 
   do n=1,mesh
     chir(n,ns)=psi(n,nwf0)
   enddo
!            if (lam.eq.3) stop
endif
!
!   compute the first and second derivative of all-electron function
!
fae=chir(ik,ns)
f1ae=deriv_7pts(chir(1,ns),ik,r(ik),dx)
f2ae=deriv2_7pts(chir(1,ns),ik,r(ik),dx)
!
!   computes the norm of the all-electron function
!
do n=1,ik+1
   gi(n)=chir(n,ns)**2  
enddo
faenor=int_0_inf_dr(gi,r,r2,dx,ik,nst)
!
!   find q_i with the correct log derivatives
!
if (lbes4)then 
   call find_qi(f1ae/fae,xc(5),ik,lam,4,1,iok)
   if (iok.ne.0) return
else 
   call find_qi(f1ae/fae,xc(4),ik,lam,3,1,iok)
   if (iok.ne.0) return
endif
!
!    compute the bessel functions
!
if (lbes4)then
  do nc=1,4
     call sph_besr(ik+5,r,xc(4+nc),lam,j1(1,nc))
     jnor=j1(ik,nc)
     fact(nc)=chir(ik,ns)/jnor
     do n=1,ik+5
        j1(n,nc)=j1(n,nc)*chir(ik,ns)/jnor
     enddo
  enddo
else
  do nc=1,3
    call sph_besr(ik+5,r,xc(3+nc),lam,j1(1,nc))
    jnor=j1(ik,nc)
    fact(nc)=chir(ik,ns)/jnor
    do n=1,ik+5
       j1(n,nc)=j1(n,nc)*chir(ik,ns)/jnor
    enddo
  enddo
endif
!
!    compute the bm functions (second derivative of the j1)
!    and the integrals for cm
!
if (lbes4)then
   ij=0
   do nc=1,4
      bm(nc)=deriv2_7pts(j1(1,nc),ik,r(ik),dx)
      do nc1=1,nc
         ij=ij+1
         do n=1,ik
            gi(n)=j1(n,nc)*j1(n,nc1)
         enddo
         cm(ij)=int_0_inf_dr(gi,r,r2,dx,ik,nst)
      enddo
   enddo
else
   ij=0
   do nc=1,3
      bm(nc)=deriv2_7pts(j1(1,nc),ik,r(ik),dx)
      do nc1=1,nc
         ij=ij+1
         do n=1,ik
            gi(n)=j1(n,nc)*j1(n,nc1)
         enddo
         cm(ij)=int_0_inf_dr(gi,r,r2,dx,ik,nst)
      enddo
   enddo
endif
!
!    solve the second order equation. See notes
!
if (lbes4) then
   wmax=0.d0
   do n=1,mesh
      if (abs(chir(n,ns)).gt.wmax.and.r(n).lt.4.d0) then
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
   if (deter.lt.0.d0) then
      call errore('compute phi','negative determinant',-1) 
      write(6,110) ns,f1ae/fae, f2ae, faenor

 110    format(/5x,i5,' ld= ',f10.6,' f2ae',f10.6,' faenor',f10.6)

!           do n=1,mesh
!              write(6,*) r(n),chir(n,ns)
!           enddo
!           stop  
      iok=1
      return
   endif  

   xc(4)=(-b+sqrt(deter))/(2.d0*a)
   xc(3)=mu0-mu4*xc(4)
   xc(2)=s0-s4*xc(4)
   xc(1)=t0-t4*xc(4)
!
!      check for the norm of the pseudo wavefunction
!
   do n=1,ik
      gi(n)=(xc(1)*j1(n,1)+xc(2)*j1(n,2) &
            +xc(3)*j1(n,3)+xc(4)*j1(n,4))**2
   enddo
   psnor=int_0_inf_dr(gi,r,r2,dx,ik,nst)
   if (iflag.eq.1) then
      write(6,120) els(ns),rcut(ns),2.d0*xc(6)**2 
120     format(/ /5x, ' Wfc  ',a3,'  rcut=',f6.3, &
                '  Estimated cut-off energy= ', f11.2,' Ry')
      write(6,130) rho0
130     format (5x,' Using 4 Bessel functions for this wfc,', &
              ' rho(0) =',f6.3)
!          write(6,'(5x,'' AE norm before r_c= '',f12.9)') faenor
!          write(6,'(5x,'' PS norm before r_c= '',f12.9)') psnor
   endif

   do n=1,ik
      phis(n,ns)=xc(1)*j1(n,1)+xc(2)*j1(n,2)  &
                +xc(3)*j1(n,3)+xc(4)*j1(n,4)
   enddo
   do n=ik+1,mesh
      phis(n,ns)= chir(n,ns)
   enddo
   nnode=0  
   do n=1,ik+1
      if ( phis(n,ns) .ne. sign(phis(n,ns),phis(n+1,ns)) ) then
         write(6,140) lam,ns,r(n)
140      format (5x,'l=',i4,' ns=',i4,' Node at ',f10.8)
         nnode=nnode+1
      endif
   enddo
   iok=nnode
   if (iflag.eq.1) &
      write(6,150) nnode,r(ik)
150   format (5x,' This function has ',i4,' nodes', &
                 ' between 0 and r=',f8.3)

      do nc=1,4
         xc(nc)=xc(nc)*fact(nc)
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
   if (deter.lt.0.d0) then
      call errore('compute phi','negative determinant',-1) 
      write(6,160) ns,f1ae/fae, f2ae, faenor
160   format (/5x,i5,' ld= ',f10.6,' f2ae',f10.6, ' faenor',f10.6)
      iok=1
      return
   endif  

   xc(3)=(-b+sqrt(deter))/(2.d0*a)
   xc(2)=-xc(3)*gamma+delta
   xc(1)=1.d0-xc(2)-xc(3)
!
!      check for the norm of the pseudo wavefunction
!
   do n=1,ik
      gi(n)=(xc(1)*j1(n,1)+xc(2)*j1(n,2)+xc(3)*j1(n,3))**2
   enddo
   psnor=int_0_inf_dr(gi,r,r2,dx,ik,nst)
   if (iflag.eq.1) then
      write(6,170) els(ns),rcut(ns),2.d0*xc(6)**2 
170   format (/ /5x, ' Wfc  ',a3,'  rcut=',f6.3, &
                '  Estimated cut-off energy= ', f11.2,' Ry')
!          write(6,'(5x,'' AE norm before r_c= '',f12.9)') faenor
!          write(6,'(5x,'' PS norm before r_c= '',f12.9)') psnor
   endif

   do n=1,ik
      phis(n,ns)= xc(1)*j1(n,1) + xc(2)*j1(n,2) &
                                + xc(3)*j1(n,3)
   enddo

   do n=ik+1,mesh
      phis(n,ns)= chir(n,ns)
   enddo

   nnode=0  
   do n=1,ik+1
      if ( phis(n,ns) .ne. sign(phis(n,ns),phis(n+1,ns)) ) then
             write(6,'(5x,''l='',i4,'' ns='',i4,'' Node at '',f10.8)') &
                 lam,ns,r(n)
         nnode=nnode+1
      endif
   enddo
   iok=nnode
   if (iflag.eq.1)  &
         write(6,180) nnode,r(ik)
180      format (5x,' This function has ',i4,' nodes', &
                      ' between 0 and r=',f8.3) 

   do nc=1,3
      xc(nc)=xc(nc)*fact(nc)
   enddo
endif

return
end
