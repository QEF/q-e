!
!--------------------------------------------------------------------------
      subroutine compute_phipot(lam,ik,nwf0,ns,xc)
!--------------------------------------------------------------------------
!
!     This routine computes the phi functions by pseusizing the
!     all_electron chi functions. In input it receives, the point
!     ik where the cut is done, the angualar momentum lam,
!     and the correspondence with the all eletron wavefunction
!
!
!      
use ld1inc
implicit none

      real(kind=dp) :: &
               fae,    & ! the value of the all-electron function
               ff,     & ! compute the second derivative
               signo,  & ! the sign of the ae wfc
               wmax,   & !
               den,    & ! denominator
               faenor   ! the norm of the function
      
      real(kind=dp), parameter :: pi=3.14159265358979d0

      integer :: &
               ik,i,  &    ! the point corresponding to rc
               ns,    &  ! the function to pseudize
               lam      ! the angualar momentum

      real(kind=dp) :: &
               jnor,psnor,fact(4), f2aep,f2aem,f3ae, &
               gi(ndm),j1(ndm,4),cm(10),bm(4),ze2,c(6),c2, &
               xc(8),int_0_inf_dr,delta,chir(ndm,nwfx),pr, &
               dpoly, d2pr, dpr,  &
               lamda0,lamda3,lamda4,mu0,mu4,s0,s4,t0,t4, rab(ndm), &
               chi_dir(ndm,2)
     
      integer :: &
               m, n, nwf0, nst, nnode, nc, nc1, ij, imax, iq

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
            call lschps(3,zed,exp(dx),dx,mesh,mesh,mesh, &
                        1,lam,enls(ns),chir(1,ns),r,vpot)
         elseif (rel.eq.2) then
            do i=1,mesh
               rab(i)=r(i)*dx
            enddo
            call dir_outward(ndm,mesh,lam,jjs(ns),enls(ns),dx, &
                        chi_dir,r,rab,vpot)
            chir(:,ns)=chi_dir(:,2)
         else
            call intref(lam,enls(ns),mesh,dx,r,r2,sqr, &
                        vpot,ze2,chir(1,ns))
         endif
!
!    fix arbitrarely the norm at the cut-off radius equal to 0.5
!
         jnor=chir(ik,ns)
         do n=1,mesh
            chir(n,ns)=chir(n,ns)*0.5d0/jnor
!            write(6,*) r(n),chir(n,ns)
         enddo
      else 
         do n=1,mesh
            chir(n,ns)=psi(n,nwf0)
         enddo
      endif

      do n=1,ik+1
         gi(n)=chir(n,ns)**2
      enddo
      faenor=int_0_inf_dr(gi,r,r2,dx,ik,nst)

      call find_coefficients &
           (ik,chir(1,ns),enls(ns),r,dx,faenor,vpot,c,c2,lam,mesh)

!
!   pseudowavefunction found
!
      signo= 1.d0 !chir(ik+1,ns)/abs(chir(ik+1,ns))
      do i=1,ik
         phis(i,ns)=signo*r(i)**(lam+1)*exp(pr(c,c2,r(i)))
      end do
      do i=ik+1,mesh
         phis(i,ns)=chir(i,ns)
      end do
      do i=1,ik
         dpoly = dpr(c,c2,r(i))    ! first derivate of polynom
         chis(i,ns) = (enls(ns) + (2*lam+2)/r(i)*dpoly + &
                     d2pr(c,c2,r(i)) + dpoly**2)*phis(i,ns)
      enddo
      do i = ik+1,mesh
         chis(i,ns) = vpot(i,1)*phis(i,ns)
      enddo

!      write(38,*) 'call for l=',lam,' jjs=',jjs(ns)
!      do i=1,mesh
!         write(38,*) r(i), phis(i,ns), chir(i,ns)
!      enddo
!      stop 

      return
      end
