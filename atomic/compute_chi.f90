!
!--------------------------------------------------------------------------
subroutine compute_chi(lam,ns,xc,lbes4)
  !--------------------------------------------------------------------------
  !
  !     This routine computes the chi functions by inversion of
  !     the schroedinger equation with the RRKJ3 method.
  !      
  use ld1inc

  integer :: &
       ns, &    ! the wavefunction
       lam      ! the angular momentum
  logical :: &
       lbes4 
  real(kind=dp) :: &
       xc(8)
  !
  real(kind=dp) :: &
       j1(ndm),aux(ndm), &
       b(4),c(4), arow(ndm),brow(ndm),crow(ndm),drow(ndm), &
       b0e, g0, g1, g2, &
       ddx12, &
       x4l6, &
       x6l12, dpoly
  real(kind=dp), external :: pr, d2pr, dpr

  integer :: &
       n, nstart

  !
  if (tm) then
     do i=1,ik
        dpoly = dpr(xc,xc(7),r(i))    ! first derivate of polynomial
        chis(i,ns) = (enls(ns) + (2*lam+2)/r(i)*dpoly + &
             d2pr(xc,xc(7),r(i)) + dpoly**2)*phis(i,ns)
     enddo
     do i = ik+1,mesh
        chis(i,ns) = vpot(i,1)*phis(i,ns)
     enddo
     return
  end if
  !
  !   compute the projectors by inverting the schroedinger equation
  !
  !   first expand in a taylor series the phis function
  !   Since we know that the phis functions are a sum of
  !   bessel function with coeeficients xc, we compute analiticaly
  !   the asymptotic expansion
  !
  !
  ddx12=dx*dx/12.d0
  x4l6=4*lam+6
  x6l12=6*lam+12

  do n=1,6
     j1(n)=phis(n,ns)/r(n)**(lam+1)
  enddo
  call seriesbes(j1,r,r2,6,c)
  !      write(6,'(''phi(serie)= '',4e15.7)') c

  if (lam.eq.0) then
     if(lbes4.or.rho0.eq.0.d0)then
        c(1)=xc(1)+  &
             xc(2)+  &
             xc(3)
        c(2)=0.d0
        c(3)=-xc(1)*(xc(4)**2/6.d0) &
             -xc(2)*(xc(5)**2/6.d0) &
             -xc(3)*(xc(6)**2/6.d0)
        c(4)=0.d0
     else
        c(1)=xc(1)+  &
             xc(2)+  &
             xc(3)+xc(4)
        c(2)=0.d0
        c(3)=-xc(1)*(xc(5)**2/6.d0)  &
             -xc(2)*(xc(6)**2/6.d0)  &
             -xc(3)*(xc(7)**2/6.d0)  &
             -xc(4)*(xc(8)**2/6.d0)
        c(4)=0.d0
     endif
  elseif (lam.eq.3) then
     c(1)=xc(1)*(48.d0*xc(4)**3/5040.d0)+   &
          xc(2)*(48.d0*xc(5)**3/5040.d0)+   &
          xc(3)*(48.d0*xc(6)**3/5040.d0)
     c(2)=0.d0
     c(3)=-xc(1)*(192.d0*xc(4)**5/362880.d0)  &
          -xc(2)*(192.d0*xc(5)**5/362880.d0)  &
          -xc(3)*(192.d0*xc(6)**5/362880.d0)
     c(4)=0.d0
  elseif (lam.eq.2) then
     c(1)=xc(1)*(xc(4)**2/15.d0)+   &
          xc(2)*(xc(5)**2/15.d0)+   &
          xc(3)*(xc(6)**2/15.d0)
     c(2)=0.d0
     c(3)=-xc(1)*(xc(4)**4/210.d0)  &
          -xc(2)*(xc(5)**4/210.d0)  &
          -xc(3)*(xc(6)**4/210.d0)
     c(4)=0.d0
  elseif (lam.eq.1) then
     c(1)=xc(1)*(xc(4)/3.d0)+  &
          xc(2)*(xc(5)/3.d0)+  &
          xc(3)*(xc(6)/3.d0)
     c(2)=0.d0
     c(3)=-xc(1)*(xc(4)**3/30.d0) &
          -xc(2)*(xc(5)**3/30.d0) &
          -xc(3)*(xc(6)**3/30.d0)
     c(4)=0.d0
  else
     call errore('compute_chi','lam not programmed',1) 
  endif
  !      write(6,'(''phi(compu)= '',4e15.7)') c
  !
  !     and the potential
  !
  do n=1,4
     j1(n)=vpsloc(n)
  enddo
  call series(j1,r,r2,b)
  !
  !   and compute the taylor expansion of the chis
  !
  b0e=(b(1)-enls(ns))

  g0=x4l6*c(3)-b0e*c(1)
  g1=x6l12*c(4)-c(1)*b(2)
  g2=-(b0e*c(3)+b(3)*c(1))
  nstart=5
  do n=1,nstart-1
     chis(n,ns)= (g0+r(n)*(g1+g2*r(n)))*r(n)**(lam+3)/sqr(n)
  enddo
  do n=1,mesh
     aux(n)= (g0+r(n)*(g1+g2*r(n)))
  enddo
  !
  !    set up the equation
  !
  do n=1,mesh
     phis(n,ns)=phis(n,ns)/sqr(n)
  enddo
  do n=1,mesh
     j1(n)=r2(n)*(vpsloc(n)-enls(ns))+(lam+0.5d0)**2
     j1(n)=1.d0-ddx12*j1(n)
  enddo

  do n=nstart,mesh-3
     drow(n)= phis(n+1,ns)*j1(n+1)   &
          + phis(n,ns)*(-12.d0+10.d0*j1(n))+ &
          phis(n-1,ns)*j1(n-1)

     brow(n)=10.d0*ddx12
     crow(n)=ddx12
     arow(n)=ddx12
  enddo
  drow(nstart)=drow(nstart)-ddx12*chis(nstart-1,ns)
  chis(mesh-2,ns)=0.d0
  chis(mesh-1,ns)=0.d0
  chis(mesh,ns)=0.d0
  !
  !    and solve it
  !
  call tridiag(arow(nstart),brow(nstart),crow(nstart), &
       drow(nstart),chis(nstart,ns),mesh-3-nstart)
  !
  !   put the correct normalization and r dependence
  !  
  do n=1,mesh
     phis(n,ns)=phis(n,ns)*sqr(n)
     chis(n,ns)=chis(n,ns)*sqr(n)/r2(n)
     !         if(lam.eq.0)
     !     +    write(*,'(5(e20.13,1x))')
     !     +          r(n),chis(n,ns),chis(n,ns)/r(n)**(lam+1),
     !     +          aux(n),aux(n)*r(n)**(lam+1)
  enddo
  !
  !    smooth close to the origin with asymptotic expansion
  !
  do n=nstart,mesh
     if (abs(chis(n,ns)/r(n)**(lam+1)-aux(n))  &
          .lt.1.d-3*abs(aux(n)) ) goto 100
     chis(n,ns)=aux(n)*r(n)**(lam+1)
  enddo

100 if (n.eq.mesh+1.or.r(n).gt.0.05d0)then
     print*,lam,ns,n,mesh,r(n)
     call errore('compute_chi','n is too large',1)
  endif
  !
  !    clean also after 9 a.u.
  !
  do n=mesh,1,-1
     if (r(n).lt.9.d0) goto 200
     chis(n,ns)=0.d0
  enddo
200 continue
  return
end subroutine compute_chi


subroutine tridiag(a,b,c,r,u,n)
  !
  !     See Numerical Recipies.
  !
  implicit none

  integer :: n
  integer, parameter :: nmax=5000, dp=kind(1.d0)

  real(kind=dp) :: a(n),b(n),c(n),r(n),u(n), gam(nmax), bet

  integer j
  if (n.gt.nmax)  &
       call errore('tridiag','nmax is too small',1)

  if (abs(b(1)).lt.1.d-10)  &
       call errore('tridiag','b(1) is too small',1)

  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     if (abs(bet).lt.1.d-10) &
          call errore('tridiag','bet is too small',1)

     u(j)=(r(j)-a(j)*u(j-1))/bet
  enddo
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  enddo
  return
end subroutine tridiag
