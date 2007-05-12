!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
subroutine compute_chi(lam,ik,ikk_in,phi_in,chi_out,xc,e,lbes4)
  !--------------------------------------------------------------------------
  !
  !     This routine computes the chi functions:
  !          |chi> = (\epsilon -T -V_{loc)) |psi>
  !      
  use kinds, only : DP
  use ld1inc

  implicit none
  integer :: &
       ik,    & ! the point corresponding to rc
       ikk_in,& ! the point after which the chi should be zero
       ns,    & ! the wavefunction
       lam      ! the angular momentum
  logical :: &
       lbes4 

  real(DP) :: &
       e,     &       ! input: the energy
       xc(8), &       ! input: the coefficients of the bessel function
       phi_in(ndm), & ! input: pseudo wavefunction
       chi_out(ndm)   ! output: the chi function
  !
  real(DP) :: &
       j1(ndm), aux(ndm), gi(ndm),&
       b(4),c(4), arow(ndm),brow(ndm),crow(ndm),drow(ndm), &
       b0e, g0, g1, g2, &
       ddx12,           &
       x4l6,            &
       x6l12,           &
       int_0_inf_dr,    &
       integral

  integer :: &
       n, nstart, nst
  !
  !   RRKJ: first expand in a taylor series the phis function
  !   Since we know that the phis functions are a sum of Bessel 
  !   functions with coefficients xc, we compute analytically
  !   the asymptotic expansion
  !
  !
  ddx12=dx*dx/12.0_dp
  x4l6=4*lam+6
  x6l12=6*lam+12

  do n=1,6
     j1(n)=phi_in(n)/r(n)**(lam+1)
  enddo
  call seriesbes(j1,r,r2,6,c)
  !
  if (lam == 0) then
     if(lbes4.or.rho0.eq.0.0_dp)then
        c(1)=xc(1)+xc(2)+xc(3)
        c(2)=0.0_dp
        c(3)=-xc(1)*(xc(4)**2/6.0_dp) &
             -xc(2)*(xc(5)**2/6.0_dp) &
             -xc(3)*(xc(6)**2/6.0_dp)
        c(4)=0.0_dp
     else
        c(1)=xc(1)+xc(2)+xc(3)+xc(4)
        c(2)=0.0_dp
        c(3)=-xc(1)*(xc(5)**2/6.0_dp)  &
             -xc(2)*(xc(6)**2/6.0_dp)  &
             -xc(3)*(xc(7)**2/6.0_dp)  &
             -xc(4)*(xc(8)**2/6.0_dp)
        c(4)=0.0_dp
     endif
  elseif (lam == 3) then
     c(1)=xc(1)*(48.0_dp*xc(4)**3/5040.0_dp)+   &
          xc(2)*(48.0_dp*xc(5)**3/5040.0_dp)+   &
          xc(3)*(48.0_dp*xc(6)**3/5040.0_dp)
     c(2)=0.0_dp
     c(3)=-xc(1)*(192.0_dp*xc(4)**5/362880.0_dp)  &
          -xc(2)*(192.0_dp*xc(5)**5/362880.0_dp)  &
          -xc(3)*(192.0_dp*xc(6)**5/362880.0_dp)
     c(4)=0.0_dp
  elseif (lam == 2) then
     c(1)=xc(1)*(xc(4)**2/15.0_dp)+   &
          xc(2)*(xc(5)**2/15.0_dp)+   &
          xc(3)*(xc(6)**2/15.0_dp)
     c(2)=0.0_dp
     c(3)=-xc(1)*(xc(4)**4/210.0_dp)  &
          -xc(2)*(xc(5)**4/210.0_dp)  &
          -xc(3)*(xc(6)**4/210.0_dp)
     c(4)=0.0_dp
  elseif (lam == 1) then
     c(1)=xc(1)*(xc(4)/3.0_dp)+  &
          xc(2)*(xc(5)/3.0_dp)+  &
          xc(3)*(xc(6)/3.0_dp)
     c(2)=0.0_dp
     c(3)=-xc(1)*(xc(4)**3/30.0_dp) &
          -xc(2)*(xc(5)**3/30.0_dp) &
          -xc(3)*(xc(6)**3/30.0_dp)
     c(4)=0.0_dp
  else
     call errore('compute_chi','lam not programmed',1) 
  endif
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
  b0e=(b(1)-e)

  g0=x4l6*c(3)-b0e*c(1)
  g1=x6l12*c(4)-c(1)*b(2)
  g2=-(b0e*c(3)+b(3)*c(1))
  nstart=5
  do n=1,nstart-1
     chi_out(n)= (g0+r(n)*(g1+g2*r(n)))*r(n)**(lam+3)/sqr(n)
  enddo
  do n=1,mesh
     aux(n)= (g0+r(n)*(g1+g2*r(n)))
  enddo
  !
  !    set up the equation
  !
  do n=1,mesh
     gi(n)=phi_in(n)/sqr(n)
  enddo
  do n=1,mesh
     j1(n)=r2(n)*(vpsloc(n)-e)+(lam+0.5_dp)**2
     j1(n)=1.0_dp-ddx12*j1(n)
  enddo

  do n=nstart,mesh-3
     drow(n)= gi(n+1)*j1(n+1)   &
            + gi(n)*(-12.0_dp+10.0_dp*j1(n))+ &
              gi(n-1)*j1(n-1)

     brow(n)=10.0_dp*ddx12
     crow(n)=ddx12
     arow(n)=ddx12
  enddo
  drow(nstart)=drow(nstart)-ddx12*chi_out(nstart-1)
  chi_out(mesh-2)=0.0_dp
  chi_out(mesh-1)=0.0_dp
  chi_out(mesh)=0.0_dp
  !
  !    and solve it
  !
  call tridiag(arow(nstart),brow(nstart),crow(nstart), &
       drow(nstart),chi_out(nstart),mesh-3-nstart)
  !
  !   put the correct normalization and r dependence
  !  
  do n=1,mesh
     chi_out(n)=chi_out(n)*sqr(n)/r2(n)
     !         if(lam.eq.0)
     !     +    write(*,'(5(e20.13,1x))')
     !     +          r(n),chi_out(n),chi_out(n)/r(n)**(lam+1),
     !     +          aux(n),aux(n)*r(n)**(lam+1)
  enddo
  !
  !    smooth close to the origin with asymptotic expansion
  !
  do n=nstart,mesh
     if (abs(chi_out(n)/r(n)**(lam+1)-aux(n))  &
          .lt.1.e-3_dp*abs(aux(n)) ) goto 100
     chi_out(n)=aux(n)*r(n)**(lam+1)
  enddo

100 if (n.eq.mesh+1.or.r(n).gt.0.05_dp)then
     print*,lam,n,mesh,r(n)
     call errore('compute_chi','n is too large',1)
  endif
  !
  !    clean also after 9 a.u.
  !
  do n=mesh,1,-1
     if (r(n).lt.9.0_dp) goto 200
     chi_out(n)=0.0_dp
  enddo
200 continue
     !    check that the chi are zero beyond ikk
  nst=0
  gi=0.0_dp
  do n=ikk_in+1,mesh
     gi(n)=chi_out(n)**2
  enddo
  do n=min(ikk_in+20,mesh),mesh
     chi_out(n)=0.0_dp
  enddo
  integral=int_0_inf_dr(gi,r,r2,dx,mesh,nst)
  if (integral > 2.e-6_dp) then
     write(6, '(5x,''ns='',i4,'' l='',i4, '' integral='',f15.9, &
          & '' r(ikk) '',f15.9)') ns, lam, integral, r(ikk_in)
     call infomsg ('gener_pseudo ','chi too large beyond r_c', -1)
     do n=ikk_in,mesh
        write(6,*) r(n),gi(n)
     enddo
     stop
  endif
  return
end subroutine compute_chi


subroutine tridiag(a,b,c,r,u,n)
  !
  !     See Numerical Recipes.
  !
  use kinds, only : DP
  implicit none

  integer :: n
  real(DP) :: a(n),b(n),c(n),r(n),u(n)
  real(DP) :: gam(n), bet

  integer j

  if (abs(b(1)).lt.1.e-10_DP)  &
       call errore('tridiag','b(1) is too small',1)

  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     if (abs(bet) < 1.e-10_DP) &
          call errore('tridiag','bet is too small',1)
     u(j)=(r(j)-a(j)*u(j-1))/bet
  enddo
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  enddo
  return
end subroutine tridiag
