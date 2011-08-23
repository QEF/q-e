!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine ascheq(nn,lam,e,mesh,grid,vpot,ze2,thresh0,y,nstop)
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation for
  !  bound states in a local potential.
  !  thresh determines the absolute accuracy for the eigenvalue
  !
  use kinds, only : DP
  use radial_grids, only: radial_grid_type, series
  implicit none
  type(radial_grid_type), intent(in) :: grid
  integer :: mesh,lam, ierr
  integer:: nn,nstop,maxter,iter,l1,i,ik,ncross,n, &
       nstart,ns,n2,nst2,ndcr
  real(DP) :: ze2,ddx12,eup,elw,b0e,ymx,rap,rstart,di,expn,  &
       c1,c2,fe,a0,a1,a2,sum0,f2,sum,sqlhf,f0,f1,dfe,de,eps, &
       yln,xp,sum1,xl1,x4l6
  real(DP):: vpot(mesh), y(mesh)
  real(DP),allocatable:: c(:), el(:), f(:)
  real(DP):: b(0:3),e,thresh0, thresh
  data maxter/50/
  !
  !  set up constants and initialize
  !
  allocate(c(mesh),stat=ierr)
  allocate(f(mesh),stat=ierr)
  allocate(el(mesh),stat=ierr)

  thresh=thresh0
  if (e<-5.e+2) thresh=thresh0*10.0_DP
  iter=0
  ddx12=grid%dx*grid%dx/12.0_dp
  l1=lam+1
  sqlhf=(DBLE(lam)+0.5_dp)**2
  ndcr=nn-lam-1
  !
  !  set initial lower and upper bounds to the eigenvalue
  !
  eup=vpot(mesh)+sqlhf/grid%r2(mesh)
  elw=eup
  do i=1,mesh
     elw=min(elw,vpot(i)+sqlhf/grid%r2(i))
  enddo
  nstop=200
  if(eup.eq.elw) go to 900
  if(e.gt.eup) e=0.9_DP*eup+0.1_DP*elw
  if(e.lt.elw) e=0.9_DP*elw+0.1_DP*eup
  !
  !  series developement of the potential near the origin
  !
  do i=1,4
     y(i)=vpot(i)-ze2/grid%r(i)
  enddo
  call series(y,grid%r,grid%r2,b)
  !
300 continue
  iter=iter+1
  nstop=300
  if(iter.gt.maxter) go to 900
  !
  !  set up the f-function and determine the position of its last
  !  change of sign
  !  f < 0 (approximatively) means classically allowed   region
  !  f > 0         "           "        "      forbidden   "
  !
  f(1)=ddx12*(grid%r2(1)*(vpot(1)-e)+sqlhf)
  do i=2,mesh
     f(i)=ddx12*(grid%r2(i)*(vpot(i)-e)+sqlhf)
     if( f(i) .ne. sign(f(i),f(i-1)) ) ik=i
  enddo
  nstop=302
  
  if(ik.ge.mesh-2) go to 900
  do i=1,mesh
     f(i)=1.0_dp-f(i)
  enddo
  !
  y(:) = 0.0_dp
  !
  !  determination of the wave-function in the first two points by
  !  series developement
  !
  xl1=lam+1.0_DP
  x4l6=4.0_dp*lam+6.0_dp
  b0e=b(0)-e
  c1=0.5_dp*ze2/xl1
  c2=(c1*ze2+b0e)/x4l6
  call start_scheq( lam, e, b, grid, ze2, y )
  !
  !  start outward integration and count number of crossings
  !
  ncross=0
  ymx=0.0_dp
  do n=2,ik-1
     y(n+1)=((12.0_dp-10.0_dp*f(n))*y(n)-f(n-1)*y(n-1))/f(n+1)
     if ( y(n) .ne. sign(y(n),y(n+1)) ) ncross=ncross+1
     ymx=max(ymx,abs(y(n+1)))
  end do
  !
  !  matching radius has been reached going out. if ncross is not
  !  equal to ndcr, modify the trial eigenvalue.
  !
  if(ndcr < ncross) then
     !
     !  too many crossings. e is an upper bound to the true eigen-
     !  value. increase abs(e)
     !
     eup=e
     rap=(DBLE(ncross+l1)/DBLE(nn))**2
     e=(e-vpot(mesh))*rap+vpot(mesh)
     if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup
     go to 300
  else if (ndcr > ncross) then
     !
     !  too few crossings. e is a lower bound to the true eigen-
     !  value. decrease abs(e)
     !
     elw=e
     rap=(DBLE(ncross+l1)/DBLE(nn))**2
     e=(e-vpot(mesh))*rap+vpot(mesh)
     if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
     go to 300
  end if
  !
  !  prepare inward integration
  !  charlotte froese can j phys 41,1895(1963)
  !
  !            start at  min( rmax, 10*rmatch )
  !
  nstart=mesh
  ns=10
  rstart=ns*grid%r(ik)
  if(rstart.lt.grid%r(mesh)) then
     do  i=ik,mesh
        nstart=i
        if(grid%r(i).ge.rstart) go to 403
     enddo
403  nstart=nstart/2
     nstart=2*nstart+1
  end if
  !
  !  set up a, l, and c vectors
  !
  n=ik+1
  el(n)=10.0_dp*f(n)-12.0_dp
  c(n)=-f(ik)*y(ik)
  n2=ik+2
  do n=n2,nstart
     di=10.0_dp*f(n)-12.0_dp
     el(n)=di-f(n)*f(n-1)/el(n-1)
     c(n)=-c(n-1)*f(n-1)/el(n-1)
  enddo
  !
  !  start inward integration by the froese's tail procedure
  !
  expn=exp(-sqrt(12.0_dp*abs(1.0_dp-f(nstart-1))))
  y(nstart-1)=c(nstart-1)/(el(nstart-1)+f(nstart)*expn)
  y(nstart)=expn*y(nstart-1)
  do n=nstart-2,ik+1,-1
    y(n)=(c(n)-f(n+1)*y(n+1))/el(n)
 end do
  !
  !  if necessary, improve the trial eigenvalue by the cooley's
  !  procedure. jw cooley math of comp 15,363(1961)
  !
  fe=(12.0_dp-10.0_dp*f(ik))*y(ik)-f(ik-1)*y(ik-1)-f(ik+1)*y(ik+1)
  !
  !  calculate the normalization
  !
  if(ymx.ge.1.0e10_dp) then
     do  i=1,mesh
        y(i)=y(i)/ymx
     enddo
  end if
  a0=1.0_dp/DBLE(2*lam+3)
  a1=c1/DBLE(lam+2)
  a2=(c1*c1+c2+c2)/DBLE(2*lam+5)
  sum0=(a0+grid%r(1)*(a1+grid%r(1)*a2))*grid%r(1)**(2*lam+3)
  nst2=nstart-2
  f2=grid%r2(1  )*y(1  )*y(1  )
  sum=grid%r(1)*f2/DBLE(2*l1+1)
  do n=1,nst2,2
     f0=f2
     f1=grid%r2(n+1)*y(n+1)*y(n+1)
     f2=grid%r2(n+2)*y(n+2)*y(n+2)
     sum=sum+f0+f2+4.0_DP*f1
  enddo
  sum=sum0+grid%dx*sum/3.0_dp
  dfe=-y(ik)*f(ik)/grid%dx/sum
  de=-fe*dfe
  eps=abs(de/e)
  if(abs(de).lt.thresh) go to 600
  if(eps.gt.0.25_dp) de=0.25_dp*de/eps
  if(de.gt.0.0_dp) elw=e
  if(de.lt.0.0_dp) eup=e
  e=e+de
  if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
  if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup
  if(iter.lt.maxter) go to 300
  nstop=50
600 continue
  !
  !  normalize the eigenfunction and exit
  !
  do n=nstart,mesh-1
     y(n+1)=0.0_dp
     if(y(n).eq.0.0_dp) go to 601
     yln=log(abs(y(n)))
     xp=-sqrt(12.0_dp*abs(1.0_dp-f(n)))
     expn=yln+xp
     if(expn.lt.-80.0_dp) go to 601
     y(n+1)=sign(exp(expn),y(n))
601 continue
  enddo
  sum1=0.0_dp
  do n=nstart,mesh-2,2
     f0=f2
     f1=grid%r2(n+1)*y(n+1)*y(n+1)
     f2=grid%r2(n+2)*y(n+2)*y(n+2)
     sum1=sum1+f0+f2+4.0_dp*f1
  enddo
  sum=sum+grid%dx*sum1/3.0_dp
  sum=sqrt(sum)
  do n=1,mesh
     y(n)=grid%sqr(n)*y(n)/sum
  enddo
  if(nstop.lt.100) go to 900
  nstop=0
  deallocate(el)
  deallocate(f )
  deallocate(c )
  return
  !
  !  error exit
  !
  ! 900  write(6,9000) nstop,nn,lam,elw,eup
  ! 9000 format(5x,'error in ascheq: nstop =',i4,'. n l =',2i3,/ &
  !     & 5x,'elw =',f15.10,' eup =',f15.10)
900 continue
  deallocate(el)
  deallocate(f )
  deallocate(c )
  return

end subroutine ascheq
