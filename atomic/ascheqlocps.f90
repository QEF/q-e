!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine ascheqlocps(nn,lam,e,mesh,ndm,dx,r,r2,sqr,vpot, &
     thresh,y)
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation for
  !  bound states in a local potential.
  !  thresh dermines the absolute accuracy for the eigenvalue
  !
  use kinds, only : DP
  implicit none
  integer :: &
       nn, &       ! main quantum number for node number
       lam,&       ! angular momentum
       mesh, &     ! size of radial mesh
       ndm         ! maximum radial mesh 

  real(kind=dp) :: &
       e,       &  ! output eigenvalue
       dx,      &  ! linear delta x for radial mesh
       r(mesh), &  ! radial mesh
       r2(mesh),&  ! square of radial mesh
       sqr(mesh),& ! square root of radial mesh
       vpot(mesh),&! the local potential 
       thresh,   & ! precision of eigenvalue
       y(mesh)     ! the output solution
  !
  !    the local variables
  !
  integer :: maxter=100       ! maximum number of iterations

  real(kind=dp) :: &
       ddx12,      &   ! dx**2/12 used for Numerov integration
       sqlhf,      &   ! the term for angular momentum in equation
       xl1, x4l6, x6l12, x8l20,& ! used for starting series expansion
       ze2,    &       ! possible coulomb term aroun the origin (set 0)
       b(0:3), &      ! coefficients of taylor expansion of potential
       c1,c2,c3,c4,b0e, & ! auxiliary for expansion of wavefunction
       rr1,rr2, &    ! values of y in the first points
       ymx,      &   ! the maximum value of the function
       eup,elw,  &   ! actual energy interval
       rap,      &   ! the ratio between the number of nodes
       fe,sum,dfe,de, &! auxiliary for numerov computation of e
       eps,        & ! the epsilon of the delta e
       yln, xp, expn,& ! used to compute the tail of the solution
       int_0_inf_dr   ! integral function

  real(kind=dp),allocatable :: &
       f(:), &       ! the f function
       el(:),c(:),&  ! auxiliary for inward integration
       yi(:)         ! the irregular solution


  integer :: &
       n, &     ! counter on mesh points
       iter, &  ! counter on iteration
       ik,   &  ! matching point
       ns,   &  ! counter on beta functions
       l1,   &  ! lam+1
       nst,  &  ! used in the integration routine
       ndcr, &  ! the required number of nodes
       ierr, &  ! used to control allocation
       ninter, & ! number of possible energy intervals
       ncross, & ! actual number of nodes
       nstart ! starting point for inward integration


  !
  !  set up constants and allocate variables the 
  !
  allocate(f(mesh),stat=ierr)
  allocate(el(mesh),stat=ierr)
  allocate(c(mesh),stat=ierr)
  allocate(yi(mesh),stat=ierr)

  ddx12=dx*dx/12.0_dp
  l1=lam+1
  nst=l1*2
  sqlhf=(dble(lam)+0.5_dp)**2
  xl1=lam+1
  x4l6=4*lam+6
  x6l12=6*lam+12
  x8l20=8*lam+20

  ndcr=nn-lam-1
  !
  !  series developement of the potential near the origin
  !
  ze2=0.0_dp
  do n=1,4
     y(n)=vpot(n)-ze2/r(n)
  enddo
  call series(y,r,r2,b)

  eup=vpot(mesh)+sqlhf/r2(mesh)
  elw=eup
  do n=1,mesh
     elw=min(elw,vpot(n)+sqlhf/r2(n))
  enddo

  if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
  if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup

  do iter=1,maxter
     !
     !  set up the f-function and determine the position of its last
     !  change of sign
     !  f < 0 (approximatively) means classically allowed   region
     !  f > 0         "           "        "      forbidden   "
     !
     ik=0
     f(1)=ddx12*(r2(1)*(vpot(1)-e)+sqlhf)
     do n=2,mesh
        f(n)=ddx12*(r2(n)*(vpot(n)-e)+sqlhf)
        if( f(n) .ne. sign(f(n),f(n-1)).and.n.ne.mesh ) ik=n
     enddo
     if (ik.eq.0) ik=mesh*3/4
     !     +      call errore('ascheqps','f(n) has constant sign',1)

     if(ik.ge.mesh-2) then
        do n=1,mesh
           write(6,*) r(n), vpot(n), f(n)
        enddo
        call errore('ascheqps','No point found for matching',1)
        stop
     endif
     !
     !     if everything is ok continue the integration and define f
     !
     do n=1,mesh
        f(n)=1.0_dp-f(n)
     enddo
     !
     !  determination of the wave-function in the first two points by
     !  series developement
     !
     b0e=b(0)-e
     c1=0.5_dp*ze2/xl1
     c2=(c1*ze2+b0e)/x4l6
     c3=(c2*ze2+c1*b0e+b(1))/x6l12
     c4=(c3*ze2+c2*b0e+c1*b(1)+b(2))/x8l20
     rr1=(1.0_dp+r(1)*(c1+r(1)*(c2+r(1)*(c3+r(1)*c4))))*r(1)**l1
     rr2=(1.0_dp+r(2)*(c1+r(2)*(c2+r(2)*(c3+r(2)*c4))))*r(2)**l1
     y(1)=rr1/sqr(1)
     y(2)=rr2/sqr(2)
     !
     !     and outward integration
     !
     do n=2,ik-1
        y(n+1)=((12.0_dp-10.0_dp*f(n))*y(n)-f(n-1)*y(n-1))/f(n+1)
     enddo
     !
     !    count the number of nodes and find the maximum of the wavefunction
     !
     ncross=0
     ymx=0.0_dp
     do n=1,ik-1
        ymx=max(ymx,abs(y(n)))
        if ( y(n) .ne. sign(y(n),y(n+1)) ) ncross=ncross+1     
     enddo
     !
     !  matching radius has been reached going out. if ncross is not
     !  equal to ndcr, modify the trial eigenvalue.
     !
     if (ndcr.lt.ncross) then
        !
        !  too many crossings. e is an upper bound to the true eigen-
        !  value. increase abs(e)
        !
        eup=e
        rap=(dble(ncross+l1)/dble(nn))**2
        e=(e-vpot(mesh))*rap+vpot(mesh)
        if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup
        if(e.gt.eup) e=0.1_dp*elw+0.9_dp*eup
        go to 300

     else if (ndcr.gt.ncross) then
        !
        !  too few crossings. e is a lower bound to the true eigen-
        !  value. decrease abs(e)
        !
        elw=e
        rap=(dble(ncross+l1)/dble(nn))**2
        e=(e-vpot(mesh))*rap+vpot(mesh)
        if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
        if(e.lt.elw) e=0.1_dp*eup+0.9_dp*elw
        go to 300
     endif

     call integrate_inward(e,mesh,ndm,dx,r,r2,sqr,f,y, &
          c,el,ik,nstart)
     !
     !  if y is too large renormalize it
     !
     if (ymx.ge.1.e10_dp) then
        do n=1,mesh
           y(n)=y(n)/ymx
        enddo
     endif
     !
     !  if necessary, improve the trial eigenvalue by the cooley's
     !  procedure. jw cooley math of comp 15,363(1961)
     !
     fe=(12.0_dp-10.0_dp*f(ik))*y(ik) &
          -f(ik-1)*y(ik-1)-f(ik+1)*y(ik+1)

     do n=1,nstart
        el(n)=r(n)*y(n)*y(n)
     enddo
     sum=int_0_inf_dr(el,r,r2,dx,nstart,nst)

     dfe=-y(ik)*f(ik)/dx/sum
     de=-fe*dfe
     !
     !    check for convergence
     !
     eps=abs(de/e)
     if(abs(de).lt.thresh) go to 900
     !        write(6,*) iter,eps
     !
     !    if convergence not achieved adjust the eigenvalue
     !
     if(eps.gt.0.25_dp) de=0.25_dp*de/eps
     if(de.gt.0.0_dp) elw=e
     if(de.lt.0.0_dp) eup=e
     !
     !    since with the semilocal potential the node theorem is not verified,
     !    the routine can stall around a false point. In this case try to avoid
     !
     e=e+de
     if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
     if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup
300  continue
  enddo
  !      write(6,'('' solution not found in ascheqps'')')
  !1000  write(6,9000) nn,lam,elw,eup,eps
  !9000  format(5x,'error in ascheq: n l =',2i3,/
  !     + 5x,'elw =',f15.10,' eup =',f15.10,'eps',f15.10)
900 continue
  !
  !   exponential tail of the solution if it was not computed
  !
  if (nstart.lt.mesh) then
     do n=nstart,mesh-1
        if (y(n).eq.0.0_dp) then
           y(n+1)=0.0_dp
        else
           yln=log(abs(y(n)))
           xp=-sqrt(12.0_dp*abs(1.0_dp-f(n)))
           expn=yln+xp
           if (expn.lt.-80.0_dp) then
              y(n+1)=0.0_dp
           else
              y(n+1)=sign(exp(expn),y(n))
           endif
        endif
     enddo
  endif
  !
  !  normalize the eigenfunction and exit
  !
  do n=1,mesh
     el(n)=r(n)*y(n)*y(n)
  enddo
  sum=int_0_inf_dr(el,r,r2,dx,mesh,nst)
  sum=sqrt(sum)
  do n=1,mesh
     y(n)=sqr(n)*y(n)/sum
  enddo

  deallocate(yi)
  deallocate(el)
  deallocate(f )
  deallocate(c )

  return
end subroutine ascheqlocps
