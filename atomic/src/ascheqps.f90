!
! Copyright (C) 2004-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine ascheqps ( nam, lam, jam, e0, mesh, ndm, grid, vpot, thresh,&
                      y, beta, ddd, qq, nbeta, nwfx, lls, jjs, ikk, nstop )
  !---------------------------------------------------------------
  !
  !  numerical integration of a generalized radial schroedinger equation,
  !  using Numerov with outward and inward integration and matching.
  !  Works for both norm-conserving nonlocal and US pseudopotentials
  !  Requires in input a good estimate "e0" of the energy
  !
  use io_global, only : stdout
  use kinds, only : DP
  use radial_grids, only: radial_grid_type, series

  implicit none

  type(radial_grid_type), intent(in) :: grid
  integer, intent(in) :: &
       nam, &
       lam, &      ! l angular momentum
       mesh,&      ! size of radial mesh
       ndm, &      ! maximum radial mesh 
       nbeta,&     ! number of beta function  
       nwfx, &     ! maximum number of beta functions
       ikk(nbeta),&! for each beta the point where it become zero
       lls(nbeta)  ! for each beta the angular momentum

  real(DP), intent(in) :: &
       jam,       & ! j angular momentum
       vpot(mesh),& ! the local potential 
       thresh,    & ! precision of eigenvalue
       jjs(nwfx), & ! the j angular momentum
       beta(ndm,nwfx), &            ! the beta functions
       ddd(nwfx,nwfx),qq(nwfx,nwfx) ! parameters for computing B_ij

  real(DP), intent(inout) :: &
       e0,      &  ! output eigenvalue
       y(mesh)     ! the output solution

  integer, intent(out) :: &
       nstop       ! error code, used to check the behavior of the routine
  !
  !    the local variables
  !
  integer :: &
       ndcr,  &    ! number of required nodes
       n1, n2, &   ! counters
       ikl         ! auxiliary variables
  real(DP) :: &
       work(nbeta),& ! auxiliary space
       e,          &  ! energy
       ddx12,      &  ! dx**2/12 used for Numerov integration
       sqlhf,      &  ! the term for angular momentum in equation
       ze2,        &  ! possible coulomb term aroun the origin (set 0)
       b(0:3),     &  ! coefficients of taylor expansion of potential
       eup,elw,    & ! actual energy interval
       ymx,        & ! the maximum value of the function
       fe,integ,dfe,de, &! auxiliary for numerov computation of e
       eps,        & ! the epsilon of the delta e
       yln, xp, expn,& ! used to compute the tail of the solution
       int_0_inf_dr  ! integral function

  real(DP), allocatable :: &
       fun(:),  &   ! integrand function
       f(:),    &   ! the f function
       el(:),c(:) ! auxiliary for inward integration

  integer, parameter :: &
       maxter=100    ! maximum number of iterations

  integer :: &
       n,  &    ! counter on mesh points
       iter,&   ! counter on iteration
       ik,  &   ! matching point
       ns,  &   ! counter on beta functions
       l1,  &   ! lam+1
       nst, &   ! used in the integration routine
       ierr, &
       ncross,& ! actual number of nodes
       nstart  ! starting point for inward integration

  logical, save :: first(0:10,0:10) = .true.


   if (mesh.ne.grid%mesh) &
        call errore('compute_solution','mesh dimension is not as expected',1)
  !
  !  set up constants and allocate variables the 
  !
  allocate(fun(mesh), stat=ierr)
  allocate(f(mesh), stat=ierr)
  allocate(el(mesh), stat=ierr)
  allocate(c(mesh), stat=ierr)
  nstop=0
  nstart=0
  e=e0
!  write(6,*) 'entering ', nam,lam, e
  eup=0.3_DP*e
  elw=1.3_dp*e
  ndcr=nam-lam-1
!  write(6,*) 'entering ascheqps', vpot(mesh-20)*grid%r(mesh-20)

  ddx12=grid%dx*grid%dx/12.0_dp
  l1=lam+1
  nst=l1*2
  sqlhf=(DBLE(lam)+0.5_dp)**2
  !
  !  series developement of the potential near the origin
  !
  do n=1,4
     y(n)=vpot(n)
  enddo
  call series(y,grid%r,grid%r2,b)

  !      write(stdout,*) 'enter lam,eup,elw,e',lam,nbeta,eup,elw,e
  !
  !  set up the f-function and determine the position of its last
  !  change of sign
  !  f < 0 (approximatively) means classically allowed   region
  !  f > 0         "           "        "      forbidden   "
  !
  do iter=1,maxter
!  write(6,*) 'starting iter', iter, elw, e, eup
  ik=1
  f(1)=ddx12*(grid%r2(1)*(vpot(1)-e)+sqlhf)
  do n=2,mesh
     f(n)=ddx12*(grid%r2(n)*(vpot(n)-e)+sqlhf)
     if( f(n).ne.sign(f(n),f(n-1)).and.n.lt.mesh-5 ) ik=n
  enddo
!  if (ik.eq.0.and.nbeta.eq.0) ik=mesh*3/4
  if (ik.eq.1.or.grid%r(ik)>4.0_DP) ik=mesh*3/4

  if(ik.ge.mesh-2) then
     do n=1,mesh
        write(stdout,*) grid%r(n), vpot(n), f(n)
     enddo
     call errore('compute_solution', 'No point found for matching',1)
  endif
  !
  !     determine if ik is sufficiently large
  !
  do ns=1,nbeta
     if (lls(ns).eq.lam .and. jjs(ns).eq.jam .and. ikk(ns).gt.ik) ik=ikk(ns)+3
  enddo
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
  ! no coulomb divergence in the origin for a pseudopotential
  ze2=0.d0 
  call start_scheq( lam, e, b, grid, ze2, y )
  !
  !    outward integration before ik
  !
  call integrate_outward (lam,jam,e,mesh,ndm,grid,f,b,y,beta,ddd,qq,&
       nbeta,nwfx,lls,jjs,ikk,ik)

  ncross=0
  ymx=0.0_dp
  do n=2,ik-1
     if ( y(n) .ne. sign(y(n),y(n+1)) ) ncross=ncross+1
     ymx=max(ymx,abs(y(n+1)))
  end do
!
!  If at this point the number of nodes is wrong it means that something
!  is probably wrong in the calling routines. A ghost might be present
!  in the pseudopotential. With a nonlocal pseudopotential there is no
!  node theorem so strictly speaking the following instructions are
!  wrong but sometimes they help so we keep them here.
!
  if( ndcr /= ncross .and. first(nam,lam) ) then
      write(stdout,"(/,7x,'Warning: n=',i1,', l=',i1,' expected ',i1,' nodes,',&
        & ' found ',i1)") nam,lam,ndcr,ncross
      write(stdout,'(7x,a,/,7x,a,/)') 'Setting wfc to zero for this iteration',&
                    '(This warning will only be printed once per wavefunction)'
      first(nam,lam) = .false.
  endif

  if(ndcr < ncross) then
     !
     !  too many crossings. e is an upper bound to the true eigenvalue.
     !  increase abs(e)
     !
     eup=e
     e=0.9_dp*elw+0.1_dp*eup
!     write(6,*) 'too many crossing', ncross, ndcr
!     call errore('aschqps','wrong number of nodes. Probably a Ghost?',1)
     y=0.0_DP
     ymx=0.0_dp
     go to 300
  else if (ndcr > ncross) then
     !
     !  too few crossings. e is a lower bound to the true eigenvalue.
     !  decrease abs(e)
     !
     elw=e
     e=0.9_dp*eup+0.1_dp*elw
!     write(6,*) 'too few crossing', ncross, ndcr
!     call errore('aschqps','wrong number of nodes. Probably a Ghost?',1)
     y=0.0_DP
     ymx=0.0_dp
     go to 300
  end if
  !
  !    inward integration up to ik
  !
  call integrate_inward(e,mesh,ndm,grid,f,y,c,el,ik,nstart)
  !
  !  if necessary, improve the trial eigenvalue by the cooley's
  !  procedure. jw cooley math of comp 15,363(1961)
  !
  fe=(12.0_dp-10.0_dp*f(ik))*y(ik)-f(ik-1)*y(ik-1)-f(ik+1)*y(ik+1)
  !
  ! audjust the normalization if needed
  !
  if(ymx.ge.1.0e10_dp) y=y/ymx
  !
  !  calculate the normalization
  !
  do n1=1,nbeta
     if (lam.eq.lls(n1).and.abs(jam-jjs(n1)).lt.1.e-7_dp) then
        ikl=ikk(n1)
        do n=1,ikl
           fun(n)=beta(n,n1)*y(n)*grid%sqr(n)
        enddo
        work(n1)=int_0_inf_dr(fun,grid,ikl,nst)
     else
        work(n1)=0.0_dp
     endif
  enddo
  do n=1,nstart
     fun(n)= y(n)*y(n)*grid%r(n)
  enddo
  integ=int_0_inf_dr(fun,grid,nstart,nst)
  do n1=1,nbeta
     do n2=1,nbeta
        integ = integ + qq(n1,n2)*work(n1)*work(n2)
     enddo
  enddo
  dfe=-y(ik)*f(ik)/grid%dx/integ
  de=-fe*dfe
  eps=abs(de/e)
!  write(6,'(i5, 3f20.12)') iter, e, de
  if(abs(de).lt.thresh) go to 600
  if(eps.gt.0.25_dp) de=0.25_dp*de/eps
  if(de.gt.0.0_dp) elw=e
  if(de.lt.0.0_dp) eup=e
  e=e+de
  if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
  if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup
300 continue
  enddo
  nstop=1
  if (nstart==0) goto 900
600 continue
  !
  !   exponential tail of the solution if it was not computed
  !
  if (nstart.lt.mesh) then
     do n=nstart,mesh-1
        if (y(n) == 0.0_dp) then
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
  !  normalize the eigenfunction as if they were norm conserving. 
  !  If this is a US PP the correct normalization is done outside this
  !  routine.
  !
  do n=1,mesh
     el(n)=grid%r(n)*y(n)*y(n)
  enddo
  integ=int_0_inf_dr(el,grid,mesh,nst)
  if (integ>0.0_DP) then
     integ=sqrt(integ)
     do n=1,mesh
        y(n)=grid%sqr(n)*y(n)/integ
     enddo
     e0=e
  else
     nstop=1
  endif

900 continue

  deallocate(el)
  deallocate(f )
  deallocate(c )
  deallocate(fun )
  return

end subroutine ascheqps
