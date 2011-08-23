!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine compute_solution(lam,jam,e,mesh,ndm,grid,vpot,y,beta,ddd,&
           qq,nbeta,nwfx,lls,jjs,ikk)
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation for
  !  bound states in a local potential.
  !
  !  This routine works at fixed e
  !  It allows a nonlocal potential and an overlap
  !  operator. Therefore it can solve a general schroedinger
  !  equation necessary to solve a US pseudopotential
  !
  use io_global, only : stdout
  use kinds, only : DP
  use radial_grids, only: radial_grid_type, series

  implicit none

  type(radial_grid_type), intent(in) :: grid
  integer :: &
       lam, &      ! l angular momentum
       mesh,&      ! size of radial mesh
       ndm, &      ! maximum radial mesh 
       nbeta,&     ! number of beta function  
       nwfx, &     ! maximum number of beta functions
       lls(nbeta),&! for each beta the angular momentum
       ikk(nbeta) ! for each beta the point where it become zero

  real(DP) :: &
       e,       &  ! output eigenvalue
       jam,     &  ! j angular momentum
       vpot(mesh),&! the local potential 
       y(mesh),  & ! the output solution
       jjs(nwfx), & ! the j angular momentum
       beta(ndm,nwfx),& ! the beta functions
       ddd(nwfx,nwfx),qq(nwfx,nwfx) ! parameters for computing B_ij
  !
  !    the local variables
  !
  real(DP) :: &
       ddx12,      &  ! dx**2/12 used for Numerov integration
       sqlhf,      &  ! the term for angular momentum in equation
       ze2,        &  ! possible coulomb term aroun the origin (set 0)
       b(0:3),     &  ! coefficients of taylor expansion of potential
       sum,       &! auxiliary 
       yln, xp, expn,& ! used to compute the tail of the solution
       int_0_inf_dr  ! integral function

  real(DP), allocatable :: &
       f(:),    &   ! the f function
       el(:),c(:) ! auxiliary for inward integration

  integer :: &
       n,  &    ! counter on mesh points
       ik,  &   ! matching point
       ns,  &   ! counter on beta functions
       l1,  &   ! lam+1
       nst, &   ! used in the integration routine
       ierr, &
       nstart  ! starting point for inward integration

   if (mesh.ne.grid%mesh) call errore('compute_solution','mesh dimension is not as expected',1)
  !
  !  set up constants and allocate variables the 
  !
  allocate(f(mesh), stat=ierr)
  allocate(el(mesh), stat=ierr)
  allocate(c(mesh), stat=ierr)

  ddx12=grid%dx*grid%dx/12.0_dp
  l1=lam+1
  nst=l1*2
  sqlhf=(DBLE(lam)+0.5_dp)**2
  !
  !  series developement of the potential near the origin
  !
  ze2=0.0_dp
  do n=1,4
     y(n)=vpot(n)-ze2/grid%r(n)
  enddo
  call series(y,grid%r,grid%r2,b)

  !      write(stdout,*) 'eneter lam,eup,elw,e',lam,nbeta,eup,elw,e
  !
  !  set up the f-function and determine the position of its last
  !  change of sign
  !  f < 0 (approximatively) means classically allowed   region
  !  f > 0         "           "        "      forbidden   "
  !
  ik=0
  f(1)=ddx12*(grid%r2(1)*(vpot(1)-e)+sqlhf)
  do n=2,mesh
     f(n)=ddx12*(grid%r2(n)*(vpot(n)-e)+sqlhf)
     if( f(n).ne.sign(f(n),f(n-1)).and.n.lt.mesh-5 ) ik=n
  enddo
  if (ik.eq.0.and.nbeta.eq.0) ik=mesh*3/4

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
     if (lls(ns).eq.lam.and.ikk(ns).gt.ik) ik=ikk(ns)
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
  call start_scheq( lam, e, b, grid, ze2, y )
  !
  !    outward integration before ik
  !
  call integrate_outward (lam,jam,e,mesh,ndm,grid,f, b,y,beta,ddd,qq,&
       nbeta,nwfx,lls,jjs,ikk,ik)
  !
  !    inward integration up to ik
  !
  call integrate_inward(e,mesh,ndm,grid,f,y,c,el,ik,nstart)
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
  !  normalize the eigenfunction and exit
  !
  do n=1,mesh
     el(n)=grid%r(n)*y(n)*y(n)
  enddo
  sum=int_0_inf_dr(el,grid,mesh,nst)
  sum=sqrt(sum)
  do n=1,mesh
     y(n)=grid%sqr(n)*y(n)/sum
  enddo

  deallocate(el)
  deallocate(f )
  deallocate(c )
  return

end subroutine compute_solution
