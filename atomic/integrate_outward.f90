!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine integrate_outward (lam,jam,e,mesh,ndm,grid,f, &
     b,y,beta,ddd,qq,nbeta,nwfx,lls,jjs,ikk,ik)
  !---------------------------------------------------------------------
  !
  !    Integrate the wavefunction from 0 to r(ik) 
  !    generalized separable or US pseudopotentials are allowed
  !    This routine assumes that y countains already the
  !    correct values in the first two points
  !
  use kinds, only : DP
  use radial_grids, only : radial_grid_type
  implicit none
  type(radial_grid_type):: grid
  integer ::   &
       lam,    &     ! l angular momentum
       mesh,   &    ! size of radial mesh
       ndm,    &   ! maximum radial mesh
       nbeta,  &   ! number of beta function
       nwfx,   &   ! maximum number of beta functions
       lls(nbeta),&! for each beta the angular momentum
       ikk(nbeta),&! for each beta the integration point
       ik         ! the last integration point

  real(DP) :: &
       e,       &  ! output eigenvalue
       jam,     &  ! j angular momentum
       f(mesh), &  ! the f function
       b(0:3), &   ! the taylor expansion of the potential
       y(mesh), &  ! the output solution
       jjs(nwfx), & ! the j angular momentum
       beta(ndm,nwfx),& ! the beta functions
       ddd(nwfx,nwfx),qq(nwfx,nwfx) ! parameters for computing B_ij

  integer ::  &
       nst, &      ! the exponential around the origin
       n,    &     ! counter on mesh points
       iib,jjb, &  ! counter on beta with correct lam
       ierr,    &  ! used to control allocation
       ib,jb,   &  ! counter on beta
       info      ! info on exit of LAPACK subroutines

  integer, allocatable :: iwork(:) ! auxiliary space  

  real(DP) :: &
       b0e,     & ! the expansion of the known part
       ddx12,   & ! the deltax enetering the equations
       x4l6,    & ! auxiliary for small r expansion
       j1(4),d(4),& ! auxiliary for starting values of chi
       delta,xc(4),& ! auxiliary for starting values of eta
       int_0_inf_dr  ! the integral function

  real(DP), allocatable :: &
       el(:), &  ! auxiliary for integration
       cm(:,:), &! the linear system
       bm(:), & ! the known part of the linear system
       c(:), &   ! the chi functions
       coef(:), & ! the solution of the linear system
       eta(:,:) ! the partial solution of the nonomogeneous


  allocate(c(ik), stat=ierr)
  allocate(el(ik), stat=ierr)
  allocate(cm(nbeta,nbeta), stat=ierr)
  allocate(bm(nbeta), stat=ierr)
  allocate(coef(nbeta), stat=ierr)
  allocate(iwork(nbeta), stat=ierr)
  allocate(eta(ik,nbeta), stat=ierr)

  if (mesh.ne.grid%mesh) call errore('integrate_outward','mesh dimension is not as expected',1)
  ddx12=grid%dx*grid%dx/12.0_DP
  b0e=b(0)-e
  x4l6=4*lam+6
  nst=(lam+1)*2
  !
  !  first solve the omogeneous equation
  !
  do n=2,ik-1
     y(n+1)=((12.0_DP-10.0_DP*f(n))*y(n)-f(n-1)*y(n-1))/f(n+1)
  enddo
  !
  !     for each beta function with correct angular momentum
  !     solve the inhomogeneous equation
  !
  iib=0
  jjb=0
  do ib=1,nbeta
     !
     !    set up the known part
     !         
     if (lls(ib).eq.lam.and.jjs(ib).eq.jam) then
        iib=iib+1
        c=0.0_DP
        do jb=1,nbeta
           if (lls(jb).eq.lam.and.jjs(jb).eq.jam) then
              do n=1,ikk(jb)
                 c(n)= c(n)+(ddd(jb,ib) -e*qq(jb,ib))*beta(n,jb)
              enddo
           endif
        enddo
        !
        !     compute the starting values of the solutions
        !
        do n=1,4
           j1(n)=c(n)/grid%r(n)**(lam+1)
        enddo
        call seriesbes(j1,grid%r,grid%r2,4,d)
        delta=b0e**2+x4l6*b(2)
        xc(1)=(-d(1)*b0e-x4l6*d(3))/delta
        xc(3)=(-b0e*d(3)+d(1)*b(2))/delta
        xc(2)=0.0_DP
        xc(4)=0.0_DP
        do n=1,3
           eta(n,iib)=grid%r(n)**(lam+1)*(xc(1)+grid%r2(n)*xc(3))/grid%sqr(n)
        enddo

        do n=1,ik
           c(n)=c(n)*grid%r2(n)/grid%sqr(n)
        enddo
        !
        !   solve the inhomogeneous equation
        !
        do n=3,ik-1
           eta(n+1,iib)=((12.0_DP-10.0_DP*f(n))*eta(n,iib) &
                -f(n-1)*eta(n-1,iib) &
                + ddx12*(10.0_DP*c(n)+c(n-1)+c(n+1)) )/f(n+1)
        enddo
        !
        !    compute the coefficents of the linear system
        !
        jjb=0
        do jb=1,nbeta
           if (lls(jb).eq.lam.and.jjs(jb).eq.jam) then
              jjb=jjb+1
              do n=1,min(ik,ikk(jb))
                 el(n)=beta(n,jb)*eta(n,iib)*grid%sqr(n)
              enddo
              cm(jjb,iib)=-int_0_inf_dr(el,grid,min(ik,ikk(jb)),nst)
           endif
        enddo

        do n=1,min(ik,ikk(ib))
           el(n)=beta(n,ib)*y(n)*grid%sqr(n)
        enddo
        bm(iib)=int_0_inf_dr(el,grid,min(ik,ikk(ib)),nst)
        cm(iib,iib)=1.0_DP+cm(iib,iib)
     endif
  enddo
  if (iib.ne.jjb) call errore('integrate_outward','jjb.ne.iib',1)

  if (iib.gt.0) then
     call dcopy(iib,bm,1,coef,1)

     call DGESV(iib,1,cm,nbeta,iwork,coef,nbeta,info)

     if (info /= 0) call errore('integrate_outward', &
                &  'problems solving the linear system',info)

     do ib=1,iib
        do n=1,ik
           y(n)= y(n)+coef(ib)*eta(n,ib)
        enddo
     enddo
  endif

  deallocate(eta)
  deallocate(iwork)
  deallocate(coef)
  deallocate(bm)
  deallocate(cm)
  deallocate(el)
  deallocate(c)

  return
end subroutine integrate_outward
