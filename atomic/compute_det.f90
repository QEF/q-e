!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine compute_det(nn,lam,jam,e,mesh,ndm,dx,r,r2,sqr,vpot, &
     beta,ddd,qq,nbeta,nwfx,lls,jjs,ikk,det)
  !---------------------------------------------------------------
  !
  !     This routine computes the function described in
  !     Khein PRB 51 16608 (1995) at the energy e.
  !     if nbeta=0 it computes the wronskian of the local
  !     potential
  !
  !
  use io_global, only : stdout
  use kinds, only : DP
  implicit none

  integer :: &
       nn, &       ! main quantum number for node number
       lam,&       ! angular momentum
       mesh,&      ! size of radial mesh
       ndm, &      ! maximum radial mesh 
       nbeta,&     ! number of beta function
       nwfx, &     ! maximum number of beta functions
       lls(nbeta),&! for each beta the angular momentum
       ikk(nbeta) ! for each beta the point where it become zero

  real(DP) :: &
       e,       &  ! output eigenvalue
       dx,      &  ! linear delta x for radial mesh
       jam,     &  ! the j of this calculation
       r(mesh), &  ! radial mesh
       r2(mesh),&  ! square of radial mesh
       sqr(mesh), &! square root of radial mesh
       vpot(mesh),&! the local potential 
       det,     & ! the value of the wronskian
       jjs(nwfx), &! the j of each beta
       beta(ndm,nwfx),& ! the beta functions
       ddd(nwfx,nwfx),qq(nwfx,nwfx)   ! parameters for computing B_ij
  !
  !    the local variables
  !
  real(DP) :: &
       ddx12,      &  ! dx**2/12 used for Numerov integration
       sqlhf,      &  ! the term for angular momentum in equation
       xl1, x4l6, x6l12, x8l20, &! used for starting series expansion
       ze2,    &      ! possible coulomb term aroun the origin (set 0)
       b(0:3), &      ! coefficients of taylor expansion of potential
       c1,c2,c3,c4,b0e, & ! auxiliary for expansion of wavefunction
       rr1,rr2,   &  ! values of y in the first points
       int_0_inf_dr, & ! the integral function
       deriv_7pts,   & ! compute the derivative
       wron1          ! the value of the wronskian

  real(DP),allocatable :: &
       bm(:,:),amat(:,:), &! matrix to compute det
       y(:),     & ! the regular solution
       yi(:),    & ! the irregular solution
       f(:),     & ! the f function
       el(:),c(:) ! auxiliary for inward integration


  integer :: &
       n, &     ! counter on mesh points
       iter, &  ! counter on iteration
       ik,   &  ! matching point
       ns,   &  ! counter on beta functions
       l1,   &  ! lam+1
       ib,jb,kb, &
       ierr,    &  ! used to control allocation
       iib,jjb,kkb, &! indeces for beta functions
       nstart ! starting point for inward integration

  !
  !  set up constants and allocate variables the 
  !
  allocate(f(mesh), stat=ierr)
  allocate(el(mesh), stat=ierr)
  allocate(c(mesh), stat=ierr)
  allocate(yi(mesh), stat=ierr)
  allocate(y(mesh), stat=ierr)
  allocate(bm(nwfx,nwfx), stat=ierr)
  allocate(amat(nwfx,nwfx), stat=ierr)


  ddx12=dx*dx/12.0_dp
  l1=lam+1
  sqlhf=(DBLE(lam)+0.5_dp)**2
  xl1=lam+1
  x4l6=4*lam+6
  x6l12=6*lam+12
  x8l20=8*lam+20
  !
  !  series developement of the potential near the origin
  !
  ze2=0.0_dp
  do n=1,4
     y(n)=vpot(n)-ze2/r(n)
  enddo
  call series(y,r,r2,b)
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
     if( f(n) .ne. sign(f(n),f(n-1)) .and.n.lt.mesh-5) ik=n
  enddo
  if (ik.eq.0.and.nbeta.eq.0) ik=mesh*3/4

  if(ik.ge.mesh-2) then
     do n=1,mesh
        write(stdout,*) r(n), vpot(n), f(n)
     enddo
     call errore('compute_det','No point found for matching',1)
     stop
  endif

  do ib=1,nbeta
     if (lls(ib).eq.lam.and.ik.lt.ikk(ib)) ik=ikk(ib)
  enddo
  if (ik.eq.0) ik=mesh*3/4
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
  !    an inward integration
  !
  yi(ik)=y(ik)
  call integrate_inward(e,mesh,ndm,dx,r,r2,sqr,f,yi, &
       c,el,ik,nstart)
  !
  !     complete the regular and irregular solutions
  !
  do n=ik+1,2,-1
     yi(n-1)=((12.0_dp-10.0_dp*f(n))*yi(n)-f(n+1)*yi(n+1))/f(n-1)
  enddo
  !
  !       compute the wronskian
  !
  do n=1,mesh
     yi(n)=yi(n)*sqr(n)
  enddo
  do n=1,ik-1
     y(n)=y(n)*sqr(n)
  enddo

  n=max(20,ik-30)
  wron1= (y(n)*deriv_7pts(yi,n,r(n),dx)-yi(n)*deriv_7pts(y,n,r(n),dx))
  !
  !    compute the vector h, and the unsymmetrized bm
  !     
  iib=0
  jjb=0
  do ib=1,nbeta
     if (lls(ib).eq.lam.and.jjs(ib).eq.jam) then
        iib=iib+1
        el=0.0_dp
        c(1)= 0.0_dp
        do n=1,ikk(ib)
           el(n)=y(n)*beta(n,ib)
        enddo
        do n=2,nstart-1,2
           c(n+1)=c(n-1)+(el(n-1)*r(n-1)+ &
                4.0_dp*el(n)*r(n) + el(n+1)*r(n+1))*dx/3.0_dp
        enddo
        c(2)=c(1)+(r(2)-r(1))*(c(3)-c(1))/(r(3)-r(1))
        do n=3,nstart-1,2
           c(n+1)=c(n-1)+(el(n-1)*r(n-1)+ &
                4.0_dp*el(n)*r(n) + el(n+1)*r(n+1))*dx/3.0_dp
        enddo
        !            do n=1,nstart
        !               write(stdout,*) r(n),c(n),c(n+1)
        !            enddo
        !            stop

        jjb=0
        do jb=1,nbeta
           if (lls(jb).eq.lam.and.jjs(jb).eq.jam) then
              jjb=jjb+1
              do n=1,ikk(jb)
                 el(n)=yi(n)*c(n)*beta(n,jb)
              enddo
              bm(jjb,iib)=0.0_dp
              do n=2,ikk(jb)-1,2
                 bm(jjb,iib)=bm(jjb,iib)+(el(n-1)*r(n-1)+ &
                      4.0_dp*el(n)*r(n) + el(n+1)*r(n+1))*dx/3.0_dp
              enddo
           endif
        enddo
     endif
  enddo
  if (iib.ne.jjb) call errore('compute_det', &
       'iib different from jjb',1)
  !
  !    symmetrize b
  !
  !      write(stdout,*) iib,2.0_DP*bm(1,1),ddd(1,1),qq(1,1)
  do ib=1,iib
     do jb=1,ib
        bm(ib,jb)=(bm(ib,jb)+bm(jb,ib))/wron1
     enddo
  enddo
  do ib=1,iib
     do jb=1,ib-1
        bm(jb,ib)=bm(ib,jb)
     enddo
  enddo
  !
  !    compute A
  !

  iib=0
  jjb=0
  kkb=0
  do ib=1,nbeta
     if (lls(ib).eq.lam.and.jjs(ib).eq.jam) then
        iib=iib+1
        jjb=0  
        do jb=1,nbeta
           if (lls(jb).eq.lam.and.jjs(jb).eq.jam) then
              jjb=jjb+1
              amat(iib,jjb)=0.0_dp
              kkb=0  
              do kb=1,nbeta
                 if (lls(kb).eq.lam.and.jjs(kb).eq.jam) then
                    kkb=kkb+1
                    amat(iib,jjb)=amat(iib,jjb)+ &
                         (ddd(kb,jb)-e*qq(kb,jb))*bm(iib,kkb)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  if (iib.ne.jjb.or.iib.ne.kkb) call errore('compute_det', &
       'iib different from jjb or kkb',1)
  !
  !    compute the determinant of A-1
  !      

  do ib=1,iib
     amat(ib,ib)=amat(ib,ib)-1.0_dp
  enddo

  if (iib == 1) then
     det=amat(1,1)
  elseif (iib == 2) then
     det=amat(1,1)*amat(2,2)-amat(2,1)*amat(1,2)
  elseif (iib == 3) then
     det=amat(1,1)*(amat(2,2)*amat(3,3) - amat(3,2)*amat(2,3)) -       &
         amat(1,2)*(amat(2,1)*amat(3,3) - amat(3,1)*amat(2,3)) +       &
         amat(1,3)*(amat(2,1)*amat(3,2) - amat(3,1)*amat(2,2)) 
  elseif (iib == 0) then
     det=wron1
  else
     call errore('compute_det','too many beta functions',1)
  endif

  deallocate(amat) 
  deallocate(bm) 
  deallocate(c) 
  deallocate(y) 
  deallocate(yi) 
  deallocate(el) 
  deallocate(f) 

  return
end subroutine compute_det
