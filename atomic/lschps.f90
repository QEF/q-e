!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine lschps(mode,z,grid,nin,mch,n,l,e,u,v)
!
  !
  ! integrates radial pauli-type scalar-relativistic equation
  ! on a logarithmic grid
  ! modified routine to be used in finding norm-conserving
  ! pseudopotential
  !
  ! mode = 1 is for full potential bound state
  ! mode = 2 is for pseudopotential bound state
  ! mode = 3 is for full potential to find log derivative
  ! mode = 4 is for pseudopotential to find energy which produces
  !            specified log derivative (fake bound state)
  ! mode = 5 is for pseudopotential to produce wavefunction beyond
  !            radius used for pseudopotential construction
  !
  use kinds, only : DP
  use radial_grids, only: radial_grid_type
  use ld1inc, only : cau_fact
  implicit none
  !
  ! I/O variables
  !
  integer, intent (in) :: mode
  real(DP), intent(in) :: z
  type (radial_grid_type), intent(in) :: grid
  integer :: nin, mch, n, l
  real(DP) :: e
  real(DP):: v(grid%mesh),u(grid%mesh)
  !
  ! local variables
  !
  real(DP),parameter:: e2=2.0_dp
  real(DP), external:: aei, aeo, aii, aio
  real(DP):: al, als, ammax,  cn
  real(DP):: de, emax, emin
  real(DP):: eps, fss, gamma, ro, sc
  real(DP):: sls, sn, tfapot, uld, uout,  upin, upout
  real(DP):: xkap, exp
  integer:: i, it, mmax, nint, node, ndm, ierr

  ! these arrays are used as work space
  real(DP),allocatable :: up(:),upp(:),cf(:),dv(:),fr(:),frp(:)

!- modified in order to pass radial grid data 
!
  ammax= exp(grid%dx)
  al   = grid%dx
  mmax = grid%mesh

  allocate(up(mmax), stat=ierr)
  allocate(upp(mmax), stat=ierr)
  allocate(cf(mmax), stat=ierr)
  allocate(dv(mmax), stat=ierr)
  allocate(fr(mmax), stat=ierr)
  allocate(frp(mmax), stat=ierr)

  uld=0.0_dp
  v=v/e2
  e=e/e2
  !
  ! convergence factor for solution of schroedinger eq.  if calculated
  ! correction to eigenvalue is smaller in magnitude than eps times
  ! the magnitude of the current guess, the current guess is not changed.
  eps=1.0e-8_dp
  !
  ! relativistic - non-relativistic switch
  !
  if(mode .eq. 1 .or. mode .eq. 3) then
!     fss=(1.0_dp/137.036_dp)**2
     fss=(1.0_dp/cau_fact)**2
     if(l == 0) gamma=sqrt(1.0_dp-fss*z**2)
     if(l .gt. 0) gamma=(l*sqrt(l**2-fss*z**2) + &
          (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
  else
     fss=1.0e-20_dp
     gamma=l+1
  end if
  !
  sls=l*(l+1)
  !
  if(mode .eq. 1 .or. mode .eq. 2) then
     emax=v(mmax)+0.5_dp*sls/grid%r(mmax)**2
     emin=0.0_dp
     do i=1,mmax
        emin=min(emin,v(i)+0.5_dp*sls/grid%r(i)**2)
        !           if (l.eq.0)  write(6,*) grid%r(i),v(i)  
     end do
     !         if (l.eq.0) stop  
     if(e .gt. emax) e=1.25_dp*emax
     if(e .lt. emin) e=0.75_dp*emin
     if(e .gt. emax) e=0.5_dp*(emax+emin)
  else if(mode .eq. 4) then
     emax=e + 10.0_dp
     emin=e - 10.0_dp
  end if
  !
  do i=1,4
     u(i)=0.0_dp
     up(i)=0.0_dp
     upp(i)=0.0_dp
  end do
  nint=0
  als=al**2
  !
  ! return point for bound state convergence
10 nint=nint+1
  if(nint .gt. 60) then
     print '('' warning: wfc '',2i2,'' not converged'')', n, l
     u=0.0_dp
     go to 999
  end if
  !
  ! coefficient array for u in differential eq.
  do i=1,mmax
     cf(i)=als*sls + 2.0_dp*als*(v(i)-e)*grid%r(i)**2
  end do
  !
  ! calculate dv/dr for darwin correction
  dv(1)=(-50.0_dp*v(1)+96.0_dp*v(2)-72.0_dp*v(3)+32.0_dp*v(4) &
       -6.0_dp*v(5))/(24.0_dp*al*grid%r(1))
  dv(2)=(-6.0_dp*v(1)-20.0_dp*v(2)+36.0_dp*v(3)-12.0_dp*v(4) &
       +2.0_dp*v(5))/(24.0_dp*al*grid%r(2))
  !
  do i=3,mmax-2
     dv(i)=(2.0_dp*v(i-2)-16.0_dp*v(i-1)+16.0_dp*v(i+1) &
          -2.0_dp*v(i+2))/(24.0_dp*al*grid%r(i))
  end do
  dv(mmax-1)=( 3.0_dp*v(mmax)+10.0_dp*v(mmax-1)-18.0_dp*v(mmax-2)+ &
       6.0_dp*v(mmax-3)-v(mmax-4))/(12.0_dp*al*grid%r(mmax-1))
  dv(mmax)=( 25.0_dp*v(mmax)-48.0_dp*v(mmax-1)+36.0_dp*v(mmax-2)-&
       16.0_dp*v(mmax-3)+3.0_dp*v(mmax-4))/(12.0_dp*al*grid%r(mmax))
  !
  !  relativistic coefficient arrays for u (fr) and up (frp).
  do i=1,mmax
     fr(i)=als*(grid%r(i)**2)*(-fss*(v(i)-e)**2 + 0.5_dp*fss*dv(i)/ &
          (grid%r(i)*(1.0_dp+0.5_dp*fss*(e-v(i)))))
     frp(i)=-al*grid%r(i)*0.5_dp*fss*dv(i)/(1.0_dp+0.5_dp*fss*(e-v(i)))
  end do
  !
  ! find classical turning point for matching
  if(mode .eq. 1 .or. mode .eq. 2) then
     do i=mmax,2,-1
        if(cf(i-1) .le. 0.0_dp .and. cf(i) .gt. 0.0_dp) then
           mch=i
           go to 40
        end if
     end do
     print '('' warning: wfc '',2i2,'' no turning point'')', n, l
     e=0.0_dp
     do i=1,mmax
        u (i)=0.0_dp
     end do
     go to 999
  else
     nin=mch
  end if
40 continue
  !
  ! start wavefunction with series
  !
  do i=1,4
     u(i)=grid%r(i)**gamma
     up(i)=al*gamma*grid%r(i)**gamma
     upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
  end do
  !
  ! outward integration using predictor once, corrector
  ! twice
  node=0
  !               
  do i=4,mch-1
     u(i+1)=u(i)+aeo(up,i)
     up(i+1)=up(i)+aeo(upp,i)
     do it=1,2
        upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
        up(i+1)=up(i)+aio(upp,i)
        u(i+1)=u(i)+aio(up,i)
     end do
     if(u(i+1)*u(i) .le. 0.0_dp) node=node+1
  end do
  !
  uout=u(mch)
  upout=up(mch)
  !
  !
  if(node-n+l+1 .eq. 0 .or. mode .eq. 3 .or. mode .eq. 5) then
     !
     if(mode .eq. 1 .or. mode .eq. 2) then
        !
        ! start inward integration at 10*classical turning
        ! point with simple exponential
        nin=mch+2.3_dp/al
        if(nin+4 .gt. mmax) nin=mmax-4
        xkap=sqrt(sls/grid%r(nin)**2 + 2.0_dp*(v(nin)-e))
        !
        do i=nin,nin+4
           u(i)=exp(-xkap*(grid%r(i)-grid%r(nin)))
           up(i)=-grid%r(i)*al*xkap*u(i)
           upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
        end do
        !
        ! integrate inward
        !
        do i=nin,mch+1,-1
           u(i-1)=u(i)+aei(up,i)
           up(i-1)=up(i)+aei(upp,i)
           do it=1,2
              upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
              up(i-1)=up(i)+aii(upp,i)
              u(i-1)=u(i)+aii(up,i)
           end do
        end do
        !
        ! scale outside wf for continuity
        sc=uout/u(mch)
        !
        do i=mch,nin
           up(i)=sc*up(i)
           u (i)=sc*u (i)
        end do
        !
        upin=up(mch)
        !
     else
        !
        upin=uld*uout
        !
     end if
     !
     ! perform normalization sum
     !
     ro=grid%r(1)/sqrt(ammax)
     sn=ro**(2.0_dp*gamma+1.0_dp)/(2.0_dp*gamma+1.0_dp)
     !
     do i=1,nin-3
        sn=sn+al*grid%r(i)*u(i)**2
     end do
     !
     sn=sn + al*(23.0_dp*grid%r(nin-2)*u(nin-2)**2 &
          + 28.0_dp*grid%r(nin-1)*u(nin-1)**2 &
          +  9.0_dp*grid%r(nin  )*u(nin  )**2)/24.0_dp
     !
     ! normalize u
     cn=1.0_dp/sqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin
     !
     do i=1,nin
        up(i)=cn*up(i)
        u(i)=cn*u(i)
     end do
     do i=nin+1,mmax
        u(i)=0.0_dp
     end do
     !
     ! exit for fixed-energy calculation
     !
     if(mode .eq. 3 .or. mode .eq. 5) go to 999

     ! perturbation theory for energy shift
     de=0.5_dp*uout*(upout-upin)/(al*grid%r(mch))
     !
     ! convergence test and possible exit
     !
     if ( abs(de) .lt. max(abs(e),0.2_dp)*eps) go to 999
     !
     if(de .gt. 0.0_dp) then 
        emin=e
     else
        emax=e
     end if
     e=e+de
     if(e .gt. emax .or. e .lt. emin) e=0.5_dp*(emax+emin)
     !
     ! loop back to converge e
     !
     go to 10
     !
  else if(node-n+l+1 .lt. 0) then
     ! too few nodes
     emin=e
     e=0.5_dp*(emin+emax)
     go to 10
     !
  else
     ! too many nodes
     emax=e
     e=0.5_dp*(emin+emax)
     go to 10
  end if
  !
  ! deallocate arrays and exit
  !
999 continue
  deallocate(frp)
  deallocate(fr)
  deallocate(dv)
  deallocate(cf)
  deallocate(upp)
  deallocate(up)
  e=e*e2
  v=v*e2
  return

end subroutine lschps
!
function aei(y,j)
  !
  use kinds, only : DP
  implicit none
  integer j
  real(DP):: y(j+3), aei
  !
  aei=-(4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j+1) &
       +37.0_dp*y(j+2)-9.0_dp*y(j+3))
  return
end function aei
!
! adams extrapolation and interpolation formulas for
! outward and inward integration, abramowitz and
! stegun, p. 896
function aeo(y,j)
  !
  use kinds, only : DP
  implicit none
  integer:: j   
  real(DP):: y(j), aeo
  !
  aeo=(4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j-1) &
       +37.0_dp*y(j-2)-9.0_dp*y(j-3))
  return
end function aeo
!
function aii(y,j)
  !
  use kinds, only : DP
  implicit none
  integer:: j
  real(DP) :: y(j+2), aii
  !
  aii=-(4.16666666667e-2_dp)*(9.0_dp*y(j-1)+19.0_dp*y(j) &
       -5.0_dp*y(j+1)+y(j+2))
  return
end function aii
!
function aio(y,j)
  !
  use kinds, only : DP
  implicit none
  integer :: j
  real(DP):: y(j+1), aio
  !
  aio=(4.16666666667e-2_dp)*(9.0_dp*y(j+1)+19.0_dp*y(j) &
       -5.0_dp*y(j-1)+y(j-2))
  return
end function aio
