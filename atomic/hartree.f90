!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine hartree(k,nst,mesh,r,r2,sqr,dx,f,vh)
  !---------------------------------------------------------------
  !
  use io_global, only : stdout
  use kinds, only : DP
  implicit none
  integer ::       & 
       k,   & ! input: the k of the equation
       nst, & ! input: at low r, f goes as r**nst
       mesh   ! input: the dimension of the mesh

  real(DP)::        &
       r(mesh),   &  ! input: the radial mesh
       r2(mesh),  &  ! input: the square of the 
                     !        radial mesh
       sqr(mesh), &  ! input: the square root of 
                     !        the radial mesh
       dx,        &  ! input: the delta x of the 
                     !        linear mesh
       f(mesh),   &  ! input: the 4\pi r2 \rho 
                     !        function
       vh(mesh)      ! output: the required solution

  integer ::        &
       k21,  &   ! 2k+1
       nk1,  &   ! nst-k-1
       ierr, &   ! integer variable for 
                 ! allocation control
       i         ! counter

  real(DP)::        &
       c0,c2,c3, & ! coefficients of the polynomial
                   ! expansion close to r=0
       ch,       & ! dx squared / 12.0
       xkh2,     & ! ch * f
       ei, di,   & ! auxiliary variables for the 
                   ! diagonal and ! off diagonal 
                   ! elements of the matrix
       f1, fn,   & ! variables used for the 
                   ! boundary condition
       vhim1, vhi  ! variables for the right hand 
                   ! side

  real(DP), allocatable:: &
       d(:), &       ! the diagonal elements of 
                     ! the tridiagonal sys.
       e(:)          ! the off diagonal elements 
                     ! of the trid. sys.
  !
  ! Allocate space for the diagonal and off diagonal elements
  !
  allocate(d(mesh),stat=ierr)
  allocate(e(mesh),stat=ierr)
  call errore('hartree', 'allocating d or e ', ierr)
  !
  !    Find the series expansion of the solution close to r=0
  !
  k21=2*k+1
  nk1=nst-k-1
  if(nk1.le.0) then
     write(stdout,100) k,nst
100  format(5x,'stop in "hartree": k=',i3,'  nst=',i3)
     call errore('hartree', 'wrong k or nst ', 1)
  else if(nk1.ge.4) then
     c2=0.0_dp
     c3=0.0_dp
  else
     e(1)=0.0_dp
     do i=1,4
        d(i)=-k21*f(i)/r(i)**nst
     end do
     call series(d,r,r2,e(nk1))
     c2=e(1)/(4.0_dp*k+6.0_dp)
     c3=e(2)/(6.0_dp*k+12.0_dp)
  end if
  !
  !  Set the main auxiliary parameters
  !
  ch=dx*dx/12.0_dp
  xkh2=ch*(DBLE(k)+0.5_dp)**2
  ei=1.0_dp-xkh2
  di=-(2.0_dp+10.0_dp*xkh2)
  !
  !  Set the diagonal and the off diagonal elements of the 
  !  linear system, compute a part of the right hand side 
  !
  do i=2,mesh
     d(i)=-di
     e(i)=-ei
     vh(i)=k21*ch*sqr(i)*f(i)
  end do
  !
  !   Use the boundary condition to eliminate the value of the 
  !   solution in the first point from the first equation. This 
  !   part for the diagonal element
  !
  f1=(sqr(1)/sqr(2))**k21
  d(2)=d(2)-ei*f1
  !
  !   Use the boundary condition to eliminate the value of the 
  !   solution in the last point from the last equation
  !
  fn=(sqr(mesh-1)/sqr(mesh))**k21
  d(mesh-1)=d(mesh-1)-ei*fn
  !
  !   In the first point vh(1) has the same definition as in 
  !   the other points
  !
  vhim1=k21*ch*sqr(1)*f(1)
  !
  !   Compute the right hand side using the auxiliary quantity
  !   vh(i).
  !
  do i=2,mesh-1
     vhi=vh(i)
     vh(i)=vhim1+10.0_dp*vhi+vh(i+1)
     vhim1=vhi
  end do
  !
  !   Use the boundary condition to eliminate the value of the
  !   solution in the first point from the first equation. This 
  !   part for the right hand side.
  !
  vh(2)=vh(2)-ei*sqr(1)**k21*(c2*(r2(2)-r2(1)) &
       +c3*(r(2)**3-r(1)**3))
  !
  !   solve the linear system with lapack routine dptsv
  !
  call dptsv(mesh-2,1,d(2),e(2),vh(2),mesh-2,ierr)
  call errore('hartree', 'problem in lapack ', ierr)
  !
  ! Set the value of the solution at the first and last point
  ! First, find c0 from the solution in the second point
  !
  c0=vh(2)/sqr(2)**k21-c2*r2(2)-c3*r(2)*r2(2)
  !
  ! and then use the series expansion at the first point
  !
  vh(1)=sqr(1)**k21*(c0+c2*r2(1)+c3*r(1)**3)
  !
  ! the solution at the last point is given  by the boundary 
  ! condition
  !
  vh(mesh)=vh(mesh-1)*fn
  !
  !     The solution must be divided by r (from the equation) 
  !     and multiplied by the square root of r (from the log 
  !     mesh transformation)
  !
  do i=1,mesh
     vh(i)= vh(i) / sqr(i)
  end do

  deallocate(e)
  deallocate(d)

  return
end subroutine hartree
