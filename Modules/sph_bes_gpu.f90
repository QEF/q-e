!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------

module sph_bes_gpum
contains
#if defined(__CUDA)
attributes(DEVICE) subroutine sph_bes_gpu (msh, r, q, l, jl)
#else
subroutine sph_bes_gpu (msh, r, q, l, jl)
#endif
  !--------------------------------------------------------------------
  !
  ! ... input:
  ! ...   msh     = number of grid points points
  ! ...   r(1:msh)= radial grid
  ! ...   q       = q
  ! ...   l       = angular momentum (-1 <= l <= 6)
  ! ... output:
  ! ...   jl(1:msh) = j_l(q*r(i))  (j_l = spherical bessel function)
  !
  use kinds, only: DP
  USE constants, ONLY : eps14
  !
  implicit none
  !
  integer, VALUE :: msh, l
  real(DP), VALUE :: q
  real(DP) :: r (msh), jl (msh)
#if defined(__CUDA)
  attributes(DEVICE) :: r, jl
#endif
  !
  ! xseries = convergence radius of the series for small x of j_l(x)
  real(DP) :: x, xl, xseries = 0.05_dp
  integer :: ir, ir0, ir_, i, n, semifact
  !
  ! 
  !  case q=0
  if (abs (q) < eps14) then
     do ir = 1, msh
     if (l == -1) then
        !call errore ('sph_bes', 'j_{-1}(0) ?!?', 1)
     elseif (l == 0) then
           jl(ir) = 1.d0
     else
           jl(ir) = 0.d0
     endif
     enddo
     return
  end if 

  !  case l=-1

  if (l == - 1) then
     ! if (abs (q * r (1) ) < eps14) call errore ('sph_bes', 'j_{-1}(0) ?!?',1)
     do ir = 1, msh
        jl(ir) = cos (q * r (ir) ) / (q * r (ir) )
     end do
     return
  end if

  ! series expansion for small values of the argument
  ! ir0 is the first grid point for which q*r(ir0) > xseries
  ! notice that for small q it may happen that q*r(msh) < xseries !

  ir0 = msh+1
  do ir = 1, msh
     if ( abs (q * r (ir) ) > xseries ) then
        ir0 = ir
        exit
     end if
  end do

  do ir = 1, ir0 - 1
     x = q * r (ir)
     if ( l == 0 ) then
        xl = 1.0_dp
     else
        xl = x**l
     end if

     semifact = 1
     n = 2*l+1
     do i = n, 1, -2
        semifact = i*semifact
     end do

     jl (ir) = xl/semifact * &
                ( 1.0_dp - x**2/1.0_dp/2.0_dp/(2.0_dp*l+3) * &
                ( 1.0_dp - x**2/2.0_dp/2.0_dp/(2.0_dp*l+5) * &
                ( 1.0_dp - x**2/3.0_dp/2.0_dp/(2.0_dp*l+7) * &
                ( 1.0_dp - x**2/4.0_dp/2.0_dp/(2.0_dp*l+9) ) ) ) )
  end do

  ! the following shouldn't be needed but do you trust compilers
  ! to do the right thing in this special case ? I don't - PG

  if ( ir0 > msh ) return

  do ir=ir0, msh
     if (l == 0) then
        jl (ir) = sin (q * r (ir) ) / (q * r (ir) )
     elseif (l == 1) then
        jl (ir) = (sin (q * r (ir) ) / (q * r (ir) ) - &
                     cos (q * r (ir) ) ) / (q * r (ir) )
     elseif (l == 2) then
     
        jl (ir) = ( (3.d0 / (q*r(ir)) - (q*r(ir)) ) * sin (q*r(ir)) - &
                       3.d0 * cos (q*r(ir)) ) / (q*r(ir))**2
     
     elseif (l == 3) then
     
     
        jl (ir) = (sin (q*r(ir)) * &
                     (15.d0 / (q*r(ir)) - 6.d0 * (q*r(ir)) ) + &
                     cos (q*r(ir)) * ( (q*r(ir))**2 - 15.d0) ) / &
                     (q*r(ir)) **3
     
     elseif (l == 4) then
     
     
        jl (ir) = (sin (q*r(ir)) * &
                     (105.d0 - 45.d0 * (q*r(ir))**2 + (q*r(ir))**4) + &
                     cos (q*r(ir)) * &
                     (10.d0 * (q*r(ir))**3 - 105.d0 * (q*r(ir))) ) / &
                        (q*r(ir))**5
     
     elseif (l == 5) then
     
        jl (ir) = (-cos(q*r(ir)) - &
                     (945.d0*cos(q*r(ir))) / (q*r(ir)) ** 4 + &
                     (105.d0*cos(q*r(ir))) / (q*r(ir)) ** 2 + &
                     (945.d0*sin(q*r(ir))) / (q*r(ir)) ** 5 - &
                     (420.d0*sin(q*r(ir))) / (q*r(ir)) ** 3 + &
                     ( 15.d0*sin(q*r(ir))) / (q*r(ir)) ) / (q*r(ir))
     
     elseif (l == 6) then
     
     
        jl (ir) = ((-10395.d0*cos(q*r(ir))) / (q*r(ir))**5 + &
                     (  1260.d0*cos(q*r(ir))) / (q*r(ir))**3 - &
                     (    21.d0*cos(q*r(ir))) / (q*r(ir))    - &
                                sin(q*r(ir))                   + &
                     ( 10395.d0*sin(q*r(ir))) / (q*r(ir))**6 - &
                     (  4725.d0*sin(q*r(ir))) / (q*r(ir))**4 + &
                     (   210.d0*sin(q*r(ir))) / (q*r(ir))**2 ) / (q*r(ir))
     
     !else
     !
     !   call errore ('sph_bes', 'not implemented', abs(l))

     endif
     !if (l > 6 ) call errore ('sph_bes', 'not implemented', abs(l))
  end do 
  !
  return
end subroutine sph_bes_gpu
end module
