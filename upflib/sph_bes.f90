!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine sph_bes (msh, r, q, l, jl)
  !$acc routine vector
  !--------------------------------------------------------------------
  !! Spherical Bessel function.
  !
  USE upf_kinds, only: DP
  USE upf_const, only: eps14
  !
  implicit none
  !
  integer :: msh
  !! number of grid points points
  integer :: l
  !! angular momentum (-1 <= l <= 6)
  real(DP) :: r (msh)
  !! radial grid
  real(DP) :: q
  !! q
  real(DP) :: jl (msh)
  !! Output: Spherical Bessel function \(j_l(q*r(i))\)
  !
  ! xseries = convergence radius of the series for small x of j_l(x)
  real(DP) :: x, xl, xseries = 0.05_dp
  integer :: i, ir, ir0
  integer :: semifact
  !
#if defined (__MASS)
  real(DP) :: qr(msh), sin_qr(msh), cos_qr(msh)
#endif
 
  !  case q=0

  if (abs (q) < eps14) then
     if (l == -1) then
        stop !call upf_error ('sph_bes', 'j_{-1}(0) ?!?', 1)
     elseif (l == 0) then
        !$acc loop vector
        do ir = 1, msh
          jl(ir) = 1.d0
        enddo
     else
        !$acc loop vector
        do ir = 1, msh
          jl(ir) = 0.d0
        enddo  
     endif
     return
  end if 

  !  case l=-1

  if (l == - 1) then
     if (abs (q * r (1) ) < eps14) stop !call upf_error ('sph_bes', 'j_{-1}(0) ?!?',1)

#if defined (__MASS)

     qr = q * r
     call vcos( cos_qr, qr, msh)
     jl = cos_qr / qr

#else
     !$acc loop vector
     do ir = 1, msh
       jl (ir) = cos (q * r (ir) ) / (q * r (ir) )
     enddo
#endif

     return

  end if

  ! series expansion for small values of the argument
  ! ir0 is the first grid point for which q*r(ir0) > xseries
  ! notice that for small q it may happen that q*r(msh) < xseries !

  ir0 = msh+1
  !$acc loop vector
  do ir = 1, msh
     if ( abs (q * r (ir) ) > xseries ) then
        ir0 = ir
        exit
     end if
  end do

  !$acc loop vector
  do ir = 1, ir0 - 1
     x = q * r (ir)
     if ( l == 0 ) then
        xl = 1.0_dp
     else
        xl = x**l
     end if
     !--
     semifact = 1
     !$acc loop seq reduction(*:semifact)
     do i = 2*l+1, 1, -2
       semifact = i*semifact
     enddo
     !---
     jl (ir) = xl/DBLE(semifact) * &
                ( 1.0_dp - x**2/1.0_dp/2.0_dp/DBLE(2*l+3) * &
                ( 1.0_dp - x**2/2.0_dp/2.0_dp/DBLE(2*l+5) * &
                ( 1.0_dp - x**2/3.0_dp/2.0_dp/DBLE(2*l+7) * &
                ( 1.0_dp - x**2/4.0_dp/2.0_dp/DBLE(2*l+9) ) ) ) )
  end do

  ! the following shouldn't be needed but do you trust compilers
  ! to do the right thing in this special case ? I don't - PG

  if ( ir0 > msh ) return

  if (l == 0) then

#if defined (__MASS)

     qr = q * r
     call vsin( sin_qr, qr, msh)
     jl (ir0:) = sin_qr(ir0:) / (q * r (ir0:) )

#else
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = sin (q * r (ir) ) / (q * r (ir) )
     enddo
#endif

  elseif (l == 1) then

#if defined (__MASS)

     qr = q * r
     call vcos( cos_qr, qr, msh)
     call vsin( sin_qr, qr, msh)
     jl (ir0:) = ( sin_qr(ir0:) / (q * r (ir0:) ) - &
                   cos_qr(ir0:) ) / (q * r (ir0:) )

#else
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (sin (q * r (ir) ) / (q * r (ir) ) - &
                  cos (q * r (ir) ) ) / (q * r (ir) )
     enddo
#endif

  elseif (l == 2) then

#if defined (__MASS)

     qr = q * r
     call vcos( cos_qr, qr, msh)
     call vsin( sin_qr, qr, msh)
     jl (ir0:) = ( (3.d0 / (q*r(ir0:)) - (q*r(ir0:)) ) * sin_qr(ir0: ) - &
                    3.d0 * cos_qr(ir0:) ) / (q*r(ir0:))**2

#else
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = ( (3.d0 / (q*r(ir)) - (q*r(ir)) ) * sin (q*r(ir)) - &
                    3.d0 * cos (q*r(ir)) ) / (q*r(ir))**2
     enddo
#endif

  elseif (l == 3) then

#if defined (__MASS)

     qr = q * r
     call vcos( cos_qr, qr, msh)
     call vsin( sin_qr, qr, msh)
     jl (ir0:) = (sin_qr (ir0:) * &
                  (15.d0 / (q*r(ir0:)) - 6.d0 * (q*r(ir0:)) ) + &
                  cos_qr (ir0:) * ( (q*r(ir0:))**2 - 15.d0) ) / &
                  (q*r(ir0:))**3

#else
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (sin (q*r(ir)) * &
                 (15.d0 / (q*r(ir)) - 6.d0 * (q*r(ir)) ) + &
                 cos (q*r(ir)) * ( (q*r(ir))**2 - 15.d0) ) / &
                 (q*r(ir)) **3
     enddo
#endif

  elseif (l == 4) then

#if defined (__MASS)

     qr = q * r
     call vcos( cos_qr, qr, msh)
     call vsin( sin_qr, qr, msh)
     jl (ir0:) = (sin_qr (ir0:) * &
                  (105.d0 - 45.d0 * (q*r(ir0:))**2 + (q*r(ir0:))**4) + &
                  cos_qr (ir0:) * &
                  (10.d0 * (q*r(ir0:))**3 - 105.d0 * (q*r(ir0:))) ) / &
                    (q*r(ir0:))**5

#else
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (sin (q*r(ir)) * &
                 (105.d0 - 45.d0 * (q*r(ir))**2 + (q*r(ir))**4) + &
                 cos (q*r(ir)) * &
                 (10.d0 * (q*r(ir))**3 - 105.d0 * (q*r(ir))) ) / &
                    (q*r(ir))**5
     enddo
#endif

  elseif (l == 5) then

#if defined (__MASS)
     qr = q * r
     call vcos( cos_qr, qr, msh)
     call vsin( sin_qr, qr, msh)
     jl (ir0:) = (-cos_qr(ir0:) - &
                  (945.d0*cos_qr(ir0:)) / (q*r(ir0:)) ** 4 + &
                  (105.d0*cos_qr(ir0:)) / (q*r(ir0:)) ** 2 + &
                  (945.d0*sin_qr(ir0:)) / (q*r(ir0:)) ** 5 - &
                  (420.d0*sin_qr(ir0:)) / (q*r(ir0:)) ** 3 + &
                  ( 15.d0*sin_qr(ir0:)) / (q*r(ir0:)) ) / (q*r(ir0:))
#else
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (-cos(q*r(ir)) - &
                  (945.d0*cos(q*r(ir))) / (q*r(ir)) ** 4 + &
                  (105.d0*cos(q*r(ir))) / (q*r(ir)) ** 2 + &
                  (945.d0*sin(q*r(ir))) / (q*r(ir)) ** 5 - &
                  (420.d0*sin(q*r(ir))) / (q*r(ir)) ** 3 + &
                  ( 15.d0*sin(q*r(ir))) / (q*r(ir)) ) / (q*r(ir))
     enddo
#endif

  elseif (l == 6) then

#if defined (__MASS)

     qr = q * r
     call vcos( cos_qr, qr, msh)
     call vsin( sin_qr, qr, msh)
     jl (ir0:) = ((-10395.d0*cos_qr(ir0:)) / (q*r(ir0:))**5 + &
                  (  1260.d0*cos_qr(ir0:)) / (q*r(ir0:))**3 - &
                  (    21.d0*cos_qr(ir0:)) / (q*r(ir0:))    - &
                             sin_qr(ir0:)                   + &
                  ( 10395.d0*sin_qr(ir0:)) / (q*r(ir0:))**6 - &
                  (  4725.d0*sin_qr(ir0:)) / (q*r(ir0:))**4 + &
                  (   210.d0*sin_qr(ir0:)) / (q*r(ir0:))**2 ) / (q*r(ir0:))
#else
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = ((-10395.d0*cos(q*r(ir))) / (q*r(ir))**5 + &
                  (  1260.d0*cos(q*r(ir))) / (q*r(ir))**3 - &
                  (    21.d0*cos(q*r(ir))) / (q*r(ir))    - &
                             sin(q*r(ir))                   + &
                  ( 10395.d0*sin(q*r(ir))) / (q*r(ir))**6 - &
                  (  4725.d0*sin(q*r(ir))) / (q*r(ir))**4 + &
                  (   210.d0*sin(q*r(ir))) / (q*r(ir))**2 ) / (q*r(ir))
     enddo
#endif

  else

     stop !call upf_error ('sph_bes', 'not implemented', abs(l))

  endif
  !
  return
end subroutine sph_bes

integer function semifact(n)
  ! semifact(n) = n!!
  implicit none
  integer :: n, i

  semifact = 1
  do i = n, 1, -2
     semifact = i*semifact
  end do
  return
end function semifact
!
SUBROUTINE sph_dbes ( nr, r, xg, l, jl, djl )
  !
  !! Calculates \(x*dj_l(x)/dx\) using the recursion formula:
  !! $$ dj_l(x)/dx = l/x j_l(x) - j_{l+1}(x) $$
  !! for \(l=0\), and for \(l>0\):
  !! $$ dj_l(x)/dx = j_{l-1}(x) - (l+1)/x j_l(x) $$
  !! Requires \(j_l(r)\) in input.
  !! Used only in CP. Note that upflib uses numerical differentiation.
  !
  USE upf_kinds, only: DP
  USE upf_const, ONLY : eps8
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l, nr
  REAL (DP), INTENT(IN) :: xg, jl(nr), r(nr)
  REAL (DP), INTENT(OUT):: djl(nr)
  !
  if ( xg < eps8 ) then
     !
     ! special case q=0
     ! note that x*dj_l(x)/dx = 0 for x = 0
     !
     djl(:) = 0.0d0
  else
     !
     if ( l > 0 ) then
        call sph_bes ( nr, r, xg, l-1, djl )
        djl(:) = djl(:) * (xg * r(:) ) - (l+1) * jl(:)
     else if ( l == 0 ) then
        call sph_bes ( nr, r, xg, l+1, djl )
        djl(:) = - djl(:) * (xg * r(:) )
     else
        call upf_error('sph_dbes','l < 0 not implemented', abs(l) )
     end if
  end if
  !
end SUBROUTINE sph_dbes
