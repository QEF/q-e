!
! Copyright (C) 2002-2003 CP group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Machine-dependent routines for:
!    erf, erfc, freq functions
!
! ================
!     for machines that do not have these routines in the math libraries
!
#if defined __INTEL || defined __PGI
!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
real(kind=8) function erf (x)  
  !---------------------------------------------------------------------
  !
  !     Error function - computed from the rational approximations of
  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
  !
  !     for abs(x) le 0.47 erf is calculated directly
  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
  !
  implicit none  
  real(kind=8) :: x, x2, p1 (4), q1 (4), erfc  
  external erfc  
  data p1 / 2.42667955230532d2, 2.19792616182942d1, &
       6.99638348861914d0, - 3.56098437018154d-2 /
  data q1 / 2.15058875869861d2, 9.11649054045149d1, &
       1.50827976304078d1, 1.00000000000000d0 /
  !
  if (abs (x) .gt.6.d0) then  
     !
     !  erf(6)=1-10^(-17) cannot be distinguished from 1 with 16-byte words
     !
     erf = sign (1.d0, x)  
  else  
     if (abs (x) .le.0.47d0) then  
        x2 = x**2  
        erf = x * (p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 ( &
             4) ) ) ) / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 ( &
             4) ) ) )
     else  
        erf = 1.d0 - erfc (x)  
     endif
  endif
  !
  return  
end function erf
!
!---------------------------------------------------------------------
real(kind=8) function erfc (x)  
  !---------------------------------------------------------------------
  !
  !     erfc(x) = 1-erf(x)  - See comments in erf
  !
  implicit none  
  real(kind=8) :: x, ax, x2, xm2, erf, p2 (8), q2 (8), p3 (5), q3 (5), &
       pim1
  external erf  
  data p2 / 3.00459261020162d2, 4.51918953711873d2, &
       3.39320816734344d2, 1.52989285046940d2, 4.31622272220567d1, &
       7.21175825088309d0, 5.64195517478974d-1, - 1.36864857382717d-7 /
  data q2 / 3.00459260956983d2, 7.90950925327898d2, &
       9.31354094850610d2, 6.38980264465631d2, 2.77585444743988d2, &
       7.70001529352295d1, 1.27827273196294d1, 1.00000000000000d0 /
  data p3 / - 2.99610707703542d-3, - 4.94730910623251d-2, - &
       2.26956593539687d-1, - 2.78661308609648d-1, - 2.23192459734185d-2 &
       /
  data q3 / 1.06209230528468d-2, 1.91308926107830d-1, &
       1.05167510706793d0, 1.98733201817135d0, 1.00000000000000d0 /

  data pim1 / 0.564189583547756d0 /  
  !        ( pim1= sqrt(1/pi) )
  ax = abs (x)  
  if (ax.gt.26.d0) then  
     !
     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
     !
     erfc = 0.d0  
  elseif (ax.gt.4.d0) then  
     x2 = x**2  
     xm2 = (1.d0 / ax) **2  
     erfc = (1.d0 / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
          (q3 (4) + xm2 * q3 (5) ) ) ) ) )
  elseif (ax.gt.0.47d0) then  
     x2 = x**2  
     erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
  else  
     erfc = 1.d0 - erf (ax)  
  endif
  !
  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
  !
  if (x.lt.0.d0) erfc = 2.d0 - erfc  
  !
  return  
end function erfc

#endif

!---------------------------------------------------------------------
real(8) function gauss_freq (x)
  !---------------------------------------------------------------------
  !
  !     gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
  !             - See comments in erf
  !
  real(kind=8) :: x, c, erf, erfc
  external erf

  data c / 0.707106781186548d0 /
  !        ( c= sqrt(1/2) )
  gauss_freq = 0.5d0 * erfc ( - x * c)
  !
  return
end function gauss_freq


