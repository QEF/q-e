!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine sph_besr (msh, r, q, l, jl)  
  !--------------------------------------------------------------------
  !
  ! input:
  !   msh     = number of grid points points
  !   r(1:msh)= radial grid
  !   q       = q
  !   l       = angular momentum (-1 <= l <= 6)
  ! output:
  !   jl(1:msh) = j_l(q*r(i))  (j_l = spherical bessel function)
  !
  implicit none  
  integer,parameter :: dp=kind(1.d0)
  !
  integer :: msh, l  
  real(kind=DP) :: r (msh), q, jl (msh)  
  !
  real(kind=DP), parameter :: eps = 1.d-8
  integer :: ir, ir0  
  !
  if (abs (q) < eps) then  
     if (l == -1) then  
        call errore ('sph_besr', 'j_{-1}(0) ?!?', 1)  
     elseif (l == 0) then  
        do ir = 1, msh  
           jl (ir) = r(ir)   
        enddo
     else  
        jl(:) = 0.d0
     endif
  else  
     if (abs (q * r (1) ) > eps) then  
        ir0 = 1  
     else  
        if (l == -1) then  
           call errore ('sph_besr', 'j_{-1}(0) ?!?', 2)  
        elseif (l == 0) then  
           jl (1) = r(1)  
        else  
           jl (1) = 0.d0  
        endif
        ir0 = 2 
     endif
     if (l == - 1) then  
        do ir = ir0, msh  
           jl (ir) = cos (q * r (ir) ) / q   
        enddo
     elseif (l == 0) then  
        do ir = ir0, msh  
           jl (ir) = sin (q * r (ir) ) / q  
        enddo
     elseif (l == 1) then  
        do ir = ir0, msh  
           if (abs(q*r(ir)).lt.6.d-2) then
              jl(ir)=(q*r(ir))/3.d0-(q*r(ir))**3/30.d0+ &
                     (q*r(ir))**5/840.d0-(q*r(ir))**7/45360.d0
              jl(ir) =  jl(ir)*r(ir)
           else
              jl(ir)=(sin(q*r(ir))/(q*r(ir))-cos(q*r(ir)))/q
           endif
        enddo
     elseif (l == 2) then  
        do ir = ir0, msh  
           if (abs(q*r(ir)).lt.6.d-2) then
              jl(ir) = (q*r(ir))**2/15.d0-(q*r(ir))**4/210.d0 &
                     +(q*r(ir))**6/7560.d0
              jl(ir) =  jl(ir)*r(ir)
           else
              jl(ir) = ( ( 3.d0-(q*r(ir))**2)*sin(q*r(ir)) &
               - 3.d0 * q*r(ir)*cos(q*r(ir)) )/(q*r(ir))**2/q
           endif
        enddo
     elseif (l == 3) then  
        do ir = ir0, msh  
           if (abs(q*r(ir)).lt.6.d-2) then
               jl(ir) = 48.d0*(q*r(ir))**3/5040.d0 &
                      - 192.d0*(q*r(ir))**5/362880.d0
               jl(ir) = jl(ir) * r(ir)
           else
              jl(ir) = (sin(q*r(ir)) * (15.d0/(q*r(ir))      &
                - 6.d0*(q*r(ir))) + cos(q*r(ir))*((q*r(ir))**2 &
                - 15.d0) ) / (q * r (ir) ) **2/q
           endif
        enddo
     elseif (l == 4) then  
        do ir = ir0, msh  
           jl (ir) = (sin (q * r (ir) ) * (105.d0 - 45.d0 * (q * r (ir) &
                ) **2 + (q * r (ir) ) **4) + cos (q * r (ir) ) * (10.d0 * &
                (q * r (ir) ) **3 - 105.d0 * (q * r (ir) ) ) ) / (q * r (ir) &
                ) **4/q
        enddo
     elseif (l == 5) then  
        do ir = ir0, msh  
           jl (ir) = (-cos(q*r(ir)) - (945.d0*cos(q*r(ir)))/(q*r(ir)) ** 4 + &
                (105.d0*cos(q*r(ir)))/ (q*r(ir)) ** 2 + &
                (945.d0*sin(q*r(ir)))/ (q*r(ir)) ** 5 - &
                (420.d0*sin(q*r(ir)))/(q*r(ir)) ** 3 +  &
                (15.d0*sin(q*r(ir)))/ (q*r(ir)))/q
        end do
     elseif (l == 6) then  
        do ir = ir0, msh  
           jl (ir) = ((-10395.d0*cos(q*r(ir)))/(q*r(ir)) ** 5 + &
                (1260.d0*cos(q*r(ir)))/(q*r(ir)) ** 3 - &
                (21.d0*cos(q*r(ir)))/ (q*r(ir)) - sin(q*r(ir)) + &
                (10395.d0*sin(q*r(ir)))/(q*r(ir)) ** 6 - &
                (4725.d0*sin(q*r(ir)))/ (q*r(ir)) ** 4 + &
                (210.d0*sin(q*r(ir)))/(q*r(ir)) ** 2)/q
        end do
     else  
        call errore ('sph_besr', 'not implemented', l)  
     endif
  endif
  !
  return  
end subroutine sph_besr
