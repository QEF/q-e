!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine rho2zeta (rho, rho_core, nrxx, nspin, iop)  
  !--------------------------------------------------------------------
  ! if (iopi.eq.1) transform the spin up spin down charge density rho(*,is
  !                into rho(*,1) = ( rho_up + rho_dw ) and
  !                     rho(*,2) = ( rho_up - rho_dw ) / rho_tot = zeta
  ! if (iopi.eq.-1) do the opposit transformation
  !
  use parameters
  implicit none  
  integer :: iop, nspin, nrxx, ir  
  ! the input option
  ! the number of spin polarizations
  ! the fft grid dimension
  ! the counter for fft grid

  real(kind=DP) :: rho (nrxx, nspin), rho_core (nrxx), rho_up, rho_dw, &
       zeta, rhox
  ! the scf charge density
  ! the core charge density
  ! auxiliary variable for rho up
  ! auxiliary variable for rho dw
  ! auxiliary variable for zeta
  ! auxiliary variable for total rho
  real(kind=DP) :: eps  
  ! a small number

  parameter (eps = 1.0d-30)  

  if (nspin.eq.1) return  
  if (iop.eq.1) then  
     do ir = 1, nrxx  
        rhox = rho (ir, 1) + rho (ir, 2) + rho_core (ir)  
        if (rhox.gt.eps) then  
           zeta = (rho (ir, 1) - rho (ir, 2) ) / rhox  
        else  
           zeta = 0.d0  
        endif
        rho (ir, 1) = rho (ir, 1) + rho (ir, 2)  
        rho (ir, 2) = zeta  
     enddo
  elseif (iop.eq. - 1) then  
     do ir = 1, nrxx  
        rhox = rho (ir, 1) + rho_core (ir)  
        if (rhox.gt.eps) then  
           rho_up = 0.5d0 * (rho (ir, 1) + rho (ir, 2) * rhox)  
           rho_dw = 0.5d0 * (rho (ir, 1) - rho (ir, 2) * rhox)  
           rho (ir, 1) = rho_up  
           rho (ir, 2) = rho_dw  
        else  
           rho (ir, 1) = 0.d0  
           rho (ir, 2) = 0.d0  

        endif
     enddo
  else  
     write ( * , * ) ' iop =', iop  
     call error ('mag2zeta', 'wrong iop', 1)  

  endif
  return  

end subroutine rho2zeta
