!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine set_rc_rv()
  !-----------------------------------------------------------------------
  !
  !      input : all-electron wavefunctions + valence states
  !      output: separated core and valence charges 
  !
  use kinds, only : dp
  use ld1_parameters, only : nwfx
  
  use ld1inc, only : grid, aeccharge, aevcharge, nwf, oc, isw, rel, psi, &
                     core_state
  implicit none

  integer :: n, ns, is
  !
  !      calculates core charge density
  !
  aevcharge=0.0_DP
  aeccharge=0.0_DP
  do n=1,grid%mesh
     do ns=1,nwf
        if (oc(ns)>0.0_DP) then
           is=isw(ns)
           if (rel==2) then
              if (core_state(ns)) then
                 aeccharge(n)=aeccharge(n) &
                              +oc(ns)*( psi(n,1,ns)**2 + psi(n,2,ns)**2 )
              else
                 aevcharge(n,is)=aevcharge(n,is)+oc(ns)*(psi(n,1,ns)**2 &
                                                       + psi(n,2,ns)**2)
              endif
           else
              if (core_state(ns)) then
                 aeccharge(n) = aeccharge(n) + oc(ns)*psi(n,1,ns)**2
              else
                 aevcharge(n,is) = aevcharge(n,is) + oc(ns)*psi(n,1,ns)**2
              endif
           endif
        endif
     enddo
  enddo
  return
end subroutine set_rc_rv
