!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*******************************************************************

subroutine magnetic_kkterm(ik,outvv, outhh)

!*******************************************************************

  use kinds,                only: DP
  use wvfct,                only: npw, npwx, nbnd
  use wavefunctions_module, only: evc
  use klist,                only: wk
  use control_ph,           only: alpha_pv

  implicit none

!  input
  integer, intent(in)      :: ik         ! kpoint

! output

  complex(DP), intent(out) :: outvv(3,3), outhh(3,3)

! local

  integer :: p0,p1,             &  !direction
       ibnd,i
  
  complex(DP), allocatable :: dev(:,:), dpsi(:,:) 
  complex(DP), external :: ZDOTC
  

  allocate (dev(npwx,nbnd))
  allocate (dpsi(npwx,nbnd))

! no shift here
  alpha_pv = 0.d0


  do p0=1,3
    
     dev=(0.d0,0.d0)
     call take_nloc_k_kq(ik, ik, evc, p0, dev )
     
     call solve_cg(dev,ik, dpsi)

     do p1=1,3
        dev=(0.d0,0.d0)
        call take_nloc_k_kq(ik, ik, dpsi, p1, dev)
        do ibnd=1,nbnd
           outvv(p0,p1) = outvv(p0,p1) + &
                ZDOTC(npw,evc(1,ibnd),1,dev(1,ibnd),1) &
                * wk(ik) / DBLE(nbnd)
        enddo

        dev=(0.d0,0.d0)
        call grad(ik, dpsi, p1, dev)
        do ibnd=1,nbnd
           outhh(p0,p1) = outhh(p0,p1) + &
                ZDOTC(npw,evc(1,ibnd),1,dev(1,ibnd),1) &
                * wk(ik) / DBLE(nbnd)
        enddo
        
     enddo
  enddo
  
  deallocate(dpsi)
  deallocate(dev)
  

end subroutine magnetic_kkterm
