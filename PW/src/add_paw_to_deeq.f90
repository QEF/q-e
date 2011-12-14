!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

SUBROUTINE add_paw_to_deeq(deeq)
     ! Add paw contributions to deeq (computed in paw_potential)
  USE kinds,                ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,    ONLY : upf, nh, nhm
  USE paw_variables, ONLY : okpaw, ddd_paw
  USE lsda_mod,      ONLY : nspin
  IMPLICIT NONE
  integer :: na, nt, ih, jh, ijh
  REAL(kind=dp), intent(inout) :: deeq( nhm, nhm, nat, nspin )

  if (okpaw) then
     do na=1,nat
        nt = ityp(na)
        IF (.not.upf(nt)%tpawp) cycle
        ijh=0
        do ih=1,nh(nt)
           do jh=ih,nh(nt)
              ijh=ijh+1
              deeq(ih,jh,na,1:nspin) = deeq(ih,jh,na,1:nspin) &
                                     + ddd_paw(ijh,na,1:nspin)
              deeq(jh,ih,na,1:nspin) = deeq(ih,jh,na,1:nspin) 
           end do
        end do
     end do
  end IF 
  RETURN
  
END SUBROUTINE add_paw_to_deeq
