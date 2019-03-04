!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE add_paw_to_deeq( deeq )
  !
  !! Add paw contributions to the integral of the perturbed potential 
  !! with the Q function (computed in paw_potential).
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,    ONLY : upf, nh, nhm
  USE paw_variables, ONLY : okpaw, ddd_paw
  USE lsda_mod,      ONLY : nspin
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), INTENT(INOUT) :: deeq( nhm, nhm, nat, nspin )
  !! integral of the perturbed potential
  !
  ! ... local variables
  !
  INTEGER :: na, nt, ih, jh, ijh
  !
  IF (okpaw) THEN
     DO na=1,nat
        nt = ityp(na)
        IF (.not.upf(nt)%tpawp) CYCLE
        ijh=0
        DO ih=1,nh(nt)
           DO jh=ih,nh(nt)
              ijh=ijh+1
              deeq(ih,jh,na,1:nspin) = deeq(ih,jh,na,1:nspin) + &
                                        ddd_paw(ijh,na,1:nspin)
              deeq(jh,ih,na,1:nspin) = deeq(ih,jh,na,1:nspin) 
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  RETURN
  !
END SUBROUTINE add_paw_to_deeq
