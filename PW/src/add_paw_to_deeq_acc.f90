!
! Copyright (C) 2010-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE add_paw_to_deeq_acc(deeq)
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
  REAL(DP), INTENT(INOUT) :: deeq(nhm,nhm,nat,nspin)
  !! integral of the perturbed potential
  !
  ! ... local variables
  !
  INTEGER :: na, nb, nab, nt, ih, jh, ijh, nhnt, is

! OPTIMIZE HERE: squeeze loop on atoms having PAW pseudo

  IF (okpaw) THEN
     !$acc data present(deeq) copyin(ddd_paw)
     DO na=1,nat
        nt = ityp(na)
        IF (.not.upf(nt)%tpawp) CYCLE
        nhnt = nh(nt)
        !$acc parallel loop collapse(3)
        DO is=1,nspin
           DO ih=1,nhnt
              DO jh=1,nhnt
                 IF (jh >= ih) then
                    ijh = jh + ((ih-1)*(2*nhnt-ih))/2
                    deeq(ih,jh,na,is) = deeq(ih,jh,na,is) &
                                           + ddd_paw(ijh,na,is)
                    deeq(jh,ih,na,is) = deeq(ih,jh,na,is) 
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$acc end data
  ENDIF
  !
  RETURN
  !
END SUBROUTINE add_paw_to_deeq_acc
