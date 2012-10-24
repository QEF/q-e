!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE add_vhub_to_deeq(deeq)
  !
  ! Add Hubbard contributions to deeq when U_projection is pseudo
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,    ONLY : nh, nhm
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : v
  USE ldaU,          ONLY : Hubbard_U, Hubbard_alpha, &
                            Hubbard_l, oatwfc, q_ae
  IMPLICIT NONE
  REAL(KIND=DP), INTENT(INOUT) :: deeq( nhm, nhm, nat, nspin )
  INTEGER :: na, nt, ih, jh, ijh, m1, m2, ow1, ow2
  !
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     ! skip atoms without Hubbard U
     IF ( Hubbard_U(nt)==0.D0 .AND. Hubbard_alpha(nt)==0.D0 ) CYCLE
     ! 
     DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
           !
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 !
                 ow1 = oatwfc(na)+m1
                 ow2 = oatwfc(na)+m2
                 deeq(ih,jh,na,1:nspin) = deeq(ih,jh,na,1:nspin) + &
                    v%ns(m1,m2,1:nspin,na)*q_ae(ow1,ih,na)*q_ae(ow2,jh,na)
                 !
              ENDDO 
           ENDDO
           !
           deeq(jh,ih,na,1:nspin) = deeq(ih,jh,na,1:nspin)
           !
        ENDDO
     ENDDO
     !
  ENDDO
  !
END SUBROUTINE add_vhub_to_deeq
