!
! Copyright (C) 2010-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE add_vhub_to_deeq_acc( deeq )
!-----------------------------------------------------------------
  !! Add Hubbard contributions to the integral of \(V_\text{eff}\)
  !! and \(Q_{nm}\) when Hubbard_projectors is pseudo.
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,    ONLY : nh, nhm
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : v
  USE ldaU,          ONLY : is_hubbard, Hubbard_l, offsetU, q_ae
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), INTENT(INOUT) :: deeq(nhm,nhm,nat,nspin)
  !! integral of \(V_\text{eff}\) and \(Q_{nm}\)
  !
  !  ... local variables
  !
  REAL(KIND=DP), ALLOCATABLE :: deeq_aux(:,:,:)
  INTEGER :: na, nt, ih, jh, ijh, m1, m2, ow1, ow2, is, nhnt
  !
  ALLOCATE( deeq_aux(nhm,nhm,nspin) )
  !$acc data create(deeq_aux) present(deeq)
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     ! skip atoms without Hubbard U
     IF ( .NOT. is_hubbard(nt) ) CYCLE
     !
     deeq_aux = 0.d0
     DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 ! This loop should be replaced by matrix multiplication:
                 !    deeq_aux = q_ae * v%ns * q_ae 
                 ow1 = offsetU(na)+m1
                 ow2 = offsetU(na)+m2
                 deeq_aux(ih,jh,1:nspin) = deeq_aux(ih,jh,1:nspin) + &
                    v%ns(m1,m2,1:nspin,na)*q_ae(ow1,ih,na)*q_ae(ow2,jh,na)
                 !
              ENDDO 
           ENDDO
           !
           deeq_aux(jh,ih,1:nspin) = deeq_aux(ih,jh,1:nspin)
           !
        ENDDO
     ENDDO
     !$acc update device(deeq_aux)
     nhnt = nh(nt)
     !$acc parallel loop collapse(3)
     DO is = 1, nspin
        DO ih = 1, nhnt
           DO jh = 1, nhnt
              deeq(jh,ih,na,is) = deeq(jh,ih,na,is) + deeq_aux(jh,ih,is)
           END DO
        END DO
     END DO
     !
  ENDDO
  !$acc end data
  DEALLOCATE(deeq_aux)
  !
END SUBROUTINE add_vhub_to_deeq_acc
