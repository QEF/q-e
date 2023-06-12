!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE init_tab_rho (omega, comm)
  !----------------------------------------------------------------------
  !
  ! Compute interpolation table for atomic charge density
  ! tab_rho(i,nt) = \rho_nt(q_i) for atom of type nt, on grid q_i
  !! Output in tab_rhc, module uspp_data
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : fpi, eps8
  USE atom,         ONLY : rgrid, msh
  USE uspp_param,   ONLY : upf, nsp
  USE uspp_data,    ONLY : nqxq, dq, tab_rho
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: comm
  !! MPI communicator, to split the workload
  REAL(dp), INTENT(IN) :: omega
  !! Unit-cell volume
  !
  INTEGER :: ndm, startq, lastq, nt, iq, ir
  !! Various indices
  REAL(dp) :: q
  REAL(dp), ALLOCATABLE :: aux (:)
  !! Work space
  !
  ndm = MAXVAL( msh(1:nsp) )
  ALLOCATE (aux(ndm))
  !
  CALL divide (comm, nqxq, startq, lastq)
  !
  DO nt = 1, nsp
     !
     tab_rho(:,nt)= 0.d0
     !
     DO iq = startq, lastq
        !
        q = (iq - 1) * dq
        DO ir = 1, msh(nt)
           IF ( iq > 1 .AND. rgrid(nt)%r(ir) > eps8 ) then
              aux (ir) = upf(nt)%rho_at(ir) * & 
                        sin(q*rgrid(nt)%r(ir)) / (q*rgrid(nt)%r(ir))
           ELSE
              aux (ir) = upf(nt)%rho_at(ir)
           ENDIF 
        ENDDO
        !
        CALL simpson ( msh(nt), aux, rgrid(nt)%rab, tab_rho(iq,nt) )
        tab_rho (iq,nt) = tab_rho (iq,nt) / omega 
        !
     ENDDO
     !
  END DO
  !
  CALL mp_sum ( tab_rho (:,1:nsp), comm )
!$acc enter data copyin (tab_rho)
  !
  DEALLOCATE (aux)
  !
END SUBROUTINE init_tab_rho
