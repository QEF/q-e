!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rho2zeta( rho, rho_core, nrxx, nspin, iop )
  !---------------------------------------------------------------------------
  !! * if ( iopi == 1 )  it transforms rho(:,2:nspin) into zeta:
  !!   \[ \text{rho}(:,2:\text{nspin}) = \text{rho}(:,2:\text{nspin})/
  !!      \text{rho_tot}(:) = \text{zeta}(:,2:\text{nspin}) \]
  !! * if ( iopi == -1)  it does the opposite transformation.
  !
  USE constants,   ONLY: eps32
  USE io_global,   ONLY: stdout
  USE kinds,       ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: iop, nspin, nrxx
  REAL(DP), INTENT(IN) :: rho_core(nrxx)
  REAL(DP), INTENT(INOUT) :: rho(nrxx,nspin)
  !
  INTEGER :: is
  !
  IF ( nspin == 1 ) RETURN  
  !
  IF ( iop == -1 ) THEN
     !
     DO is = 2, nspin
        WHERE( rho(:,1)+rho_core(:) > eps32 )
           !
           rho(:,is) = rho(:,is)*(rho(:,1)+rho_core(:))
           !
        ELSEWHERE
           !
           rho(:,is) = 0.0_DP
           !
        END WHERE
     END DO
     !
  ELSE IF ( iop == 1 ) THEN
        !
     DO is = 2, nspin
        WHERE( rho(:,1)+rho_core(:) > eps32 )
           !
           rho(:,is) = rho(:,is) / ( rho(:,1)+rho_core(:) )
           !
        ELSEWHERE
           !
           rho(:,is) = 0.0_DP
           !
        END WHERE
     END DO
     !
  ELSE
     !
     WRITE ( stdout , '(5X,"iop = ",I5)' ) iop
     !
     CALL errore( 'rho2zeta', 'wrong iop', 1 )
     !
  END IF
  !
END SUBROUTINE rho2zeta
