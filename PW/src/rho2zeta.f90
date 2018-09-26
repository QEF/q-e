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
  !
  ! ... if ( iopi == 1 )  transform the spin up spin down charge density 
  ! ...                   rho(*,is) into :
  !
  ! ...                      rho(*,1) = ( rho_up + rho_dw ) and
  ! ...                      rho(*,2) = ( rho_up - rho_dw ) / rho_tot = zeta
  !
  ! ... if ( iopi == -1)  do the opposit transformation
  !
  USE constants, ONLY : eps32
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER :: iop, nspin, nrxx, ir
    ! the input option
    ! the number of spin polarizations
    ! the fft grid dimension
    ! the counter for fft grid
  REAL(DP) :: rho(nrxx,nspin), rho_core(nrxx), &
                   rho_up, rho_dw, zeta, rhox
    ! the scf charge density
    ! the core charge density
    ! auxiliary variable for rho up
    ! auxiliary variable for rho dw
    ! auxiliary variable for zeta
    ! auxiliary variable for total rho
  !
  !
  IF ( nspin == 1 ) RETURN
  !
  IF ( iop == 1 ) THEN
     !
     DO ir = 1, nrxx
        !
        rhox = rho(ir,1) + rho(ir,2) + rho_core(ir)
        !
        IF ( rhox > eps32 ) THEN
           !
           zeta = ( rho(ir,1) - rho(ir,2) ) / rhox
           !
        ELSE
           !
           zeta = 0.D0
           !
        END IF
        !
        rho(ir,1) = rho(ir,1) + rho(ir,2)
        rho(ir,2) = zeta
        !
     END DO
     !
  ELSE IF ( iop == - 1 ) THEN
     !
     DO ir = 1, nrxx
        !
        rhox = rho(ir,1) + rho_core(ir)
        !
        rho_up = 0.5D0 * ( rho(ir,1) + rho(ir,2) * rhox )
        rho_dw = 0.5D0 * ( rho(ir,1) - rho(ir,2) * rhox )
        !
        rho(ir,1) = rho_up
        rho(ir,2) = rho_dw
        !
     END DO
     !
  ELSE
     !
     WRITE( stdout , '(5X,"iop = ",I5)' ) iop
     !
     CALL errore( 'rho2zeta', 'wrong iop', 1 )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE rho2zeta
!
!----------------------------------------------------------------------------
SUBROUTINE rho2mag( rho, rho_core, nrxx, nspin, iop )
  !---------------------------------------------------------------------------
  ! analogous to rho2zeta, for noncolinear magnetization
  !
  USE constants, ONLY : eps32
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in)    :: iop, nspin, nrxx
  REAL(DP), intent(in)   :: rho_core(nrxx)
  REAL(DP), intent(inout):: rho(nrxx,nspin)
  !
  INTEGER :: is
  !
  IF ( iop == -1 ) THEN
     !
     DO is = 2, nspin
        !
        WHERE( rho(:,1)+rho_core(:) > eps32 )
           !
           rho(:,is) = rho(:,is)*(rho(:,1)+rho_core(:))
           !
        ELSEWHERE
           !
           rho(:,is) = 0.D0
              !
        END WHERE
        !
     END DO
     !
     ELSE IF ( iop == 1 ) THEN
        !
     DO is = 2, nspin
        !
        WHERE( rho(:,1)+rho_core(:) > eps32 )
           !
           rho(:,is) = rho(:,is) / ( rho(:,1)+rho_core(:) )
           !
        ELSEWHERE
           !
           rho(:,is) = 0.D0
           !
        END WHERE
        !
     END DO
     !
  ELSE
     !
     WRITE( stdout , '(5X,"iop = ",I5)' ) iop
     !
     CALL errore( 'rho2mag', 'wrong iop', 1 )
     !
  END IF
  !
END SUBROUTINE rho2mag
