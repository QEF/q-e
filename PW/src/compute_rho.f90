!
! Copyright (C) 2005-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
   SUBROUTINE compute_rho(rho,rhoout,segni,nrxx)
!
!   This subroutine diagonalizes the spin density matrix and gives as output
!   the spin up and spin down compotents of the charge
!   If lsign is true up and dw are with respect to the fixed quantization
!   axis ux, otherwise rho+|m| is always rhoup and rho-|m| is always
!   rhodw
!
   USE kinds, ONLY : dp
   USE noncollin_module, ONLY : lsign, ux
   IMPLICIT NONE
   INTEGER :: nrxx  ! input: the dimension of the mesh
         
   REAL(DP), INTENT(IN) ::  rho(nrxx,4)
   REAL(DP), INTENT(OUT) :: rhoout(nrxx,2)
   REAL(DP), INTENT(OUT) :: segni(nrxx)
                 ! input: the four components of the charge 
                 ! output: the spin up and spin down charge
                 ! output: the orientation when needed
   REAL(DP) :: amag
   INTEGER :: ir           ! counter on mesh points

   segni=1.0_DP
   IF (lsign) THEN
      DO ir=1,nrxx
         segni(ir)=SIGN(1.0_DP,rho(ir,2)*ux(1)+rho(ir,3)*ux(2)+rho(ir,4)*ux(3))
      ENDDO
   ENDIF

   DO ir=1,nrxx
      amag=SQRT(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
      rhoout(ir,1)=0.5d0*(rho(ir,1)+segni(ir)*amag)
      rhoout(ir,2)=0.5d0*(rho(ir,1)-segni(ir)*amag)
   ENDDO

   RETURN
END SUBROUTINE compute_rho
