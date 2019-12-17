!
! Copyright (C) 2005-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
   SUBROUTINE compute_rho( rho, rhoout, segni, nrxx )
   !--------------------------------------------------------------------------
   !! This subroutine diagonalizes the spin density matrix and gives as output
   !! the spin-up and spin-down components of the charge.  
   !! If \(\text{lsign}\) is .TRUE. up and down are with respect to the fixed 
   !! quantization axis ux, otherwise \(\text{rho}+|m|\) is always rho-up and
   !! \(\text{rho}-|m|\) is always rho-down.
   !
   USE kinds,            ONLY : dp
   USE noncollin_module, ONLY : lsign, ux
   !
   IMPLICIT NONE
   !
   INTEGER :: nrxx
   !! input: the dimension of the mesh
   REAL(DP), INTENT(IN) :: rho(nrxx,4)
   !! the four components of the charge 
   REAL(DP), INTENT(OUT) :: rhoout(nrxx,2)
   !! the spin up and spin down charge
   REAL(DP), INTENT(OUT) :: segni(nrxx)
   !! the orientation when needed
   !
   ! ... local variables
   !
   REAL(DP) :: amag
   INTEGER :: ir           ! counter on mesh points
   IF (lsign) THEN 
      !$omp parallel do default(shared) private(ir, amag) 
      DO ir = 1, nrxx 
         segni(ir)=SIGN(1.0_DP,rho(ir,2)*ux(1)+rho(ir,3)*ux(2)+rho(ir,4)*ux(3))
         amag=SQRT(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
         rhoout(ir,1)=0.5d0*(rho(ir,1)+segni(ir)*amag)
         rhoout(ir,2)=0.5d0*(rho(ir,1)-segni(ir)*amag)
      ENDDO
      !$omp end parallel do
    ELSE 
      !$omp parallel do default(shared) private(ir, amag)
      DO ir =1, nrxx 
         segni(ir) = 1.0_DP
         amag=SQRT(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
         rhoout(ir,1)=0.5d0*(rho(ir,1) + amag)
         rhoout(ir,2)=0.5d0*(rho(ir,1) - amag)
      END DO
      !$omp end parallel do 
    END IF
   RETURN
END SUBROUTINE compute_rho
