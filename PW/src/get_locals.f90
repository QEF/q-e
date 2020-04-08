!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE get_locals( rholoc, magloc, rho )
  !----------------------------------------------------------------------------
      !! Here local integrations are carried out around atoms.
      !! The points and weights for these integrations are determined in the
      !! subroutine make_pointlists, the result may be printed in the
      !! subroutine report_mag. If constraints are present, the results of this
      !! calculation are used in v_of_rho for determining the penalty functional.
      !
      USE kinds,             ONLY: DP
      USE ions_base,         ONLY: nat
      USE cell_base,         ONLY: omega
      USE lsda_mod,          ONLY: nspin
      USE mp_bands,          ONLY: intra_bgrp_comm
      USE mp,                ONLY: mp_sum
      USE fft_base,          ONLY: dfftp
      USE noncollin_module,  ONLY: pointlist, factlist, noncolin
      !
      IMPLICIT NONE
      !
      REAL(DP) :: rholoc(nat)
      !! integrated charge arount the atoms
      REAL(DP) :: magloc(nspin-1,nat) 
      !! integrated magnetic moment around the atom
      REAL(DP) :: rho(dfftp%nnr,nspin)
      !! charge density
      !
      ! ... local variables
      !
      INTEGER :: i,ipol
      REAL(DP) :: fact
      REAL(DP), ALLOCATABLE :: auxrholoc(:,:)
      !
      ALLOCATE( auxrholoc(0:nat,nspin) )
      auxrholoc(:,:) = 0.0_DP
      !
      DO i = 1, dfftp%nnr
         auxrholoc(pointlist(i),1:nspin) = auxrholoc(pointlist(i),1:nspin) + &
                                           rho(i,1:nspin) * factlist(i)
      ENDDO
      !
      CALL mp_sum( auxrholoc(0:nat,1:nspin), intra_bgrp_comm )
      !     
      fact = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
      !
      rholoc(1:nat) = auxrholoc(1:nat,1) * fact
      !
      DO ipol = 1, nspin-1
         magloc(ipol,1:nat) = auxrholoc(1:nat,ipol+1) * fact
      ENDDO
      !
      DEALLOCATE( auxrholoc )
      !
END SUBROUTINE get_locals
