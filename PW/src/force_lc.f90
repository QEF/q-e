!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE force_lc( nat, tau, ityp, alat, omega, ngm, ngl, &
                     igtongl, g, rho, nl, gstart, gamma_only, vloc, forcelc )
  !----------------------------------------------------------------------
  !! It calculates the local-potential contribution to forces on atoms.
  !
  USE kinds
  USE constants,       ONLY : tpi
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft
  USE esm,             ONLY : esm_force_lc, do_comp_esm, esm_bc
  USE Coul_cut_2D,     ONLY : do_cutoff_2D, cutoff_force_lc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nat
  !! number of atoms in the cell
  INTEGER, INTENT(IN) :: ngm
  !! number of G vectors
  INTEGER, INTENT(IN) :: ngl
  !! number of shells
  INTEGER, INTENT(IN) :: gstart
  !! index of the first G vector whose module is > 0 (see
  !! module 'gvect' in Modules/recvec.f90)
  INTEGER, INTENT(IN) :: igtongl(ngm)
  !! correspondence G <-> shell of G
  INTEGER, INTENT(IN) :: nl(ngm)
  !! correspondence fft mesh <-> G vec
  INTEGER, INTENT(IN) :: ityp(nat)
  !! types of atoms
  LOGICAL, INTENT(IN) :: gamma_only
  !! gamma only
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! coordinates of the atoms
  REAL(DP), INTENT(IN) :: g(3,ngm)
  !! coordinates of G vectors
  REAL(DP), INTENT(IN) :: vloc(ngl,*)
  !! local potential
  REAL(DP), INTENT(IN) :: rho(dfftp%nnr)
  !! valence charge
  REAL(DP), INTENT(IN) :: alat
  !! lattice parameter
  REAL(DP), INTENT(IN) :: omega
  !! unit cell volume
  REAL(DP), INTENT(OUT) :: forcelc(3, nat)
  !! the local-potential contribution to forces on atoms
  !
  ! ... local variables
  !
  INTEGER :: ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  ! auxiliary space for FFT
  REAL(DP) :: arg, fact
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  ALLOCATE( aux(dfftp%nnr) )
  !
  aux(:) = CMPLX( rho(:), 0.0_DP, KIND=DP )
  !
  CALL fwfft( 'Rho', aux, dfftp )
  !
  ! aux contains now  n(G)
  !
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  !
!$omp parallel do private(arg)
  DO na = 1, nat
     !
     forcelc (1:3, na) = 0.d0
     ! contribution from G=0 is zero
     DO ig = gstart, ngm
        arg = (g(1,ig) * tau(1,na) + g(2,ig) * tau(2,na) + &
               g(3,ig) * tau(3,na) ) * tpi
        forcelc(1:3,na) = forcelc(1:3, na) + &
                          g(1:3,ig) * vloc(igtongl(ig),ityp(na) ) * &
                          (SIN(arg)*DBLE(aux(nl(ig))) + COS(arg)*AIMAG(aux(nl(ig))) )
     ENDDO
     !
     forcelc(1:3,na) = fact * forcelc(1:3,na) * omega * tpi / alat
  ENDDO
!$omp end parallel do
  !
  IF ( do_comp_esm .AND. (esm_bc .NE. 'pbc') ) THEN
     !
     ! ... Perform corrections for ESM method (add long-range part)
     CALL esm_force_lc( aux, forcelc )
     !
  ENDIF
  !
  ! IN 2D calculations: re-add the erf/r contribution to the forces. It was substracted from
  ! vloc (in vloc_of_g) and readded to vltot only (in setlocal)
  IF ( do_cutoff_2D ) CALL cutoff_force_lc( aux, forcelc )
  !
  CALL mp_sum( forcelc, intra_bgrp_comm )
  !
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE force_lc
