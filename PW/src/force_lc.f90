!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE force_lc( nat, tau, ityp, ntyp, alat, omega, ngm, ngl, &
                     igtongl, g, rho, gstart, gamma_only, vloc, forcelc )
  !----------------------------------------------------------------------
  !! It calculates the local potential contribution to forces on atoms.
  !
  USE kinds
  USE constants,       ONLY : tpi
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  USE fft_base,        ONLY : dfftp
  USE fft_rho,         ONLY : rho_r2g
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
  INTEGER, INTENT(IN) :: ityp(nat)
  !! types of atoms
  INTEGER, INTENT(IN) :: ntyp
  !! number of types of atoms
  LOGICAL, INTENT(IN) :: gamma_only
  !! gamma only
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! coordinates of the atoms
  REAL(DP), INTENT(IN) :: g(3,ngm)
  !! coordinates of G vectors
  REAL(DP), INTENT(IN) :: vloc(ngl,ntyp)
  !! local potential
  REAL(DP), INTENT(IN) :: rho(dfftp%nnr)
  !! valence charge
  REAL(DP), INTENT(IN) :: alat
  !! lattice parameter
  REAL(DP), INTENT(IN) :: omega
  !! unit cell volume
  REAL(DP), INTENT(OUT) :: forcelc(3,nat)
  !! the local potential contribution to forces on atoms
  !
  ! ... local variables
  !
  INTEGER :: ig, na, ityp_na
  ! counter on G vectors
  ! counter on atoms
  ! atom type index
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  ! auxiliary space for FFT
  REAL(DP) :: arg, arg0, fact
  REAL(DP) :: forcelc_x, forcelc_y, forcelc_z, tau1, tau2, tau3
  !
  ! ... contribution to the force from the local part of the bare potential
  !     F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  !$acc data present_or_copyin( igtongl, g, rho, vloc )
  !
  ALLOCATE( aux(dfftp%nnr,1) )
  !$acc data create( aux )
  !
  CALL rho_r2g( dfftp, rho, aux )
  !
  IF ( ( do_comp_esm .AND. (esm_bc .NE. 'pbc') ) .OR. do_cutoff_2D ) THEN
    !$acc update self(aux)
  ENDIF
  !
  ! ... aux contains now n(G)
  !
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  !
#if !defined(_OPENACC)
!$omp parallel do private( ityp_na,tau1,tau2,tau3,forcelc_x,forcelc_y,forcelc_z,&
!$omp                      ig,arg,arg0 )
#endif
  DO na = 1, nat
     !
     ityp_na = ityp(na)
     tau1 = tau(1,na) ; tau2 = tau(2,na) ; tau3 = tau(3,na)
     forcelc_x = 0.d0 ; forcelc_y = 0.d0 ; forcelc_z = 0.d0
     !
     ! ... contribution from G=0 is zero
     !
     !$acc parallel loop reduction(+:forcelc_x,forcelc_y,forcelc_z)
     DO ig = gstart, ngm
        arg0 = (g(1,ig)*tau1 + g(2,ig)*tau2 + g(3,ig)*tau3) * tpi
        !
        arg = vloc(igtongl(ig),ityp_na) * &
              ( SIN(arg0)*DBLE(aux(ig,1)) + COS(arg0)*AIMAG(aux(ig,1)) )
        !
        forcelc_x = forcelc_x + g(1,ig) * arg
        forcelc_y = forcelc_y + g(2,ig) * arg
        forcelc_z = forcelc_z + g(3,ig) * arg
     ENDDO
     !
     forcelc(1,na) = fact * forcelc_x * omega * tpi / alat
     forcelc(2,na) = fact * forcelc_y * omega * tpi / alat
     forcelc(3,na) = fact * forcelc_z * omega * tpi / alat
  ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
  !
  IF ( do_comp_esm .AND. (esm_bc .NE. 'pbc') ) THEN
     !
     ! ... Perform corrections for ESM method (add long-range part)
     CALL esm_force_lc( aux(:,1), forcelc )
     !
  ENDIF
  !
  ! 2D calculations: re-add the erf/r contribution to the forces that was 
  ! substracted from vloc (in vloc_of_g) and re-added to vltot only (in setlocal)
  IF ( do_cutoff_2D ) CALL cutoff_force_lc( gamma_only, aux(:,1), forcelc )
  !
  CALL mp_sum( forcelc, intra_bgrp_comm )
  !
  !$acc end data
  DEALLOCATE( aux )
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE force_lc
