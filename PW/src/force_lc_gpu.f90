!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE force_lc_gpu( nat, tau, ityp, alat, omega, ngm, ngl, &
                         igtongl_d, g_d, rho, nl_d, gstart, gamma_only, vloc_d, forcelc )
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
  USE device_fbuff_m,        ONLY : dev_buf
  USE device_memcpy_m,   ONLY : dev_memcpy
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
  INTEGER, INTENT(IN) :: igtongl_d(ngm)
  !! correspondence G <-> shell of G
  INTEGER, INTENT(IN) :: nl_d(ngm)
#if defined(__CUDA)
  attributes(DEVICE) :: igtongl_d, nl_d
#endif
  !! correspondence fft mesh <-> G vec
  INTEGER, INTENT(IN) :: ityp(nat)
  !! types of atoms
  LOGICAL, INTENT(IN) :: gamma_only
  !! gamma only
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! coordinates of the atoms
  REAL(DP), INTENT(IN) :: g_d(3,ngm)
  !! coordinates of G vectors
  REAL(DP), INTENT(IN) :: vloc_d(ngl,*)
#if defined(__CUDA)
  attributes(DEVICE) :: g_d, vloc_d
#endif
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
  INTEGER :: ipol, ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  COMPLEX(DP), ALLOCATABLE :: aux(:)
  COMPLEX(DP), POINTER     :: aux_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: aux_d
#endif
  ! auxiliary space for FFT
  REAL(DP) :: arg, fact
  !
  ! auxiliary for cuf kernels and error checking
  INTEGER :: ityp_na, ierr
  REAL(DP):: forcelc_x, forcelc_y, forcelc_z, tau1, tau2, tau3
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  ALLOCATE( aux(dfftp%nnr) )
  CALL dev_buf%lock_buffer(aux_d, dfftp%nnr, ierr)
  IF (ierr /= 0) CALL errore( 'force_lc_gpu', 'cannot allocate buffers', -1 )
  !
  aux(:) = CMPLX( rho(:), 0.0_DP, KIND=DP )
  CALL dev_memcpy( aux_d, aux )
  !
  CALL fwfft( 'Rho', aux_d, dfftp )
  IF ( ( do_comp_esm .AND. (esm_bc .NE. 'pbc') ) .or. do_cutoff_2D ) &
     CALL dev_memcpy( aux, aux_d )
  !
  ! aux contains now  n(G)
  !
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  !
  DO na = 1, nat
     ! This is no longher needed
     !do ipol = 1, 3
     !   forcelc (ipol, na) = 0.d0
     !enddo
     !
     ityp_na = ityp (na)
     tau1 = tau (1, na); tau2 = tau (2, na); tau3 = tau (3, na);
     forcelc_x = 0.d0 ; forcelc_y = 0.d0 ; forcelc_z = 0.d0
     ! contribution from G=0 is zero
     !$cuf kernel do(1)
     DO ig = gstart, ngm
        arg = (g_d (1, ig) * tau1 + g_d (2, ig) * tau2 + &
               g_d (3, ig) * tau3 ) * tpi
        !
        ! arg is used as auxiliary variable here
        arg = vloc_d (igtongl_d (ig), ityp_na ) * &
                (sin(arg)*DBLE(aux_d(nl_d(ig))) + cos(arg)*AIMAG(aux_d(nl_d(ig))) )
        !
        forcelc_x = forcelc_x + g_d (1, ig) * arg
        forcelc_y = forcelc_y + g_d (2, ig) * arg
        forcelc_z = forcelc_z + g_d (3, ig) * arg
        !
     ENDDO
     forcelc(1, na) = forcelc_x; forcelc(2, na) = forcelc_y; forcelc(3, na) = forcelc_z;
     forcelc(1:3,na) = fact * forcelc(1:3,na) * omega * tpi / alat
  ENDDO
  IF ( do_comp_esm .AND. (esm_bc .NE. 'pbc') ) THEN
     !
     ! ... Perform corrections for ESM method (add long-range part)
     !
     CALL esm_force_lc ( aux, forcelc )
  ENDIF
  !
  ! IN 2D calculations: re-add the erf/r contribution to the forces. It was substracted from
  ! vloc (in vloc_of_g) and readded to vltot only (in setlocal)
  IF ( do_cutoff_2D ) CALL cutoff_force_lc( aux, forcelc )
  !
  CALL mp_sum( forcelc, intra_bgrp_comm )
  !
  DEALLOCATE( aux )
  CALL dev_buf%release_buffer(aux_d, ierr)
  !
  RETURN
  !
END SUBROUTINE force_lc_gpu
