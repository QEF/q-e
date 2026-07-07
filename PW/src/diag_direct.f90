!
! Copyright (C) 2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE diag_direct
  !----------------------------------------------------------------------------
  !!
  !! Direct diagonalization of the dense Hamiltonian matrix H(G,G').
  !! First build the matrix in the plane-wave basis, then diagonalize it using
  !! LAXlib routines.
  !!
  !! H(G,G') = δ_GG' × (k+G)² + V_loc(G-G') + Σ_ij |β_i(G)⟩ D_ij ⟨β_j(G')|
  !!
  !! RESTRICTIONS: NCPP only, no special features (see diag_direct_check_compat)
  !!
  !! This module is inspired by the ParaBands code, part of the BerkeleyGW package,
  !! distributed under a 3-Clause BSD license:
  !! BerkeleyGW, Copyright (c) 2011, The Regents of the University of California.
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: diag_direct_check_compat
  PUBLIC :: diag_direct_run_k
  PUBLIC :: diag_direct_test_hamiltonian
  !
  INTEGER :: ngk_g
  !! global total number of k+G vectors at this k point
  INTEGER, ALLOCATABLE :: igk_l2g_kdip(:)
  !! Local to global mappings
  INTEGER, ALLOCATABLE :: mill_k_global(:, :)
  !! Miller indices in the global k+G ordering
  !
  ! Derived type for Miller index lookup table
  !
  TYPE mill_lookup_type
     INTEGER, ALLOCATABLE :: map(:,:,:)
     !! 3D lookup array: map(nx,ny,nz) = G-vector index or 0 if not present
     INTEGER :: imin, imax
     !! Bounds in first dimension
     INTEGER :: jmin, jmax
     !! Bounds in second dimension
     INTEGER :: kmin, kmax
     !! Bounds in third dimension
  END TYPE mill_lookup_type
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  INTEGER FUNCTION lookup_mill_index( lookup, mill_vec ) RESULT(ig)
    !-----------------------------------------------------------------------
    !! Safe lookup of G-vector index from Miller indices
    !! Returns 0 if Miller indices are out of bounds or G-vector not present
    !
    IMPLICIT NONE
    !
    TYPE(mill_lookup_type), INTENT(IN) :: lookup
    !! Lookup table structure
    INTEGER, INTENT(IN) :: mill_vec(3)
    !! Miller indices to look up
    !
    ! Check bounds
    IF ( mill_vec(1) < lookup%imin .OR. mill_vec(1) > lookup%imax .OR. &
         mill_vec(2) < lookup%jmin .OR. mill_vec(2) > lookup%jmax .OR. &
         mill_vec(3) < lookup%kmin .OR. mill_vec(3) > lookup%kmax ) THEN
       ig = 0
    ELSE
       ! Lookup in table
       ig = lookup%map(mill_vec(1), mill_vec(2), mill_vec(3))
    ENDIF
    !
  END FUNCTION lookup_mill_index
  !
  !-----------------------------------------------------------------------
  SUBROUTINE build_mill_lookup_table(lookup)
    !-----------------------------------------------------------------------
    !! Build 3D lookup table: Miller indices (Gx,Gy,Gz) → global G-vector index
    !!
    !! Uses GLOBAL Miller indices to ensure all G-vectors are indexed,
    !! including those from other processors in parallel execution.
    !!
    !! Returns global G-vector indices (ig_l2g) for use across processors.
    !
    USE mp,            ONLY : mp_max, mp_min, mp_sum
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE gvect,         ONLY : mill, ngm, ig_l2g
    !
    IMPLICIT NONE
    !
    TYPE(mill_lookup_type), INTENT(OUT) :: lookup
    !! Lookup table structure (output)
    !
    INTEGER :: ig
    INTEGER :: mill_max(3), mill_min(3)
    !
    ! Determine bounds from actual Miller indices using intrinsic functions
    !
    mill_max(1) = MAXVAL(mill(1, 1:ngm))
    mill_max(2) = MAXVAL(mill(2, 1:ngm))
    mill_max(3) = MAXVAL(mill(3, 1:ngm))
    !
    mill_min(1) = MINVAL(mill(1, 1:ngm))
    mill_min(2) = MINVAL(mill(2, 1:ngm))
    mill_min(3) = MINVAL(mill(3, 1:ngm))
    !
    CALL mp_max(mill_max, intra_bgrp_comm)
    CALL mp_min(mill_min, intra_bgrp_comm)
    !
    ! Set bounds (no padding needed)
    lookup%imin = mill_min(1)
    lookup%imax = mill_max(1)
    lookup%jmin = mill_min(2)
    lookup%jmax = mill_max(2)
    lookup%kmin = mill_min(3)
    lookup%kmax = mill_max(3)
    !
    ! Allocate 3D array
    ALLOCATE(lookup%map(lookup%imin:lookup%imax, &
                        lookup%jmin:lookup%jmax, &
                        lookup%kmin:lookup%kmax))
    lookup%map(:,:,:) = 0
    !
    ! Fill lookup table for the global G-vector indices.
    ! Each processor fills its own G-vectors, and then we sum across all.
    !
    DO ig = 1, ngm
       lookup%map(mill(1,ig), mill(2,ig), mill(3,ig)) = ig_l2g(ig)
    ENDDO
    !
    CALL mp_sum(lookup%map, intra_bgrp_comm)
    !
  END SUBROUTINE build_mill_lookup_table
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_check_compat()
    !-----------------------------------------------------------------------
    !! Check if calculation settings are compatible with full-band diagonalization
    !!
    !! Calls errore and stops if incompatible
    !
    USE kinds,            ONLY : DP
    USE control_flags,    ONLY : gamma_only
    USE realus,           ONLY : real_space
    USE noncollin_module, ONLY : noncolin
    USE uspp,             ONLY : okvan
    USE ldaU,             ONLY : lda_plus_u
    USE bp,               ONLY : lelfield
    USE xc_lib,           ONLY : exx_is_active, xclib_dft_is
    USE fft_base,         ONLY : dffts
    !
    IMPLICIT NONE
    !
    ! Check for incompatible features
    !
    IF ( okvan ) CALL errore('diag_direct_check_compat', &
       'Full-band diagonalization requires NCPP (no USPP/PAW)', 1)
    !
    IF ( gamma_only ) CALL errore('diag_direct_check_compat', &
       'gamma_only not supported in full-band diagonalization', 1)
    !
    IF ( real_space ) CALL errore('diag_direct_check_compat', &
       'real_space algorithms not supported in full-band diagonalization', 1)
    !
    IF ( noncolin ) CALL errore('diag_direct_check_compat', &
       'non-collinear calculation not supported in full-band diagonalization', 1)
    !
    IF ( xclib_dft_is('meta') ) CALL errore('diag_direct_check_compat', &
       'meta-GGA functionals not supported in full-band diagonalization', 1)
    !
    IF ( lda_plus_u ) CALL errore('diag_direct_check_compat', &
       'DFT+U not supported in full-band diagonalization', 1)
    !
    IF ( exx_is_active() ) CALL errore('diag_direct_check_compat', &
       'Exact exchange/hybrid functionals not supported in full-band diagonalization', 1)
    !
    IF ( lelfield ) CALL errore('diag_direct_check_compat', &
       'Electric field (Berry phase) not supported in full-band diagonalization', 1)
    !
    IF ( dffts%has_task_groups ) CALL errore('diag_direct_check_compat', &
       'Task groups not supported in full-band diagonalization', 1)
    !
  END SUBROUTINE diag_direct_check_compat
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_hamiltonian_kinetic(hmat, npw)
    !-----------------------------------------------------------------------
    !! Build kinetic energy contribution (diagonal) for local columns
    !! Matrix is hmat(ngk_g, npw) - full rows, local columns
    !
    USE kinds,  ONLY : DP
    USE wvfct,  ONLY : g2kin
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: hmat(ngk_g, npw)
    !! Hamiltonian matrix (full rows, local columns)
    INTEGER, INTENT(IN) :: npw
    !! number of local plane waves
    !
    INTEGER :: ig, ig_global
    !
    CALL start_clock('direct_kin')
    !
    ! Build kinetic energy diagonal for local columns
    ! For column ig (local), diagonal element is at row ig_global
    DO ig = 1, npw
       ig_global = igk_l2g_kdip(ig)
       hmat(ig_global, ig) = CMPLX(g2kin(ig), 0.0_DP, KIND=DP)
    ENDDO
    !
    CALL stop_clock('direct_kin')
    !
  END SUBROUTINE diag_direct_hamiltonian_kinetic
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_hamiltonian_vloc(hmat, npw)
    !-----------------------------------------------------------------------
    !! Build local potential contribution V_loc(G-G') for local columns
    !! Matrix is hmat(ngk_g, npw) - full rows, local columns
    !
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, ngm_g, ig_l2g
    USE lsda_mod,      ONLY : current_spin
    USE scf,           ONLY : vrs
    USE fft_base,      ONLY : dffts
    USE fft_interfaces,ONLY : fwfft
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: hmat(ngk_g, npw)
    !! Hamiltonian matrix (full rows, local columns)
    INTEGER, INTENT(IN) :: npw
    !! number of local plane waves
    !
    ! Local variables
    COMPLEX(DP), ALLOCATABLE :: vrs_g(:)
    !! Local potential in G-space
    COMPLEX(DP), ALLOCATABLE :: vrs_aux(:)
    !! FFT workspace for V_loc
    TYPE(mill_lookup_type) :: mill_lookup
    !! Lookup table structure for Miller indices → G-vector index
    INTEGER :: ig, jg_local, jg_global, ig_global, igg
    INTEGER :: mill_diff(3)
    !
    CALL start_clock('direct_vloc')
    !
    ALLOCATE(vrs_g(ngm_g))
    ALLOCATE(vrs_aux(dffts%nnr))
    vrs_g(:) = (0.0_DP, 0.0_DP)
    !
    ! Transform V_loc to G-space
    vrs_aux(:) = CMPLX(vrs(:,current_spin), 0.0_DP, KIND=DP)
    CALL fwfft('Rho', vrs_aux, dffts)
    !
    ! Extract G-space components using global G-vector indices
    DO ig = 1, ngm
       vrs_g(ig_l2g(ig)) = vrs_aux(dffts%nl(ig))
    ENDDO
    CALL mp_sum(vrs_g, intra_bgrp_comm)
    !
    ! Build Miller index lookup table
    CALL build_mill_lookup_table(mill_lookup)
    !
    ! Build V_loc contributions for local columns (full rows)
    ! For each local column jg_local, build full row using mill_k_global
    DO jg_local = 1, npw
       jg_global = igk_l2g_kdip(jg_local)
       !
       ! Build full row for this column
       DO ig_global = 1, ngk_g
          ! Compute Miller indices difference G - G' using global arrays
          mill_diff(1:3) = mill_k_global(1:3, ig_global) - mill_k_global(1:3, jg_global)
          !
          ! Lookup the index of G - G' in the global list of G vectors
          igg = lookup_mill_index(mill_lookup, mill_diff)
          !
          ! If found, add V_loc(G-G') to H(G,G')
          IF (igg > 0) THEN
             hmat(ig_global, jg_local) = hmat(ig_global, jg_local) + vrs_g(igg)
          ENDIF
       ENDDO
    ENDDO
    !
    DEALLOCATE(vrs_g)
    DEALLOCATE(vrs_aux)
    DEALLOCATE(mill_lookup%map)
    !
    CALL stop_clock('direct_vloc')
    !
  END SUBROUTINE diag_direct_hamiltonian_vloc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_hamiltonian_vnl(hmat, npw)
    !-----------------------------------------------------------------------
    !! Build nonlocal pseudopotential contribution
    !! H(G,G') += Σ_ij |β_i(G)⟩ D_ij ⟨β_j(G')|
    !
    USE kinds,      ONLY : DP
    USE mp,         ONLY : mp_sum
    USE mp_bands,   ONLY : intra_bgrp_comm
    USE lsda_mod,   ONLY : current_spin
    USE wvfct,      ONLY : npwx
    USE uspp,       ONLY : nkb, vkb, deeq, ofsbeta
    USE uspp_param, ONLY : nh
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: hmat(ngk_g, npw)
    !! Hamiltonian matrix (full rows, local columns)
    INTEGER, INTENT(IN) :: npw
    !! number of local plane waves
    !
    ! Local variables
    INTEGER :: nt, na, ib, jb, ig
    !! Counters
    COMPLEX(DP), ALLOCATABLE :: vkb_global(:, :)
    !! Global vkb array collected from all processors in intra_bgrp_comm
    COMPLEX(DP), ALLOCATABLE :: dmat(:, :)
    !! Pseudopotential D matrix for all atoms
    COMPLEX(DP), ALLOCATABLE :: tmp(:, :)
    !! Temporary matrix for ZGEMM
    !
    IF (nkb == 0) RETURN
    !
    CALL start_clock('direct_vnl')
    !
    ALLOCATE(vkb_global(ngk_g, nkb))
    ALLOCATE(dmat(nkb, nkb))
    ALLOCATE(tmp(nkb, npw))
    !
    ! Collect vkb globally
    !
    vkb_global = (0.0_DP, 0.0_DP)
    DO ig = 1, npw
       vkb_global(igk_l2g_kdip(ig), 1:nkb) = vkb(ig, 1:nkb)
    ENDDO
    CALL mp_sum(vkb_global, intra_bgrp_comm)
    !
    ! Construct D matrix
    !
    dmat(:, :) = (0.0_DP, 0.0_DP)
    !
    ! Loop over atom types and atoms
    DO nt = 1, ntyp
       IF ( nh(nt) == 0 ) CYCLE
       !
       DO na = 1, nat
          IF ( ityp(na) /= nt ) CYCLE
          !
          ! Copy real deeq into complex dmat
          !
          DO ib = 1, nh(nt)
             DO jb = 1, nh(nt)
                dmat(ofsbeta(na)+ib, ofsbeta(na)+jb) = CMPLX(deeq(ib, jb, na, current_spin), 0.0_DP, KIND=DP)
             ENDDO
          ENDDO
          !
       ENDDO  ! na
       !
    ENDDO  ! nt
    !
    ! H(ig_global, jg) += Σ_{b, b'} vkb_global(ig_global, b) * D(b, b') * vkb^\dagger(jg, b')
    !
    CALL ZGEMM('N', 'C', nkb, npw, nkb, (1.0_DP, 0.0_DP), &
               dmat, nkb, vkb, npwx, (0.0_DP, 0.0_DP), tmp, nkb)
    !
    CALL ZGEMM('N', 'N', ngk_g, npw, nkb, (1.0_DP, 0.0_DP), &
               vkb_global, ngk_g, tmp, nkb, (1.0_DP, 0.0_DP), hmat, ngk_g)
    !
    DEALLOCATE(vkb_global)
    DEALLOCATE(dmat)
    DEALLOCATE(tmp)
    !
    CALL stop_clock('direct_vnl')
    !
  END SUBROUTINE diag_direct_hamiltonian_vnl
  !-----------------------------------------------------------------------
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_run_k(ik, npw, nbnd, evc, et_k, notconv)
    !-----------------------------------------------------------------------
    !! Construct and diagonalize explicit Hamiltonian matrix for k-point ik
    !!
    !! This is the main coordinator that:
    !! 1. Sets up the local-to-global k+G mapping
    !! 2. Builds H column-block by column-block via hamiltonian_kinetic,
    !!    hamiltonian_vloc, hamiltonian_vnl
    !! 3. Gathers the column blocks into a global H matrix
    !! 4. Diagonalizes via direct_diag_k
    !! 5. Scatters the global eigenvectors back to the local PW layout
    !!
    !! Optional: diag_direct_test_hamiltonian verifies H*v against h_psi(v)
    !! (call is commented out by default; uncomment for debugging).
    !
    USE kinds,         ONLY : DP
    USE wvfct,         ONLY : npwx
    USE io_global,     ONLY : stdout
    USE mp,            ONLY : mp_sum, mp_max
    USE control_flags, ONLY : use_para_diag
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    INTEGER, INTENT(IN) :: nbnd
    !! number of bands
    COMPLEX(DP), INTENT(OUT) :: evc(npwx,nbnd)
    !! eigenvectors (output)
    REAL(DP), INTENT(OUT) :: et_k(nbnd)
    !! eigenvalues (output)
    INTEGER, INTENT(OUT) :: notconv
    !! number of non-converged bands (always 0 for direct method)
    !
    ! ... local variables
    !
    INTEGER :: ierr
    !! allocation error code
    INTEGER :: ibnd, ig
    !! loop indices for the global -> local eigenvector scatter
    REAL(DP) :: mem_gb
    !! memory estimate
    COMPLEX(DP), ALLOCATABLE :: hmat(:,:)
    !! Hamiltonian matrix (ngk_g, npw) - full rows, local columns
    COMPLEX(DP), ALLOCATABLE :: hmat_global(:,:)
    !! Global Hamiltonian matrix (ngk_g, ngk_g)
    COMPLEX(DP), ALLOCATABLE :: evc_global(:,:)
    !! Global eigenvectors (ngk_g, nbnd) returned by direct_diag_k
    !
    CALL start_clock('diag_direct')
    !
    ! Setup global G-vector mapping for this k-point
    !
    CALL diag_direct_setup_gmap(npw, ik)
    !
    ! Report memory requirements (full rows, local columns)
    !
    mem_gb = REAL(ngk_g,DP)**2 * 16.0_DP / 1024.0_DP**3
    WRITE(stdout,'(/,5X,"Dense H construction for k-point",I5)') ik
    WRITE(stdout,'(5X,"npw (local) =",I8,", ngk_g (global) =",I8)') npw, ngk_g
    WRITE(stdout,'(5X,"Matrix size: ngk_g =",I8,", Memory:",F8.3," GB")') ngk_g, mem_gb
    !
    ! Allocate matrices with full rows, local columns
    ! Each processor builds H(ngk_g, npw) - full rows for local columns
    !
    ALLOCATE( hmat(ngk_g, npw), STAT=ierr )
    IF ( ierr /= 0 ) CALL errore('diag_direct_run_k', 'H matrix allocation failed', ierr)
    !
    ! Initialize matrices
    !
    hmat(:,:) = (0.0_DP, 0.0_DP)
    !
    ! Build Hamiltonian terms (full rows for local columns)
    !
    CALL diag_direct_hamiltonian_kinetic(hmat, npw)
    CALL diag_direct_hamiltonian_vloc(hmat, npw)
    CALL diag_direct_hamiltonian_vnl(hmat, npw)
    !
    ! Collect Hamiltonian columns from all processors to global matrix
    ! Each processor has hmat(ngk_g, npw) - full rows, local columns
    ! Need to gather to hmat_global(ngk_g, ngk_g) - full matrix
    !
    ALLOCATE(hmat_global(ngk_g, ngk_g), STAT=ierr)
    IF ( ierr /= 0 ) CALL errore('diag_direct_run_k', 'hmat_global allocation failed', ierr)
    !
    CALL diag_direct_collect_matrices(hmat, npw, hmat_global)
    !
    ! Verification against h_psi (compares H*v from the explicit matrix
    ! against h_psi(v) for a random vector). Disabled by default because
    ! it adds an O(ngk_g^2) matvec + an h_psi call at every k-point of
    ! every SCF iteration. Uncomment when developing or debugging the
    ! Hamiltonian construction.
    !
    ! CALL diag_direct_test_hamiltonian(hmat_global, npw)
    !
    ! Diagonalize global matrix using the KS_Solvers/Direct routine,
    ! then scatter the global eigenvectors to local PW layout (PW-side concern).
    !
    ALLOCATE( evc_global(ngk_g, nbnd), STAT=ierr )
    IF ( ierr /= 0 ) CALL errore('diag_direct_run_k', 'evc_global allocation failed', ierr)
    evc_global(:,:) = (0.0_DP, 0.0_DP)
    !
    CALL direct_diag_k( ngk_g, nbnd, hmat_global, evc_global, et_k, &
                        use_para_diag, notconv )
    !
    CALL start_clock('direct_distrib')
    evc(:,:) = (0.0_DP, 0.0_DP)
    DO ibnd = 1, nbnd
       DO ig = 1, npw
          evc(ig, ibnd) = evc_global(igk_l2g_kdip(ig), ibnd)
       END DO
    END DO
    CALL stop_clock('direct_distrib')
    !
    ! Deallocate
    !
    DEALLOCATE(hmat)
    DEALLOCATE(hmat_global)
    DEALLOCATE(evc_global)
    DEALLOCATE(igk_l2g_kdip)
    DEALLOCATE(mill_k_global)
    !
    CALL stop_clock('diag_direct')
    !
    ! Uncomment for profiling:
    !
    ! WRITE(stdout,'(/,5X,"Timing breakdown for dense H construction:")')
    ! CALL print_clock('direct_kin')
    ! CALL print_clock('direct_vloc')
    ! CALL print_clock('direct_vnl')
    ! CALL print_clock('direct_collect')
    ! CALL print_clock('direct_check')
    ! CALL print_clock('direct_diag')
    ! CALL print_clock('diag_direct')
    !
  END SUBROUTINE diag_direct_run_k
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_test_hamiltonian(h, npw)
    !-----------------------------------------------------------------------
    !! Verify constructed H matrix by comparing H*v with h_psi(v)
    !! for a random test vector v
    !
    USE kinds,            ONLY : DP
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE wvfct,            ONLY : npwx
    USE noncollin_module, ONLY : npol
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    COMPLEX(DP), INTENT(IN) :: h(ngk_g, ngk_g)
    !! Hamiltonian matrix to verify
    INTEGER, INTENT(IN) :: npw
    !! number of local plane waves
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: v(:), v_global(:)
    !! random test vector
    COMPLEX(DP), ALLOCATABLE :: hv_mat_global(:)
    !! H*v from matrix multiply
    COMPLEX(DP), ALLOCATABLE :: hv_hpsi(:), hv_hpsi_global(:)
    !! H*v from h_psi
    REAL(DP), ALLOCATABLE :: rand_vec(:)
    !! random numbers
    REAL(DP) :: vnorm, max_error, rms_error
    INTEGER :: ig
    !
    EXTERNAL :: h_psi
    ! subroutine h_psi(lda,n,m,psi,hpsi)
    !
    CALL start_clock('direct_check')
    !
    ! Allocate arrays
    !
    ALLOCATE(v(npw))
    ALLOCATE(hv_hpsi(npwx))
    ALLOCATE(rand_vec(npw))
    ALLOCATE(v_global(ngk_g))
    ALLOCATE(hv_mat_global(ngk_g))
    ALLOCATE(hv_hpsi_global(ngk_g))
    !
    ! Generate random normalized vector
    !
    CALL RANDOM_NUMBER(rand_vec)
    v(:) = CMPLX(rand_vec(:), 0.0_DP, KIND=DP)
    vnorm = SQRT(SUM(ABS(v)**2))
    v(:) = v(:) / vnorm
    !
    ! Collect to a global vector
    !
    v_global(:) = (0.0_DP, 0.0_DP)
    DO ig = 1, npw
       v_global(igk_l2g_kdip(ig)) = v(ig)
    ENDDO
    CALL mp_sum(v_global, intra_bgrp_comm)
    !
    ! Compute H*v using matrix multiply
    !
    hv_mat_global(:) = (0.0_DP, 0.0_DP)
    DO ig = 1, ngk_g
       hv_mat_global(ig) = SUM(h(ig,:) * v_global(:))
    ENDDO
    !
    ! Compute H*v using h_psi
    !
    hv_hpsi(:) = (0.0_DP, 0.0_DP)
    CALL h_psi( npwx, npw, 1, v, hv_hpsi )
    !
    ! Collect to a global vector
    !
    hv_hpsi_global(:) = (0.0_DP, 0.0_DP)
    DO ig = 1, npw
       hv_hpsi_global(igk_l2g_kdip(ig)) = hv_hpsi(ig)
    ENDDO
    CALL mp_sum(hv_hpsi_global, intra_bgrp_comm)
    !
    ! Compare results
    !
    max_error = 0.0_DP
    rms_error = 0.0_DP
    DO ig = 1, ngk_g
       max_error = MAX(max_error, ABS(hv_mat_global(ig) - hv_hpsi_global(ig)))
       rms_error = rms_error + ABS(hv_mat_global(ig) - hv_hpsi_global(ig))**2
    ENDDO
    rms_error = SQRT(rms_error / REAL(npw,DP))
    !
    ! Report results: stop the run if H*v disagrees with h_psi(v).
    ! Continuing past a failed verification means SCF would proceed with
    ! a wrong Hamiltonian, so this must be an errore, not a stdout note.
    !
    IF (max_error > 1.d-10) THEN
      WRITE(stdout,'(5X,"max |H_mat*v - h_psi(v)| =",E12.4)') max_error
      WRITE(stdout,'(5X,"RMS |H_mat*v - h_psi(v)| =",E12.4)') rms_error
      CALL errore('diag_direct_test_hamiltonian', &
                  'explicit H disagrees with h_psi', 1)
    ENDIF
    !
    DEALLOCATE( v, v_global, hv_mat_global, hv_hpsi_global, rand_vec, hv_hpsi )
    !
    CALL stop_clock('direct_check')
    !
  END SUBROUTINE diag_direct_test_hamiltonian
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_collect_matrices(hmat, npw, hmat_global)
    !-----------------------------------------------------------------------
    !! Collect local matrix columns to global matrix
    !!
    !! Each processor has hmat(ngk_g, npw) - full rows, local columns
    !! Gather to hmat_global(ngk_g, ngk_g) - full matrix on all processors
    !!
    !! NCPP-only: no S matrix is needed (S=I is handled inside direct_diag_k).
    !
    USE kinds,     ONLY : DP
    USE mp_bands,  ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN)  :: hmat(ngk_g, npw)
    !! Local Hamiltonian matrix (full rows, local columns)
    INTEGER, INTENT(IN) :: npw
    !! Number of local plane waves
    COMPLEX(DP), INTENT(OUT) :: hmat_global(ngk_g, ngk_g)
    !! Global Hamiltonian matrix
    !
    INTEGER :: ig, ig_global
    !
    CALL start_clock('direct_collect')
    !
    ! Initialize global matrix to zero
    hmat_global(:,:) = (0.0_DP, 0.0_DP)
    !
    ! Collect H matrix: each processor contributes its local columns
    DO ig = 1, npw
       ig_global = igk_l2g_kdip(ig)
       hmat_global(:, ig_global) = hmat(:, ig)
    END DO
    !
    ! Sum H matrix across all processors (non-overlapping contributions)
    CALL mp_sum(hmat_global, intra_bgrp_comm)
    !
    ! Note: No S matrix needed for NCPP (S=I handled by direct_diag_k)
    !
    CALL stop_clock('direct_collect')
    !
  END SUBROUTINE diag_direct_collect_matrices
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_direct_setup_gmap(npw, ik)
     !-----------------------------------------------------------------------
     !! Setup the mapping of the local to global list of k+G vectors.
     !! ngk_g is the size of the Hamiltonian at the given k point.
     !! igk_l2g_kdip maps local k+G indices (1 to npw = ngk(ik)) to global ones (1 to ngk_k).
     !!
     !! See comment in SUBROUTINE gk_l2gmap_kdip in pw_restart_new.f90 for details.
     !-----------------------------------------------------------------------
     !
     USE mp,             ONLY : mp_sum, mp_max
     USE mp_bands,       ONLY : intra_bgrp_comm
     USE gvect,          ONLY : ig_l2g, mill
     USE klist,          ONLY : igk_k
     USE pw_restart_new, ONLY : gk_l2gmap_kdip
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: npw
     !! number of local plane waves
     INTEGER, INTENT(IN) :: ik
     !! k-point index
     !
     INTEGER :: ig
     !! G vector index
     INTEGER :: npw_g
     !! maximum global G-vector index
     INTEGER, ALLOCATABLE :: igk_l2g(:)
     !! Local to global mappings
     !
     ! Setup global plane wave mapping (following write_collected_wfc pattern)
     !
     ! Build local-to-global mapping for k+G vectors
     !
     ALLOCATE(igk_l2g(npw))
     igk_l2g = 0
     DO ig = 1, npw
        igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
     END DO
     !
     ! Find maximum global G-vector index
     !
     npw_g = MAXVAL(igk_l2g(1:npw))
     CALL mp_max(npw_g, intra_bgrp_comm)
     !
     ! Get global number of plane waves
     !
     ngk_g = npw
     CALL mp_sum(ngk_g, intra_bgrp_comm)
     !
     ! Build local-to-global-kdip mapping for collection/distribution
     !
     ALLOCATE(igk_l2g_kdip(npw))
     igk_l2g_kdip = 0
     CALL gk_l2gmap_kdip(npw_g, ngk_g, npw, igk_l2g, igk_l2g_kdip)
     !
     ! Collect Miller indices from all processors
     ! mill uses G-list type 2 (see comment in SUBROUTINE gk_l2gmap_kdip)
     ! We convert to type 4 in two steps:
     !   - (2) -> (1) : inverse of igk_k(:, ik)
     !   - (1) -> (4) : ig_l2g_kdip
     !
     ALLOCATE(mill_k_global(3, ngk_g))
     !
     mill_k_global(:, :) = 0
     DO ig = 1, npw
        mill_k_global(:, igk_l2g_kdip(ig)) = mill(:, igk_k(ig, ik))
     ENDDO
     CALL mp_sum(mill_k_global, intra_bgrp_comm)
     !
  END SUBROUTINE diag_direct_setup_gmap
  !-----------------------------------------------------------------------
  !
END MODULE diag_direct
