!
! Copyright (C) 2003-2024 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module wannier
   !
   ! Wannier module
   !
   USE kinds,      ONLY : DP
   USE fft_types,  ONLY : fft_type_descriptor
   !
   !integer, allocatable :: nnb(:)       ! #b  (ik)
   integer              :: nnb          ! #b
   integer, allocatable :: kpb(:,:)     ! k+b (ik,ib)
   integer, allocatable :: g_kpb(:,:,:) ! G_k+b (ipol,ik,ib)
   integer, allocatable :: ig_(:,:)     ! G_k+b (ipol,ik,ib)
   integer, allocatable :: lw(:,:), mw(:,:) ! l and m of wannier (16,n_wannier)
   integer, allocatable :: num_sph(:)   ! num. func. in lin. comb., (n_wannier)
   logical, allocatable :: excluded_band(:)
   ! begin change Lopez, Thonhauser, Souza
   integer  :: iun_nnkp,iun_mmn,iun_amn,iun_band,iun_spn,iun_plot,iun_parity,&
        nnbx,nexband,iun_uHu,&
        iun_uIu,& !ivo
   ! end change Lopez, Thonhauser, Souza
        iun_sHu, iun_sIu ! shc
   integer  :: n_wannier !number of WF
   integer  :: n_proj    !number of projection
   complex(DP), allocatable :: gf(:,:)  ! guding_function(npwx,n_wannier)
   complex(DP), allocatable :: gf_spinor(:,:)
   complex(DP), allocatable :: sgf_spinor(:,:)
   integer               :: ispinw, ikstart, ikstop, iknum, reduce_unk_factor
   character(LEN=15)     :: wan_mode    ! running mode
   logical               :: logwann, wvfn_formatted, write_unk, write_eig, &
   ! begin change Lopez, Thonhauser, Souza
                            write_amn,write_mmn,reduce_unk,write_spn,&
                            write_unkg,write_uhu,&
                            write_dmn,read_sym, & !YN
                            write_uIu, spn_formatted, uHu_formatted, uIu_formatted, & !ivo
   ! end change Lopez, Thonhauser, Souza
   ! shc
                            write_sHu, write_sIu, sHu_formatted, sIu_formatted, &
   ! end shc
   ! irreducible BZ
                            irr_bz, &
   ! vv: Begin SCDM keywords
                            scdm_proj
   character(LEN=15)     :: scdm_entanglement
   real(DP)              :: scdm_mu, scdm_sigma
   ! vv: End SCDM keywords
   ! input data from nnkp file
   real(DP), allocatable :: center_w(:,:)     ! center_w(3,n_wannier)
   integer,  allocatable :: spin_eig(:)
   real(DP), allocatable :: spin_qaxis(:,:)
   integer, allocatable  :: l_w(:), mr_w(:) ! l and mr of wannier (n_wannier) as from table 3.1,3.2 of spec.
   integer, allocatable  :: r_w(:)      ! index of radial function (n_wannier) as from table 3.3 of spec.
   real(DP), allocatable :: xaxis(:,:),zaxis(:,:) ! xaxis and zaxis(3,n_wannier)
   real(DP), allocatable :: alpha_w(:)  ! alpha_w(n_wannier) ( called zona in wannier spec)
   !
   real(DP), allocatable :: csph(:,:)    ! expansion coefficients of gf on QE ylm function (16,n_wannier)
   CHARACTER(len=256) :: seedname  = 'wannier'  ! prepended to file names in wannier90
   ! For implementation of wannier_lib
   integer               :: mp_grid(3)            ! dimensions of MP k-point grid
   real(DP)              :: rlatt(3,3),glatt(3,3) ! real and recip lattices (Cartesian co-ords, units of Angstrom)
   real(DP), allocatable :: kpt_latt(:,:)  ! k-points in crystal co-ords. kpt_latt(3,iknum)
   real(DP), allocatable :: atcart(:,:)    ! atom centres in Cartesian co-ords and Angstrom units. atcart(3,nat)
   integer               :: num_bands      ! number of bands left after exclusions
   character(len=3), allocatable :: atsym(:) ! atomic symbols. atsym(nat)
   integer               :: num_nnmax=12
   complex(DP), allocatable :: m_mat(:,:,:,:), a_mat(:,:,:)
   complex(DP), allocatable :: u_mat(:,:,:), u_mat_opt(:,:,:)
   logical, allocatable     :: lwindow(:,:)
   real(DP), allocatable    :: wann_centers(:,:),wann_spreads(:)
   real(DP)                 :: spreads(3)
   real(DP), allocatable    :: eigval(:,:)
   logical                  :: old_spinor_proj  ! for compatability for nnkp files prior to W90v2.0
   integer,allocatable :: rir(:,:)
   logical,allocatable :: zerophase(:,:)
   real(DP), allocatable :: bvec(:,:), xbvec(:,:)    ! bvectors
   integer, PARAMETER :: header_len = 60
   !! The length of header in amn/mmn/eig/... files
   !! For unformatted stream IO, this must be the same as wannier90 when reading those files.
   !
   REAL(DP), ALLOCATABLE :: xk_all(:, :)
   !! xk vector at all k points
   ! DFT+U variables
   logical               :: hubbard          ! If .true. then write WF projectors for DFT+U
   integer               :: exclude_ks_bands ! how many lowest-lying KS bands were excluded from the 
                                             ! wannierization (those which are below the energy where
                                             ! the wannierization starts)
   logical               :: wan2hub(500)     ! select WFs which will be used as the basis set
                                             ! for the Hubbard projectors
   TYPE(fft_type_descriptor) :: dfftp_c, dffts_c ! custom FFTs
   !
   CONTAINS
   !
   !----------------------------------------------------------------------------
   SUBROUTINE utility_setup_wfc_and_pw(ik_global, evc, npw, igk_k_ik)
      !-------------------------------------------------------------------------
      !! Read wavefunction evc from file for global k point ik_global.
      !! Setup the plane wave variables npw and igk_k_ik
      !-------------------------------------------------------------------------
      USE kinds,           ONLY : DP, i8b
      USE mp_pools,        ONLY : my_pool_id, nproc_pool, me_pool
      USE io_files,        ONLY : prefix, tmp_dir, nwordwfc, iunwfc
      USE wvfct,           ONLY : npwx
      USE pwcom,           ONLY : nbnd
      USE noncollin_module,ONLY : npol
      USE klist,           ONLY : nkstot, igk_k, ngk
      USE gvect,           ONLY : g, ngm
      USE gvecw,           ONLY : gcutw
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ik_global
      !! Global k point index
      COMPLEX(DP), INTENT(OUT) :: evc(npwx * npol, nbnd)
      !! wavefunction at ik_global, read from file
      INTEGER, INTENT(OUT) :: npw
      !! number of plane waves at ik_global
      INTEGER, INTENT(OUT) :: igk_k_ik(npwx)
      !! index of G vectors at ik_global
      !
      CHARACTER(LEN = 256) :: wfcfile
      !! Temp file
      INTEGER :: ipool
      !! Pool index where ik_global belongs to
      INTEGER :: ik_local
      !! Local k index of ik_global in ipool
      INTEGER :: iproc
      !! Processer index, which is the postfix of the filename.
      INTEGER :: iun
      !! File unit
      INTEGER(KIND=i8b) :: unf_recl
      !! Record length. Double precision to prevent integer overflow
      INTEGER :: direct_io_factor
      !! Factor for record length
      INTEGER :: ios
      !! Error number
      REAL(KIND = DP) :: dummy
      !! Dummy variable
      REAL(DP) :: g2kin_(npwx)
      !! Dummy g2kin_ to call gk_sort
      CHARACTER(len=6), EXTERNAL :: int_to_char
      !
      evc = (0.0_DP, 0.0_DP)
      !
      ! For a global k point index, find the pool and local k point index.
      CALL pool_and_local_kpoint_index(nkstot, ik_global, ipool, ik_local)
      !
      IF (ipool == my_pool_id) THEN
! print*, "pool ", my_pool_id, " Reading locally  pool, ik_local ", ipool, ik_local
         !
         ! If ipool is this pool, just read from the buffer
         !
         CALL davcio(evc, 2*nwordwfc, iunwfc, ik_local, -1)
         !
      ELSE
         !
         ! If ipool is not this pool, open wfc file, read, and close.
         !
         ! Processer id
         iproc = ipool * nproc_pool + me_pool
         !
         wfcfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc' // TRIM(int_to_char( iproc+1 ))
! print*, "pool ", my_pool_id, " Reading globally pool, ik_local ", ipool, ik_local, trim(wfcfile)
         !
         INQUIRE(IOLENGTH = direct_io_factor) dummy
         unf_recl = direct_io_factor * INT(2*nwordwfc, KIND=KIND(unf_recl))
         !
         OPEN(NEWUNIT = iun, FILE = TRIM(wfcfile), FORM = 'unformatted', &
            ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl)
         IF (ios /= 0) CALL errore('utility_setup_wfc_and_pw', &
            'error opening wfc file', 1)
         READ(iun, REC = ik_local) evc
         CLOSE(iun, STATUS = 'KEEP')
      ENDIF
      !
      ! Setup npw and igk_k_ik, the G vector ordering at ik_global.
      !
      IF (ipool == my_pool_id) THEN
         ! Use local G vector ordering
         npw = ngk(ik_local)
         igk_k_ik = igk_k(:, ik_local)
      ELSE
         ! k point from different pool. Calculate G vector ordering.
         CALL gk_sort(xk_all(1, ik_global), ngm, g, gcutw, npw, igk_k_ik, g2kin_)
      ENDIF
      !
   !----------------------------------------------------------------------------
   END SUBROUTINE utility_setup_wfc_and_pw
   !----------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------
   SUBROUTINE utility_compute_u_kb(ik, i_b, evc_kb, evc_kb_m, becp_kb)
      !-------------------------------------------------------------------------
      !!
      !! For a given k point ik and b vector i_b, compute |u_{n, k+b}>.
      !! Optionally compute becp for k+b (for USPP case).
      !! Do not calculate the excluded bands, using the variable excluded_band.
      !!
      !! Uses pre-computed information about b vector: kpb, zerophase, and ig_.
      !!
      !! TODO: Better name for evc_b
      !!
      !-------------------------------------------------------------------------
      !
      USE kinds,           ONLY : DP
      USE fft_base,        ONLY : dffts
      USE fft_interfaces,  ONLY : fwfft, invfft
      USE mp_pools,        ONLY : my_pool_id
      USE io_files,        ONLY : nwordwfc, iunwfc
      USE control_flags,   ONLY : gamma_only
      USE gvect,           ONLY : gstart, g, ngm
      USE gvecw,           ONLY : gcutw
      USE wvfct,           ONLY : nbnd, npwx
      USE wavefunctions,   ONLY : psic
      USE klist,           ONLY : xk, ngk, igk_k, nkstot
      USE noncollin_module,ONLY : npol
      USE becmod,          ONLY : bec_type, calbec
      USE uspp,            ONLY : vkb
      USE uspp_init,       ONLY : init_us_2
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ik
      !! Local index of the k point
      INTEGER, INTENT(IN) :: i_b
      !! Index of the b vector
      COMPLEX(DP), INTENT(OUT) :: evc_kb(npol*npwx, num_bands)
      !! Wavefunction u_{n,k+b}(G), with the G-vector ordering at k.
      COMPLEX(DP), OPTIONAL, INTENT(OUT) :: evc_kb_m(npwx, num_bands)
      !! Optional. For Gamma-only case. Complex conjugate of u_{n,k+b}(-G).
      TYPE(bec_type), OPTIONAL, INTENT(INOUT) :: becp_kb
      !! Optional. Product of beta functions with |psi_k+b>
      !
      LOGICAL :: zerophase_kb
      !! zerophase at (k, b)
      INTEGER :: ig_kb
      !! G vector index ig_ for (k, b)
      INTEGER :: ik_g_w90
      !! Global index of k+b. For lsda with spin down, subtract nkstot/2 to
      !! start from 1.
      INTEGER :: ikp_b
      !! Index of k+b
      INTEGER :: npw
      !! Number of plane waves at k
      INTEGER :: npw_b
      !! Number of plane waves at k+b
      INTEGER :: istart
      !! Starting G vector index
      INTEGER :: iend
      !! Ending G vector index
      INTEGER :: ipol
      !! Polarization index
      INTEGER :: n
      !! Band index
      INTEGER :: n_incl
      !! Band index after exclusion
      INTEGER :: ipool_b
      !! Pool index where ikp_b belongs to
      INTEGER :: ik_b_local
      !! Local k index of ikp_b in ipool_b
      INTEGER :: igk_kb(npwx)
      !! G vector index at k+b
      REAL(DP) :: g2kin_(npwx)
      !! Dummy g2kin_ to call gk_sort
      !
      COMPLEX(DP), ALLOCATABLE :: phase(:)
      !! Temporary storage for phase e&{i * G_kpb * r}
      COMPLEX(DP), ALLOCATABLE :: evc_b(:, :)
      !! Temporary storage for wavefunction at k+b
      INTEGER :: ierr
      !! Error number
      !
      INTEGER, EXTERNAL :: global_kpoint_index
      !
      IF (PRESENT(evc_kb_m) .AND. (.NOT. gamma_only)) CALL errore("utility_compute_u_kb", &
         "evc_kb_m can be used only in the Gamma-only case", 1)
      !
      CALL start_clock("compute_u_kb")
      !
      ALLOCATE(phase(dffts%nnr), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phase', 1)
      ALLOCATE(evc_b(npol*npwx, nbnd), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_b', 1)
      !
      evc_kb = (0.0_DP, 0.0_DP)
      !
      ik_g_w90 = global_kpoint_index(nkstot, ik) - ikstart + 1
      !
      npw = ngk(ik)
      !
      ikp_b = kpb(ik_g_w90, i_b) + ikstart - 1
      ig_kb = ig_(ik_g_w90, i_b)
      zerophase_kb = zerophase(ik_g_w90, i_b)
      !
      ! Read wavefunction and setup ngk, igk_k for ikp_b
      !
      CALL utility_setup_wfc_and_pw(ikp_b, evc_b, npw_b, igk_kb)
      !
      ! Exclude the excluded bands.
      !
      n_incl = 0
      DO n = 1, nbnd
         IF (excluded_band(n)) CYCLE
         n_incl = n_incl + 1
         evc_b(:, n_incl) = evc_b(:, n)
      ENDDO
      !
      ! Compute the phase e^{i * g_kpb * r} if phase is not 1.
      ! Computed phase is used inside the loop over bands.
      !
      IF (.NOT. zerophase_kb) THEN
         phase(:) = (0.0_DP, 0.0_DP)
         IF (ig_kb > 0) phase( dffts%nl(ig_kb) ) = (1.0_DP, 0.0_DP)
         CALL invfft('Wave', phase, dffts)
      ENDIF
      !
      ! Loop over bands
      !
      DO n = 1, num_bands
         !
         DO ipol = 1, npol
            psic = (0.0_DP, 0.0_DP)
            !
            ! Copy evc_b to psic using the G-vector ordering at k+b.
            !
            istart = 1 + (ipol-1) * npwx
            iend = istart + npw_b - 1
            psic(dffts%nl( igk_kb(1:npw_b) )) = evc_b(istart:iend, n)
            !
            IF (gamma_only) psic(dffts%nlm( igk_kb(1:npw_b) )) = CONJG(evc_b(istart:iend, n))
            !
            ! Multiply e^{i * g_kpb * r} phase if necessary.
            ! TODO: Add explanation on the phase factor
            !
            IF (.NOT. zerophase_kb) THEN
               CALL invfft('Wave', psic, dffts)
               psic(1:dffts%nnr) = psic(1:dffts%nnr) * CONJG(phase(1:dffts%nnr))
               CALL fwfft('Wave', psic, dffts)
            ENDIF
            !
            ! Save |u_{n, k+b}> in evc_kb using the G-vector ordering at k.
            !
            istart = 1 + (ipol-1) * npwx
            iend = istart + npw - 1
            evc_kb(istart:iend, n) = psic(dffts%nl( igk_k(1:npw, ik) ))
            !
            IF (gamma_only) THEN
               IF (gstart == 2) psic(dffts%nlm(1)) = (0.0_DP, 0.0_DP)
               evc_kb_m(1:npw, n) = CONJG(psic(dffts%nlm( igk_k(1:npw, ik) )))
            ENDIF
            !
         ENDDO ! npol
         !
      ENDDO ! nbnd
      !
      IF (PRESENT(becp_kb)) THEN
         !
         ! Compute the product of beta functions with |psi_k+b>
         !
         CALL init_us_2(npw_b, igk_kb, xk_all(1, ikp_b), vkb)
         CALL calbec(npw_b, vkb, evc_b, becp_kb, num_bands)
         !
      ENDIF
      !
      DEALLOCATE(phase)
      DEALLOCATE(evc_b)
      !
      CALL stop_clock("compute_u_kb")
      !
   !----------------------------------------------------------------------------
   END SUBROUTINE utility_compute_u_kb
   !----------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------
   SUBROUTINE utility_merge_files(postfix, formatted, ndata)
      !-------------------------------------------------------------------------
      !! For pool parallelization, each root_pool writes to different files.
      !! Here, concatenate all prefix.postfix files
      !-------------------------------------------------------------------------
      !
      USE kinds,           ONLY : DP
      USE io_global,       ONLY : ionode
      USE mp_pools,        ONLY : npool
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: postfix
      !! postfix for filename
      LOGICAL, INTENT(IN) :: formatted
      !! True if formatted file, false if unformatted file.
      INTEGER, INTENT(IN), OPTIONAL :: ndata
      !! Used only if formatted = .FALSE. Length of data on each line.
      !
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=256) :: line
      INTEGER :: ipool, iun, iun2, i, ierr
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      COMPLEX(DP), ALLOCATABLE :: arr(:)
      !
      IF ((.NOT. formatted) .AND. (.NOT. PRESENT(ndata))) THEN
         CALL errore("utility_merge_files", "If formatted is false, ndata must be provided", 1)
      ENDIF
      !
      IF (npool == 1) RETURN
      IF (.NOT. ionode) RETURN
      !
      IF (formatted) THEN
         !
         filename = TRIM(seedname) // "." // TRIM(postfix)
         OPEN (NEWUNIT=iun, file=TRIM(filename), form='formatted', STATUS="OLD", POSITION="APPEND")
         !
         DO ipool = 2, npool
            filename = TRIM(seedname) // "." // TRIM(postfix) // TRIM(int_to_char(ipool))
            OPEN(NEWUNIT=iun2, FILE=TRIM(filename), FORM='formatted')
            DO WHILE (.TRUE.)
               READ(iun2, '(A)', END=200) line
               WRITE(iun, '(A)') TRIM(line)
            ENDDO
200         CLOSE(iun2, STATUS="DELETE")
         ENDDO
         !
         CLOSE(iun, STATUS="KEEP")
         !
      ELSE ! .NOT. formatted
         !
         ALLOCATE(arr(ndata), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating arr', 1)
         !
         filename = TRIM(seedname) // "." // TRIM(postfix)
         OPEN (NEWUNIT=iun, file=TRIM(filename), form='unformatted', STATUS="OLD", POSITION="APPEND")
         !
         DO ipool = 2, npool
            filename = TRIM(seedname) // "." // TRIM(postfix) // TRIM(int_to_char(ipool))
            OPEN(NEWUNIT=iun2, FILE=TRIM(filename), FORM='unformatted')
            DO WHILE (.TRUE.)
               READ(iun2, END=201) (arr(i), i=1, ndata)
               WRITE(iun) (arr(i), i=1, ndata)
            ENDDO
201      CLOSE(iun2, STATUS="DELETE")
         ENDDO
         !
         CLOSE(iun, STATUS="KEEP")
         !
         DEALLOCATE(arr)
         !
      ENDIF ! formatted
   !----------------------------------------------------------------------------
   END SUBROUTINE utility_merge_files
   !----------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------
   SUBROUTINE print_progress(i, n)
      !-------------------------------------------------------------------------
      !! Print progress of iterations
      !-------------------------------------------------------------------------
      !
      USE io_global,  ONLY : stdout
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: i
      !! Current iteration
      INTEGER, INTENT(IN) :: n
      !! Endpoint of the iteration
      !
      WRITE(stdout, '(I8)', ADVANCE='NO') i
      !! Newline every 10 iterations
      IF (MOD(i, 10) == 0) WRITE(stdout, *)
      !! Newline at the end of the iteration if not already done
      IF ((i == n) .AND. (MOD(n, 10) /= 0)) WRITE(stdout, *)
      FLUSH(stdout)
      !
   !----------------------------------------------------------------------------
   END SUBROUTINE print_progress
   !----------------------------------------------------------------------------
   !
end module wannier
