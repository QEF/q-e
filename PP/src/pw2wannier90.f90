!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! pw2wannier was written by Stefano de Gironcoli
! with later additions by
! Jonathan Yates - spinors
! Arash Mostofi - gamma point and transport things
! Timo Thonhauser, Graham Lopez, Ivo Souza
!         uHu, uIu terms for orbital magnetisation
! please send bugs and comments to
! Jonathan Yates and Arash Mostofi
! Takashi Koretsune and Florian Thoele -- noncollinear and USPPs
! Valerio Vitale - Selected columns of density matrix (SCDM)
! Jae-Mo Lihm - SCDM with noncollinear
! Ji Hoon Ryoo, Minsu Ghim - sHu, sIu terms for spin Hall conductivity
!
! NOTE: old_spinor_proj is deprecated.
!
!
module wannier
   USE kinds, only : DP
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
   !
   CONTAINS
   !
   !----------------------------------------------------------------------------
   SUBROUTINE utility_setup_wfc_and_pw(ik_global, evc, npw, igk_k_ik)
      !-------------------------------------------------------------------------
      !! Read wavefunction evc from file for global k point ik_global.
      !! Setup the plane wave variables npw and igk_k_ik
      !-------------------------------------------------------------------------
      USE kinds,           ONLY : DP
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
      INTEGER*8 :: unf_recl
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
!

module atproj
   ! module for atomic projectors
   USE kinds, ONLY: DP

   ! Atomic projectors for species
   TYPE atproj_type
      CHARACTER(len=3) :: atsym ! atomic symbol
      INTEGER :: ngrid ! number of grid points
      REAL(DP), ALLOCATABLE :: xgrid(:)
      REAL(DP), ALLOCATABLE :: rgrid(:) ! exp(xgrid(:))
      INTEGER :: nproj ! total number of projectors for this species
      INTEGER, ALLOCATABLE :: l(:) ! angular momentum of each wfc
      REAL(DP), ALLOCATABLE :: radial(:, :) ! magnitude at each point, ngrid x nproj
   END TYPE atproj_type

   ! atomic proj input variables, start with atom_proj*
   LOGICAL :: atom_proj
   CHARACTER(LEN=256) :: atom_proj_dir ! directory of external projectors
   LOGICAL :: atom_proj_ext ! switch for using external files instead of orbitals from UPF
   INTEGER, PARAMETER :: nexatproj_max = 2000 ! max allowed number of projectors to be excluded
   INTEGER :: atom_proj_exclude(nexatproj_max) ! index starts from 1
   LOGICAL :: atom_proj_ortho ! whether perform Lowdin orthonormalization
   LOGICAL :: atom_proj_sym ! whether perform symmetrization

   ! atomic proj internal variables, using *atproj*
   INTEGER :: nexatproj ! actual number of excluded projectors
   INTEGER :: natproj ! total number of projectors = n_proj + nexatproj
   LOGICAL, ALLOCATABLE :: atproj_excl(:) ! size = total num of projectors
   INTEGER :: iun_atproj
   TYPE(atproj_type), ALLOCATABLE :: atproj_typs(:) ! all atom proj types

   REAL(DP), ALLOCATABLE :: tab_at(:, :, :)
   !! interpolation table for atomic projectors

CONTAINS

   SUBROUTINE skip_comments(file_unit)
      ! read a file_unit and skip lines starting with #
      !
      USE, INTRINSIC :: ISO_FORTRAN_ENV
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: file_unit
      CHARACTER(len=256) :: line
      INTEGER :: ret_code
      !
      DO
      READ (file_unit, '(A)', iostat=ret_code) line
      IF (ret_code == IOSTAT_END) EXIT
      IF (ret_code /= 0) THEN
         ! read error
         EXIT
      END IF
      IF (INDEX(ADJUSTL(line), "#") == 1) CYCLE
      EXIT
      END DO

      BACKSPACE file_unit

      RETURN
   END SUBROUTINE skip_comments

   SUBROUTINE allocate_atproj_type(typ, ngrid, nproj)
      !
      IMPLICIT NONE
      !
      TYPE(atproj_type), INTENT(INOUT) :: typ
      INTEGER, INTENT(IN) :: ngrid, nproj
      INTEGER :: ierr
      !
      typ%ngrid = ngrid
      typ%nproj = nproj
      ALLOCATE (typ%xgrid(ngrid), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating typ%xgrid', 1)
      ALLOCATE (typ%rgrid(ngrid), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating typ%rgrid', 1)
      ALLOCATE (typ%l(nproj), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating typ%l', 1)
      ALLOCATE (typ%radial(ngrid, nproj), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating typ%radial', 1)

      RETURN
   END SUBROUTINE allocate_atproj_type

   SUBROUTINE read_atomproj(typs)
      !
      ! read data files for atom proj
      ! should be called only by one node
      !
      USE kinds, ONLY: dp
      USE io_global, ONLY: stdout
      USE ions_base, ONLY: nsp, atm
      !
      IMPLICIT NONE
      !
      TYPE(atproj_type), INTENT(INOUT) :: typs(nsp)
      !
      INTEGER :: i, j, it
      LOGICAL :: file_exists
      CHARACTER(len=256) :: filename
      INTEGER :: ngrid, nproj

      DO it = 1, nsp
         filename = TRIM(atom_proj_dir)//'/'//TRIM(atm(it))//".dat"
         INQUIRE (FILE=TRIM(filename), EXIST=file_exists)
         IF (.NOT. file_exists) &
            CALL errore('pw2wannier90', 'file not exists: '//TRIM(filename), 1)
         OPEN (NEWUNIT=iun_atproj, file=TRIM(filename), form='formatted')
         CALL skip_comments(iun_atproj)

         READ (iun_atproj, *) ngrid, nproj
         WRITE (stdout, *) " Read from "//TRIM(filename)
         WRITE (stdout, '((A),(I4))') "   number of grid points   = ", ngrid
         WRITE (stdout, '((A),(I4))') "   number of projectors    = ", nproj

         CALL allocate_atproj_type(typs(it), ngrid, nproj)
         typs(it)%atsym = atm(it)

         READ (iun_atproj, *) (typs(it)%l(i), i=1, nproj)
         WRITE (stdout, '((A))', advance='no') "   ang. mom. of projectors = "
         DO i = 1, nproj
            WRITE (stdout, '(I4)', advance='no') typs(it)%l(i)
         END DO
         WRITE (stdout, *)
         WRITE (stdout, *)

         DO i = 1, ngrid
            READ (iun_atproj, *) typs(it)%xgrid(i), typs(it)%rgrid(i), &
            (typs(it)%radial(i, j), j=1, nproj)
            ! PRINT *, i, typs(it)%xgrid(i), typs(it)%rgrid(i)
         END DO

         CLOSE (iun_atproj)
      END DO

      RETURN
   END SUBROUTINE read_atomproj

   !-----------------------------------------------------------------------
   SUBROUTINE atomproj_wfc(ik, wfcatom)
      !-----------------------------------------------------------------------
      !! This routine computes the superposition of atomic wavefunctions
      !! for k-point "ik" - output in "wfcatom".
      !
      !  adapted from PW/src/atomic_wfc.f90, PP/src/atomic_wfc_nc_proj.f90
      !  to use external atomic wavefunctions other than UPF ones.
      !
      USE kinds, ONLY: DP
      USE constants, ONLY: tpi, fpi, pi
      USE cell_base, ONLY: omega, tpiba
      USE ions_base, ONLY: nat, ntyp => nsp, ityp, tau
      USE gvect, ONLY: mill, eigts1, eigts2, eigts3, g
      USE klist, ONLY: xk, igk_k, ngk
      USE wvfct, ONLY: npwx
      USE noncollin_module, ONLY: noncolin, npol, angle1, angle2, lspinorb, &
                                  domag, starting_spin_angle
      USE upf_spinorb,ONLY : rot_ylm, lmaxx, fcoef, lmaxx
      USE wannier, ONLY: n_proj
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ik
      !! k-point index
      COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx, npol, natproj)
      !! Superposition of atomic wavefunctions
      !
      ! ... local variables
      !
      INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
                 i0, i1, i2, i3, npw
      REAL(DP), ALLOCATABLE :: qg(:), ylm(:, :), chiq(:, :, :), gk(:, :)
      COMPLEX(DP), ALLOCATABLE :: sk(:), aux(:)
      COMPLEX(DP) :: kphase, lphase
      REAL(DP)    :: arg
      INTEGER :: nwfcm
      !! max number of radial atomic projectors across atoms
      INTEGER :: ierr

      CALL start_clock('atomproj_wfc')

      ! calculate max angular momentum required in wavefunctions
      lmax_wfc = 0
      nwfcm = 0
      DO nt = 1, ntyp
         lmax_wfc = MAX(lmax_wfc, MAXVAL(atproj_typs(nt)%l))
         nwfcm = MAX(nwfcm, atproj_typs(nt)%nproj)
      END DO
      !
      npw = ngk(ik)
      !
      ALLOCATE (ylm(npw, (lmax_wfc + 1)**2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ylm', 1)
      ALLOCATE (chiq(npw, nwfcm, ntyp), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating chiq', 1)
      ALLOCATE (gk(3, npw), qg(npw), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating gk/qg', 1)
      !
      DO ig = 1, npw
         iig = igk_k(ig, ik)
         gk(1, ig) = xk(1, ik) + g(1, iig)
         gk(2, ig) = xk(2, ik) + g(2, iig)
         gk(3, ig) = xk(3, ik) + g(3, iig)
         qg(ig) = gk(1, ig)**2 + gk(2, ig)**2 + gk(3, ig)**2
      END DO
      !
      !  ylm = spherical harmonics
      !
      CALL ylmr2((lmax_wfc + 1)**2, npw, gk, qg, ylm)
      !
      ! set now q=|k+G| in atomic units
      !
      DO ig = 1, npw
         qg(ig) = SQRT(qg(ig))*tpiba
      END DO
      !
      CALL interp_atproj(npw, qg, nwfcm, chiq)
      !
      DEALLOCATE (qg, gk)
      ALLOCATE (aux(npw), sk(npw), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating aux/sk', 1)
      !
      wfcatom(:, :, :) = (0.0_DP, 0.0_DP)
      n_starting_wfc = 0
      !
      DO na = 1, nat
        arg = (xk(1, ik)*tau(1, na) + xk(2, ik)*tau(2, na) + xk(3, ik)*tau(3, na))*tpi
        kphase = CMPLX(COS(arg), -SIN(arg), KIND=DP)
        !
        !     sk is the structure factor
        !
        DO ig = 1, npw
          iig = igk_k(ig, ik)
          sk(ig) = kphase*eigts1(mill(1, iig), na)* &
                   eigts2(mill(2, iig), na)* &
                   eigts3(mill(3, iig), na)
        END DO
        !
        nt = ityp(na)
        DO nb = 1, atproj_typs(nt)%nproj
          l = atproj_typs(nt)%l(nb)
          lphase = (0.D0, 1.D0)**l
          !
          !  only support without spin-orbit coupling
          !  if needed in the future, add code according to
          !  PP/src/atomic_wfc_nc_proj.f90
          !
          CALL atomic_wfc___()
          !
          ! END IF
          !
        END DO
        !
      END DO

      DEALLOCATE (aux, sk, chiq, ylm)

      CALL stop_clock('atomproj_wfc')

      RETURN

   CONTAINS

      SUBROUTINE atomic_wfc___()
        !
        ! ... LSDA or nonmagnetic case
        !
        DO m = 1, 2*l + 1
          lm = l**2 + m
          n_starting_wfc = n_starting_wfc + 1
          IF (n_starting_wfc > natproj) CALL errore &
            ('atomic_wfc___', 'internal error: too many wfcs', 1)
          !
          DO ig = 1, npw
            wfcatom(ig, 1, n_starting_wfc) = lphase* &
                                             sk(ig)*ylm(ig, lm)*chiq(ig, nb, nt)
          END DO
          !
        END DO
        !
      END SUBROUTINE atomic_wfc___
      !
   END SUBROUTINE atomproj_wfc

   SUBROUTINE interp_atproj(npw, qg, nwfcm, chiq)
      !-----------------------------------------------------------------------
      !
      ! computes chiq: radial fourier transform of atomic projector chi
      !
      !! adapted from upflib/interp_atwfc.f90
      !  to support external projectors
      !
      USE kinds, ONLY: dp
      USE ions_base, ONLY: nsp
      USE uspp_data, ONLY: dq
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: npw
      INTEGER, INTENT(IN)  :: nwfcm
      REAL(dp), INTENT(IN) :: qg(npw)
      REAL(dp), INTENT(OUT):: chiq(npw, nwfcm, nsp)
      !
      INTEGER :: nt, nb, ig
      INTEGER :: i0, i1, i2, i3
      REAL(dp):: qgr, px, ux, vx, wx
      !
      DO nt = 1, nsp
        DO nb = 1, atproj_typs(nt)%nproj
          DO ig = 1, npw
            qgr = qg(ig)
            px = qgr/dq - INT(qgr/dq)
            ux = 1.D0 - px
            vx = 2.D0 - px
            wx = 3.D0 - px
            i0 = INT(qgr/dq) + 1
            i1 = i0 + 1
            i2 = i0 + 2
            i3 = i0 + 3
            chiq(ig, nb, nt) = &
              tab_at(i0, nb, nt)*ux*vx*wx/6.D0 + &
              tab_at(i1, nb, nt)*px*vx*wx/2.D0 - &
              tab_at(i2, nb, nt)*px*ux*wx/2.D0 + &
              tab_at(i3, nb, nt)*px*ux*vx/6.D0
          END DO
        END DO
      END DO

   END SUBROUTINE interp_atproj

   SUBROUTINE init_tab_atproj(intra_bgrp_comm)
      !-----------------------------------------------------------------------
      !! This routine computes a table with the radial Fourier transform
      !! of the atomic wavefunctions.
      !!
      !! adapted from upflib/init_tab_atwfc.f90
      !! to support external projectors
      !
      USE kinds, ONLY: DP
      USE upf_const, ONLY: fpi
      USE uspp_data, ONLY: nqx, dq
      USE ions_base, ONLY: nsp
      USE cell_base, ONLY: omega
      USE mp, ONLY: mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: intra_bgrp_comm
      !
      INTEGER :: nt, nb, iq, ir, l, startq, lastq, ndm, nwfcm, ierr
      !
      REAL(DP), ALLOCATABLE :: aux(:), vchi(:), rab(:)
      REAL(DP) :: vqint, pref, q
      !
      ndm = 0
      nwfcm = 0
      DO nt = 1, nsp
        ndm = MAX(ndm, atproj_typs(nt)%ngrid)
        nwfcm = MAX(nwfcm, atproj_typs(nt)%nproj)
      END DO
      ALLOCATE (aux(ndm), vchi(ndm), rab(ndm), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating aux/vchi/rab', 1)
      !
      ! chiq = radial fourier transform of atomic orbitals chi
      !
      pref = fpi/SQRT(omega)
      ! needed to normalize atomic wfcs (not a bad idea in general)
      CALL divide(intra_bgrp_comm, nqx, startq, lastq)
      !
      ! nqx = INT( (SQRT(ecutwfc) / dq + 4) )
      ALLOCATE (tab_at(nqx, nwfcm, nsp), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating tab_at', 1)
      tab_at(:, :, :) = 0.0_DP
      !
      DO nt = 1, nsp
        rab = (atproj_typs(nt)%xgrid(2) - atproj_typs(nt)%xgrid(1))* &
              atproj_typs(nt)%rgrid
        DO nb = 1, atproj_typs(nt)%nproj
          !
          l = atproj_typs(nt)%l(nb)
          !
          DO iq = startq, lastq
            q = dq*(iq - 1)
            CALL sph_bes(atproj_typs(nt)%ngrid, atproj_typs(nt)%rgrid, q, l, aux)
            DO ir = 1, atproj_typs(nt)%ngrid
              vchi(ir) = atproj_typs(nt)%radial(ir, nb)*aux(ir)*atproj_typs(nt)%rgrid(ir)
            END DO
            CALL simpson(atproj_typs(nt)%ngrid, vchi, rab, vqint)
            tab_at(iq, nb, nt) = vqint*pref
          END DO
          !
        END DO
      END DO
      !
      CALL mp_sum(tab_at, intra_bgrp_comm)
      !
      DEALLOCATE (aux, vchi, rab)
      !
      RETURN
      !
   END SUBROUTINE init_tab_atproj

   SUBROUTINE deallocate_atproj
      IMPLICIT NONE

      IF (ALLOCATED(tab_at)) DEALLOCATE (tab_at)

      RETURN
   END SUBROUTINE deallocate_atproj

end module atproj
!
!------------------------------------------------------------------------
PROGRAM pw2wannier90
  ! This is the interface to the Wannier90 code: see http://www.wannier.org
  !------------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup
  USE mp_pools,   ONLY : npool
  USE mp_bands,   ONLY : nbgrp
  USE mp,         ONLY : mp_bcast, mp_barrier
  USE mp_world,   ONLY : world_comm
  USE cell_base,  ONLY : at, bg
  USE lsda_mod,   ONLY : nspin, isk
  USE klist,      ONLY : nkstot, nks, xk
  USE io_files,   ONLY : prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE control_flags,    ONLY : gamma_only
  USE environment,ONLY : environment_start, environment_end
  USE wannier
  use atproj, only : atom_proj, atom_proj_dir, atom_proj_ext, &
                     atom_proj_exclude, atom_proj_ortho, atom_proj_sym
  USE read_namelists_module, only : check_namelist_read
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, ierr
  CHARACTER(len=4) :: spin_component
  CHARACTER(len=256) :: outdir

  ! these are in wannier module.....-> integer :: ispinw, ikstart, ikstop, iknum
  NAMELIST / inputpp / outdir, prefix, spin_component, wan_mode, &
       seedname, write_unk, write_amn, write_mmn, write_spn, write_eig,&
   ! begin change Lopez, Thonhauser, Souza
       wvfn_formatted, reduce_unk, reduce_unk_factor, &
       write_unkg, write_uhu, write_dmn, read_sym, & !YN:
       write_uIu, spn_formatted, uHu_formatted, uIu_formatted,& !ivo
   ! end change Lopez, Thonhauser, Souza
   ! shc
       write_sHu, write_sIu, sHu_formatted, sIu_formatted,&
   ! end shc
       irr_bz,& ! Koretsune
   ! begin change Vitale
       scdm_proj, scdm_entanglement, scdm_mu, scdm_sigma, &
   ! end change Vitale
       atom_proj, atom_proj_dir, atom_proj_ext, atom_proj_exclude, &
       atom_proj_ortho
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  !! not sure if this should be called also in 'library' mode or not !!
  CALL environment_start ( 'PW2WANNIER' )
  !
  CALL start_clock( 'init_pw2wan' )
  !
  ! Read input on i/o node and broadcast to the rest
  !
  ios = 0
  IF(ionode) THEN
     !
     ! Check to see if we are reading from a file
     !
     CALL input_from_file()
     !
     !   set default values for variables in namelist
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     prefix = ' '
     seedname = 'wannier'
     spin_component = 'none'
     wan_mode = 'standalone'
     wvfn_formatted = .false.
     spn_formatted=.false.
     uHu_formatted=.false.
     uIu_formatted=.false.
     write_unk = .false.
     write_amn = .true.
     write_mmn = .true.
     write_spn = .false.
     write_eig = .true.
     ! begin change Lopez, Thonhauser, Souza
     write_uhu = .false.
     write_uIu = .false. !ivo
     ! end change Lopez, Thonhauser, Souza
     ! shc
     write_sHu = .false.
     write_sIu = .false.
     sHu_formatted=.false.
     sIu_formatted=.false.
     ! end shc
     reduce_unk= .false.
     reduce_unk_factor = -9999
     write_unkg= .false.
     write_dmn = .false. !YN:
     read_sym  = .false. !YN:
     irr_bz = .false.
     scdm_proj = .false.
     scdm_entanglement = 'isolated'
     scdm_mu = 0.0_dp
     scdm_sigma = 1.0_dp
     atom_proj = .false.
     atom_proj_dir = './'
     atom_proj_ext = .false.
     atom_proj_ortho = .true.
     atom_proj_exclude = -1
     ! Haven't tested symmetrization with external projectors, disable it for now
     atom_proj_sym = .false.
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat=ios)
     !
     !     Check of namelist variables
     !
     tmp_dir = trimcheck(outdir)
     ! back to all nodes
  ENDIF
  !
  CALL mp_bcast(ios,ionode_id, world_comm)
  CALL check_namelist_read(ios, 5, "inputpp")
  !
  ! broadcast input variable to all nodes
  !
  CALL mp_bcast(outdir,ionode_id, world_comm)
  CALL mp_bcast(tmp_dir,ionode_id, world_comm)
  CALL mp_bcast(prefix,ionode_id, world_comm)
  CALL mp_bcast(seedname,ionode_id, world_comm)
  CALL mp_bcast(spin_component,ionode_id, world_comm)
  CALL mp_bcast(wan_mode,ionode_id, world_comm)
  CALL mp_bcast(wvfn_formatted,ionode_id, world_comm)
  CALL mp_bcast(write_unk,ionode_id, world_comm)
  CALL mp_bcast(write_amn,ionode_id, world_comm)
  CALL mp_bcast(write_mmn,ionode_id, world_comm)
  CALL mp_bcast(write_eig,ionode_id, world_comm)
  ! begin change Lopez, Thonhauser, Souza
  CALL mp_bcast(write_uhu,ionode_id, world_comm)
  CALL mp_bcast(write_uIu,ionode_id, world_comm) !ivo
  ! end change Lopez, Thonhauser, Souza
  ! shc
  CALL mp_bcast(write_sHu,ionode_id, world_comm)
  CALL mp_bcast(write_sIu,ionode_id, world_comm)
  ! end shc
  CALL mp_bcast(write_spn,ionode_id, world_comm)
  CALL mp_bcast(reduce_unk,ionode_id, world_comm)
  CALL mp_bcast(reduce_unk_factor,ionode_id, world_comm)
  CALL mp_bcast(write_unkg,ionode_id, world_comm)
  CALL mp_bcast(write_dmn,ionode_id, world_comm)
  CALL mp_bcast(read_sym,ionode_id, world_comm)
  CALL mp_bcast(irr_bz,ionode_id, world_comm)
  CALL mp_bcast(scdm_proj,ionode_id, world_comm)
  CALL mp_bcast(scdm_entanglement,ionode_id, world_comm)
  CALL mp_bcast(scdm_mu,ionode_id, world_comm)
  CALL mp_bcast(scdm_sigma,ionode_id, world_comm)
  CALL mp_bcast(spn_formatted, ionode_id, world_comm)
  CALL mp_bcast(uHu_formatted, ionode_id, world_comm)
  CALL mp_bcast(uIu_formatted, ionode_id, world_comm)
  CALL mp_bcast(sHu_formatted, ionode_id, world_comm)
  CALL mp_bcast(sIu_formatted, ionode_id, world_comm)
  CALL mp_bcast(atom_proj, ionode_id, world_comm)
  CALL mp_bcast(atom_proj_dir, ionode_id, world_comm)
  CALL mp_bcast(atom_proj_ext,ionode_id, world_comm)
  CALL mp_bcast(atom_proj_ortho, ionode_id, world_comm)
  CALL mp_bcast(atom_proj_sym, ionode_id, world_comm)
  CALL mp_bcast(atom_proj_exclude, ionode_id, world_comm)
  !
  ! Check: kpoint distribution with pools in library mode not implemented
  !
  IF (npool > 1 .and. wan_mode == 'library') CALL errore('pw2wannier90', &
      'pools not implemented for library mode', 1)
  !
  ! Check: bands distribution not implemented
  IF (nbgrp > 1) CALL errore('pw2wannier90', 'bands (-nb) not implemented', nbgrp)
  !
  IF (reduce_unk_factor == -9999) THEN
     ! If reduce_unk_factor is not provided, set it based on reduce_unk.
     IF (reduce_unk) THEN
        reduce_unk_factor = 2
     ELSE
        reduce_unk_factor = 1
     ENDIF
  ENDIF
  IF (reduce_unk_factor < 1) CALL errore('pw2wannier90', 'reduce_unk_factor < 1', 1)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  logwann = .true.
  WRITE(stdout,*)
  WRITE(stdout,'(5x,A)') 'Reading nscf_save data'
  CALL read_file
  WRITE(stdout,*)
  !
  !
  IF (npool > 1 .and. gamma_only) CALL errore('pw2wannier90', &
      'pools not compatible with gamma_only', 1)
  !
  IF (noncolin.and.gamma_only) CALL errore('pw2wannier90',&
       'Non-collinear and gamma_only not implemented',1)
  IF (gamma_only.and.scdm_proj) CALL errore('pw2wannier90',&
       'Gamma_only and SCDM not implemented',1)
  IF (scdm_proj) then
    IF ((trim(scdm_entanglement) /= 'isolated') .AND. &
        (trim(scdm_entanglement) /= 'erfc') .AND. &
        (trim(scdm_entanglement) /= 'gaussian')) then
        call errore('pw2wannier90', &
             'Can not recognize the choice for scdm_entanglement. ' &
             //'Valid options are: isolated, erfc and gaussian', 1)
    ENDIF
  ENDIF
  IF (scdm_sigma <= 0._dp) &
    call errore('pw2wannier90','Sigma in the SCDM method must be positive.', 1)
  IF (irr_bz) THEN
     IF (gamma_only) CALL errore('pw2wannier90', "irr_bz and gamma_only are not compatible", 1)
     IF (write_spn) CALL errore('pw2wannier90', "irr_bz and write_spn not implemented", 1)
     IF (write_unk) CALL errore('pw2wannier90', "irr_bz and write_unk not implemented", 1)
     IF (write_uHu) CALL errore('pw2wannier90', "irr_bz and write_uHu not implemented", 1)
     IF (write_uIu) CALL errore('pw2wannier90', "irr_bz and write_uIu not implemented", 1)
     IF (write_sHu) CALL errore('pw2wannier90', "irr_bz and write_sHu not implemented", 1)
     IF (write_sIu) CALL errore('pw2wannier90', "irr_bz and write_sIu not implemented", 1)
     IF (write_dmn) CALL errore('pw2wannier90', "irr_bz and write_dmn not implemented", 1)
     IF (scdm_proj) CALL errore('pw2wannier90', "irr_bz and SCDM not implemented", 1)
     IF (write_unkg) CALL errore('pw2wannier90', "irr_bz and write_unkg not implemented", 1)
  ENDIF
  IF (atom_proj) then
    IF (atom_proj_ext .and. noncolin) CALL errore('pw2wannier90', &
         "atom_proj_ext and noncolin not implemented", 1)
  ENDIF
  !
  SELECT CASE ( trim( spin_component ) )
  CASE ( 'up' )
     WRITE(stdout,*) ' Spin CASE ( up )'
     ispinw  = 1
     ikstart = 1
     ikstop  = nkstot/2
     iknum   = nkstot/2
  CASE ( 'down' )
     WRITE(stdout,*) ' Spin CASE ( down )'
     ispinw = 2
     ikstart = nkstot/2 + 1
     ikstop  = nkstot
     iknum   = nkstot/2
  CASE DEFAULT
     IF(noncolin) THEN
        WRITE(stdout,*) ' Spin CASE ( non-collinear )'
     ELSE
        WRITE(stdout,*) ' Spin CASE ( default = unpolarized )'
     ENDIF
     ispinw = 0
     ikstart = 1
     ikstop  = nkstot
     iknum   = nkstot
  END SELECT
  !
  ! Setup for pool parallelization
  ALLOCATE(xk_all(3, nkstot), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating xk_all', 1)
  CALL poolcollect(3, nks, xk, nkstot, xk_all)
  !
  CALL stop_clock( 'init_pw2wan' )
  !
  WRITE(stdout,*)
  WRITE(stdout,*) ' Wannier mode is: ',wan_mode
  WRITE(stdout,*)
  !
  IF(wan_mode=='standalone') THEN
     !
     WRITE(stdout,*) ' -----------------'
     WRITE(stdout,*) ' *** Reading nnkp '
     WRITE(stdout,*) ' -----------------'
     WRITE(stdout,*)
     CALL read_nnkp
     WRITE(stdout,*) ' Opening pp-files '
     CALL openfil_pp
     CALL ylm_expansion
     WRITE(stdout,*)
     WRITE(stdout,*)
     if(write_dmn)then
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*) ' *** Compute DMN '
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*)
        CALL compute_dmn !YN:
        WRITE(stdout,*)
     end if
     IF (write_amn) THEN
        IF (scdm_proj) THEN
            WRITE(stdout,*) ' --------------------------'
            WRITE(stdout,*) ' *** Compute  A with SCDM-k'
            WRITE(stdout,*) ' --------------------------'
            WRITE(stdout,*)
            CALL compute_amn_with_scdm
         ELSEIF (atom_proj) THEN
            WRITE(stdout,*) ' -------------------------------------'
            WRITE(stdout,*) ' *** Compute  A with atomic projectors'
            WRITE(stdout,*) ' -------------------------------------'
            WRITE(stdout,*)
            CALL compute_amn_with_atomproj
        ELSE
            WRITE(stdout,*) ' --------------------------'
            WRITE(stdout,*) ' *** Compute  A projections'
            WRITE(stdout,*) ' --------------------------'
            WRITE(stdout,*)
            CALL compute_amn
        ENDIF
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** A matrix is not computed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
     IF(write_mmn) THEN
        WRITE(stdout,*) ' ---------------'
        WRITE(stdout,*) ' *** Compute  M '
        WRITE(stdout,*) ' ---------------'
        WRITE(stdout,*)
        IF(irr_bz) THEN
           CALL compute_mmn_ibz
        ELSE
           CALL compute_mmn
        ENDIF
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** M matrix is not computed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
     if(noncolin) then
        IF(write_spn) THEN
           WRITE(stdout,*) ' ------------------'
           WRITE(stdout,*) ' *** Compute  Spin '
           WRITE(stdout,*) ' ------------------'
           WRITE(stdout,*)
           CALL compute_spin
           WRITE(stdout,*)
        ELSE
           WRITE(stdout,*) ' --------------------------------'
           WRITE(stdout,*) ' *** Spin matrix is not computed '
           WRITE(stdout,*) ' --------------------------------'
           WRITE(stdout,*)
        ENDIF
     elseif(write_spn) then
        write(stdout,*) ' -----------------------------------'
        write(stdout,*) ' *** Non-collinear calculation is   '
        write(stdout,*) '     required for spin              '
        write(stdout,*) '     term  to be computed           '
        write(stdout,*) ' -----------------------------------'
     endif
     IF(write_uHu.or.write_uIu) THEN
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*) ' *** Compute Orb '
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*)
        CALL compute_orb
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------------'
        WRITE(stdout,*) ' *** Orbital terms are not computed '
        WRITE(stdout,*) ' -----------------------------------'
        WRITE(stdout,*)
     ENDIF
     IF(write_sHu.or.write_sIu) THEN
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*) ' *** Compute shc '
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*)
        CALL compute_shc
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------------'
        WRITE(stdout,*) ' *** SHC terms are not computed '
        WRITE(stdout,*) ' -----------------------------------'
        WRITE(stdout,*)
     ENDIF
     IF(write_eig) THEN
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*) ' *** Write bands '
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*)
        CALL write_band
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' --------------------------'
        WRITE(stdout,*) ' *** Bands are not written '
        WRITE(stdout,*) ' --------------------------'
        WRITE(stdout,*)
     ENDIF
     IF(write_unk) THEN
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*) ' *** Write plot info '
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*)
        CALL write_plot
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** Plot info is not printed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
     IF(write_unkg) THEN
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*) ' *** Write parity info '
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*)
        CALL write_parity
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** Parity info is not printed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
     WRITE(stdout,*) ' ------------'
     WRITE(stdout,*) ' *** Stop pp '
     WRITE(stdout,*) ' ------------'
     WRITE(stdout,*)
     !
     WRITE(stdout, *)
     CALL print_clock('init_pw2wan')
     CALL print_clock('compute_dmn')
     CALL print_clock('compute_amn')
     CALL print_clock('compute_mmn')
     CALL print_clock('compute_spin')
     CALL print_clock('compute_immn')
     CALL print_clock('compute_shc')
     CALL print_clock('compute_orb')
     CALL print_clock('write_unk')
     CALL print_clock('write_parity')
     !
     WRITE(stdout, '(/5x, "Internal routines:")')
     CALL print_clock('scdm_QRCP')
     CALL print_clock('compute_u_kb')
     CALL print_clock('h_psi')
     !
     CALL mp_barrier(world_comm)
     !
     ! not sure if this should be called also in 'library' mode or not !!
     CALL environment_end ( 'PW2WANNIER' )
     IF ( ionode ) WRITE( stdout, *  )
     CALL stop_pp
     !
  ENDIF
  !
  IF(wan_mode=='library') THEN
     !
!     seedname='wannier'
     WRITE(stdout,*) ' Setting up...'
     CALL setup_nnkp
     WRITE(stdout,*)
     WRITE(stdout,*) ' Opening pp-files '
     CALL openfil_pp
     WRITE(stdout,*)
     WRITE(stdout,*) ' Ylm expansion'
     CALL ylm_expansion
     WRITE(stdout,*)
     CALL compute_amn
     CALL compute_mmn
     if(noncolin) then
        IF(write_spn) THEN
           CALL compute_spin
        ENDIF
     ENDIF
     IF(write_uHu.or.write_uIu) THEN
        CALL compute_orb
     ENDIF
     CALL write_band
     IF(write_unk) CALL write_plot
     IF(write_unkg) THEN
        CALL write_parity
     ENDIF
     CALL run_wannier
     CALL lib_dealloc
     CALL stop_pp
     !
  ENDIF
  !
  IF(wan_mode=='wannier2sic') THEN
     !
     CALL read_nnkp
     CALL wan2sic
     !
  ENDIF
  !
  STOP
END PROGRAM pw2wannier90
!
!-----------------------------------------------------------------------
SUBROUTINE lib_dealloc
  !-----------------------------------------------------------------------
  !
  USE wannier

  IMPLICIT NONE

  DEALLOCATE(m_mat,u_mat,u_mat_opt,a_mat,eigval)

  RETURN
END SUBROUTINE lib_dealloc
!
!-----------------------------------------------------------------------
SUBROUTINE setup_nnkp
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps6, tpi, bohr => BOHR_RADIUS_ANGS
  USE cell_base, ONLY : at, bg, alat
  USE gvect,     ONLY : g, gg
  USE ions_base, ONLY : nat, tau, ityp, atm
  USE klist,     ONLY : xk
  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_pools,  ONLY : intra_pool_comm
  USE mp_world,  ONLY : world_comm
  USE wvfct,     ONLY : nbnd,npwx
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin
  USE wannier

  IMPLICIT NONE
  real(DP) :: g_(3), gg_
  INTEGER  :: ik, ib, ig, iw, ia, indexb, TYPE, ierr
  INTEGER, ALLOCATABLE :: ig_check(:,:)
  real(DP) :: xnorm, znorm, coseno
  INTEGER  :: exclude_bands(nbnd)

  ! aam: translations between PW2Wannier90 and Wannier90
  ! pw2wannier90   <==>   Wannier90
  !    nbnd                num_bands_tot
  !    n_wannier           num_wann
  !    num_bands           num_bands
  !    nat                 num_atoms
  !    iknum               num_kpts
  !    rlatt               transpose(real_lattice)
  !    glatt               transpose(recip_lattice)
  !    kpt_latt            kpt_latt
  !    nnb                 nntot
  !    kpb                 nnlist
  !    g_kpb               nncell
  !    mp_grid             mp_grid
  !    center_w            proj_site
  !    l_w,mr_w,r_w        proj_l,proj_m,proj_radial
  !    xaxis,zaxis         proj_x,proj_z
  !    alpha_w             proj_zona
  !    exclude_bands       exclude_bands
  !    atcart              atoms_cart
  !    atsym               atom_symbols

  ALLOCATE( kpt_latt(3,iknum), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating kpt_latt', 1)
  ALLOCATE( atcart(3,nat), atsym(nat), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating atcart/atsym', 1)
  ALLOCATE( kpb(iknum,num_nnmax), g_kpb(3,iknum,num_nnmax), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating kpb/g_kpb', 1)
  ALLOCATE( center_w(3,nbnd), alpha_w(nbnd), l_w(nbnd), &
       mr_w(nbnd), r_w(nbnd), zaxis(3,nbnd), xaxis(3,nbnd), stat=ierr)
       IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating center_w/alpha_w/l_w/...', 1)
  ALLOCATE( excluded_band(nbnd), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating excluded_band', 1)

  ! real lattice (Cartesians, Angstrom)
  rlatt(:,:) = transpose(at(:,:))*alat*bohr
  ! reciprocal lattice (Cartesians, Angstrom)
  glatt(:,:) = transpose(bg(:,:))*tpi/(alat*bohr)
  ! convert Cartesian k-points to crystallographic co-ordinates
  kpt_latt(:,1:iknum)=xk(:,1:iknum)
  CALL cryst_to_cart(iknum,kpt_latt,at,-1)
  ! atom co-ordinates in Cartesian co-ords and Angstrom units
  atcart(:,:) = tau(:,:)*bohr*alat
  ! atom symbols
  DO ia=1,nat
     TYPE=ityp(ia)
     atsym(ia)=atm(TYPE)
  ENDDO

  ! MP grid dimensions
  CALL find_mp_grid()

  WRITE(stdout,'("  - Number of atoms is (",i3,")")') nat

#if defined(__WANLIB)
  IF (ionode) THEN
     CALL wannier_setup(seedname,mp_grid,iknum,rlatt, &               ! input
          glatt,kpt_latt,nbnd,nat,atsym,atcart,gamma_only,noncolin, & ! input
          nnb,kpb,g_kpb,num_bands,n_wannier,center_w, &               ! output
          l_w,mr_w,r_w,zaxis,xaxis,alpha_w,exclude_bands)             ! output
  ENDIF
#endif

  CALL mp_bcast(nnb,ionode_id, world_comm)
  CALL mp_bcast(kpb,ionode_id, world_comm)
  CALL mp_bcast(g_kpb,ionode_id, world_comm)
  CALL mp_bcast(num_bands,ionode_id, world_comm)
  CALL mp_bcast(n_wannier,ionode_id, world_comm)
  CALL mp_bcast(center_w,ionode_id, world_comm)
  CALL mp_bcast(l_w,ionode_id, world_comm)
  CALL mp_bcast(mr_w,ionode_id, world_comm)
  CALL mp_bcast(r_w,ionode_id, world_comm)
  CALL mp_bcast(zaxis,ionode_id, world_comm)
  CALL mp_bcast(xaxis,ionode_id, world_comm)
  CALL mp_bcast(alpha_w,ionode_id, world_comm)
  CALL mp_bcast(exclude_bands,ionode_id, world_comm)

  IF(noncolin) THEN
     n_proj=n_wannier/2
  ELSE
     n_proj=n_wannier
  ENDIF

  ALLOCATE( gf(npwx,n_proj), csph(16,n_proj), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating gf/csph', 1)

  WRITE(stdout,'("  - Number of wannier functions is (",i3,")")') n_wannier

  excluded_band(1:nbnd)=.false.
  nexband=0
  band_loop: DO ib=1,nbnd
     indexb=exclude_bands(ib)
     IF (indexb>nbnd .or. indexb<0) THEN
        CALL errore('setup_nnkp',' wrong excluded band index ', 1)
     ELSEIF (indexb==0) THEN
        exit band_loop
     ELSE
        nexband=nexband+1
        excluded_band(indexb)=.true.
     ENDIF
  ENDDO band_loop

  IF ( (nbnd-nexband)/=num_bands ) &
       CALL errore('setup_nnkp',' something wrong with num_bands',1)

  DO iw=1,n_proj
     xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
          xaxis(3,iw)*xaxis(3,iw))
     IF (xnorm < eps6) CALL errore ('setup_nnkp',' |xaxis| < eps ',1)
     znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
          zaxis(3,iw)*zaxis(3,iw))
     IF (znorm < eps6) CALL errore ('setup_nnkp',' |zaxis| < eps ',1)
     coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
          xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
     IF (abs(coseno) > eps6) &
          CALL errore('setup_nnkp',' xaxis and zaxis are not orthogonal !',1)
     IF (alpha_w(iw) < eps6) &
          CALL errore('setup_nnkp',' zona value must be positive', 1)
     ! convert wannier center in cartesian coordinates (in unit of alat)
     CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
  ENDDO
  WRITE(stdout,*) ' - All guiding functions are given '

  nnbx=0
  nnb=max(nnbx,nnb)

  ALLOCATE( ig_(iknum,nnb), ig_check(iknum,nnb), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ig_/ig_check', 1)
  ALLOCATE( zerophase(iknum,nnb), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating zerophase', 1)
  zerophase = .false.

  DO ik=1, iknum
     DO ib = 1, nnb
        IF ( (g_kpb(1,ik,ib).eq.0) .and.  &
             (g_kpb(2,ik,ib).eq.0) .and.  &
             (g_kpb(3,ik,ib).eq.0) ) zerophase(ik,ib) = .true.
        g_(:) = REAL( g_kpb(:,ik,ib) )
        CALL cryst_to_cart (1, g_, bg, 1)
        gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
        ig_(ik,ib) = 0
        ig = 1
        DO WHILE  (gg(ig) <= gg_ + eps6)
           IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
           ig= ig +1
        ENDDO
     ENDDO
  ENDDO

  ig_check(:,:) = ig_(:,:)
  CALL mp_sum( ig_check, intra_pool_comm )
  DO ik=1, iknum
     DO ib = 1, nnb
        IF (ig_check(ik,ib) ==0) &
          CALL errore('setup_nnkp', &
                      ' g_kpb vector is not in the list of Gs', 100*ik+ib )
     ENDDO
  ENDDO
  DEALLOCATE (ig_check)

  WRITE(stdout,*) ' - All neighbours are found '
  WRITE(stdout,*)

  RETURN
END SUBROUTINE setup_nnkp
 !
 !-----------------------------------------------------------------------
SUBROUTINE run_wannier
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE control_flags, ONLY : gamma_only
  USE wannier

  IMPLICIT NONE

  INTEGER :: ierr

  ALLOCATE(u_mat(n_wannier,n_wannier,iknum), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating u_mat', 1)
  ALLOCATE(u_mat_opt(num_bands,n_wannier,iknum), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating u_mat_opt', 1)
  ALLOCATE(lwindow(num_bands,iknum), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating lwindow', 1)
  ALLOCATE(wann_centers(3,n_wannier), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wann_centers', 1)
  ALLOCATE(wann_spreads(n_wannier), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wann_spreads', 1)

#if defined(__WANLIB)
  IF (ionode) THEN
     CALL wannier_run(seedname,mp_grid,iknum,rlatt, &                ! input
          glatt,kpt_latt,num_bands,n_wannier,nnb,nat, &              ! input
          atsym,atcart,gamma_only,m_mat,a_mat,eigval, &              ! input
          u_mat,u_mat_opt,lwindow,wann_centers,wann_spreads,spreads) ! output
  ENDIF
#endif

  CALL mp_bcast(u_mat,ionode_id, world_comm)
  CALL mp_bcast(u_mat_opt,ionode_id, world_comm)
  CALL mp_bcast(lwindow,ionode_id, world_comm)
  CALL mp_bcast(wann_centers,ionode_id, world_comm)
  CALL mp_bcast(wann_spreads,ionode_id, world_comm)
  CALL mp_bcast(spreads,ionode_id, world_comm)

  RETURN
END SUBROUTINE run_wannier
!-----------------------------------------------------------------------
!
SUBROUTINE find_mp_grid()
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY: DP
  USE wannier

  IMPLICIT NONE

  ! <<<local variables>>>
  INTEGER  :: ik,ntemp,ii
  real(DP) :: min_k,temp(3,iknum),mpg1

  min_k=minval(kpt_latt(1,:))
  ii=0
  DO ik=1,iknum
     IF (kpt_latt(1,ik)==min_k) THEN
        ii=ii+1
        temp(:,ii)=kpt_latt(:,ik)
     ENDIF
  ENDDO
  ntemp=ii

  min_k=minval(temp(2,1:ntemp))
  ii=0
  DO ik=1,ntemp
     IF (temp(2,ik)==min_k) THEN
        ii=ii+1
     ENDIF
  ENDDO
  mp_grid(3)=ii

  min_k=minval(temp(3,1:ntemp))
  ii=0
  DO ik=1,ntemp
     IF (temp(3,ik)==min_k) THEN
        ii=ii+1
     ENDIF
  ENDDO
  mp_grid(2)=ii

  IF ( (mp_grid(2)==0) .or. (mp_grid(3)==0) ) &
       CALL errore('find_mp_grid',' one or more mp_grid dimensions is zero', 1)

  mpg1=iknum/(mp_grid(2)*mp_grid(3))

  mp_grid(1) = nint(mpg1)

  WRITE(stdout,*)
  WRITE(stdout,'(3(a,i3))') '  MP grid is ',mp_grid(1),' x',mp_grid(2),' x',mp_grid(3)

  IF (real(mp_grid(1),kind=DP)/=mpg1) &
       CALL errore('find_mp_grid',' determining mp_grid failed', 1)

  RETURN
END SUBROUTINE find_mp_grid
!-----------------------------------------------------------------------
!
SUBROUTINE read_nnkp
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE kinds,     ONLY: DP
  USE constants, ONLY : eps6, tpi, bohr => BOHR_RADIUS_ANGS
  USE cell_base, ONLY : at, bg, alat
  USE gvect,     ONLY : g, gg
  USE klist,     ONLY : nkstot, xk
  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_pools,  ONLY : intra_pool_comm
  USE mp_world,  ONLY : world_comm
  USE wvfct,     ONLY : npwx, nbnd
  USE noncollin_module, ONLY : noncolin
  USE wannier
  USE atproj,    ONLY : atom_proj

  IMPLICIT NONE
  !
  real(DP) :: g_(3), gg_
  INTEGER :: ik, ib, ig, ipol, iw, idum, indexb
  INTEGER :: numk, i, j, ierr
  INTEGER, ALLOCATABLE :: ig_check(:,:)
  real(DP) :: xx(3), xnorm, znorm, coseno
  LOGICAL :: have_nnkp,found
  INTEGER :: tmp_auto ! vv: Needed for the selection of projections with SCDM
  REAL(DP), ALLOCATABLE :: xkc_full(:,:)

  IF (ionode) THEN  ! Read nnkp file on ionode only

     INQUIRE(file=trim(seedname)//".nnkp",exist=have_nnkp)
     IF(.not. have_nnkp) THEN
        CALL errore( 'pw2wannier90', 'Could not find the file '&
           &//trim(seedname)//'.nnkp', 1 )
     ENDIF

     OPEN (NEWUNIT=iun_nnkp, file=trim(seedname)//".nnkp",form='formatted', status="old")

  ENDIF

  nnbx=0

  !   check the information from *.nnkp with the nscf_save data
  WRITE(stdout,*) ' Checking info from wannier.nnkp file'
  WRITE(stdout,*)

  IF (ionode) THEN   ! read from ionode only

     CALL scan_file_to('real_lattice',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find real_lattice block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     DO j=1,3
        READ(iun_nnkp,*) (rlatt(i,j),i=1,3)
        DO i = 1,3
           rlatt(i,j) = rlatt(i,j)/(alat*bohr)
        ENDDO
     ENDDO
     DO j=1,3
        DO i=1,3
           IF(abs(rlatt(i,j)-at(i,j))>eps6) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*)  ' rlatt(i,j) =',rlatt(i,j),  ' at(i,j)=',at(i,j)
              CALL errore( 'pw2wannier90', 'Direct lattice mismatch', 3*j+i )
           ENDIF
        ENDDO
     ENDDO
     WRITE(stdout,*) ' - Real lattice is ok'

     CALL scan_file_to('recip_lattice',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find recip_lattice block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     DO j=1,3
        READ(iun_nnkp,*) (glatt(i,j),i=1,3)
        DO i = 1,3
           glatt(i,j) = (alat*bohr)*glatt(i,j)/tpi
        ENDDO
     ENDDO
     DO j=1,3
        DO i=1,3
           IF(abs(glatt(i,j)-bg(i,j))>eps6) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*)  ' glatt(i,j)=',glatt(i,j), ' bg(i,j)=',bg(i,j)
              CALL errore( 'pw2wannier90', 'Reciprocal lattice mismatch', 3*j+i )
           ENDIF
        ENDDO
     ENDDO
     WRITE(stdout,*) ' - Reciprocal lattice is ok'

     CALL scan_file_to('kpoints',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find kpoints block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     READ(iun_nnkp,*) numk
     IF(irr_bz) THEN
        ALLOCATE(xkc_full(3,numk), stat=ierr)
        IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating xkc_full', 1)
        DO i=1,numk
           READ(iun_nnkp,*) xkc_full(:,i)
        END DO
        !IF(any(abs(xkc_full(:,:)) > 0.5)) THEN
        !   CALL errore( 'pw2wannier90', 'kpoints should be -0.5 <= k < 0.5', 1 )
        !ENDIF
     ELSE
        IF(numk/=iknum) THEN
           WRITE(stdout,*)  ' Something wrong! '
           WRITE(stdout,*)  ' numk=',numk, ' iknum=',iknum
           CALL errore( 'pw2wannier90', 'Wrong number of k-points', numk)
        ENDIF
        ! Check that the k-points agree with those in the nnkp file.
        ! We do this check only at the ionode which has all the k points even
        ! with pool parallelization, so calling poolcollect is not necessary.
        DO i=1,numk
           READ(iun_nnkp,*) xx(1), xx(2), xx(3)
           CALL cryst_to_cart( 1, xx, bg, 1 )
           IF (ANY(ABS(xx - xk(:,i))>eps6)) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*) ' k-point ',i,' is wrong'
              WRITE(stdout,*) xx(1), xx(2), xx(3)
              WRITE(stdout,*) xk(1,i), xk(2,i), xk(3,i)
              CALL errore( 'pw2wannier90', 'problems with k-points', i )
           ENDIF
        ENDDO
     ENDIF
     WRITE(stdout,*) ' - K-points are ok'

  ENDIF ! ionode

  ! Broadcast
  CALL mp_bcast(rlatt,ionode_id, world_comm)
  CALL mp_bcast(glatt,ionode_id, world_comm)

  IF (ionode) THEN   ! read from ionode only
     if(noncolin) then
        old_spinor_proj=.false.
        CALL scan_file_to('spinor_projections',found)
        if(.not.found) then
           !try old style projections
           CALL scan_file_to('projections',found)
           if(found) then
              old_spinor_proj=.true.
           else
              CALL errore( 'pw2wannier90', 'Could not find projections block in '&
                 &//trim(seedname)//'.nnkp', 1 )
           endif
        end if
     else
        old_spinor_proj=.false.
        CALL scan_file_to('projections',found)
        if(.not.found) then
           CALL errore( 'pw2wannier90', 'Could not find projections block in '&
              &//trim(seedname)//'.nnkp', 1 )
        endif
     endif
     READ(iun_nnkp,*) n_proj
  ENDIF

  ! Broadcast
  CALL mp_bcast(n_proj,ionode_id, world_comm)
  CALL mp_bcast(old_spinor_proj,ionode_id, world_comm)

  IF(old_spinor_proj)THEN
     WRITE(stdout,'(//," ****** begin WARNING ****** ",/)')
     WRITE(stdout,'(" The pw.x calculation was done with non-collinear spin ")')
     WRITE(stdout,'(" but spinor = T was not specified in the wannier90 .win file!")')
     WRITE(stdout,'(" Please set spinor = T and rerun wannier90.x -pp  ")')
     WRITE(stdout,'(/," ******  end WARNING  ****** ",//)')
     CALL errore("pw2wannier90", "Spinorbit without spinor=T",1)
  ENDIF

  ! It is not clear if the next instruction is required or not, it probably depend
  ! on the version of wannier90 that was used to generate the nnkp file:
  n_wannier = n_proj

  ALLOCATE( center_w(3,n_proj), alpha_w(n_proj), gf(npwx,n_proj), &
       l_w(n_proj), mr_w(n_proj), r_w(n_proj), &
       zaxis(3,n_proj), xaxis(3,n_proj), csph(16,n_proj), stat=ierr)
       IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating center_w/alpha_w/...', 1)
  if (noncolin) then
     ALLOCATE( spin_eig(n_proj),spin_qaxis(3,n_proj), stat=ierr)
     IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating spin_eig/spin_qaxis', 1)
  endif

  IF (ionode) THEN   ! read from ionode only
     DO iw=1,n_proj
        READ(iun_nnkp,*) (center_w(i,iw), i=1,3), l_w(iw), mr_w(iw), r_w(iw)
        READ(iun_nnkp,*) (zaxis(i,iw),i=1,3),(xaxis(i,iw),i=1,3),alpha_w(iw)
        xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
             xaxis(3,iw)*xaxis(3,iw))
        IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
        znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
             zaxis(3,iw)*zaxis(3,iw))
        IF (znorm < eps6) CALL errore ('read_nnkp',' |zaxis| < eps ',1)
        coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
             xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
        IF (abs(coseno) > eps6) &
             CALL errore('read_nnkp',' xaxis and zaxis are not orthogonal !',1)
        IF (alpha_w(iw) < eps6) &
             CALL errore('read_nnkp',' zona value must be positive', 1)
        ! convert wannier center in cartesian coordinates (in unit of alat)
        CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
        if (noncolin) then
           READ(iun_nnkp,*) spin_eig(iw),(spin_qaxis(i,iw),i=1,3)
           xnorm = sqrt(spin_qaxis(1,iw)*spin_qaxis(1,iw) + spin_qaxis(2,iw)*spin_qaxis(2,iw) + &
             spin_qaxis(3,iw)*spin_qaxis(3,iw))
           IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
           spin_qaxis(:,iw)=spin_qaxis(:,iw)/xnorm
        endif
     ENDDO
  ENDIF

  ! automatic projections
  IF (ionode) THEN
     CALL scan_file_to('auto_projections',found)
     IF (found) THEN
        READ (iun_nnkp, *) n_wannier
        READ (iun_nnkp, *) tmp_auto

        IF (scdm_proj) THEN
           IF (n_proj > 0) THEN
              WRITE(stdout,'(//, " ****** begin Error message ******",/)')
              WRITE(stdout,'(/," Found a projection block, an auto_projections block",/)')
              WRITE(stdout,'(/," and scdm_proj = T in the input file. These three options are inconsistent.",/)')
              WRITE(stdout,'(/," Please refer to the Wannier90 User guide for correct use of these flags.",/)')
              WRITE(stdout,'(/, " ****** end Error message ******",//)')
              CALL errore( 'pw2wannier90', 'Inconsistent options for projections.', 1 )
           ELSE
              IF (tmp_auto /= 0) CALL errore( 'pw2wannier90', 'Second entry in auto_projections block is not 0. ' // &
              'See Wannier90 User Guide in the auto_projections section for clarifications.', 1 )
           ENDIF
        ELSE IF (atom_proj) THEN
           continue
        ELSE
           ! Fire an error whether or not a projections block is found
           CALL errore( 'pw2wannier90', 'scdm_proj = F and atom_proj = F '&
                &'but found an auto_projections block in '&
                &//trim(seedname)//'.nnkp', 1 )
        ENDIF
     ELSE
        IF (scdm_proj) THEN
           ! Fire an error whether or not a projections block is found
           CALL errore( 'pw2wannier90', 'scdm_proj = T but cannot find an auto_projections block in '&
                &//trim(seedname)//'.nnkp', 1 )
        ENDIF
        IF (atom_proj) THEN
          CALL errore( 'pw2wannier90', 'atom_proj = T but cannot find an auto_projections block in '&
                      &//trim(seedname)//'.nnkp', 1 )
        ENDIF
     ENDIF
  ENDIF

  ! Broadcast
  CALL mp_bcast(n_wannier,ionode_id, world_comm)
  CALL mp_bcast(center_w,ionode_id, world_comm)
  CALL mp_bcast(l_w,ionode_id, world_comm)
  CALL mp_bcast(mr_w,ionode_id, world_comm)
  CALL mp_bcast(r_w,ionode_id, world_comm)
  CALL mp_bcast(zaxis,ionode_id, world_comm)
  CALL mp_bcast(xaxis,ionode_id, world_comm)
  CALL mp_bcast(alpha_w,ionode_id, world_comm)
  if (noncolin) then
     CALL mp_bcast(spin_eig,ionode_id, world_comm)
     CALL mp_bcast(spin_qaxis,ionode_id, world_comm)
  end if

  WRITE(stdout,'("  - Number of wannier functions is ok (",i3,")")') n_wannier

  IF (.not. scdm_proj) WRITE(stdout,*) ' - All guiding functions are given '
  !
  WRITE(stdout,*)
  WRITE(stdout,*) ' Projections:'
  DO iw=1,n_proj
     WRITE(stdout,'(3f12.6,3i3,f12.6)') &
          center_w(1:3,iw),l_w(iw),mr_w(iw),r_w(iw),alpha_w(iw)
  ENDDO

  IF (ionode) THEN   ! read from ionode only
     CALL scan_file_to('nnkpts',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find nnkpts block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     READ (iun_nnkp,*) nnb
  ENDIF

  ! Broadcast
  CALL mp_bcast(nnb,ionode_id, world_comm)
  !
  nnbx = max (nnbx, nnb )
  !
  ALLOCATE ( kpb(iknum,nnbx), g_kpb(3,iknum,nnbx),&
             ig_(iknum,nnbx), ig_check(iknum,nnbx), stat=ierr)
             IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating kpb/g_kpb/...', 1)
  ALLOCATE( zerophase(iknum,nnbx), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating zerophase', 1)
  zerophase = .false.

  !  read data about neighbours
  WRITE(stdout,*)
  WRITE(stdout,*) ' Reading data about k-point neighbours '
  WRITE(stdout,*)

  IF (ionode) THEN
     DO ik=1, iknum
        DO ib = 1, nnb
           READ(iun_nnkp,*) idum, kpb(ik,ib), (g_kpb(ipol,ik,ib), ipol =1,3)
        ENDDO
     ENDDO
  ENDIF

  ! Broadcast
  CALL mp_bcast(kpb,ionode_id, world_comm)
  CALL mp_bcast(g_kpb,ionode_id, world_comm)

  DO ik=1, iknum
     DO ib = 1, nnb
        IF ( (g_kpb(1,ik,ib).eq.0) .and.  &
             (g_kpb(2,ik,ib).eq.0) .and.  &
             (g_kpb(3,ik,ib).eq.0) ) zerophase(ik,ib) = .true.
        g_(:) = REAL( g_kpb(:,ik,ib), KIND=DP )
        CALL cryst_to_cart (1, g_, bg, 1)
        gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
        ig_(ik,ib) = 0
        ig = 1
        DO WHILE  (gg(ig) <= gg_ + eps6)
           IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
           ig= ig +1
        ENDDO
     ENDDO
  ENDDO
  ig_check(:,:) = ig_(:,:)
  CALL mp_sum( ig_check, intra_pool_comm )
  DO ik=1, iknum
     DO ib = 1, nnb
        IF (ig_check(ik,ib) ==0) &
          CALL errore('read_nnkp', &
                      ' g_kpb vector is not in the list of Gs', 100*ik+ib )
     ENDDO
  ENDDO
  DEALLOCATE (ig_check)

  WRITE(stdout,*) ' All neighbours are found '
  WRITE(stdout,*)

  ALLOCATE( excluded_band(nbnd), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating excluded_band', 1)

  IF (ionode) THEN     ! read from ionode only
     CALL scan_file_to('exclude_bands',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find exclude_bands block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     READ (iun_nnkp,*) nexband
     excluded_band(1:nbnd)=.false.
     DO i=1,nexband
        READ(iun_nnkp,*) indexb
        IF (indexb<1 .or. indexb>nbnd) &
             CALL errore('read_nnkp',' wrong excluded band index ', 1)
        excluded_band(indexb)=.true.
     ENDDO
  ENDIF
  num_bands=nbnd-nexband

  ! Broadcast
  CALL mp_bcast(nexband,ionode_id, world_comm)
  CALL mp_bcast(excluded_band,ionode_id, world_comm)
  CALL mp_bcast(num_bands,ionode_id, world_comm)

  !
  ALLOCATE( bvec(3,nnb), xbvec(3,nnb), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating bvec/xbvec', 1)
  IF (ionode) THEN
     xbvec = 0
     IF (irr_bz) THEN
        DO i=1, nnb
           xbvec(:,i) = xkc_full(:,kpb(1,i)) - xkc_full(:,1) + g_kpb(:,1,i)
        END DO
        DEALLOCATE(xkc_full)
     ENDIF
  ENDIF
  CALL mp_bcast(xbvec, ionode_id, world_comm)
  bvec = xbvec
  CALL cryst_to_cart(nnb, bvec, bg, +1)

  IF (ionode) CLOSE (iun_nnkp)   ! ionode only

  RETURN
END SUBROUTINE read_nnkp
!
!-----------------------------------------------------------------------
SUBROUTINE scan_file_to (keyword,found)
   !-----------------------------------------------------------------------
   !
   USE wannier, ONLY :iun_nnkp
   USE io_global,  ONLY : stdout
   IMPLICIT NONE
   CHARACTER(len=*), intent(in) :: keyword
   logical, intent(out) :: found
   CHARACTER(len=80) :: line1, line2
!
! by uncommenting the following line the file scan restarts every time
! from the beginning thus making the reading independent on the order
! of data-blocks
!   rewind (iun_nnkp)
!
10 CONTINUE
   READ(iun_nnkp,*,end=20) line1, line2
   IF(line1/='begin')  GOTO 10
   IF(line2/=keyword) GOTO 10
   found=.true.
   RETURN
20 found=.false.
   rewind (iun_nnkp)

END SUBROUTINE scan_file_to
!
!-----------------------------------------------------------------------
SUBROUTINE pw2wan_set_symm (nsym, sr, tvec)
   !-----------------------------------------------------------------------
   !
   ! Uses nkqs and index_sym from module pw2wan, computes rir
   !
   USE symm_base,       ONLY : s, ft, allfrac
   USE fft_base,        ONLY : dffts
   USE cell_base,       ONLY : at, bg
   USE wannier,         ONLY : rir, read_sym
   USE kinds,           ONLY : DP
   USE io_global,       ONLY : stdout
   !
   IMPLICIT NONE
   !
   INTEGER  , intent(in) :: nsym
   REAL(DP) , intent(in) :: sr(3,3,nsym), tvec(3,nsym)
   REAL(DP) :: st(3,3), v(3)
   INTEGER, allocatable :: s_in(:,:,:)
   REAL(DP), allocatable:: ft_in(:,:)
   INTEGER :: nxxs, nr1,nr2,nr3, nr1x,nr2x,nr3x
   INTEGER :: ikq, isym, i,j,k, ri,rj,rk, ir, ierr
   LOGICAL :: ispresent(nsym)
   !
   nr1 = dffts%nr1
   nr2 = dffts%nr2
   nr3 = dffts%nr3
   nr1x= dffts%nr1x
   nr2x= dffts%nr2x
   nr3x= dffts%nr3x
   nxxs = nr1x*nr2x*nr3x
   !
   !  sr -> s
   ALLOCATE(s_in(3,3,nsym), ft_in(3,nsym), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating s_in/ft_in', 1)
   IF(read_sym ) THEN
      IF(allfrac) THEN
         call errore("pw2wan_set_symm", "use_all_frac = .true. + read_sym = .true. not supported", 1)
      END IF
      DO isym = 1, nsym
         !st = transpose( MATMUL(transpose(bg), sr(:,:,isym)) )
         st = transpose( MATMUL(transpose(bg), transpose(sr(:,:,isym))) )
         s_in(:,:,isym) = nint( MATMUL(transpose(at), st) )
         v = MATMUL(transpose(bg), tvec(:,isym))
         ft_in(1,isym) = v(1)
         ft_in(2,isym) = v(2)
         ft_in(3,isym) = v(3)
      END DO
      IF( any(s(:,:,1:nsym) /= s_in(:,:,1:nsym)) .or. any(ft_in(:,1:nsym) /= ft(:,1:nsym)) ) THEN
         write(stdout,*) " Input symmetry is different from crystal symmetry"
         write(stdout,*)
      END IF
   ELSE
      s_in = s(:,:,1:nsym)
      ft_in = ft(:,1:nsym)
   END IF
   !
   IF(.not. allocated(rir)) ALLOCATE(rir(nxxs,nsym))
   rir = 0
   ispresent(1:nsym) = .false.

   DO isym = 1, nsym
      ! scale sym.ops. with FFT dimensions, check consistency
      ! FIXME: what happens with fractional translations?
      IF ( mod(s_in(2, 1, isym) * nr1, nr2) /= 0 .or. &
           mod(s_in(3, 1, isym) * nr1, nr3) /= 0 .or. &
           mod(s_in(1, 2, isym) * nr2, nr1) /= 0 .or. &
           mod(s_in(3, 2, isym) * nr2, nr3) /= 0 .or. &
           mod(s_in(1, 3, isym) * nr3, nr1) /= 0 .or. &
           mod(s_in(2, 3, isym) * nr3, nr2) /= 0 ) THEN
         CALL errore ('pw2waninit',' smooth grid is not compatible with &
                                   & symmetry: change cutoff',isym)
      ENDIF
      s_in (2,1,isym) = s_in (2,1,isym) * nr1 / nr2
      s_in (3,1,isym) = s_in (3,1,isym) * nr1 / nr3
      s_in (1,2,isym) = s_in (1,2,isym) * nr2 / nr1
      s_in (2,2,isym) = s_in (2,2,isym)
      s_in (3,2,isym) = s_in (3,2,isym) * nr2 / nr3
      s_in (1,3,isym) = s_in (1,3,isym) * nr3 / nr1
      s_in (2,3,isym) = s_in (2,3,isym) * nr3 / nr2
      s_in (3,3,isym) = s_in (3,3,isym)

      DO ir=1, nxxs
         rir(ir,isym) = ir
      ENDDO
      DO k = 1, nr3
         DO j = 1, nr2
            DO i = 1, nr1
               CALL rotate_grid_point (s_in(:,:,isym), (/0,0,0/), i,j,k, &
                    nr1,nr2,nr3, ri,rj,rk)
               !
               ir =   i + ( j-1)*nr1x + ( k-1)*nr1x*nr2x
               rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(s_in, ft_in)
END SUBROUTINE pw2wan_set_symm

!-----------------------------------------------------------------------
SUBROUTINE compute_dmn
   !Calculate d_matrix_wann/band for site-symmetry mode given by Rei Sakuma.
   !Contributions for this subroutine:
   !  Yoshiro Nohara (June to July, 2016)
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout, ionode, ionode_id
   USE kinds,           ONLY: DP
   USE wvfct,           ONLY : nbnd, npwx
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions, ONLY : evc, psic, psic_nc
   USE fft_base,        ONLY : dffts, dfftp
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE klist,           ONLY : nkstot, xk, igk_k, ngk
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g, ngm, gstart
   USE cell_base,       ONLY : omega, alat, tpiba, at, bg
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,       ONLY : tpi, bohr => BOHR_RADIUS_ANGS
   USE uspp,            ONLY : nkb, vkb, okvan
   USE uspp_param,      ONLY : upf, nh, lmaxq, nhm
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE mp_pools,        ONLY : intra_pool_comm, inter_pool_comm, me_pool, root_pool, my_pool_id
   USE mp,              ONLY : mp_sum, mp_bcast, mp_barrier
   USE mp_world,        ONLY : world_comm
   USE noncollin_module,ONLY : noncolin, npol
   USE gvecw,           ONLY : gcutw
   USE wannier
   USE symm_base,       ONLY : nsymin=>nsym,srin=>sr,ftin=>ft,invsin=>invs
   USE fft_base,        ONLY : dffts
   USE scatter_mod,     ONLY : gather_grid, scatter_grid
   USE uspp_init,       ONLY : init_us_2
   IMPLICIT NONE
   !
   complex(DP), parameter :: cmplx_i=(0.0_DP,1.0_DP)
   !
   real(DP), parameter :: p12(3,12)=reshape(                            &
      (/0d0, 0d0, 1.00000000000000d0,                                   &
        0.894427190999916d0, 0d0, 0.447213595499958d0,                  &
        0.276393202250021d0, 0.850650808352040d0, 0.447213595499958d0,  &
       -0.723606797749979d0, 0.525731112119134d0, 0.447213595499958d0,  &
       -0.723606797749979d0, -0.525731112119134d0, 0.447213595499958d0, &
        0.276393202250021d0, -0.850650808352040d0, 0.447213595499958d0, &
        0.723606797749979d0, 0.525731112119134d0, -0.447213595499958d0, &
       -0.276393202250021d0, 0.850650808352040d0, -0.447213595499958d0, &
       -0.894427190999916d0, 0d0, -0.447213595499958d0,                 &
       -0.276393202250021d0, -0.850650808352040d0, -0.447213595499958d0,&
        0.723606797749979d0, -0.525731112119134d0, -0.447213595499958d0,&
        0d0, 0d0, -1.00000000000000d0/),(/3,12/))
   real(DP), parameter :: p20(3,20)=reshape(                            &
      (/0.525731112119134d0, 0.381966011250105d0, 0.850650808352040d0,  &
       -0.200811415886227d0, 0.618033988749895d0, 0.850650808352040d0,  &
       -0.649839392465813d0, 0d0, 0.850650808352040d0,                  &
       -0.200811415886227d0, -0.618033988749895d0, 0.850650808352040d0, &
        0.525731112119134d0, -0.381966011250105d0, 0.850650808352040d0, &
        0.850650808352040d0, 0.618033988749895d0, 0.200811415886227d0,  &
       -0.324919696232906d0, 1.00000000000000d0, 0.200811415886227d0,   &
       -1.05146222423827d0, 0d0, 0.200811415886227d0,                   &
      -0.324919696232906d0, -1.00000000000000d0, 0.200811415886227d0,   &
       0.850650808352040d0, -0.618033988749895d0, 0.200811415886227d0,  &
       0.324919696232906d0, 1.00000000000000d0, -0.200811415886227d0,   &
      -0.850650808352040d0, 0.618033988749895d0, -0.200811415886227d0,  &
      -0.850650808352040d0, -0.618033988749895d0, -0.200811415886227d0, &
       0.324919696232906d0, -1.00000000000000d0, -0.200811415886227d0,  &
       1.05146222423827d0, 0d0, -0.200811415886227d0,                   &
       0.200811415886227d0, 0.618033988749895d0, -0.850650808352040d0,  &
      -0.525731112119134d0, 0.381966011250105d0, -0.850650808352040d0,  &
      -0.525731112119134d0, -0.381966011250105d0, -0.850650808352040d0, &
       0.200811415886227d0, -0.618033988749895d0, -0.850650808352040d0, &
      0.649839392465813d0, 0d0, -0.850650808352040d0/),(/3,20/))
   real(DP), parameter :: pwg(2)=(/2.976190476190479d-2,3.214285714285711d-2/)
   !
   INTEGER :: npw, mmn_tot, ik, ikp, ipol, isym, npwq, i, m, n, ir, jsym
   INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ind, nir
   INTEGER :: ikevc, ikpevcq, s, counter, iun_dmn, iun_sym, ig, igp, ip, jp, np, iw, jw
   INTEGER :: ir_start, ir_end
   INTEGER :: ik_global, isk_global, ipool, ik_local
   COMPLEX(DP), ALLOCATABLE :: phase(:), aux(:), aux2(:), &
                               becp2(:,:), Mkb(:,:), aux_nc(:,:)
   real(DP), ALLOCATABLE    :: rbecp2(:,:),sr(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: qb(:,:,:), phs(:,:)
   COMPLEX(DP)              :: qgm
   REAL(DP)                 :: qg
   real(DP), ALLOCATABLE    :: workg(:)
   real(DP), ALLOCATABLE    :: ylm(:), dxk(:), tvec(:,:), dylm(:,:), wws(:,:,:), vps2t(:,:,:), vaxis(:,:,:)
   INTEGER, ALLOCATABLE     :: iks2k(:,:),iks2g(:,:),ik2ir(:),ir2ik(:)
   INTEGER, ALLOCATABLE     :: iw2ip(:),ip2iw(:),ips2p(:,:),invs(:)
   logical, ALLOCATABLE     :: lfound(:)
   COMPLEX(DP)              :: mmn, zdotc, phase1
   real(DP)                 :: arg, g_(3),v1(3),v2(3),v3(3),v4(3),v5(3),err,ermx,dvec(3,32),dwgt(32),dvec2(3,32),dmat(3,3)
   INTEGER                  :: nn,inn,loop,loop2,ierr
   LOGICAL                  :: nn_found
   INTEGER                  :: istart,iend
   INTEGER                  :: ibnd_n, ibnd_m,nsym, nxxs
   COMPLEX(DP), ALLOCATABLE :: psic_all(:), temppsic_all(:)
   LOGICAL                  :: have_sym
   COMPLEX(DP), ALLOCATABLE :: evc_k(:, :), evc_sk(:, :)
   INTEGER :: igk_k_ik(npwx)
   !! G vector index at irreducible k point
   INTEGER :: igk_k_sk(npwx)
   !! G vector index at S*k point
   REAL(DP) :: g2kin_(npwx)
   !! Dummy g2kin_ to call gk_sort
   !
   IF (noncolin) CALL errore('compute_dmn', 'Non-collinear not implemented', 1)
   IF (gamma_only) CALL errore('compute_dmn', 'gamma-only not implemented', 1)
   IF (wan_mode=='library') CALL errore('compute_dmn', 'library mode not implemented', 1)
   !
   CALL start_clock( 'compute_dmn' )
   !
   dmat=0d0
   dmat(1,1)=1d0
   dmat(2,2)=1d0
   dmat(3,3)=1d0
   if(read_sym)then
      write(stdout,*) ' Reading symmetry from file '//trim(seedname)//'.sym'
      write(stdout,*) ' '
      if(ionode) then
         inquire(file=trim(seedname)//".sym",exist=have_sym)
         if(.not. have_sym) then
            call errore( 'pw2wannier90', 'Could not find the file '&
               &//trim(seedname)//'.sym', 1 )
         endif
         open(NEWUNIT=iun_sym, file=trim(seedname)//".sym",form='formatted')
         read(iun_sym,*) nsym
      end if
      call mp_bcast(nsym,ionode_id, world_comm)
      allocate(invs(nsym),sr(3,3,nsym),tvec(3,nsym), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating invs/sr/tvec', 1)
      invs=-999
      if(ionode) then
         DO isym = 1, nsym
            read(iun_sym,*)
            read(iun_sym,*) sr(:,:,isym), tvec(:,isym)
         end do
         close(iun_sym)
      end if
      call mp_bcast(sr, ionode_id, world_comm)
      call mp_bcast(tvec, ionode_id, world_comm)
      DO isym = 1, nsym
         do jsym=1,nsym
            if(invs(jsym).ge.1) cycle
            v1=MATMUL(MATMUL(tvec(:,isym),sr(:,:,jsym))+tvec(:,jsym),bg)
            if(sum(abs(MATMUL(sr(:,:,isym),sr(:,:,jsym))-dmat))+sum(abs(v1-dble(nint(v1))))<1d-3) then
               invs(isym)=jsym
               invs(jsym)=isym
            end if
         end do
      end do
   else
      nsym=nsymin
      allocate(sr(3,3,nsym),invs(nsym),tvec(3,nsym), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating sr/invs/tvec', 1)
      ! original sr corresponds to transpose(s)
      ! so here we use sr = transpose(original sr)
      DO isym = 1, nsym
        sr(:,:,isym)=transpose(srin(:,:,isym))
      end do
      invs=invsin(1:nsym)
      tvec=MATMUL(at(:,:),ftin(:,1:nsym))
      if(ionode)then
         open(NEWUNIT=iun_sym, file=trim(seedname)//".sym",form='formatted')
         write(iun_sym,"(i5)") nsym
         DO isym = 1, nsym
            write(iun_sym,*)
            write(iun_sym,"(1p,3e23.15)") sr(:,:,isym), tvec(:,isym)
         end do
         close(iun_sym)
      end if
   end if
   DO isym = 1, nsym
      IF (invs(isym) <= 0 .OR. invs(isym) >= nsym + 1) THEN
         CALL errore("compute_dmn", "out of range in invs", 1)
      ENDIF
      v1 = MATMUL(MATMUL(tvec(:,isym),sr(:,:,invs(isym)))+tvec(:,invs(isym)),bg)
      IF (ANY( ABS(MATMUL(sr(:,:,isym), sr(:,:,invs(isym))) - dmat) > 1.d-3) .OR. &
          ANY( ABS(v1 - REAL(NINT(v1), DP)) > 1.d-3 )) THEN
         CALL errore("compute_dmn", "inconsistent invs", 1)
      ENDIF
   ENDDO
   !
   CALL pw2wan_set_symm ( nsym, sr, tvec )
   !
   ! Write rotation matrix and translation vector in Cartesian coordinates.
   WRITE(stdout, '(a,i5)') "  Number of symmetry operators = ", nsym
   DO isym = 1, nsym
      WRITE(stdout, '(2x,i5,a)') isym, "-th symmetry operators is"
      WRITE(stdout, '(3f15.7)') sr(:,:,isym), tvec(:,isym)
      IF (isym == 1) THEN
         dmat = sr(:, :, isym)
         dmat(1, 1) = dmat(1, 1) - 1.d0
         dmat(2, 2) = dmat(2, 2) - 1.d0
         dmat(3, 3) = dmat(3, 3) - 1.d0
         IF (ANY(ABS(dmat) > 1.d-5) .OR. ANY(ABS(tvec(:, isym)) > 1.d-5)) THEN
            CALL errore("compute_dmn", "1st-symmetry operator is the identity operation.", 1)
         ENDIF
      ENDIF
   ENDDO
   !
   ! Setup iks2k and iks2g, such that symmetry operation isym moves
   ! k(iks2k(ik,isym)) to k(ik) + G(iks2g(ik,isym)).
   !
   ALLOCATE(iks2k(iknum, nsym), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating iks2k', 1)
   ALLOCATE(iks2g(iknum, nsym), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating iks2g', 1)
   iks2k(:, :) = -999
   iks2g(:, :) = -999
   !
   ALLOCATE(lfound(MAX(iknum, ngm)), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating lfound', 1)
   DO isym = 1, nsym
      lfound(:) = .FALSE.
      DO ik = 1, iknum
         v1 = xk_all(:, ik)
         v2 = MATMUL(sr(:,:,isym), v1)
         DO ikp = 1, iknum
            IF (lfound(ikp)) CYCLE
            v3 = xk_all(:, ikp)
            v4 = MATMUL(v2-v3, at)
            IF (ALL( ABS( REAL(NINT(v4), DP) - v4 ) < 1.d-5 )) THEN
               iks2k(ik, isym) = ikp
               lfound(ikp) = .TRUE.
            ENDIF
            IF (iks2k(ik, isym) >= 1) EXIT
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(lfound)
   DO isym = 1, nsym
      DO ik = 1, iknum
         ikp = iks2k(ik, isym)
         v1 = xk_all(:, ikp)
         v2 = MATMUL(v1, sr(:,:,isym))
         v3 = xk_all(:, ik)
         DO ig = 1, ngm
            v4 = g(:, ig)
            IF ( ALL( ABS(v3 + v4 - v2) < 1.d-5) ) THEN
               iks2g(ik, isym) = ig
               EXIT
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   IF (COUNT(iks2k <= 0) /= 0) CALL errore("compute_dmn", "inconsistency in iks2k", COUNT(iks2k <= 0))
   !
   ! Setup ik2ir and ir2ik
   ! ik2ir: Gives irreducible-k points from regular-k points. (Global k index)
   ! ir2ik: Gives regular-k points from irreducible-k points. (Global k index)
   !
   ALLOCATE(ik2ir(iknum), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ik2ir', 1)
   ALLOCATE(ir2ik(iknum), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ir2ik', 1)
   ik2ir(:) = -999
   ir2ik(:) = -999
   !
   nir = 0
   DO ik = 1, iknum
      IF (ik2ir(ik) > 0) CYCLE
      nir = nir + 1
      ir2ik(nir) = ik
      ik2ir(ik) = nir
      DO isym = 1, nsym
         ikp = iks2k(ik, isym)
         IF (ik2ir(ikp) > 0) CYCLE
         ik2ir(ikp) = nir
      ENDDO
   ENDDO
   !
   ! Setup iw2ip and ip2iw, which are the conversion table between Wannier
   ! functions and position indexes.
   ! The position indexes are for the list of non-identical Wannier centers.
   ! The Wannier center of iw-th Wannier function is equal to that of
   ! the ip2iw(iw2ip(iw))-th Wannier function.
   !
   ALLOCATE(iw2ip(n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating iw2ip', 1)
   ALLOCATE(ip2iw(n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ip2iw', 1)
   !
   np = 0
   DO iw = 1, n_wannier
      v1 = center_w(:, iw)
      jp = 0
      DO ip = 1, np
         IF ( ALL( ABS( v1 - center_w(:, ip2iw(ip)) ) < 1.d-2 ) ) THEN
            jp = ip
            EXIT
         ENDIF
      ENDDO
      IF (jp == 0) THEN
         np = np + 1
         iw2ip(iw) = np
         ip2iw(np) = iw
      ELSE
         iw2ip(iw) = jp
      ENDIF
   ENDDO
   !
   ALLOCATE(ips2p(np, nsym), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ips2p', 1)
   ALLOCATE(lfound(np), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating lfound', 1)
   ips2p(:, :) = -999
   DO isym = 1, nsym
      lfound=.false.
      do ip=1,np
         v1=center_w(:,ip2iw(ip))
         v2=MATMUL(sr(:,:,isym),(v1+tvec(:,isym)))
         do jp=1,np
            if(lfound(jp)) cycle
            v3=center_w(:,ip2iw(jp))
            v4=MATMUL(v3-v2,bg)
            if(sum(abs(dble(nint(v4))-v4))<1d-2) then
               lfound(jp)=.true.
               ips2p(ip,isym)=jp
               exit !Sym.op.(isym) moves position(ips2p(ip,isym)) to position(ip) + T, where
            end if                                       !T is given by vps2t(:,ip,isym).
         end do
         if(ips2p(ip,isym) <= 0) then
            write(stdout,"(a,3f18.10,a,3f18.10,a)")"  Could not find ",v2,"(",MATMUL(v2,bg),")"
            write(stdout,"(a,3f18.10,a,3f18.10,a)")"  coming from    ",v1,"(",MATMUL(v1,bg),")"
            write(stdout,"(a,i5,a               )")"  of Wannier site",ip,"."
            call errore("compute_dmn", "Error: missing Wannier sites, see the output.", 1)
         end if
      end do
   end do
   allocate(vps2t(3,np,nsym), stat=ierr) !See above.
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating vps2t', 1)
   DO isym = 1, nsym
      do ip=1,np
         v1=center_w(:,ip2iw(ip))
         jp=ips2p(ip,isym)
         v2=center_w(:,ip2iw(jp))
         v3=MATMUL(v2,sr(:,:,isym))-tvec(:,isym)
         vps2t(:,ip,isym)=v3-v1
      end do
   end do
   dvec(:,1:12)=p12
   dvec(:,13:32)=p20
   do ip=1,32
      dvec(:,ip)=dvec(:,ip)/sqrt(sum(dvec(:,ip)**2))
   end do
   dwgt(1:12)=pwg(1)
   dwgt(13:32)=pwg(2)
   !write(stdout,*) sum(dwgt) !Checking the weight sum to be 1.
   allocate(dylm(32,5),vaxis(3,3,n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating dylm', 1)
   dylm=0d0
   vaxis=0d0
   do ip=1,5
      CALL ylm_wannier(dylm(1,ip),2,ip,dvec,32)
   end do
   !do ip=1,5
   !   write(stdout,"(5f25.15)") (sum(dylm(:,ip)*dylm(:,jp)*dwgt)*2d0*tpi,jp=1,5)
   !end do !Checking spherical integral.
   allocate(wws(n_wannier,n_wannier,nsym), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wws', 1)
   wws=0d0
   DO iw = 1, n_wannier
      call set_u_matrix (xaxis(:,iw),zaxis(:,iw),vaxis(:,:,iw))
   end do
   DO isym = 1, nsym
      DO iw = 1, n_wannier
         ip=iw2ip(iw)
         jp=ips2p(ip,isym)
         CALL ylm_wannier(dylm(1,1),l_w(iw),mr_w(iw),MATMUL(vaxis(:,:,iw),dvec),32)
         do jw=1,n_wannier
            if(iw2ip(jw).ne.jp) cycle
            do ir=1,32
               dvec2(:,ir)=MATMUL(sr(:,:,isym),dvec(:,ir))
            end do
            CALL ylm_wannier(dylm(1,2),l_w(jw),mr_w(jw),MATMUL(vaxis(:,:,jw),dvec2),32)
            wws(jw,iw,isym)=sum(dylm(:,1)*dylm(:,2)*dwgt)*2d0*tpi !<Rotated Y(jw)|Not rotated Y(iw)> for sym.op.(isym).
         end do
      end do
   end do
   deallocate(dylm,vaxis)
   DO isym = 1, nsym
      DO iw = 1, n_wannier
         err=abs((sum(wws(:,iw,isym)**2)+sum(wws(iw,:,isym)**2))*.5d0-1d0)
         if(err.gt.1d-3) then
            write(stdout, '(a,i5,a,i5,a)') "compute_dmn: Symmetry operator (", isym, &
                    ") could not transform Wannier function (", iw, ")."
            write(stdout, '(a,f15.7,a  )') "compute_dmn: The error is ", err, "."
            call errore("compute_dmn", "Error: missing Wannier functions, see the output.", 1)
         end if
      end do
   end do

   CALL utility_open_output_file("dmn", .TRUE., iun_dmn)
   IF (ionode) THEN
      WRITE(iun_dmn, '(4i9)') num_bands, nsym, nir, iknum
      WRITE(iun_dmn, *)
      WRITE(iun_dmn, '(10i9)') ik2ir(1:iknum)
      WRITE(iun_dmn, *)
      WRITE(iun_dmn, '(10i9)') ir2ik(1:nir)
      DO ir = 1, nir
         WRITE(iun_dmn, *)
         WRITE(iun_dmn, '(10i9)') iks2k(ir2ik(ir), :)
      ENDDO
   ENDIF ! ionode
   !
   ! Compute and write d matrix for Wannier functions to file
   !
   IF (ionode) THEN
      ALLOCATE(phs(n_wannier, n_wannier), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phs', 1)
      phs = (0.d0, 0.d0)
      WRITE(stdout,'(/)')
      WRITE(stdout,'(a,i8)') '  DMN(d_matrix_wann): nir = ', nir
      DO ir = 1, nir
         CALL print_progress(ir, nir)
         ik = ir2ik(ir)
         DO isym = 1, nsym
            DO iw = 1, n_wannier
               ip = iw2ip(iw)
               jp = ips2p(ip, invs(isym))
               jw = ip2iw(jp)
               v1 = xk(:,iks2k(ik,isym)) - MATMUL(sr(:,:,isym), xk_all(:,ik))
               v2 = MATMUL(v1, sr(:,:,isym))
               phs(iw,iw) = EXP(CMPLX(0.d0, SUM(vps2t(:,jp,isym)*xk_all(:,ik))*tpi, KIND=DP)) & ! Phase of T.k with lattice vectors T of above.
                          * EXP(CMPLX(0.d0, SUM(tvec(:,isym)*v2)*tpi, KIND=DP)) ! Phase of t.G with translation vector t(isym).
            ENDDO
            !
            WRITE(iun_dmn, *)
            WRITE(iun_dmn, '( " (", ES18.10, ",", ES18.10, ")" )') MATMUL(phs, CMPLX(wws(:,:,isym), 0.d0, KIND=DP))
         ENDDO
      ENDDO
      WRITE(stdout, *) ' DMN(d_matrix_wann) calculated'
      DEALLOCATE(phs)
   ENDIF ! ionode
   !
   ! Compute d matrix for Kohn-Sham wavefunctions
   !
   ALLOCATE(phase(dffts%nnr), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phase', 1)
   ALLOCATE(aux(npwx), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating aux', 1)
   ALLOCATE(evc_k(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_k', 1)
   ALLOCATE(evc_sk(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_sk', 1)
   !
   !
   !   USPP
   !
   !
   IF(okvan) THEN
      CALL allocate_bec_type ( nkb, nbnd, becp )
      IF (gamma_only) THEN
         call errore("compute_dmn", "gamma-only mode not implemented", 1)
      ELSE
         ALLOCATE ( becp2(nkb,nbnd), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating becp2', 1)
      ENDIF
   ENDIF
   !
   !     qb is  FT of Q(r)
   !
   IF(okvan) THEN
      ALLOCATE(dxk(3), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating dxk', 1)
      ALLOCATE(ylm(lmaxq*lmaxq), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ylm', 1)
      ALLOCATE(qb(nhm, nhm, ntyp), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating qb', 1)
      !
      ! Unlike in compute_mmn, here dxk is always 0 because we compute overlap between
      ! |psi(k)> and S^{-1}*|psi(Sk)>, which both have periodicity exp(ikr).
      qg = 0.d0
      dxk(:) = 0.d0
      CALL ylmr2(lmaxq*lmaxq, 1, dxk, qg, ylm)
      !
      qg = SQRT(qg) * tpiba
      !
      DO nt = 1, ntyp
         IF (upf(nt)%tvanp) THEN
            DO ih = 1, nh(nt)
               DO jh = 1, nh(nt)
                  CALL qvan2(1, ih, jh, nt, qg, qgm, ylm)
                  qb(ih, jh, nt) = omega * qgm
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      !
      DEALLOCATE(dxk)
      DEALLOCATE(ylm)
      !
   ENDIF
   !
   ALLOCATE(Mkb(num_bands, nbnd), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating Mkb', 1)
   ALLOCATE( workg(npwx), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating workg', 1)
   !
   ! Set up variables and stuff needed to rotate wavefunctions
   nxxs = dffts%nr1x *dffts%nr2x *dffts%nr3x
   ALLOCATE(psic_all(nxxs), temppsic_all(nxxs), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psic_all', 1)
   !
   ! Pool parallelization: divide irreducible k points, not regular k points
   CALL divide(inter_pool_comm, nir, ir_start, ir_end)
   !
   WRITE(stdout,'(/)')
   WRITE(stdout,'(a,i8)') '  DMN(d_matrix_band): nir = ', nir
   WRITE(stdout,'(a,i8)') '  DMN(d_matrix_band): nir in this pool = ', ir_end - ir_start + 1
   !
   DO ir = ir_start, ir_end
      CALL print_progress(ir, ir_end)
      ik_global = ir2ik(ir) ! global index of the ir-th irreducible k point
      !
      ! Read wavefunction and setup npw and igk_k_ik at k
      ikevc = ik_global + ikstart - 1
      CALL utility_setup_wfc_and_pw(ikevc, evc, npw, igk_k_ik)
      !
      ! Trim excluded bands from evc
      evc_k(:, :) = (0.d0, 0.d0)
      ibnd_m = 0
      DO m = 1, nbnd
         IF (excluded_band(m)) CYCLE
         ibnd_m = ibnd_m + 1
         evc_k(:, ibnd_m) = evc(:, m)
      ENDDO
      !
      !  USPP
      !
      IF (okvan) THEN
         CALL init_us_2(npw, igk_k_ik, xk_all(1, ik_global), vkb)
         CALL calbec(npw, vkb, evc_k, becp, num_bands)
      ENDIF
      !
      !
      DO isym = 1, nsym
         isk_global = iks2k(ik_global, isym)
         !
         ! Read wavefunction and setup npw and igk_k_ik at S*k
         ikpevcq = isk_global + ikstart - 1
         CALL utility_setup_wfc_and_pw(ikpevcq, evc, npwq, igk_k_sk)
         !
         ! Trim excluded bands from evc
         evc_sk(:, :) = (0.d0, 0.d0)
         ibnd_m = 0
         DO m = 1, nbnd
            IF (excluded_band(m)) CYCLE
            ibnd_m = ibnd_m + 1
            evc_sk(:, ibnd_m) = evc(:, m)
         ENDDO
         !
         ! apply translation vector t.
         DO ig = 1, npwq
            arg = SUM( ( MATMUL(g(:,igk_k_sk(ig)), sr(:,:,isym)) + xk_all(:, ik_global) ) * tvec(:,isym) ) * tpi
            phase1 = CMPLX(COS(arg), SIN(arg), KIND=DP)
            DO n = 1, num_bands
               evc_sk(ig, n) = evc_sk(ig, n) * phase1
            ENDDO
         ENDDO
         ! compute the phase
         phase(:) = (0.d0,0.d0)
         ! missing phase G of above is given here and below.
         IF(iks2g(ik_global, isym) >= 0) phase(dffts%nl(iks2g(ik_global, isym)))=(1d0,0d0)
         CALL invfft ('Wave', phase, dffts)
         DO n = 1, num_bands
            psic(:) = (0.d0, 0.d0)
            psic(dffts%nl(igk_k_sk(1:npwq))) = evc_sk(1:npwq, n)
            ! go to real space
            CALL invfft ('Wave', psic, dffts)
#if defined(__MPI)
            ! gather among all the CPUs
            CALL gather_grid(dffts, psic, temppsic_all)
            ! apply rotation
            !psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
            psic_all(rir(1:nxxs,isym)) = temppsic_all(1:nxxs)
            ! scatter back a piece to each CPU
            CALL scatter_grid(dffts, psic_all, psic)
#else
            psic(rir(1:nxxs, isym)) = psic(1:nxxs)
#endif
            ! apply phase k -> k+G
            psic(1:dffts%nnr) = psic(1:dffts%nnr) * phase(1:dffts%nnr)
            ! go back to G space
            CALL fwfft ('Wave', psic, dffts)
            evc_sk(1:npw, n)  = psic(dffts%nl (igk_k_ik(1:npw) ) )
         ENDDO
         !
         !  USPP
         !
         IF (okvan) THEN
            IF (gamma_only) THEN
               CALL errore("compute_dmn", "gamma-only mode not implemented", 1)
            ELSE
               CALL calbec(npw, vkb, evc_sk, becp2, num_bands)
            ENDIF
         ENDIF
         !
         Mkb(:,:) = (0.0d0, 0.0d0)
         !
         ! Compute Mkb
         ! Mkb(m,n) = < psi_m,k1 | psi_n,k2 >
         !          + \sum_{ijI} qb_{ij}^I * e^-i(0*tau_I)
         !            * <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 >
         !
         IF (gamma_only) THEN
            CALL errore("compute_dmn", "gamma-only mode not implemented", 1)
         ELSEIF(noncolin) THEN
            CALL errore("compute_dmn", "Non-collinear not implemented", 1)
         ELSE
            CALL ZGEMM('C', 'N', num_bands, num_bands, npw, &
               (1.d0, 0.d0), evc_k, npwx, evc_sk, npwx, &
               (0.d0, 0.d0), Mkb, num_bands)
         ENDIF
         CALL mp_sum(Mkb, intra_pool_comm)
         !
         ! USPP contribution to Mkb
         !
         IF (okvan) THEN
            ijkb0 = 0
            DO nt = 1, ntyp
               IF ( upf(nt)%tvanp ) THEN
                  DO na = 1, nat
                     !
                     IF ( ityp(na) == nt ) THEN
                        DO jh = 1, nh(nt)
                           jkb = ijkb0 + jh
                           DO ih = 1, nh(nt)
                              ikb = ijkb0 + ih
                              !
                              DO m = 1, num_bands
                                 IF (gamma_only) THEN
                                    CALL errore("compute_dmn", "gamma-only mode not implemented", 1)
                                 ELSE
                                    DO n = 1, num_bands
                                       Mkb(m,n) = Mkb(m,n) + &
                                       qb(ih,jh,nt) * CONJG(becp%k(ikb, m)) * becp2(jkb, n)
                                    ENDDO
                                 ENDIF
                              ENDDO ! m
                           ENDDO !ih
                        ENDDO !jh
                        ijkb0 = ijkb0 + nh(nt)
                     ENDIF  !ityp
                  ENDDO  !nat
               ELSE  !tvanp
                  DO na = 1, nat
                     IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                  ENDDO
               ENDIF !tvanp
            ENDDO !ntyp
         ENDIF ! okvan
         !
         ! Write Mkb to file
         !
         IF (me_pool == root_pool) THEN
            WRITE (iun_dmn, *)
            DO n = 1, num_bands
               DO m = 1, num_bands
                  WRITE (iun_dmn, '( " (", ES18.10, ",", ES18.10, ")" )') CONJG(Mkb(n,m))
               ENDDO
            ENDDO
         ENDIF
      ENDDO !isym
   ENDDO  !ik
   !
   WRITE(stdout, *) ' DMN(d_matrix_band) calculated'
   !
   IF (me_pool == root_pool .AND. wan_mode=='standalone') CLOSE (iun_dmn, STATUS="KEEP")
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   CALL utility_merge_files("dmn", .TRUE.)
   !
   DEALLOCATE (Mkb, phase)
   DEALLOCATE(temppsic_all, psic_all)
   DEALLOCATE(aux)
   DEALLOCATE(evc_k)
   DEALLOCATE(evc_sk)
   !
   IF(okvan) THEN
      DEALLOCATE (qb)
      CALL deallocate_bec_type (becp)
      IF (gamma_only) THEN
         CALL errore('compute_dmn','gamma-only not implemented',1)
      ELSE
         DEALLOCATE (becp2)
      ENDIF
   ENDIF
   !
   CALL stop_clock( 'compute_dmn' )
   !
END SUBROUTINE compute_dmn
!
!-----------------------------------------------------------------------
SUBROUTINE compute_mmn
   !-----------------------------------------------------------------------
   !
   USE kinds,           ONLY : DP
   USE mp_world,        ONLY : world_comm
   USE io_global,       ONLY : stdout, ionode
   USE wvfct,           ONLY : nbnd, npwx
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions,   ONLY : evc, psic, psic_nc
   USE fft_base,        ONLY : dffts, dfftp
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE klist,           ONLY : nks, nkstot, xk, igk_k, ngk
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g, ngm, gstart
   USE cell_base,       ONLY : omega, alat, tpiba, at, bg
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,       ONLY : tpi
   USE uspp,            ONLY : nkb, vkb, okvan
   USE uspp_param,      ONLY : upf, nh, lmaxq, nhm
   USE upf_spinorb,     ONLY : transform_qq_so
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE mp_pools,        ONLY : intra_pool_comm, root_pool, my_pool_id, me_pool, npool
   USE mp,              ONLY : mp_sum, mp_barrier
   USE noncollin_module,ONLY : noncolin, npol, lspinorb
   USE lsda_mod,        ONLY : lsda, isk
   USE gvecw,           ONLY : gcutw
   USE wannier
   USE uspp_init,       ONLY : init_us_2

   IMPLICIT NONE
   !
   complex(DP), parameter :: cmplx_i=(0.0_DP,1.0_DP)
   !
   INTEGER :: npw, mmn_tot, ik, ikp, ipol, ib, npwq, i, m, n
   INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0
   INTEGER :: ikevc, ikpevcq, s, counter, ik_g_w90
   COMPLEX(DP), ALLOCATABLE :: phase(:), aux(:), evc_kb_m(:,:), evc_kb(:,:), &
                               Mkb(:,:), aux_nc(:,:)
   COMPLEX(DP), ALLOCATABLE :: evc_k(:, :)
   !! Wavefunction at k. Contains only the included bands.
   COMPLEX(DP), ALLOCATABLE :: qb(:,:,:,:,:), qgm(:), qq_so(:,:,:,:)
   real(DP), ALLOCATABLE    :: qg(:), ylm(:,:), dxk(:,:,:)
   COMPLEX(DP)              :: mmn, phase1
   real(DP)                 :: arg, g_(3)
   INTEGER                  :: nn,inn,loop,loop2
   LOGICAL                  :: nn_found
   INTEGER                  :: istart,iend
   INTEGER                  :: ibnd_n, ibnd_m
   INTEGER :: ierr
   TYPE(bec_type) :: becp2
   !
   INTEGER, EXTERNAL :: global_kpoint_index
   !
   CALL start_clock( 'compute_mmn' )
   !
   ALLOCATE(phase(dffts%nnr), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phase', 1)
   ALLOCATE(evc_kb(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_kb', 1)
   ALLOCATE(evc_k(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_k', 1)
   ALLOCATE(Mkb(num_bands, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating Mkb', 1)
   !
   IF(noncolin) THEN
      ALLOCATE( aux_nc(npwx,npol), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating aux_nc', 1)
   ELSE
      ALLOCATE( aux(npwx), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating aux', 1)
   ENDIF

   IF (gamma_only) THEN
      ALLOCATE(evc_kb_m(npwx, num_bands), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_kb_m', 1)
   ENDIF

   IF (wan_mode=='library') THEN
      ALLOCATE(m_mat(num_bands, num_bands, nnb, iknum))
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating m_mat', 1)
   ENDIF

   IF (wan_mode=='standalone') THEN
      CALL utility_open_output_file("mmn", .TRUE., iun_mmn)
      IF (ionode) WRITE(iun_mmn, *) num_bands, iknum, nnb
   ENDIF
   !
   !   USPP
   !
   IF(okvan) THEN
      CALL allocate_bec_type(nkb, num_bands, becp)
      CALL allocate_bec_type(nkb, num_bands, becp2)
      !
      !     qb is  FT of Q(r)
      !
      ALLOCATE(qg(nnb), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating qg', 1)
      ALLOCATE(dxk(3, nnb, nks), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating dxk', 1)
      ALLOCATE(ylm(nnb, lmaxq*lmaxq), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ylm', 1)
      ALLOCATE(qgm(nnb), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating qgm', 1)
      ALLOCATE(qb(nhm, nhm, ntyp, nnb, nks), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating qb', 1)
      qb = (0.0_DP, 0.0_DP)
      ALLOCATE(qq_so(nhm, nhm, 4, ntyp), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating qq_so', 1)
      !
      DO ik = 1, nks
         IF (lsda .AND. isk(ik) /= ispinw) CYCLE
         ik_g_w90 = global_kpoint_index(nkstot, ik) - ikstart + 1
         !
         DO ib = 1, nnb
            ikp = kpb(ik_g_w90, ib)
            !
            g_(:) = REAL(g_kpb(:, ik_g_w90, ib), KIND=DP)
            CALL cryst_to_cart (1, g_, bg, 1)
            dxk(:, ib, ik) = xk_all(:, ikp) + g_(:) - xk(:, ik)
            qg(ib) = SUM(dxk(:, ib, ik) * dxk(:, ib, ik))
         ENDDO
!         write (stdout,'(i3,12f8.4)')  ik, qg((ik-1)*nnb+1:ik*nnb)
         !
         CALL ylmr2 (lmaxq*lmaxq, nnb, dxk(:, :, ik), qg, ylm)
         qg(:) = sqrt(qg(:)) * tpiba
         !
         DO nt = 1, ntyp
            IF (upf(nt)%tvanp ) THEN
               DO ih = 1, nh (nt)
                  DO jh = 1, nh (nt)
                     CALL qvan2 (nnb, ih, jh, nt, qg, qgm, ylm)
                     qb(ih, jh, nt, 1:nnb, ik) = omega * qgm(1:nnb)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      !
      DEALLOCATE(qg, qgm, ylm)
      !
   ENDIF
   !
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      !
      ik_g_w90 = global_kpoint_index(nkstot, ik) - ikstart + 1
      npw = ngk(ik)
      !
      ! Read wavefunctions at k, exclude the excluded bands
      !
      CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1 )
      !
      ! Trim excluded bands from evc
      ibnd_m = 0
      DO m = 1, nbnd
         IF (excluded_band(m)) CYCLE
         ibnd_m = ibnd_m + 1
         evc_k(:, ibnd_m) = evc(:, m)
      ENDDO
      !
      !  USPP
      !
      IF(okvan) THEN
         !
         ! Compute the product of beta functions with |psi_k>
         CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
         CALL calbec (npw, vkb, evc_k, becp, num_bands)
         !
      ENDIF
      !
      !
      !do ib=1,nnb(ik)
      DO ib=1,nnb
         !
         ! Read wavefunction at k+b. If okvan, also compute the product
         ! of beta functions with the wavefunctions.
         !
         IF (okvan) THEN
            IF (gamma_only) THEN
               CALL utility_compute_u_kb(ik, ib, evc_kb, becp_kb=becp2, evc_kb_m=evc_kb_m)
            ELSE
               CALL utility_compute_u_kb(ik, ib, evc_kb, becp_kb=becp2)
            ENDIF
         ELSE
            IF (gamma_only) THEN
               CALL utility_compute_u_kb(ik, ib, evc_kb, evc_kb_m=evc_kb_m)
            ELSE
               CALL utility_compute_u_kb(ik, ib, evc_kb)
            ENDIF
         ENDIF
         !
         IF(okvan .AND. lspinorb) CALL transform_qq_so(qb(:,:,:,ib,ik), qq_so)
         !
         !
         Mkb(:,:) = (0.0d0,0.0d0)
         !
         ! loops on bands
         !
         IF (wan_mode=='standalone') THEN
            IF (me_pool == root_pool) WRITE (iun_mmn, '(5i10)') &
               ik_g_w90, kpb(ik_g_w90, ib), (g_kpb(ipol, ik_g_w90, ib), ipol=1,3)
         ENDIF
         !
         !  Mkb(m,n) = Mkb(m,n) + \sum_{ijI} qb_{ij}^I * e^-i(b*tau_I)
         !             <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 >
         !
         IF (gamma_only) THEN
            DO m = 1, num_bands
               DO n = 1, m ! Mkb(m,n) is symmetric in m and n for gamma_only case
                  ! The upper half is filled later by symmetry
                  mmn = dot_product(evc_k(1:npw,m), evc_kb(1:npw, n)) &
                      + CONJG(dot_product(evc_k(1:npw, m), evc_kb_m(1:npw, n)))
                  Mkb(m, n) = mmn + Mkb(m, n)
               ENDDO
            ENDDO
         ELSEIF(noncolin) THEN
            CALL ZGEMM('C', 'N', num_bands, num_bands, npwx*npol, &
               (1.d0, 0.d0), evc_k, npwx*npol, evc_kb, npwx*npol, &
               (0.d0, 0.d0), Mkb, num_bands)
         ELSE
            CALL ZGEMM('C', 'N', num_bands, num_bands, npw, &
               (1.d0, 0.d0), evc_k, npwx, evc_kb, npwx, &
               (0.d0, 0.d0), Mkb, num_bands)
         ENDIF
         !
         ! updating of the elements of the matrix Mkb
         !
         CALL mp_sum(Mkb, intra_pool_comm)
         !
         IF (okvan) THEN
            ijkb0 = 0
            DO nt = 1, ntyp
               !
               IF (.NOT. upf(nt)%tvanp) THEN
                  DO na = 1, nat
                     IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                  ENDDO
                  CYCLE
               ENDIF
               !
               DO na = 1, nat
                  !
                  IF ( ityp(na) /= nt ) CYCLE
                  !
                  arg = dot_product( dxk(:,ib,ik), tau(:,na) ) * tpi
                  phase1 = cmplx( cos(arg), -sin(arg) ,kind=DP)
                  !
                  DO jh = 1, nh(nt)
                     jkb = ijkb0 + jh
                     DO ih = 1, nh(nt)
                        ikb = ijkb0 + ih
                        !
                        DO m = 1, num_bands
                           IF (gamma_only) THEN
                              DO n = 1, m ! Mkb(m,n) is symmetric in m and n for gamma_only case
                                 ! The upper half is filled later by symmetry
                                 Mkb(m,n) = Mkb(m,n) + &
                                      phase1 * qb(ih,jh,nt,ib,ik) * &
                                      becp%r(ikb,m) * becp2%r(jkb,n)
                              ENDDO
                           ELSEIF (noncolin) then
                              DO n = 1, num_bands
                                 IF (lspinorb) THEN
                                    Mkb(m,n) = Mkb(m,n) + &
                                      phase1 * ( &
                                         qq_so(ih,jh,1,nt) * conjg( becp%nc(ikb, 1, m) ) * becp2%nc(jkb, 1, n) &
                                       + qq_so(ih,jh,2,nt) * conjg( becp%nc(ikb, 1, m) ) * becp2%nc(jkb, 2, n) &
                                       + qq_so(ih,jh,3,nt) * conjg( becp%nc(ikb, 2, m) ) * becp2%nc(jkb, 1, n) &
                                       + qq_so(ih,jh,4,nt) * conjg( becp%nc(ikb, 2, m) ) * becp2%nc(jkb, 2, n) &
                                       )
                                 ELSE
                                    Mkb(m,n) = Mkb(m,n) + &
                                      phase1 * qb(ih,jh,nt,ib,ik) * &
                                      (conjg( becp%nc(ikb, 1, m) ) * becp2%nc(jkb, 1, n) &
                                        + conjg( becp%nc(ikb, 2, m) ) * becp2%nc(jkb, 2, n) )
                                 ENDIF
                              ENDDO
                           ELSE
                              DO n = 1, num_bands
                                 Mkb(m,n) = Mkb(m,n) + &
                                      phase1 * qb(ih,jh,nt,ib,ik) * &
                                      conjg( becp%k(ikb,m) ) * becp2%k(jkb,n)
                              ENDDO
                           ENDIF
                        ENDDO ! m
                     ENDDO !ih
                  ENDDO !jh
                  ijkb0 = ijkb0 + nh(nt)
                  !
               ENDDO  !nat
            ENDDO !ntyp
         ENDIF ! okvan
         !
         ! Mkb(m,n) is symmetric in m and n for gamma_only case. Fill the upper half.
         IF (gamma_only) THEN
            DO m = 1, num_bands
               DO n = m + 1, num_bands
                  Mkb(m, n) = Mkb(n, m)
               ENDDO
            ENDDO
         ENDIF
         !
         IF (wan_mode=='standalone') THEN
            DO n = 1, num_bands
               DO m = 1, num_bands
                  IF (me_pool == root_pool) WRITE (iun_mmn,'(2f18.12)') Mkb(m,n)
               ENDDO
            ENDDO
         ELSEIF (wan_mode=='library') THEN
            m_mat(:, :, ib, ik_g_w90) = Mkb(:, :)
         ELSE
            CALL errore('compute_mmn',' value of wan_mode not recognised',1)
         ENDIF
         !
      ENDDO !ib
   ENDDO  !ik
   !
   IF (me_pool == root_pool .AND. wan_mode=='standalone') CLOSE (iun_mmn, STATUS="KEEP")
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   !
   CALL utility_merge_files("mmn", .TRUE.)
   !
   IF (gamma_only) DEALLOCATE(evc_kb_m)
   DEALLOCATE(Mkb)
   DEALLOCATE(phase)
   DEALLOCATE(evc_k)
   DEALLOCATE(evc_kb)
   IF (noncolin) THEN
      DEALLOCATE(aux_nc)
   ELSE
      DEALLOCATE(aux)
   ENDIF
   !
   IF(okvan) THEN
      DEALLOCATE(dxk)
      DEALLOCATE(qb)
      DEALLOCATE(qq_so)
      CALL deallocate_bec_type(becp)
      CALL deallocate_bec_type(becp2)
    ENDIF
   !
   WRITE(stdout,*) ' MMN calculated'
   !
   CALL stop_clock( 'compute_mmn' )
   !
END SUBROUTINE compute_mmn

!--------------------------------------------------------------------------
SUBROUTINE compute_mmn_ibz
   !-----------------------------------------------------------------------
   !
   USE io_global,       ONLY : stdout, ionode
   USE kinds,           ONLY : DP
   USE wvfct,           ONLY : nbnd, npwx
   USE control_flags,   ONLY : gamma_only
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE wavefunctions,   ONLY : evc, psic
   USE fft_base,        ONLY : dffts
   !USE klist,           ONLY : nkstot, xk, ngk
   USE klist,           ONLY : ngk, igk_k, xk
   USE noncollin_module,ONLY : noncolin, npol, lspinorb
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,       ONLY : tpi
   USE cell_base,       ONLY : at
   USE uspp,            ONLY : nkb, vkb, okvan
   USE uspp_param,      ONLY : upf, nh, nhm
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE symm_base,       ONLY : invs
   USE mp,              ONLY : mp_sum
   USE mp_pools,        ONLY : intra_pool_comm
   USE wannier,         ONLY : iknum, ikstart, nnb, nexband, excluded_band, &
                               seedname, bvec, header_len
   USE uspp_init,       ONLY : init_us_2
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), EXTERNAL  :: zdotc
   !
   complex(DP), parameter :: cmplx_i=(0.0_DP,1.0_DP)
   !
   INTEGER                  :: ik, ib, ikp, ikevc1, ikevc2, isym, m, n, iun_mmn
   INTEGER                  :: ih, jh, nt, na, ikb, jkb, ijkb0
   INTEGER                  :: npw, npwq, ierr
   INTEGER, ALLOCATABLE     :: rir(:,:)
   CHARACTER (len=9)        :: cdate, ctime
   CHARACTER (len=header_len) :: header
   REAL(DP)                 :: kdiff(3), skp(3), arg, kpb(3)
   REAL(DP), ALLOCATABLE    :: xkc(:,:)
   COMPLEX(DP)              :: mmn, phase1
   COMPLEX(DP), ALLOCATABLE :: Mkb(:,:), evc2(:,:), evc_kb(:,:)
   COMPLEX(DP), ALLOCATABLE :: qb(:,:,:,:), qq_so(:,:,:,:,:)
   TYPE(bec_type)           :: becp1, becp2
   ! symmetry operation with time reversal symmetry
   INTEGER                  :: nsym2, s2(3,3,96), invs2(96), t_rev2(96), t_rev_spin(2,2)
   REAL(DP)                 :: sr2(3,3,96), ft2(3,96)
   !
   CALL start_clock( 'compute_immn' )
   !
   CALL setup_symm_time_reversal()
   !
   CALL pw2wan_set_symm_with_time_reversal( dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )
   !
   !
   CALL save_sym_info()
   !
   CALL date_and_tim( cdate, ctime )
   header='IBZ Mmn created on '//cdate//' at '//ctime
   IF (ionode) THEN
      OPEN(NEWUNIT=iun_mmn, file=trim(seedname)//".immn",form='formatted')
      WRITE (iun_mmn,*) header
      WRITE (iun_mmn,*) nbnd-nexband, iknum, nnb
   ENDIF
   !
   WRITE(stdout,'(a,i8)') '  MMN: iknum = ',iknum
   !
   IF (okvan) THEN
      ALLOCATE(qb(nhm, nhm, ntyp, nnb), qq_so(nhm, nhm, 4, ntyp, nnb), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating qb/qq_so', 1)
      CALL init_qb_so(qb, qq_so)
      CALL allocate_bec_type(nkb, nbnd, becp)
      CALL allocate_bec_type(nkb, nbnd, becp1)
      CALL allocate_bec_type(nkb, nbnd, becp2)
   END IF
   !
   ALLOCATE( evc2(npwx*npol,nbnd), evc_kb(npwx*npol,nbnd), Mkb(nbnd,nbnd), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc2/evc_kb/Mkb', 1)
   !
   ! calculate <psi_k | e^{-ibr} | psi_k+b>
   ! ik : <psi_k|
   DO ik=1, iknum
      WRITE (stdout,'(i8)',advance='no') ik
      IF( MOD(ik,10) == 0 ) WRITE (stdout,*)
      FLUSH(stdout)
      ikevc1 = ik + ikstart - 1
      CALL davcio (evc, 2*nwordwfc, iunwfc, ikevc1, -1 )
      npw = ngk(ik)
      !
      !  USPP
      !
      IF (okvan) THEN
         ! becp for k
         CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
         CALL calbec (npw, vkb, evc, becp)
      END IF
      !
      ! |psi_k+b>,  k+b = xk(:,ik) + xbvec(:,ib) = s.k + G = s.xk(:,ikp) + kdiff
      DO ib=1, nnb
         CALL kpb_search(ik, ib, ikp, isym, kdiff, skp)
         IF (ionode) WRITE (iun_mmn,'(7i5)') ik, ikp, (nint(kdiff(n)), n=1,3)
         ikevc2 = ikp + ikstart - 1
         CALL davcio (evc2, 2*nwordwfc, iunwfc, ikevc2, -1 )
         CALL rotate_evc(isym, ikp, ik, kdiff, evc2, evc_kb)
         !
         ! USPP
         !
         IF (okvan) THEN
            ! calculate becp1 = <beta|psi_kp> and rotate => becp2
            npwq = ngk(ikp)
            kpb = skp + kdiff
            call cryst_to_cart(1, kpb, at, 1)
            call cryst_to_cart(1, skp, at, 1)
            CALL init_us_2 (npwq, igk_k(1,ikp), xk(1,ikp), vkb)
            CALL calbec ( npwq, vkb, evc2, becp1 )
            ! t_rev2(isym) = 0 or 1  ==>  1-2*t_rev2(isym) = 1 or -1
            CALL rotate_becp(isym, 1-2*t_rev2(isym), xk(1,ikp), skp, becp1, becp2)
         END IF
         !
         ! plane wave part
         !
         Mkb = (0.0d0,0.0d0)
         DO n=1, nbnd
            IF (excluded_band(n)) CYCLE
            DO m=1, nbnd
               IF (excluded_band(m)) CYCLE
               IF (noncolin) THEN
                  mmn = zdotc (npw, evc(1,m),1,evc_kb(1,n),1) &
                      + zdotc (npw, evc(npwx+1,m),1,evc_kb(npwx+1,n),1)
               ELSE
                  mmn = zdotc (npw, evc(1,m),1,evc_kb(1,n),1)
               END IF
               !CALL mp_sum(mmn, intra_pool_comm)
               Mkb(m,n) = mmn
            END DO
         END DO
         CALL mp_sum(Mkb, intra_pool_comm)
         !
         ! USPP/PAW part
         !
         IF (okvan) THEN
            ijkb0 = 0
            DO nt = 1, ntyp
               IF ( upf(nt)%tvanp ) THEN
                  DO na = 1, nat
                     !
                     ! e^{-ibr}
                     arg = dot_product( bvec(:,ib), tau(:,na) ) * tpi
                     phase1 = cmplx( cos(arg), -sin(arg) ,kind=DP)
                     !
                     IF ( ityp(na) == nt ) THEN
                        DO jh = 1, nh(nt)
                           jkb = ijkb0 + jh
                           DO ih = 1, nh(nt)
                              ikb = ijkb0 + ih
                              !
                              DO m = 1,nbnd
                                 IF (excluded_band(m)) CYCLE
                                 DO n = 1,nbnd
                                    IF (excluded_band(n)) CYCLE
                                    IF (gamma_only) THEN
                                       mmn = phase1 * qb(ih,jh,nt,ib) * becp%r(ikb,m) * becp2%r(jkb,n)
                                    ELSE IF (noncolin) THEN
                                       IF (lspinorb) THEN
                                          mmn = phase1 * ( &
                                               qq_so(ih,jh,1,nt,ib) * conjg( becp%nc(ikb, 1, m) ) * becp2%nc(jkb, 1, n) &
                                             + qq_so(ih,jh,2,nt,ib) * conjg( becp%nc(ikb, 1, m) ) * becp2%nc(jkb, 2, n) &
                                             + qq_so(ih,jh,3,nt,ib) * conjg( becp%nc(ikb, 2, m) ) * becp2%nc(jkb, 1, n) &
                                             + qq_so(ih,jh,4,nt,ib) * conjg( becp%nc(ikb, 2, m) ) * becp2%nc(jkb, 2, n) &
                                             )
                                       ELSE
                                          mmn =  phase1 * qb(ih,jh,nt,ib) * &
                                            (conjg( becp%nc(ikb, 1, m) ) * becp2%nc(jkb, 1, n) &
                                           + conjg( becp%nc(ikb, 2, m) ) * becp2%nc(jkb, 2, n) )
                                       END IF
                                    ELSE
                                       mmn = phase1 * qb(ih,jh,nt,ib) * conjg( becp%k(ikb,m) ) * becp2%k(jkb,n)
                                    END IF
                                    Mkb(m,n) = Mkb(m,n) + mmn
                                 ENDDO
                              ENDDO ! m
                           ENDDO !ih
                        ENDDO !jh
                        ijkb0 = ijkb0 + nh(nt)
                     ENDIF  !ityp
                  ENDDO  !nat
               ELSE  !tvanp or tpanp
                  DO na = 1, nat
                     IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                  ENDDO
               ENDIF !tvanp
            ENDDO !ntyp
         END IF !okvan
         !
         ! write Mkb(m,n) => prefix.mmn
         !
         DO n=1, nbnd
            IF (excluded_band(n)) CYCLE
            DO m=1, nbnd
               IF (excluded_band(m)) CYCLE
               IF (ionode) WRITE (iun_mmn,'(2f18.12)') Mkb(m,n)
            END DO
         END DO
         !
      END DO !ib
      !
   END DO !ik
   !
   IF (okvan) THEN
      CALL deallocate_bec_type(becp)
      CALL deallocate_bec_type(becp1)
      CALL deallocate_bec_type(becp2)
   END IF
   !
   IF (ionode) CLOSE(iun_mmn)
   WRITE(stdout,'(/)')
   WRITE(stdout,*) ' IBZ MMN calculated'
   !
   CALL stop_clock( 'compute_immn' )
   !
   CONTAINS
   !
   SUBROUTINE setup_symm_time_reversal()
      ! generate symmetry operation with time reversal symmetry
      ! 1..nsym => 1..2*nsym  if time reversal
      ! nsym, s, ft, t_rev, invs, sr ==> nsym2, s2, ft2, t_rev2, invs2, sr2
      USE symm_base,       ONLY : nsym, s, ft, time_reversal, t_rev, invs, sr
      s2(:,:,1:nsym) = s(:,:,1:nsym)
      sr2(:,:,1:nsym) = sr(:,:,1:nsym)
      ft2(:,1:nsym) = ft(:,1:nsym)
      t_rev2(1:nsym) = t_rev(1:nsym)
      invs2(1:nsym) = invs(1:nsym)
      nsym2 = nsym
      ! t_rev_spin = -i sigma_y (for time reversal)
      t_rev_spin = 0
      t_rev_spin(1,2) = -1
      t_rev_spin(2,1) = 1
      IF (time_reversal) THEN
         nsym2 = 2*nsym
         s2(:,:,nsym+1:2*nsym) = s(:,:,1:nsym)
         sr2(:,:,nsym+1:2*nsym) = sr(:,:,1:nsym)
         ft2(:,nsym+1:2*nsym) = ft(:,1:nsym)
         t_rev2(nsym+1:2*nsym) = 1
         invs2(nsym+1:2*nsym) = invs(1:nsym) + nsym
      END IF
   END SUBROUTINE
   !
   SUBROUTINE save_sym_info()
      !
      ! save symmetry information, representation matrix, repmat(nbnd, nbnd, nsym) and rotation matrix, rotmat(nwan, nwan, nsym)
      ! repmat: representation of little group G_k
      ! rotmat: rotation matrix for initial Wannier function
      !
      USE klist,           ONLY : nkstot
      USE wvfct,           ONLY : nbnd
      USE symm_base,       ONLY : nsym, time_reversal, sname
      USE start_k,         ONLY : k1, k2, k3, nk1, nk2, nk3
      USE wannier,         ONLY : seedname, n_wannier, nexband
      !
      INTEGER                  :: i, m, n, iun_sym, mi, ni, a, ncount
      CHARACTER (len=9)        :: cdate,ctime
      REAL(DP)                 :: srt(3,3)
      COMPLEX(DP)              :: tr, u_spin(2,2)
      COMPLEX(DP), ALLOCATABLE :: rotmat(:,:,:), repmat(:,:,:,:)
      !
      ALLOCATE( repmat(nbnd, nbnd, nsym2, nkstot), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating repmat', 1)
      CALL get_representation_matrix(repmat)
      !
      ALLOCATE( rotmat(n_wannier, n_wannier, nsym2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating rotmat', 1)
      CALL get_rotation_matrix(rotmat)
      !
      ! save repmat
      !
      CALL date_and_tim( cdate, ctime )
      IF (ionode) THEN
         OPEN(NEWUNIT=iun_sym, file=trim(seedname)//".isym", form='formatted')
         !
         ! save symmetry operation
         !
         WRITE (iun_sym,*) 'Symmetry information created on  '//cdate//' at '//ctime
         IF (noncolin) THEN
            WRITE (iun_sym,*) nsym2, 1
         ELSE
            WRITE (iun_sym,*) nsym2, 0
         END IF
         DO i=1, nsym2
            IF (i > nsym) THEN
               WRITE (iun_sym,'(i3,":  ",a, "+T")') i, sname(i-nsym)
            ELSE
               WRITE (iun_sym,'(i3,":  ",a)') i, sname(i)
            END IF
            DO a=1, 3
               write(iun_sym,'(3i10)') s2(a,:,i)
            END DO
            write(iun_sym,'(3f20.12)') ft2(:,i)
            write(iun_sym,*) t_rev2(i)
            !
            IF (noncolin) THEN
               srt = transpose(sr2(:,:,i))
               CALL find_u(srt, u_spin)
               write(iun_sym,'(2f20.12)') (u_spin(1,a), a=1,2)
               write(iun_sym,'(2f20.12)') (u_spin(2,a), a=1,2)
            END IF
            !
            write(iun_sym,*) invs2(i)
         END DO

         write(iun_sym, *)
         WRITE (iun_sym,*) 'K points'
         write(iun_sym,'(3i10)') iknum
         DO i=1, iknum
           write(iun_sym, '(3f16.10)') xkc(:,i)
         END DO
         !
         ! save representation matrix
         !
         write(iun_sym, *)
         WRITE (iun_sym,*) 'Representation matrix of G_k'
         ncount = 0
         DO ik=1, nkstot
            DO isym=1, nsym2
               IF (.not. all(repmat(:,:,isym,ik) == 0)) ncount = ncount + 1
            END DO
         END DO
         WRITE (iun_sym,'(3i10)') nbnd - nexband, ncount
         DO ik=1, nkstot
            DO isym=1, nsym2
               IF (all(repmat(:,:,isym,ik) == 0)) CYCLE
               ncount = count(repmat(:,:,isym,ik) /= 0)
               WRITE(iun_sym, *) ik, isym, ncount
               !tr = 0
               !DO m=1, nbnd
               !   tr = tr + repmat(m,m,isym,ik)
               !END DO
               !WRITE(iun_sym, '(2f20.10)') real(tr), aimag(tr)
               mi = 0
               DO m=1, nbnd
                  IF (excluded_band(m)) CYCLE
                  mi = mi + 1
                  ni = 0
                  DO n=1, nbnd
                     IF (excluded_band(n)) CYCLE
                     ni = ni + 1
                     IF (repmat(m,n,isym,ik) == 0) CYCLE
                     WRITE(iun_sym, '(2i5,2f22.15)') mi, ni, real(repmat(m,n,isym,ik)), aimag(repmat(m,n,isym,ik))
                  END DO
               END DO
            END DO
         END DO
         !
         ! save rotation matrix
         !
         write(iun_sym, *)
         WRITE(iun_sym,*) 'Rotation matrix of Wannier functions'
         WRITE(iun_sym,*) n_wannier
         DO isym=1, nsym2
            IF (all(abs(rotmat(:,:,isym)) < 1e-10)) CYCLE
            ncount = count(abs(rotmat(:,:,isym)) >= 1e-10)
            WRITE(iun_sym, *) isym, ncount
            DO m=1, n_wannier
               DO n=1, n_wannier
                  IF (abs(rotmat(m,n,isym)) < 1e-10) CYCLE
                  WRITE(iun_sym, '(2i5,2f22.15)') m, n, real(rotmat(m,n,isym)), aimag(rotmat(m,n,isym))
               END DO
            END DO
         END DO
         CLOSE(iun_sym)
      END IF
   END SUBROUTINE
   !
   SUBROUTINE get_representation_matrix(repmat)
      USE klist,           ONLY : nkstot, xk, ngk, igk_k
      USE wavefunctions,   ONLY : evc
      USE io_files,        ONLY : nwordwfc, iunwfc
      USE wvfct,           ONLY : nbnd, npwx, et
      USE uspp,            ONLY : nkb, vkb, okvan
      USE becmod,          ONLY : bec_type, becp, calbec, &
                                  allocate_bec_type, deallocate_bec_type
      USE cell_base,       ONLY : at
      USE noncollin_module,ONLY : noncolin, npol
      USE mp,              ONLY : mp_sum
      USE mp_pools,        ONLY : intra_pool_comm
      !
      COMPLEX(DP), INTENT(OUT) :: repmat(nbnd, nbnd, nsym2, nkstot)
      INTEGER                  :: i, ik, isym, m, n, npw
      REAL(DP)                 :: k(3), sk(3), kdiff(3)
      COMPLEX(DP), ALLOCATABLE :: spsi(:,:), gpsi(:,:)
      !
      repmat = 0
      !
      ALLOCATE(xkc(3,nkstot), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating xkc', 1)
      xkc(:,1:nkstot)=xk(:,1:nkstot)
      CALL cryst_to_cart(nkstot,xkc,at,-1)
      !
      IF (okvan) CALL allocate_bec_type(nkb, nbnd, becp)
      !
      ALLOCATE (spsi(npwx*npol,nbnd), gpsi(npwx*npol,nbnd), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating spsi/gpsi', 1)
      !
      DO ik=1, nkstot
         k = xkc(:,ik)
         npw = ngk(ik)
         CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
         IF (okvan) THEN
            CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
            CALL calbec(npw, vkb, evc, becp)
            CALL s_psi(npwx, npw, nbnd, evc, spsi)
         ELSE
            spsi = evc
         END IF
         DO isym=1, nsym2
            ! check k = s.k + G(kdiff)
            sk = matmul(s2(:,:,isym), k)
            if (t_rev2(isym) == 1) sk = -sk
            kdiff = k - sk
            IF (any(abs(kdiff - nint(kdiff)) > 1e-5)) CYCLE
            !
            CALL rotate_evc(isym, ik, ik, kdiff, evc, gpsi)
            !
            DO m=1, nbnd
               IF (excluded_band(m)) CYCLE
               DO n=1, nbnd
                 IF (excluded_band(n)) CYCLE
                 IF (abs(et(m,ik) - et(n,ik)) < 1e-5) THEN
                    IF (noncolin) THEN
                       repmat(m,n,isym,ik) = zdotc(npw, spsi(1,m), 1, gpsi(1,n), 1) &
                                           + zdotc(npw, spsi(npwx+1,m), 1, gpsi(npwx+1,n), 1)
                    ELSE
                       repmat(m,n,isym,ik) = zdotc(npw, spsi(1,m), 1, gpsi(1,n), 1)
                    END IF
                 END IF
               END DO
            END DO
         END DO
      END DO
      !
      CALL mp_sum(repmat, intra_pool_comm)
      !
      IF (okvan) CALL deallocate_bec_type(becp)
      !
   END SUBROUTINE
   !
   SUBROUTINE get_rotation_matrix(rotmat)
      !
      !  rotmat(m,n) = <g_m| S^-1 |g_n>
      !
      USE cell_base,       ONLY : at, bg
      USE constants,       ONLY : tpi
      USE wannier,         ONLY : n_wannier, l_w, mr_w, xaxis, zaxis, center_w, &
                                  spin_eig, spin_qaxis
      USE symm_base,       ONLY : d1, d2, d3, nsym
      !
      COMPLEX(DP), INTENT(OUT) :: rotmat(n_wannier, n_wannier, nsym2)
      !
      REAL(DP), parameter :: p12(3,12)=reshape(                            &
         (/0d0, 0d0, 1.00000000000000d0,                                   &
           0.894427190999916d0, 0d0, 0.447213595499958d0,                  &
           0.276393202250021d0, 0.850650808352040d0, 0.447213595499958d0,  &
          -0.723606797749979d0, 0.525731112119134d0, 0.447213595499958d0,  &
          -0.723606797749979d0, -0.525731112119134d0, 0.447213595499958d0, &
           0.276393202250021d0, -0.850650808352040d0, 0.447213595499958d0, &
           0.723606797749979d0, 0.525731112119134d0, -0.447213595499958d0, &
          -0.276393202250021d0, 0.850650808352040d0, -0.447213595499958d0, &
          -0.894427190999916d0, 0d0, -0.447213595499958d0,                 &
          -0.276393202250021d0, -0.850650808352040d0, -0.447213595499958d0,&
           0.723606797749979d0, -0.525731112119134d0, -0.447213595499958d0,&
           0d0, 0d0, -1.00000000000000d0/),(/3,12/))
      REAL(DP), parameter :: p20(3,20)=reshape(                            &
         (/0.525731112119134d0, 0.381966011250105d0, 0.850650808352040d0,  &
          -0.200811415886227d0, 0.618033988749895d0, 0.850650808352040d0,  &
          -0.649839392465813d0, 0d0, 0.850650808352040d0,                  &
          -0.200811415886227d0, -0.618033988749895d0, 0.850650808352040d0, &
           0.525731112119134d0, -0.381966011250105d0, 0.850650808352040d0, &
           0.850650808352040d0, 0.618033988749895d0, 0.200811415886227d0,  &
          -0.324919696232906d0, 1.00000000000000d0, 0.200811415886227d0,   &
          -1.05146222423827d0, 0d0, 0.200811415886227d0,                   &
         -0.324919696232906d0, -1.00000000000000d0, 0.200811415886227d0,   &
          0.850650808352040d0, -0.618033988749895d0, 0.200811415886227d0,  &
          0.324919696232906d0, 1.00000000000000d0, -0.200811415886227d0,   &
         -0.850650808352040d0, 0.618033988749895d0, -0.200811415886227d0,  &
         -0.850650808352040d0, -0.618033988749895d0, -0.200811415886227d0, &
          0.324919696232906d0, -1.00000000000000d0, -0.200811415886227d0,  &
          1.05146222423827d0, 0d0, -0.200811415886227d0,                   &
          0.200811415886227d0, 0.618033988749895d0, -0.850650808352040d0,  &
         -0.525731112119134d0, 0.381966011250105d0, -0.850650808352040d0,  &
         -0.525731112119134d0, -0.381966011250105d0, -0.850650808352040d0, &
          0.200811415886227d0, -0.618033988749895d0, -0.850650808352040d0, &
         0.649839392465813d0, 0d0, -0.850650808352040d0/),(/3,20/))
      REAL(DP), parameter :: pwg(2)=(/2.976190476190479d-2,3.214285714285711d-2/)
      !
      INTEGER               :: ip, jp, isym, iw, jw, np
      REAL(DP)              :: v1(3), v2(3), v3(3), v4(3), err, srt(3,3), tvec(3,96)
      REAL(DP)              :: dvec(3,32), dvec_in(3,32), dwgt(32), dylm1(32), dylm2(32)
      COMPLEX(DP)           :: spin1(2), spin2(2), u_spin(2,2)
      INTEGER, ALLOCATABLE  :: ip2iw(:), iw2ip(:), ips2p(:,:)
      REAL(DP), ALLOCATABLE :: vaxis(:,:,:)
      logical, ALLOCATABLE  :: lfound(:)
      COMPLEX(DP), ALLOCATABLE :: check_mat(:,:)
      INTEGER               :: l, m1, m2
      REAL(DP)              :: mat
      TYPE symmetrization_tensor
          REAL(DP), POINTER :: d(:,:,:)
      END TYPE symmetrization_tensor
      TYPE(symmetrization_tensor) :: D(0:3)
      !
      tvec=matmul(at(:,:),ft2(:,1:nsym2))
      !
      ! set dvec and dwgt for integration of spherical harmonics
      dvec(:,1:12) = p12
      dvec(:,13:32) = p20
      DO ip=1,32
         dvec(:,ip)=dvec(:,ip)/sqrt(sum(dvec(:,ip)**2))
      END DO
      dwgt(1:12)=pwg(1)
      dwgt(13:32)=pwg(2)
      !
      !Conversion table between Wannier and position indexes.
      allocate(iw2ip(n_wannier),ip2iw(n_wannier), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating iw2ip/ip2iw', 1)
      np=0
      do iw=1,n_wannier
         v1=center_w(:,iw)
         jp=0
         do ip=1,np
            if(sum(abs(v1-center_w(:,ip2iw(ip)))).lt.1d-2) then
               jp=ip
               exit
            end if
         end do
         if(jp.eq.0) then
            np=np+1
            iw2ip(iw)=np
            ip2iw(np)=iw
         else
            iw2ip(iw)=jp
         end if
      end do
      !
      allocate(ips2p(np,nsym2),lfound(np), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ips2p/lfound', 1)
      ips2p=-999 ! < 0
      do isym=1,nsym2
         lfound=.false.
         do ip=1,np
            v1=center_w(:,ip2iw(ip))
            v2=matmul(v1+tvec(:,isym), sr2(:,:,isym))
            do jp=1,np
               if(lfound(jp)) cycle
               v3=center_w(:,ip2iw(jp))
               v4=matmul(v3-v2,bg)
               if(sum(abs(dble(nint(v4))-v4)).lt.1d-2) then
                  lfound(jp)=.true.
                  ips2p(ip,isym)=jp
                  exit !Sym.op.(isym) moves position(ips2p(ip,isym)) to position(ip) + T, where
               end if                                       !T is given by vps2t(:,ip,isym).
            end do
            if(ips2p(ip,isym).le.0) then
               write(stdout,"(a,3f18.10,a,3f18.10,a)")"  Could not find ",v2,"(",matmul(v2,bg),")"
               write(stdout,"(a,3f18.10,a,3f18.10,a)")"  coming from    ",v1,"(",matmul(v1,bg),")"
               write(stdout,"(a,i5,a               )")"  of Wannier site",ip,"."
               call errore("compute_mmn_sym", "Error: missing Wannier sites, see the output.", 1)
            end if
         end do
      end do
      !
      allocate( vaxis(3,3,n_wannier), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating vaxis', 1)
      rotmat=0.0d0
      do iw=1,n_wannier
         call set_u_matrix (xaxis(:,iw),zaxis(:,iw),vaxis(:,:,iw))
      end do
      do isym=1,nsym2
         srt = transpose(sr2(:,:,isym))
         CALL find_u(srt, u_spin)
         !IF (t_rev2(isym) == 1) u_spin = matmul(t_rev_spin, conjg(u_spin))
         do iw=1,n_wannier
            ip=iw2ip(iw)
            jp=ips2p(ip,isym)
            !
            dvec_in = matmul(vaxis(:,:,iw),dvec)
            CALL ylm_wannier(dylm1,l_w(iw),mr_w(iw),dvec_in,32)
            !
            do jw=1,n_wannier
               if(iw2ip(jw).ne.jp) cycle
               !
               dvec_in = matmul( vaxis(:,:,jw), matmul( srt, dvec ) )
               CALL ylm_wannier(dylm2,l_w(jw),mr_w(jw),dvec_in,32)
               !
               ! <Y(iw) | S(isym)^-1 Y(jw)>
               rotmat(iw,jw,isym)=sum(dylm1(:)*dylm2(:)*dwgt)*2d0*tpi
               IF (noncolin) THEN
                  spin1(:) = spinor( spin_qaxis(:,iw), spin_eig(iw) )
                  spin2(:) = spinor( spin_qaxis(:,jw), spin_eig(jw) )
                  IF (t_rev2(isym) == 1) THEN
                     ! t_rev_spin^-1 = tr(t_rev_spin)
                     spin2(:) = matmul(transpose(t_rev_spin), conjg(spin2))
                  END IF
                  spin2(:) = matmul(transpose(conjg(u_spin)), spin2)
                  rotmat(iw,jw,isym)=rotmat(iw,jw,isym) * dot_product(spin1,spin2)
               END IF
            end do
         end do
      end do
      deallocate(vaxis)
      deallocate(ips2p, lfound, iw2ip, ip2iw)
      allocate(check_mat(n_wannier, n_wannier), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating check_mat', 1)
      do isym=1,nsym2
         if(t_rev2(isym) == 0) then
            ! rotmat(:,:,isym) is unitary
            check_mat = matmul(rotmat(:,:,isym), transpose(conjg(rotmat(:,:,isym))))
         else
            ! rotmat(:,:,isym) is anti-unitary
            check_mat = (0.d0, 1.d0) * rotmat(:,:,isym)
            check_mat = matmul(check_mat(:,:), transpose(conjg(check_mat(:,:))))
         end if
         do iw=1,n_wannier
            do jw=1, n_wannier
               if( (iw == jw .and. abs(check_mat(iw,jw) - 1) > 1d-3) .or. &
                   (iw /= jw .and. abs(check_mat(iw,jw)) > 1d-3) ) then
                  write(stdout,"(a,i5,a,i5,a)") "compute_mmn_ibz: Symmetry operator (", isym, &
                          ") could not transform Wannier function (", iw, ")."
                  write(stdout,"(a,f15.7,a  )") "compute_mmn_ibz: The error is ", check_mat(iw,jw) , "."
                  write(stdout,*) "rotmat(:,:,isym)"
                  write(stdout,*) rotmat(:,:,isym)
                  write(stdout,*) "check_mat"
                  write(stdout,*) check_mat
                  call errore("compute_mmn_ibz", "Error: missing Wannier functions, see the output.", 1)
               end if
            end do
         end do
      end do
      deallocate(check_mat)
      !
      ! check d_matrix and this calculation
      !
      !write(stdout, *) '  checking d_matrix... \n'
      !CALL d_matrix(d1, d2, d3)
      !D(1)%d => d1 ! d1(3,3,48)
      !D(2)%d => d2 ! d2(5,5,48)
      !D(3)%d => d3 ! d3(7,7,48)
      !do isym=1, nsym
      !   srt = transpose(sr2(:,:,isym))
      !   l = 1
      !   do l=1, 3
      !      do m1 = 1, 2*l+1
      !         dvec_in = dvec
      !         CALL ylm_wannier(dylm1,l,m1,dvec_in,32)
      !         do m2 = 1, 2*l+1
      !            dvec_in = matmul( srt, dvec )
      !            CALL ylm_wannier(dylm2,l,m2,dvec_in,32)
      !            mat = sum(dylm1(:)*dylm2(:)*dwgt)*2d0*tpi  ! rotmat
      !            if(abs(mat - D(l)%d(m1,m2,isym)) > 1e-8) write(stdout, '(2i5, 2f15.7)') m1, m2, mat, D(l)%d(m1,m2,isym)
      !         end do
      !      end do
      !      do m1 = 1, 2*l+1
      !         write(stdout, *) dot_product(D(l)%d(m1,1:2*l+1,isym), D(l)%d(m1,1:2*l+1,isym))
      !      end do
      !   end do
      !end do
   END SUBROUTINE get_rotation_matrix
   !
   SUBROUTINE rotate_evc(isym, ik1, ik2, kdiff, psi, gpsi)
      !-----------------------------------------------------------------------
      ! g psi_k = e^{ik1.rS} u_k(rS) = e^{iSk1.r} u_k(rS)
      !         = e^{ik2.r} [ e^{-iGr} u_k(rS) ]
      ! k2 = s.k1 + G
      ! S=s(:,:,isym)   G=kdiff
      ! gvector order is k(:,ik1) for input and k(:,ik2) for output
      !
      ! with T  rotation + T (apply rotation and then apply T)
      ! k2 = -s.k1 + G
      ! g psi_k = [ e^{iSk1.r} u_k(rS) ]^*
      !         = e^{ik2.r} e^{-iGr} [ u_k(rS) ]^*
      !
      USE wvfct,           ONLY : nbnd, npwx
      USE wavefunctions,   ONLY : evc, psic, psic_nc
      USE fft_base,        ONLY : dffts, dfftp
      USE fft_interfaces,  ONLY : fwfft, invfft
      USE cell_base,       ONLY : bg
      USE constants,       ONLY : tpi
      USE gvect,           ONLY : g, ngm
      USE klist,           ONLY : igk_k, ngk
      USE mp,              ONLY : mp_sum
      USE mp_pools,        ONLY : intra_pool_comm
      USE fft_interfaces,  ONLY : invfft
      USE scatter_mod,     ONLY : gather_grid, scatter_grid
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN):: isym, ik1, ik2
      REAL(DP), INTENT(IN):: kdiff(3)
      COMPLEX(DP), INTENT(IN):: psi(npwx*npol,nbnd)
      COMPLEX(DP), INTENT(OUT):: gpsi(npwx*npol,nbnd)
      !
      INTEGER:: ig, igk_local, igk_global, npw1, npw2, n, nxxs, ipol, istart, isym0
      REAL(DP)                 :: kdiff_cart(3), srt(3,3)
      REAL(DP)                 :: phase_arg
      COMPLEX(DP)              :: phase_factor, u_spin(2,2), u_spin2(2,2)
      COMPLEX(DP), ALLOCATABLE :: phase(:), gpsi_tmp(:,:)
      COMPLEX(DP), ALLOCATABLE :: psic_all(:), temppsic_all(:)
      !
      nxxs = dffts%nr1x *dffts%nr2x *dffts%nr3x
      ALLOCATE( psic_all(nxxs), temppsic_all(nxxs), gpsi_tmp(npwx,2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psic_all/temppsic_all/gpsi_tmp', 1)
      ALLOCATE( phase(dffts%nnr), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phase', 1)
      !
      ! for spin space rotation (only used for noncollinear case)
      !
      srt = transpose(sr2(:,:,isym))
      CALL find_u(srt, u_spin)
      !
      ! kdiff = g(:,igk_global)
      !
      kdiff_cart = kdiff
      CALL cryst_to_cart(1, kdiff_cart, bg, 1)
      igk_local = 0
      DO ig = 1, ngm
         IF ( all (abs(kdiff_cart - g(:,ig)) < 1d-5) ) THEN
            igk_local = ig
            exit
         END IF
      END DO
      igk_global = igk_local
      CALL mp_sum(igk_global, intra_pool_comm)
      !
      ! phase = e^{iGr}
      !
      IF (igk_global > 0) THEN
         phase(:) = (0.d0, 0.d0)
         IF (igk_local > 0) phase( dffts%nl(igk_local) ) = (1.d0, 0.d0)
         CALL invfft ('Wave', phase, dffts)
      END IF
      !
      gpsi = 0
      !
      npw1 = ngk(ik1)
      npw2 = ngk(ik2)
      DO n=1, nbnd
         IF (excluded_band(n)) CYCLE
         !
         ! real space rotation
         !
         gpsi_tmp = 0
         DO ipol=1,npol
            istart = npwx*(ipol-1)
            psic(:) = (0.d0, 0.d0)
            psic(dffts%nl (igk_k (1:npw1,ik1) ) ) = psi (istart+1:istart+npw1, n)
            !
            IF (igk_global > 0 .or. isym > 1) THEN
               CALL invfft ('Wave', psic, dffts)
               IF (isym > 1) THEN
#if defined(__MPI)
                  ! gather among all the CPUs
                  CALL gather_grid(dffts, psic, temppsic_all)
                  ! apply rotation
                  psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                  ! scatter back a piece to each CPU
                  CALL scatter_grid(dffts, psic_all, psic)
#else
                  psic(1:nxxs) = psic(rir(1:nxxs,isym))
#endif
               ENDIF
               IF(t_rev2(isym) == 1) psic = conjg(psic)
               ! apply phase e^{-iGr}
               IF(igk_global > 0) psic(1:dffts%nnr) = psic(1:dffts%nnr) * conjg(phase(1:dffts%nnr))

               CALL fwfft ('Wave', psic, dffts)
            END IF
            !
            gpsi_tmp(1:npw2,ipol)  = psic(dffts%nl (igk_k(1:npw2,ik2) ) )
         END DO
         !
         ! spin space rotation
         !
         DO ipol=1,npol
            istart = npwx*(ipol-1)
            IF (noncolin) THEN
               u_spin2 = u_spin
               IF (t_rev2(isym) == 1) u_spin2 = matmul(t_rev_spin, conjg(u_spin))
               gpsi(istart+1:istart+npw2,n) = matmul(gpsi_tmp(1:npw2,:), u_spin2(ipol,:))
            ELSE
               gpsi(istart+1:istart+npw2,n) = gpsi_tmp(1:npw2,ipol)
            END IF
         END DO
      END DO
      !
      phase_arg = -tpi * dot_product(xkc(:,ik1), ft2(:,isym))
      IF (t_rev2(isym) == 1) phase_arg = -phase_arg
      phase_factor = CMPLX(COS(phase_arg), SIN(phase_arg), KIND=DP)
      gpsi = gpsi * phase_factor
      !
      DEALLOCATE( phase, psic_all, temppsic_all, gpsi_tmp )
   END SUBROUTINE rotate_evc
   !
   SUBROUTINE kpb_search(ik, ib, ikp, isym, kdiff, skp)
      ! input: ik, ib
      ! output: ikp, isym, kdiff, skp = s.xkc
      ! xkc(:,ik) + xbvec(:,ib) = s(:,:,isym) . xkc(:,ikp) + kdiff
      ! if with T
      ! xkc(:,ik) + xbvec(:,ib) = - s(:,:,isym) . xkc(:,ikp) + kdiff
      !
      USE kinds,           ONLY : DP
      USE klist,           ONLY : nkstot
      USE wannier,         ONLY : xbvec
      USE io_global,       ONLY : stdout
      IMPLICIT NONE
      INTEGER, INTENT(in):: ik, ib
      INTEGER, INTENT(out):: ikp, isym
      REAL(DP), INTENT(out):: kdiff(3), skp(3)
      !
      REAL(DP):: k(3), kb(3)
      INTEGER:: i

      k = xkc(:,ik)
      kb = xkc(:,ik) + xbvec(:,ib)
      do isym = 1, nsym2
        do ikp = 1, nkstot
          skp = matmul(s2(:,:,isym), xkc(:,ikp))
          IF (t_rev2(isym) == 1) skp = -skp
          kdiff = kb - skp
          if ( all( abs(nint(kdiff) - kdiff) < 1d-5 ) ) then
            return
          end if
        end do
      end do
      WRITE(stdout, *) "ikp, isym not found for kb = ", kb
      CALL errore('kpb_search', 'ikp, isym not found', 1)
   END SUBROUTINE
   !
   SUBROUTINE init_qb_so(qb, qq_so)
      ! TODO: Use this also in compute_mmn
      USE uspp_param,      ONLY : upf, nh, lmaxq, nhm
      USE upf_spinorb,     ONLY : transform_qq_so
      USE cell_base,       ONLY : omega, tpiba
      USE ions_base,       ONLY : ntyp => nsp
      USE noncollin_module,ONLY : lspinorb
      USE wannier,         ONLY : nnb
      !
      COMPLEX(DP), INTENT(out) :: qb(nhm, nhm, ntyp, nnb), qq_so(nhm, nhm, 4, ntyp, nnb)
      !
      INTEGER :: ih, jh, nt, ib
      REAL(DP), ALLOCATABLE    :: qg(:), ylm(:,:)
      COMPLEX(DP), ALLOCATABLE :: qgm(:), qq_so_tmp(:,:,:,:)
      !
      ALLOCATE( ylm(nnb, lmaxq*lmaxq), qg(nnb), qgm(nnb), qq_so_tmp(nhm, nhm, 4, ntyp), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ylm/qg/qgm/qq_so_tmp', 1)
      !
      DO ib=1, nnb
         qg(ib) = dot_product(bvec(:,ib), bvec(:,ib))
      END DO
      CALL ylmr2(lmaxq*lmaxq, nnb, bvec, qg, ylm)
      qg = sqrt(qg) * tpiba
      !
      DO nt =1, ntyp
         IF( .not. upf(nt)%tvanp ) CYCLE
         DO ih = 1, nh(nt)
            DO jh = 1, nh(nt)
               CALL qvan2(nnb, ih, jh, nt, qg, qgm, ylm)
               qb(ih, jh, nt, :) = omega * qgm(:)
            END DO
         END DO
      END DO
      !
      IF (lspinorb) THEN
         qq_so = (0.d0, 0.d0)
         DO ib=1, nnb
            qq_so_tmp = 0
            CALL transform_qq_so(qb(:,:,:,ib), qq_so_tmp)
            qq_so(:,:,:,:,ib) = qq_so_tmp
         END DO
      END IF
      !
      DEALLOCATE( qg, qgm, ylm, qq_so_tmp )
      !
   END SUBROUTINE
   !
   SUBROUTINE rotate_becp(isym, sgn_sym, xk0, xk, becp1, becp2)
      !
      ! based on becp_rotate_k in us_exx.f90
      !
      USE kinds,           ONLY : DP
      USE cell_base,       ONLY : at, bg
      USE ions_base,       ONLY : tau, nat, ityp
      USE symm_base,       ONLY : irt, d1, d2, d3, nsym, invs
      USE uspp,            ONLY : nkb, ofsbeta, nhtolm, nhtol
      USE uspp_param,      ONLY : nh, upf
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)            :: isym, sgn_sym
      REAL(DP), INTENT(IN)           :: xk0(3), xk(3)
      TYPE (bec_type), INTENT(IN)    :: becp1
      TYPE (bec_type), INTENT(INOUT) :: becp2
      !
      INTEGER     :: ia, nt, ma, ih, m_o, lm_i, l_i, m_i, ikb, oh, okb, isym0, isym_inv
      REAL(DP)    :: xau(3,nat), rau(3,nat), tau_phase
      COMPLEX(DP) :: tau_fact
      COMPLEX(DP) :: u_spin(2,2)
      !
      REAL(DP), TARGET :: d0(1,1,48)
      TYPE symmetrization_tensor
          REAL(DP), POINTER :: d(:,:,:)
      END TYPE symmetrization_tensor
      TYPE(symmetrization_tensor) :: D(0:3)
      !
      !
      IF ( isym == 1 ) THEN
         IF (sgn_sym > 0) THEN
            IF (noncolin) THEN
               becp2%nc = becp1%nc
            ELSE
               becp2%k = becp1%k
            END IF
         ELSE
            IF (noncolin) THEN
               becp2%nc = CONJG(becp1%nc)
            ELSE
               becp2%k = CONJG(becp1%k)
            END IF
         ENDIF
         RETURN
      ENDIF
      !
      CALL d_matrix(d1,d2,d3)
      d0(1,1,:) = 1._dp
      D(0)%d => d0 ! d0(1,1,48)
      D(1)%d => d1 ! d1(3,3,48)
      D(2)%d => d2 ! d2(5,5,48)
      D(3)%d => d3 ! d3(7,7,48)
      ! set isym0 for d0~d3 and irt
      IF (isym > nsym) THEN
        isym0 = isym - nsym
      ELSE
        isym0 = isym
      END IF
      !
      IF (ABS(sgn_sym) /= 1) CALL errore( "becp_rotate", "sign must be +1 or -1", 1 )
      !
      xau = tau
      CALL cryst_to_cart( nat, xau, bg, -1 )
      !
      DO ia = 1,nat
         rau(:,ia) = matmul(xau(:,ia), s2(:,:,isym)) - ft2(:,isym)
      ENDDO
      !
      CALL cryst_to_cart( nat, rau, at, +1 )
      !
      IF (noncolin) THEN
         becp2%nc = 0._dp
         call find_u( transpose(sr2(:,:,isym)), u_spin )
         IF (t_rev2(isym) == 1) u_spin = matmul(t_rev_spin, conjg(u_spin))
      ELSE
         becp2%k = 0._dp
      END IF
      DO ia = 1,nat
         nt = ityp(ia)
         ma = irt(isym0,ia)
         tau_phase = -tpi*( sgn_sym*SUM((tau(:,ma)-rau(:,ia))*xk0) )

         tau_fact = CMPLX(COS(tau_phase), SIN(tau_phase), KIND=DP)
         !
         DO ih = 1, nh(nt)
            !
            lm_i  = nhtolm(ih,nt)
            l_i   = nhtol(ih,nt)
            m_i   = lm_i - l_i**2
            ikb = ofsbeta(ma) + ih
            !
            DO m_o = 1, 2*l_i +1
               oh = ih - m_i + m_o
               okb = ofsbeta(ia) + oh
               IF (noncolin) THEN
                  IF (sgn_sym > 0) THEN
                     becp2%nc(okb, :, :) = becp2%nc(okb, :, :) &
                       + matmul(u_spin, D(l_i)%d(m_i,m_o, isym0) * tau_fact*becp1%nc(ikb, :, :))
                  ELSE
                     becp2%nc(okb, :, :) = becp2%nc(okb, :, :) &
                       + matmul(u_spin, D(l_i)%d(m_i,m_o, isym0) * tau_fact*CONJG(becp1%nc(ikb, :, :)))
                  ENDIF
               ELSE
                  IF (sgn_sym > 0) THEN
                     becp2%k(okb, :) = becp2%k(okb, :) &
                       + D(l_i)%d(m_i,m_o, isym0) * tau_fact*becp1%k(ikb, :)
                  ELSE
                     becp2%k(okb, :) = becp2%k(okb, :) &
                       + D(l_i)%d(m_i,m_o, isym0) * tau_fact*becp1%k(ikb, :)
                  ENDIF
               END IF
            ENDDO ! m_o
         ENDDO ! ih
      ENDDO ! nat
   END SUBROUTINE
   !
   FUNCTION spinor(r, eig)
      ! Compute a spin-1/2 state polarized along r.
      ! Return spin up state if eig > 0, spin down state otherwise.
      COMPLEX(DP):: spinor(2), ci
      INTEGER:: eig
      REAL(DP):: r(3), theta, phi
      ci = (0.d0, 1.d0)
      theta = acos(r(3)/sqrt(r(1)**2+r(2)**2+r(3)**2))
      phi = atan2(r(2), r(1))
      if(eig > 0) then
         spinor(1) = cos(theta/2)
         spinor(2) = exp( ci * phi ) * sin(theta/2)
      else
         spinor(1) = -exp( -ci * phi ) * sin(theta/2)
         spinor(2) = cos(theta/2)
      end if
      RETURN
   END FUNCTION
   !
   SUBROUTINE pw2wan_set_symm_with_time_reversal(nr1, nr2, nr3, nr1x, nr2x, nr3x)
      ! based on exx_set_symm in exx_base.f90
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x
      !
      ! ... local variables
      !
      INTEGER :: ikq, isym, i,j,k, ri,rj,rk, ir, nxxs
      INTEGER, allocatable :: ftau(:,:), s_scaled(:,:,:)
      !
      nxxs = nr1x*nr2x*nr3x
      !
      ALLOCATE( rir(nxxs,nsym2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating rir', 1)
      !
      rir = 0
      ALLOCATE ( ftau(3,nsym2), s_scaled(3,3,nsym2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ftau/s_scaled', 1)
      CALL scale_sym_ops (nsym2, s2, ft2, nr1, nr2, nr3, s_scaled, ftau)

      DO isym = 1, nsym2
         DO k = 1, nr3
            DO j = 1, nr2
               DO i = 1, nr1
                  CALL rotate_grid_point( s_scaled(1,1,isym), ftau(1,isym), &
                       i, j, k, nr1, nr2, nr3, ri, rj, rk )
                  ir = i + (j-1)*nr1x + (k-1)*nr1x*nr2x
                  rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      DEALLOCATE ( s_scaled, ftau )
      !
   END SUBROUTINE pw2wan_set_symm_with_time_reversal
   !
END SUBROUTINE compute_mmn_ibz

!-----------------------------------------------------------------------
SUBROUTINE compute_spin
   !-----------------------------------------------------------------------
   !
   USE kinds,           ONLY : DP
   USE mp,              ONLY : mp_sum, mp_barrier
   USE mp_world,        ONLY : world_comm
   USE mp_pools,        ONLY : intra_pool_comm, me_pool, root_pool, my_pool_id
   USE io_global,       ONLY : stdout, ionode
   USE wvfct,           ONLY : nbnd, npwx
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions,   ONLY : evc, psic, psic_nc
   USE fft_base,        ONLY : dffts
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE klist,           ONLY : nkstot, xk, ngk, igk_k, nks
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g
   USE cell_base,       ONLY : alat, at, bg
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,       ONLY : tpi
   USE uspp,            ONLY : nkb, vkb, okvan
   USE uspp_param,      ONLY : upf, nh, lmaxq
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE noncollin_module,ONLY : noncolin, npol
   ! begin change Lopez, Thonhauser, Souza
   USE mp,              ONLY : mp_barrier
   USE scf,             ONLY : vrs, vltot, v, kedtau
   USE gvecs,           ONLY : doublegrid
   USE lsda_mod,        ONLY : lsda, nspin, isk
   USE constants,       ONLY : rytoev
   USE uspp_param,      ONLY : upf, nh, nhm
   USE uspp,            ONLY : qq_nt, nhtol, nhtoj, indv
   USE upf_spinorb,     ONLY : fcoef
   USE uspp_init,       ONLY : init_us_2
   USE wannier

   IMPLICIT NONE
   !
   complex(DP), parameter :: cmplx_i=(0.0_DP,1.0_DP)
   !
   INTEGER :: npw, ik, ikp, ipol, ib, i, m, n
   INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, nbt, ibnd_m
   INTEGER :: ikevc, s, counter
   real(DP)                 :: arg, g_(3)
   INTEGER                  :: nn,inn,loop,loop2
   LOGICAL                  :: nn_found
   INTEGER                  :: istart,iend, ierr
   COMPLEX(DP)              :: sigma_x,sigma_y,sigma_z,cdum1,cdum2
   complex(DP), allocatable :: spn(:, :), spn_aug(:, :)
   !
   integer  :: np, is1, is2, kh, kkb
   complex(dp) :: sigma_x_aug, sigma_y_aug, sigma_z_aug
   COMPLEX(DP), ALLOCATABLE :: be_n(:,:), be_m(:,:)
   COMPLEX(DP), ALLOCATABLE :: evc_k(:, :)
   !! Wavefunction at k. Contains only the included bands.
   !
   INTEGER, EXTERNAL :: global_kpoint_index
   !
   IF (.NOT. noncolin) THEN
      WRITE(stdout, *)
      WRITE(stdout, *) " compute_spin is meaningful only for noncollinear cases"
      WRITE(stdout, *) " Doing nothing and return."
      WRITE(stdout, *)
      RETURN
   ENDIF
   !
   CALL start_clock("compute_spin")
   !
   IF (okvan) THEN
      CALL allocate_bec_type ( nkb, nbnd, becp )
      ALLOCATE(be_n(nhm,2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating be_n', 1)
      ALLOCATE(be_m(nhm,2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating be_m', 1)
      ALLOCATE(spn_aug(3, (num_bands*(num_bands+1))/2), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating spn_aug', 1)
   ENDIF
   !
   ALLOCATE(evc_k(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_k', 1)
   ALLOCATE(spn(3, (num_bands*(num_bands+1))/2), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating spn', 1)
   !
   !ivo
   ! not sure this is really needed
   IF (wan_mode=='library') CALL errore('pw2wannier90',&
        'write_spn not meant to work library mode', 1)
   !endivo
   !
   CALL utility_open_output_file("spn", spn_formatted, iun_spn)
   !
   IF (ionode) THEN
      IF (spn_formatted) THEN
         WRITE (iun_spn, *) num_bands, iknum
      ELSE
         WRITE (iun_spn) num_bands, iknum
      ENDIF
   ENDIF
   !
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      !
      CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
      npw = ngk(ik)
      !
      ! Trim excluded bands from evc
      ibnd_m = 0
      DO m = 1, nbnd
         IF (excluded_band(m)) CYCLE
         ibnd_m = ibnd_m + 1
         evc_k(:, ibnd_m) = evc(:, m)
      ENDDO
      !
      !  USPP
      !
      IF(okvan) THEN
         CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
         ! below we compute the product of beta functions with |psi>
         CALL calbec (npw, vkb, evc_k, becp, num_bands)
      ENDIF
      !
      counter = 0
      DO m = 1, num_bands
         DO n = 1, m
            counter = counter + 1
            !
            cdum1 = dot_product(evc_k(1:npw, n), evc_k(npwx+1:npwx+npw, m))
            cdum2 = dot_product(evc_k(npwx+1:npwx+npw, n), evc_k(1:npw, m))
            sigma_x = cdum1 + cdum2
            sigma_y = cmplx_i * (cdum2 - cdum1)
            sigma_z = dot_product(evc_k(1:npw, n), evc_k(1:npw, m)) &
                    - dot_product(evc_k(npwx+1:npwx+npw, n), evc_k(npwx+1:npwx+npw, m))
            spn(1, counter) = sigma_x
            spn(2, counter) = sigma_y
            spn(3, counter) = sigma_z
            !
            IF (okvan) THEN
               !
               sigma_x_aug = (0.0d0, 0.0d0)
               sigma_y_aug = (0.0d0, 0.0d0)
               sigma_z_aug = (0.0d0, 0.0d0)
               ijkb0 = 0
               !
               DO np = 1, ntyp
                  IF (.NOT. upf(np)%tvanp) THEN
                     DO na = 1, nat
                        IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                     ENDDO
                     CYCLE
                  ENDIF
                  !
                  DO na = 1, nat
                     IF (ityp(na) /= np) CYCLE
                     !
                     be_m = 0.d0
                     be_n = 0.d0
                     DO ih = 1, nh(np)
                        ikb = ijkb0 + ih
                        IF (upf(np)%has_so) THEN
                           DO kh = 1, nh(np)
                              IF ((nhtol(kh,np)==nhtol(ih,np)).AND. &
                                   (nhtoj(kh,np)==nhtoj(ih,np)).AND.     &
                                   (indv(kh,np)==indv(ih,np))) THEN
                                 kkb=ijkb0 + kh
                                 DO is1=1,2
                                    DO is2=1,2
                                       be_n(ih,is1)=be_n(ih,is1)+  &
                                          fcoef(ih,kh,is1,is2,np)*  &
                                          becp%nc(kkb,is2,n)
                                       !
                                       be_m(ih,is1)=be_m(ih,is1)+  &
                                          fcoef(ih,kh,is1,is2,np)*  &
                                          becp%nc(kkb,is2,m)
                                    ENDDO
                                 ENDDO
                              ENDIF
                           ENDDO
                        ELSE
                           DO is1 = 1, 2
                              be_n(ih, is1) = becp%nc(ikb, is1, n)
                              be_m(ih, is1) = becp%nc(ikb, is1, m)
                           ENDDO
                        ENDIF
                     ENDDO
                     !
                     DO ih = 1, nh(np)
                        DO jh = 1, nh(np)
                           sigma_x_aug = sigma_x_aug &
                             + qq_nt(ih,jh,np) * ( be_m(jh,2)*conjg(be_n(ih,1))+ be_m(jh,1)*conjg(be_n(ih,2)) )
                           !
                           sigma_y_aug = sigma_y_aug &
                              + qq_nt(ih,jh,np) * (  &
                              be_m(jh,1) * conjg(be_n(ih,2)) &
                              - be_m(jh,2) * conjg(be_n(ih,1)) &
                            ) * (0.0d0, 1.0d0)
                           !
                           sigma_z_aug = sigma_z_aug &
                              + qq_nt(ih,jh,np) * ( be_m(jh,1)*conjg(be_n(ih,1)) - be_m(jh,2)*conjg(be_n(ih,2)) )
                        ENDDO
                     ENDDO
                     !
                     ijkb0 = ijkb0 + nh(np)
                  ENDDO
               ENDDO
               spn_aug(1, counter) = sigma_x_aug
               spn_aug(2, counter) = sigma_y_aug
               spn_aug(3, counter) = sigma_z_aug
            ENDIF ! okvan
         ENDDO ! n
      ENDDO ! m
      !
      CALL mp_sum(spn, intra_pool_comm)
      IF (okvan) spn = spn + spn_aug
      !
      ! Write to file
      !
      IF (me_pool == root_pool) THEN
         CALL utility_write_array(iun_spn, spn_formatted, 3, ((num_bands*(num_bands+1))/2), spn)
      ENDIF
      !
   ENDDO ! ik
   !
   IF (me_pool == root_pool) CLOSE (iun_spn, STATUS="KEEP")
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   !
   CALL utility_merge_files("spn", spn_formatted, 3*((num_bands*(num_bands+1))/2))
   !
   DEALLOCATE(evc_k)
   DEALLOCATE(spn)
   IF (okvan) THEN
      DEALLOCATE(spn_aug)
      DEALLOCATE(be_n, be_m)
      CALL deallocate_bec_type(becp)
   ENDIF
   !
   WRITE(stdout,*) ' SPIN calculated'
   !
   CALL stop_clock("compute_spin")
   !
END SUBROUTINE compute_spin

!--------------------------------------------------------------------------
SUBROUTINE compute_orb
   !-----------------------------------------------------------------------
   !!
   !! Calculate and write uHu and uIu matrix, which are defined as follows.
   !! uHu(n, m, ib1, ib2, ik) = <u_{m,k+b1}|H(k)|u_{n,k+b2}>
   !! uIu(n, m, ib1, ib2, ik) = <u_{m,k+b1}|u_{n,k+b2}>
   !!
   !! Note that index n, not m, is the fastest index. This non-natural
   !! indexing is used just to be consistent with the previous versions.
   !!
   !-----------------------------------------------------------------------
   !
   USE kinds,           ONLY : DP
   USE io_global,       ONLY : stdout, ionode
   USE mp,              ONLY : mp_sum, mp_barrier
   USE mp_world,        ONLY : world_comm
   USE mp_pools,        ONLY : intra_pool_comm, my_pool_id, me_pool, root_pool
   USE wvfct,           ONLY : npwx, current_k
   USE control_flags,   ONLY : gamma_only
   USE fft_base,        ONLY : dfftp
   USE klist,           ONLY : xk, ngk, igk_k, nks
   USE gvect,           ONLY : g, gstart
   USE ions_base,       ONLY : ntyp => nsp
   USE uspp,            ONLY : nkb, vkb, okvan
   USE uspp_param,      ONLY : upf
   USE becmod,          ONLY : becp, allocate_bec_type, deallocate_bec_type, calbec
   USE noncollin_module,ONLY : noncolin, npol
   USE wannier
   ! begin change Lopez, Thonhauser, Souza
   USE scf,             ONLY : vrs, vltot, v, kedtau
   USE gvecs,           ONLY : doublegrid
   USE lsda_mod,        ONLY : lsda, nspin, isk, current_spin
   USE constants,       ONLY : rytoev
   USE uspp_init,       ONLY : init_us_2
   USE basis,           ONLY : natomwfc, swfcatom
   USE ldaU,            ONLY : lda_plus_u, copy_U_wfc
   !
   IMPLICIT NONE
   !
   INTEGER                  :: ik, npw, m, n, ierr
   INTEGER                  :: i_b, i_b1, i_b2
   COMPLEX(DP), ALLOCATABLE :: evc_b(:, :, :)
   !! Wavefunction at k+b for all b vectors
   COMPLEX(DP), ALLOCATABLE :: evc_b1(:,:)
   !! Wavefunction at k+b1
   COMPLEX(DP), ALLOCATABLE :: evc_b2(:,:)
   !! Wavefunction at k+b2
   COMPLEX(DP), ALLOCATABLE :: H_evc_b2(:, :)
   !! H times Wavefunction at k+b2
   COMPLEX(DP), ALLOCATABLE :: uHu(:, :, :, :)
   !! Computed uHu matrix
   COMPLEX(DP), ALLOCATABLE :: uIu(:, :, :, :)
   !! Computed uIu matrix
   COMPLEX(DP), ALLOCATABLE :: arr(:, :)
   !! Temporary storage
   COMPLEX(DP), ALLOCATABLE :: wfcatom(:, :)
   !! For DFT+U calculations
   !
   INTEGER, EXTERNAL :: global_kpoint_index
   !
   IF (.NOT. (write_uHu .OR. write_uIu)) THEN
      WRITE(stdout, *)
      WRITE(stdout, *) ' ----------------------------------------'
      WRITE(stdout, *) ' *** uHu and uIu matrices are not computed '
      WRITE(stdout, *) ' ----------------------------------------'
      WRITE(stdout, *)
      RETURN
   ENDIF
   !
   CALL start_clock('compute_orb')
   !
   IF(gamma_only) CALL errore('pw2wannier90',&
        'write_uHu and write_uIu not yet implemented for gamma_only case',1) !ivo
   IF(okvan) CALL errore('pw2wannier90',&
        'write_uHu and write_uIu not yet implemented with USPP',1) !ivo
!ivo
! not sure this is really needed
   if((write_uhu.or.write_uIu).and.wan_mode=='library')&
        call errore('pw2wannier90',&
        'write_uhu, and write_uIu not meant to work library mode',1)
!endivo
   !
   ALLOCATE(evc_b(npol*npwx, num_bands, nnb), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_b', 1)
   ALLOCATE(evc_b1(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_b1', 1)
   ALLOCATE(evc_b2(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_b2', 1)
   ALLOCATE(arr(num_bands, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating arr', 1)
   IF (lda_plus_u) THEN
      ALLOCATE(wfcatom(npwx*npol, natomwfc), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wfcatom', 1)
      ALLOCATE(swfcatom(npwx*npol, natomwfc), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating swfcatom', 1)
   ENDIF
   !
   IF (write_uHu) THEN
      ALLOCATE(H_evc_b2(npol*npwx, num_bands), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating H_evc_b2', 1)
      ALLOCATE(uHu(num_bands, num_bands, nnb, nnb), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating uHu', 1)
      uHu = (0.d0, 0.d0)
   ENDIF
   !
   IF (write_uIu) THEN
      ALLOCATE(uIu(num_bands, num_bands, nnb, nnb), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating uIu', 1)
      uIu = (0.d0, 0.d0)
   ENDIF
   !
   !====================================================================
   !
   ! The following code was inserted by Timo Thonhauser, Ivo Souza, and
   ! Graham Lopez in order to calculate the matrix elements
   ! <u_n(q+b1)|H(q)|u_m(q+b2)> necessary for the Wannier interpolation
   ! of the orbital magnetization
   !
   !====================================================================
   !
   IF (write_uHu) THEN
      WRITE(stdout, *) ' *** Compute  uHu '
      CALL utility_open_output_file("uHu", uHu_formatted, iun_uHu)
      IF (ionode) THEN
         IF (uHu_formatted) THEN
            WRITE(iun_uHu, *) num_bands, iknum, nnb
         ELSE
            WRITE(iun_uHu) num_bands, iknum, nnb
         ENDIF
      ENDIF
   ENDIF ! write_uHu
   !
   IF (write_uIu) THEN
      WRITE(stdout, *) ' *** Compute  uIu '
      CALL utility_open_output_file("uIu", uIu_formatted, iun_uIu)
      IF (ionode) THEN
         IF (uIu_formatted) THEN
            WRITE(iun_uIu, *) num_bands, iknum, nnb
         ELSE
            WRITE(iun_uIu) num_bands, iknum, nnb
         ENDIF
      ENDIF
   ENDIF ! write_uIu
   !
   CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
   CALL allocate_bec_type(nkb, num_bands, becp)
   !
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks ! loop over k points
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      !
      npw = ngk(ik)
      !
      ! sort the wfc at k and set up stuff for h_psi
      current_k = ik
      IF (lsda) current_spin = isk(ik)
      CALL init_us_2(npw, igk_k(1, ik), xk(1, ik), vkb)
      CALL g2_kin(ik)
      !
      IF (lda_plus_u) THEN
         IF (noncolin) THEN
            CALL atomic_wfc_nc_updown(ik, wfcatom)
         ELSE
            CALL atomic_wfc(ik, wfcatom)
         END IF
         CALL calbec(npw, vkb, wfcatom, becp)
         CALL s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)
         CALL copy_U_wfc(swfcatom, noncolin)
      END IF
      !
      ! Loop over the nearest neighbor b vectors and compute |u_{n,k+b}> and
      ! save them to evc_b.
      !
      DO i_b = 1, nnb
         !
         CALL utility_compute_u_kb(ik, i_b, evc_b(:, :, i_b))
         !
      ENDDO ! nnb
      !
      ! Main loop to compute uHu and uIu.
      !
      DO i_b2 = 1, nnb
         !
         ! Copy |u_{n,k+b2}> from evc_b
         evc_b2 = evc_b(:, :, i_b2)
         !
         ! Compute H(k) * |u_{n,k+b2}>
         IF (write_uHu) CALL h_psi(npwx, npw, num_bands, evc_b2, H_evc_b2)
         !
         ! Compute uHu and uIu only if i_b1 >= i_b2. For the other cases, use
         ! <u_{m,k+b1}|u_{n,k+b2}> = [<u_{n,k+b2}|u_{m,k+b1}>]*
         !
         DO i_b1 = i_b2, nnb
            !
            ! Copy |u_{n,k+b1}> from evc_b
            evc_b1 = evc_b(:, :, i_b1)
            !
            ! To be consistent with previous versions, we compute
            ! uHu(n, m) = <evc_b1(m)|H_evc_b2(n)>
            ! uIu(n, m) = <evc_b1(m)|evc_b2(n)>
            !
            IF (write_uHu) THEN
               IF (noncolin) THEN
                  CALL ZGEMM('C', 'N', num_bands, num_bands, npwx*npol, &
                     (1.d0, 0.d0), evc_b1, npwx*npol, H_evc_b2, npwx*npol, &
                     (0.d0, 0.d0), arr, num_bands)
               ELSE
                  CALL ZGEMM('C', 'N', num_bands, num_bands, npw, &
                     (1.d0, 0.d0), evc_b1, npwx, H_evc_b2, npwx, &
                     (0.d0, 0.d0), arr, num_bands)
               ENDIF
               !
               ! we need uHu(n,m) = <evc_b1(m)|H_evc_b2(n)>.
               uHu(:, :, i_b1, i_b2) = TRANSPOSE(arr)
               !
               ! Fill in the case i_b1 < i_b2
               IF (i_b1 > i_b2) uHu(:, :, i_b2, i_b1) = CONJG(arr)
               !
            ENDIF
            !
            IF (write_uIu) THEN
               IF (noncolin) THEN
                  CALL ZGEMM('C', 'N', num_bands, num_bands, npwx*npol, &
                     (1.d0, 0.d0), evc_b1, npwx*npol, evc_b2, npwx*npol, &
                     (0.d0, 0.d0), arr, num_bands)
               ELSE
                  CALL ZGEMM('C', 'N', num_bands, num_bands, npw, &
                     (1.d0, 0.d0), evc_b1, npwx, evc_b2, npwx, &
                     (0.d0, 0.d0), arr, num_bands)
               ENDIF
               !
               ! we need uIu(n, m) = <evc_b1(m)|evc_b2(n)>.
               uIu(:, :, i_b1, i_b2) = TRANSPOSE(arr)
               !
               ! Fill in the case i_b1 < i_b2
               IF (i_b1 > i_b2) uIu(:, :, i_b2, i_b1) = CONJG(arr)
               !
            ENDIF
            !
         ENDDO ! i_b1
      ENDDO ! i_b2
      !
      IF (write_uHu) uHu = uHu * rytoev ! because wannier90 works in eV
      IF (write_uHu) CALL mp_sum(uHu, intra_pool_comm)
      IF (write_uIu) CALL mp_sum(uIu, intra_pool_comm)
      !
      ! write the files out to disk
      !
      DO i_b2 = 1, nnb
         DO i_b1 = 1, nnb
            IF (me_pool == root_pool) THEN
               IF (write_uHu) THEN
                  CALL utility_write_array(iun_uHu, uHu_formatted, num_bands, num_bands, uHu(:, :, i_b1, i_b2))
               ENDIF
               IF (write_uIu) THEN
                  CALL utility_write_array(iun_uIu, uIu_formatted, num_bands, num_bands, uIu(:, :, i_b1, i_b2))
               ENDIF
            ENDIF ! root_pool
         ENDDO ! i_b1
      ENDDO ! i_b2
   ENDDO ! ik
   !
   IF (me_pool == root_pool) THEN
      IF (write_uHu) CLOSE(iun_uHu, STATUS="KEEP")
      IF (write_uIu) CLOSE(iun_uIu, STATUS="KEEP")
   ENDIF
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   !
   IF (write_uHu) CALL utility_merge_files("uHu", uHu_formatted, num_bands**2)
   IF (write_uIu) CALL utility_merge_files("uIu", uIu_formatted, num_bands**2)
   !
   ! Deallocate variables
   !
   CALL deallocate_bec_type(becp)
   DEALLOCATE(evc_b1)
   DEALLOCATE(evc_b2)
   DEALLOCATE(evc_b)
   IF (lda_plus_u) THEN
      DEALLOCATE(wfcatom)
      DEALLOCATE(swfcatom)
   ENDIF
   !
   IF (write_uHu) THEN
      DEALLOCATE(H_evc_b2)
      DEALLOCATE(uHu)
   ENDIF ! write_uHu
   !
   IF (write_uIu) THEN
      DEALLOCATE(uIu)
   ENDIF ! write_uIu
   !
   IF (write_uHu) WRITE(stdout, *) ' uHu calculated'
   IF (write_uIu) WRITE(stdout, *) ' uIu calculated'
   !
   CALL stop_clock('compute_orb')
   !
END SUBROUTINE compute_orb
!
!-------------------------------------------------------------------------
SUBROUTINE compute_shc
   !-----------------------------------------------------------------------
   !!
   !! Calculate and write sHu and sIu matrix, which are defined as follows.
   !! sHu(n, m, ispol, ib, ik) = <u_{m,k}|S(ispol) H(k)|u_{n,k+b}>
   !! sIu(n, m, ispol, ib, ik) = <u_{m,k}|S(ispol)|u_{n,k+b}>
   !!
   !! Note that index n, not m, is the fastest index. This non-natural
   !! indexing is used just to be consistent with the previous versions.
   !!
   !-----------------------------------------------------------------------
   !
   USE kinds,           ONLY : DP
   USE mp,              ONLY : mp_sum, mp_barrier
   USE mp_world,        ONLY : world_comm
   USE mp_pools,        ONLY : intra_pool_comm, my_pool_id, me_pool, root_pool
   USE io_global,       ONLY : stdout, ionode
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE constants,       ONLY : rytoev
   USE fft_base,        ONLY : dfftp
   USE control_flags,   ONLY : gamma_only
   USE wvfct,           ONLY : nbnd, npwx, current_k
   USE wavefunctions,   ONLY : evc
   USE klist,           ONLY : xk, ngk, igk_k, nks
   USE uspp,            ONLY : nkb, vkb, okvan
   USE uspp_param,      ONLY : upf
   USE becmod,          ONLY : bec_type, becp, allocate_bec_type, &
                               deallocate_bec_type, calbec
   USE gvecs,           ONLY : doublegrid
   USE noncollin_module,ONLY : noncolin, npol
   USE lsda_mod,        ONLY : nspin, lsda, isk
   USE scf,             ONLY : vrs, vltot, v, kedtau
   USE uspp_init,       ONLY : init_us_2
   USE basis,           ONLY : natomwfc, swfcatom
   USE ldaU,            ONLY : lda_plus_u, copy_U_wfc
   USE wannier
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), parameter :: cmplx_i = (0.0_DP, 1.0_DP)
   !
   INTEGER :: ik, npw, m, n, ibnd_m, ispol, i_b, ierr
   COMPLEX(DP) :: sigma_x, sigma_y, sigma_z, cdum1, cdum2
   !
   COMPLEX(DP), ALLOCATABLE :: evc_kb(:, :), H_evc_kb(:, :)
   COMPLEX(DP), ALLOCATABLE :: sHu(:, : ,:), sIu(:, :, :)
   COMPLEX(DP), ALLOCATABLE :: evc_k(:, :)
   !! Wavefunction at k. Contains only the included bands.
   COMPLEX(DP), ALLOCATABLE :: wfcatom(:, :)
   !! For DFT+U calculations
   !
   IF (.NOT. (write_sHu .OR. write_sIu)) THEN
      WRITE(stdout, *)
      WRITE(stdout, *) ' ----------------------------------------'
      WRITE(stdout, *) ' *** sHu and sIu matrix are not computed '
      WRITE(stdout, *) ' ----------------------------------------'
      WRITE(stdout, *)
      !
      RETURN
      !
   ENDIF
   !
   CALL start_clock('compute_shc')
   !
   IF (wan_mode == 'library') CALL errore('pw2wannier90', &
      'write_sHu and write_sIu not meant to work with library mode', 1)
   IF (gamma_only) CALL errore('pw2wannier90', &
      'write_sHu and write_sIu not yet implemented for gamma_only case', 1)
   IF (okvan) CALL errore('pw2wannier90', &
      'write_sHu and write_sIu not yet implemented with USPP', 1)
   IF (.NOT. noncolin) CALL errore('pw2wannier90', &
      'write_sHu and write_sIu only works with noncolin == .true.', 1)
   !
   ALLOCATE(evc_k(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_k', 1)
   ALLOCATE(evc_kb(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_kb', 1)
   !
   IF (write_sHu) THEN
      ALLOCATE(sHu(num_bands, num_bands, 3), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating sHu', 1)
      ALLOCATE(H_evc_kb(npol*npwx, num_bands), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating H_evc_kb', 1)
      IF (lda_plus_u) THEN
         ALLOCATE(wfcatom(npwx*npol, natomwfc), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wfcatom', 1)
         ALLOCATE(swfcatom(npwx*npol, natomwfc), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating swfcatom', 1)
      ENDIF
   ENDIF
   IF (write_sIu) THEN
      ALLOCATE(sIu(num_bands, num_bands, 3), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating sIu', 1)
   ENDIF
   !
   !
   IF (write_sHu) THEN
      write(stdout,*) ' *** Computing  sHu '
      CALL utility_open_output_file("sHu", sHu_formatted, iun_sHu)
      IF (ionode) THEN
         IF (sHu_formatted) THEN
            WRITE(iun_sHu, *) num_bands, iknum, nnb
         ELSE
            WRITE(iun_sHu) num_bands, iknum, nnb
         ENDIF
      ENDIF
   ENDIF
   !
   IF (write_sIu) THEN
      WRITE(stdout,*) ' *** Computing  sIu '
      CALL utility_open_output_file("sIu", sIu_formatted, iun_sIu)
      IF (ionode) THEN
         IF (sIu_formatted) THEN
            WRITE(iun_sIu, *) num_bands, iknum, nnb
         ELSE
            WRITE(iun_sIu) num_bands, iknum, nnb
         ENDIF
      ENDIF
   ENDIF
   !
   CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
   CALL allocate_bec_type(nkb, num_bands, becp)
   !
   WRITE(stdout, *)
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks ! loop over k points
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      !
      npw = ngk(ik)
      !
      ! Read wavefunctions at k, exclude the excluded bands
      !
      CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1 )
      !
      ibnd_m = 0
      DO m = 1, nbnd
         IF (excluded_band(m)) CYCLE
         ibnd_m = ibnd_m + 1
         evc_k(:, ibnd_m) = evc(:, m)
      ENDDO
      !
      ! sort the wfc at k and set up stuff for h_psi
      current_k = ik
      CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
      CALL g2_kin(ik)
      !
      IF (lda_plus_u) THEN
         IF (noncolin) THEN
            CALL atomic_wfc_nc_updown(ik, wfcatom)
         ELSE
            CALL atomic_wfc(ik, wfcatom)
         END IF
         CALL calbec(npw, vkb, wfcatom, becp)
         CALL s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)
         CALL copy_U_wfc(swfcatom, noncolin)
      END IF
      !
      DO i_b = 1, nnb ! nnb = # of nearest neighbors
         !
         ! compute  |u_{n,k+b2}> and H(k) * |u_{n,k+b2}>
         !
         CALL utility_compute_u_kb(ik, i_b, evc_kb)
         !
         IF (write_sHu) CALL h_psi(npwx, npw, num_bands, evc_kb, H_evc_kb)
         !
         sHu = (0.D0, 0.D0)
         sIu = (0.D0, 0.D0)
         !
         ! loop on bands
         DO m = 1, num_bands
            !
            DO n = 1, num_bands
               !
               ! <a|sx|b> = (a2, b1) + (a1, b2)
               ! <a|sy|b> = I (a2, b1) - I (a1, b2)
               ! <a|sz|b> = (a1, b1) - (a2, b2)
               !
               IF (write_sHu) THEN !ivo
                  cdum1 = dot_product(evc_k(1:npw, m), H_evc_kb(npwx+1:npwx+npw, n))
                  cdum2 = dot_product(evc_k(npwx+1:npwx+npw, m), H_evc_kb(1:npw, n))
                  sigma_x = cdum1 + cdum2
                  sigma_y = cmplx_i * (cdum2 - cdum1)
                  sigma_z = dot_product(evc_k(1:npw, m), H_evc_kb(1:npw, n)) &
                          - dot_product(evc_k(npwx+1:npwx+npw, m), H_evc_kb(npwx+1:npwx+npw, n))
                  !
                  sHu(n, m, 1) = sigma_x
                  sHu(n, m, 2) = sigma_y
                  sHu(n, m, 3) = sigma_z
               ENDIF ! write_sHu
               !
               IF (write_sIu) THEN !ivo
                  cdum1 = dot_product(evc_k(1:npw, m), evc_kb(npwx+1:npwx+npw, n))
                  cdum2 = dot_product(evc_k(npwx+1:npwx+npw, m), evc_kb(1:npw, n))
                  sigma_x = cdum1 + cdum2
                  sigma_y = cmplx_i * (cdum2 - cdum1)
                  sigma_z = dot_product(evc_k(1:npw, m), evc_kb(1:npw, n)) &
                          - dot_product(evc_k(npwx+1:npwx+npw, m), evc_kb(npwx+1:npwx+npw, n))
                  !
                  sIu(n, m, 1) = sigma_x
                  sIu(n, m, 2) = sigma_y
                  sIu(n, m, 3) = sigma_z
               ENDIF ! write_sIu
            ENDDO ! n
            !
         ENDDO ! m
         !
         IF (write_sHu) sHu = sHu * rytoev
         IF (write_sHu) CALL mp_sum(sHu, intra_pool_comm)
         IF (write_sIu) CALL mp_sum(sIu, intra_pool_comm)
         !
         IF (me_pool == root_pool) THEN  ! write the files out to disk
            DO ispol = 1, 3
               IF (write_sHu) THEN
                  CALL utility_write_array(iun_sHu, sHu_formatted, num_bands, num_bands, sHu(:, :, ispol))
               ENDIF
               IF (write_sIu) THEN
                  CALL utility_write_array(iun_sIu, sIu_formatted, num_bands, num_bands, sIu(:, :, ispol))
               ENDIF
            ENDDO
         ENDIF ! end of io
         !
      ENDDO ! i_b
   ENDDO ! ik
   !
   IF (me_pool == root_pool) THEN
      IF (write_sHu) CLOSE(iun_sHu, STATUS="KEEP")
      IF (write_sIu) CLOSE(iun_sIu, STATUS="KEEP")
   ENDIF
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   !
   IF (write_sHu) CALL utility_merge_files("sHu", sHu_formatted, num_bands**2)
   IF (write_sIu) CALL utility_merge_files("sIu", sIu_formatted, num_bands**2)
   !
   DEALLOCATE(evc_k)
   DEALLOCATE(evc_kb)
   IF (write_sHu) THEN
      DEALLOCATE(H_evc_kb)
      DEALLOCATE(sHu)
      IF (lda_plus_u) THEN
         DEALLOCATE(wfcatom)
         DEALLOCATE(swfcatom)
      ENDIF
   ENDIF
   IF (write_sIu) DEALLOCATE(sIu)
   !
   WRITE(stdout,*) ' shc calculated'
   !
   CALL stop_clock('compute_shc')
   !
   RETURN
   !
END SUBROUTINE

!--------------------------------------------------------------------------
SUBROUTINE utility_write_array(iun, formatted, ndim1, ndim2, arr)
   !------------------------------------------------------------------------
   !!
   !! Write a (ndim1, ndim2) array arr to file unit iun.
   !!
   !-----------------------------------------------------------------------
   !
   USE kinds, ONLY : DP
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iun
   !! Unit of the file to write. Must be already opened.
   LOGICAL, INTENT(IN) :: formatted
   !! If true, write formatted. If false, write unformatted.
   INTEGER, INTENT(IN) :: ndim1
   !! First dimension of the array.
   INTEGER, INTENT(IN) :: ndim2
   !! Second dimension of the array.
   COMPLEX(DP), INTENT(IN) :: arr(ndim1, ndim2)
   !! Array to write to file.
   !
   INTEGER :: m
   !! Array index
   INTEGER :: n
   !! Array index
   !
   IF (formatted) THEN
      ! Write formatted file. Slow bulky way for transferable files
      DO n = 1, ndim2
         DO m = 1, ndim1
            WRITE(iun, '(2ES20.10)') arr(m, n)
         ENDDO
      ENDDO
   ELSE
      ! Write unformatted file. The fast way
      WRITE(iun) ((arr(m, n), m=1,ndim1), n=1,ndim2)
   ENDIF
   !
!--------------------------------------------------------------------------
END SUBROUTINE utility_write_array
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
SUBROUTINE utility_open_output_file(postfix, formatted, iun)
   !-----------------------------------------------------------------------
   !! In ionode, open file seedname.postfix for output and write header line.
   !! For pool parallelization, the root of each pool opens file
   !! seedname.postfixID, where ID is the 1-based id of the pool.
   !-----------------------------------------------------------------------
   USE kinds,           ONLY : DP
   USE io_global,       ONLY : ionode
   USE mp_pools,        ONLY : my_pool_id, me_pool, root_pool
   USE wannier,         ONLY : seedname
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN) :: postfix
   !! postfix for filename
   LOGICAL, INTENT(IN) :: formatted
   !! True if formatted file, false if unformatted file.
   INTEGER, INTENT(OUT) :: iun
   !! Unit of the file
   !
   CHARACTER(LEN=9) :: cdate, ctime
   CHARACTER(LEN=60) :: header
   CHARACTER(LEN=256) :: filename
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   !
   IF (me_pool == root_pool) THEN
      filename = TRIM(seedname) // "." // TRIM(postfix)
      IF (.NOT. ionode) filename = TRIM(filename) // TRIM(int_to_char(my_pool_id+1))
      IF (formatted) THEN
         OPEN(NEWUNIT=iun, FILE=TRIM(filename), FORM='formatted', STATUS="REPLACE")
      ELSE
         OPEN(NEWUNIT=iun, FILE=TRIM(filename), FORM='unformatted', STATUS="REPLACE")
      ENDIF
   ENDIF
   !
   IF (ionode) THEN
      CALL date_and_tim(cdate, ctime)
      header = 'Created on ' // cdate // ' at ' // ctime
      IF (formatted) THEN
         WRITE(iun, *) header
      ELSE
         WRITE(iun) header
      ENDIF
   ENDIF
!--------------------------------------------------------------------------
END SUBROUTINE utility_open_output_file
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
SUBROUTINE compute_amn
   !-----------------------------------------------------------------------
   !!
   !! In the collinear case, n_proj and n_wannier are always the same.
   !! In the noncollinear case,
   !! 1) standalone mode: n_proj is the number of spinor projections.
   !!                     n_wannier = n_proj, but is never used here.
   !! 2) library mode: n_proj is the number of scaler projections.
   !!                  n_wannier = 2 * n_proj
   !! (For the standalone mode, n_wannier used only in SCDM projections.)
   !!
   !! library mode for noncollinear spinor projection is not working.
   !!
   !! nocolin: we have half as many projections g(r) defined as wannier
   !!          functions. We project onto (1,0) (ie up spin) and then onto
   !!          (0,1) to obtain num_wann projections. jry
   !!
   !-----------------------------------------------------------------------
   !
   USE kinds,           ONLY : DP
   USE io_global,       ONLY : stdout, ionode
   USE mp,              ONLY : mp_sum, mp_barrier
   USE mp_world,        ONLY : world_comm
   USE mp_pools,        ONLY : intra_pool_comm, me_pool, root_pool, my_pool_id
   USE klist,           ONLY : nkstot, xk, ngk, igk_k, nks
   USE wvfct,           ONLY : nbnd, npwx
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions,   ONLY : evc
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g, ngm, gstart
   USE uspp,            ONLY : nkb, vkb, okvan
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE uspp_param,      ONLY : upf
   USE noncollin_module,ONLY : noncolin, npol
   USE lsda_mod,        ONLY : lsda, isk
   USE constants,       ONLY : eps6
   USE uspp_init,       ONLY : init_us_2
   USE wannier
   !
   IMPLICIT NONE
   !
   COMPLEX(DP) :: amn_tmp, fac(2)
   COMPLEX(DP), ALLOCATABLE :: amn(:, :)
   COMPLEX(DP), ALLOCATABLE :: sgf(:,:)
   COMPLEX(DP), ALLOCATABLE :: evc_k(:, :)
   !! Wavefunction at k. Contains only the included bands.
   INTEGER :: ik, npw, ibnd, ibnd1, iw, i, nt, ipol, ik_g_w90
   LOGICAL            :: opnd, exst,spin_z_pos, spin_z_neg
   INTEGER            :: istart, ierr
   !
   INTEGER, EXTERNAL :: global_kpoint_index
   REAL(DP), EXTERNAL :: ddot
   !
   CALL start_clock( 'compute_amn' )
   !
   ALLOCATE(amn(num_bands, n_proj), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating amn', 1)
   ALLOCATE(sgf(npwx,n_proj), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating sgf', 1)
   IF (noncolin) THEN
      ALLOCATE(gf_spinor(2*npwx, n_proj), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating gf_spinor', 1)
      ALLOCATE(sgf_spinor(2*npwx, n_proj), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating sgf_spinor', 1)
   ENDIF
   ALLOCATE(evc_k(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_k', 1)
   !
   IF (wan_mode=='library') THEN
      ALLOCATE(a_mat(num_bands, n_wannier, iknum), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating a_mat', 1)
   ENDIF
   !
   IF (wan_mode=='standalone') THEN
      IF (irr_bz) THEN
         CALL utility_open_output_file("iamn", .TRUE., iun_amn)
      ELSE
         CALL utility_open_output_file("amn", .TRUE., iun_amn)
      ENDIF
      IF (ionode) WRITE(iun_amn, *) num_bands, iknum, n_proj
   ENDIF
   !
   IF (okvan) CALL allocate_bec_type(nkb, n_proj, becp)
   !
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      !
      ik_g_w90 = global_kpoint_index(nkstot, ik) - ikstart + 1
      npw = ngk(ik)
      !
      ! Read wavefunctions at k, exclude the excluded bands
      !
      CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1 )
      !
      ibnd1 = 0
      DO ibnd = 1, nbnd
         IF (excluded_band(ibnd)) CYCLE
         ibnd1 = ibnd1 + 1
         evc_k(:, ibnd1) = evc(:, ibnd)
      ENDDO
      !
      CALL generate_guiding_functions(ik)   ! they are called gf(npw,n_proj)
      !
      IF (noncolin) THEN
        sgf_spinor = (0.d0,0.d0)
        CALL orient_gf_spinor(npw)
      ENDIF
      !
      !  USPP
      !
      IF (okvan) THEN
         CALL init_us_2 (npw, igk_k(1, ik), xk(1, ik), vkb)
         !
         ! below we compute the product of beta functions with trial func.
         IF (gamma_only) THEN
            CALL calbec ( npw, vkb, gf, becp, n_proj )
         ELSEIF (noncolin) THEN
            CALL calbec ( npw, vkb, gf_spinor, becp, n_proj )
         ELSE
            CALL calbec ( npw, vkb, gf, becp, n_proj )
         ENDIF
         !
         ! and we use it for the product S|trial_func>
         IF (noncolin) THEN
           CALL s_psi (npwx, npw, n_proj, gf_spinor, sgf_spinor)
         ELSE
           CALL s_psi (npwx, npw, n_proj, gf, sgf)
         ENDIF
      ELSE
         sgf(:,:) = gf(:,:)
      ENDIF
      !
      amn = (0.0_dp, 0.0_dp)
      !
      IF (noncolin) THEN
         DO iw = 1, n_proj
            spin_z_pos = .FALSE.
            spin_z_neg = .FALSE.
            !
            ! detect if spin quantisation axis is along z
            IF((abs(spin_qaxis(1,iw)-0.0d0)<eps6).and.(abs(spin_qaxis(2,iw)-0.0d0)<eps6) &
                 .and.(abs(spin_qaxis(3,iw)-1.0d0)<eps6)) then
               spin_z_pos = .TRUE.
            ELSEIF(abs(spin_qaxis(1,iw)-0.0d0)<eps6.and.abs(spin_qaxis(2,iw)-0.0d0)<eps6 &
                 .and.abs(spin_qaxis(3,iw)+1.0d0)<eps6) THEN
               spin_z_neg = .TRUE.
            ENDIF
            !
            IF (spin_z_pos .OR. spin_z_neg) THEN
               !
               IF (spin_z_pos) THEN
                  ipol = (3-spin_eig(iw))/2
               ELSE
                  ipol = (3+spin_eig(iw))/2
               ENDIF
               istart = (ipol-1)*npwx + 1
               !
               IF (okvan) THEN
                  CALL ZGEMV('C', npw, num_bands, (1.d0, 0.d0), evc_k, &
                     npol*npwx, sgf_spinor(1, iw), 1, (0.d0, 0.d0), amn(1, iw), 1)
                  CALL ZGEMV('C', npw, num_bands, (1.d0, 0.d0), evc_k(npwx+1, 1), &
                     npol*npwx, sgf_spinor(npwx+1, iw), 1, (1.d0, 0.d0), amn(1, iw), 1)
               ELSE
                  CALL ZGEMV('C', npw, num_bands, (1.d0, 0.d0), evc_k(istart, 1), &
                     npol*npwx, sgf(1, iw), 1, (0.d0, 0.d0), amn(1, iw), 1)
               ENDIF
            ELSE ! .NOT. (spin_z_pos .OR. spin_z_neg)
               ! general routine
               ! for quantisation axis (a,b,c)
               ! 'up'    eigenvector is 1/sqrt(1+c) [c+1,a+ib]
               ! 'down'  eigenvector is 1/sqrt(1-c) [c-1,a+ib]
               IF (spin_eig(iw)==1) THEN
                  fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw)+1)*cmplx(1.0d0,0.0d0,dp)
                  fac(2)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
               ELSE
                  fac(1)=(1.0_dp/sqrt(1-spin_qaxis(3,iw)))*(spin_qaxis(3,iw)-1)*cmplx(1.0d0,0.0d0,dp)
                  fac(2)=(1.0_dp/sqrt(1-spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
               ENDIF
               !
               DO ipol = 1, npol
                  istart = (ipol-1)*npwx + 1
                  IF (okvan) THEN
                     CALL ZGEMV('C', npw, num_bands, (1.d0, 0.d0), evc_k(istart, 1), &
                        npol*npwx, sgf_spinor(istart, iw), 1, (1.d0, 0.d0), amn(1, iw), 1)
                  ELSE
                     CALL ZGEMV('C', npw, num_bands, fac(ipol), evc_k(istart, 1), &
                        npol*npwx, sgf(1, iw), 1, (1.d0, 0.d0), amn(1, iw), 1)
                  ENDIF
               ENDDO
            ENDIF ! spin_z_pos .OR. spin_z_neg
         ENDDO ! iw
      ELSE ! .NOT. noncolin
         IF (gamma_only) THEN
            DO iw = 1, n_proj
               DO ibnd = 1, num_bands
                  amn_tmp = 2.0_dp*ddot(2*npw,evc_k(1,ibnd),1,sgf(1,iw),1)
                  IF (gstart==2) amn_tmp = amn_tmp - real(conjg(evc_k(1,ibnd))*sgf(1,iw))
                  amn(ibnd, iw) = amn_tmp
               ENDDO ! ibnd
            ENDDO ! iw
         ELSE
            CALL ZGEMM('C', 'N', num_bands, n_proj, npw, &
               (1.d0, 0.d0), evc_k, npwx, sgf, npwx, &
               (0.d0, 0.d0), amn, num_bands)
         ENDIF ! gamma_only
      ENDIF ! noncolin
      !
      CALL mp_sum(amn, intra_pool_comm)
      !
      ! Write amn to file or save to a_mat
      !
      IF (wan_mode == 'standalone') THEN
         DO iw = 1, n_proj
            DO ibnd = 1, num_bands
               IF (me_pool == root_pool) WRITE(iun_amn,'(3i10,2f18.12)') &
                  ibnd, iw, ik_g_w90, amn(ibnd, iw)
            ENDDO
         ENDDO
      ELSEIF (wan_mode=='library') THEN
         DO iw = 1, n_proj
            DO ibnd = 1, num_bands
               a_mat(ibnd, iw, ik_g_w90) = amn(ibnd, iw)
            ENDDO
         ENDDO
      ENDIF
      !
   ENDDO  ! k-points
   !
   IF (me_pool == root_pool .AND. wan_mode=='standalone') CLOSE (iun_amn, STATUS="KEEP")
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   !
   IF (irr_bz) THEN
      CALL utility_merge_files("iamn", .TRUE.)
   ELSE
      CALL utility_merge_files("amn", .TRUE.)
   ENDIF
   !
   DEALLOCATE(sgf)
   DEALLOCATE(csph)
   IF (noncolin) THEN
      DEALLOCATE(sgf_spinor)
      DEALLOCATE(gf_spinor)
   ENDIF
   DEALLOCATE(amn)
   !
   IF(okvan) THEN
     CALL deallocate_bec_type (becp)
   ENDIF
   !
   WRITE(stdout,*) ' AMN calculated'
   !
   CALL stop_clock( 'compute_amn' )
   !
END SUBROUTINE compute_amn

SUBROUTINE compute_amn_with_scdm
   USE constants,       ONLY : rytoev, pi
   USE io_global,       ONLY : stdout, ionode, ionode_id
   USE mp,              ONLY : mp_bcast, mp_barrier, mp_sum
   USE mp_world,        ONLY : world_comm
   USE mp_pools,        ONLY : intra_pool_comm, inter_pool_comm, my_pool_id, &
                               me_pool, root_pool, nproc_pool
   USE wvfct,           ONLY : nbnd, et, npwx
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions,   ONLY : evc, psic, psic_nc
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE wannier
   USE klist,           ONLY : nkstot, xk, ngk, igk_k, nks
   USE gvect,           ONLY : g, mill
   USE fft_base,        ONLY : dffts
   USE scatter_mod,     ONLY : gather_grid
   USE fft_interfaces,  ONLY : invfft
   USE noncollin_module,ONLY : noncolin, npol
   USE cell_base,       ONLY : at
   USE ions_base,       ONLY : ntyp => nsp, tau
   USE uspp,            ONLY : okvan
   USE lsda_mod,        ONLY : isk, lsda
   !
   IMPLICIT NONE
   !
   INTEGER :: ik, npw, ibnd, iw, nrtot, info, lcwork, ib, gamma_idx, &
              minmn, minmn2, ig, ipool_gamma, ik_gamma_loc, i, j, k, ik_g_w90, &
              nxxs, count_piv_spin_up, m, ibnd_m, ipol, ierr
   REAL(DP):: norm_psi, focc, arg, tpi_r_dot_g, xk_cry(3), rpos_cart(3)
   COMPLEX(DP) :: tmp_cwork(2)
   COMPLEX(DP) :: nowfc_tmp
   INTEGER, ALLOCATABLE :: piv(:)
   !! vv: Pivot array in the QR factorization
   INTEGER, ALLOCATABLE :: piv_pos(:)
   !! Position of the pivot points
   INTEGER, ALLOCATABLE :: piv_spin(:)
   !! Spin index of the pivot points. 1 for spin up, 2 for spin down.
   REAL(DP), ALLOCATABLE :: rwork(:), rwork2(:), singval(:), rpos(:,:)
   !! vv: Real array for the QR factorization and SVD
   REAL(DP), ALLOCATABLE :: et_k(:)
   !! Energy eigenvalues at k. Contains only the included bands.
   COMPLEX(DP), ALLOCATABLE :: phase(:), nowfc(:,:), psi_gamma(:,:), &
       qr_tau(:), cwork(:), Umat(:,:), VTmat(:,:), Amat(:,:)
   !! vv: complex arrays for the SVD factorization
   COMPLEX(DP), ALLOCATABLE :: phase_g(:,:)
   !! exp(iGr) phase for pivot positions. Used for slow Fourier transformation.
   COMPLEX(DP), ALLOCATABLE :: psic_all(:, :)
   !! wavefunction in real space
   COMPLEX(DP), ALLOCATABLE :: psic_1d(:)
   !! wavefunction in real space, flattened if noncollinear.
   COMPLEX(DP), ALLOCATABLE :: evc_k(:, :)
   !! Wavefunction at k. Contains only the included bands.
   !
   INTEGER, EXTERNAL :: global_kpoint_index
   !
#if defined(__SCALAPACK_QRCP)
   REAL(DP) :: tmp_rwork(2)
   INTEGER :: lrwork, context, nprow, npcol, myrow, mycol, descG(9)
   INTEGER :: nblocks, rem, nblocks_loc, rem_loc=0, ibl
   INTEGER, ALLOCATABLE :: piv_p(:)
#endif
   !
   ! vv: Write info about SCDM in output
   IF (TRIM(scdm_entanglement) == 'isolated') THEN
      WRITE(stdout,'(2x,a,a/)') 'Case  : ',trim(scdm_entanglement)
   ELSEIF (TRIM(scdm_entanglement) == 'erfc' .OR. &
        TRIM(scdm_entanglement) == 'gaussian') THEN
      WRITE(stdout,'(2x,a,a)') 'Case  : ',trim(scdm_entanglement)
      WRITE(stdout,'(2x,a,f10.3,a/,1x,a,f10.3,a/)') 'mu    = ', scdm_mu, ' eV', 'sigma =', scdm_sigma, ' eV'
   ENDIF
   !
   CALL start_clock( 'compute_amn' )
   !
   ALLOCATE(et_k(num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating et_k', 1)
   ALLOCATE(evc_k(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_k', 1)
   !
   ! vv: Error for using SCDM with gamma_only
   IF (gamma_only) THEN
      call errore('pw2wannier90', 'The SCDM method does not work with gamma_only calculations.', 1)
   ENDIF
   ! vv: Allocate all the variables for the SCDM method:
   !     1)For the QR decomposition
   !     2)For the unk's on the real grid
   !     3)For the SVD
   nrtot = dffts%nr1*dffts%nr2*dffts%nr3
   nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
   info = 0
   minmn = MIN(num_bands, npol*nrtot)
   ALLOCATE(qr_tau(2*minmn), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating qr_tau', 1)
#if defined(__SCALAPACK_QRCP)
   ! Dimensions of the process grid
   nprow = 1
   npcol = nproc_pool
   ! Initialization of a default BLACS context and the processes grid
   call blacs_get( -1, 0, context )
   call blacs_gridinit( context, 'Row-major', nprow, npcol )
   call blacs_gridinfo( context, nprow, npcol, myrow, mycol )
   call descinit(descG, num_bands, npol*nrtot, minmn, minmn, 0, 0, context, max(1,minmn), info)
   ! Global blocks
   nblocks = npol*nrtot / minmn
   rem = MOD(npol*nrtot, minmn)
   IF (rem > 0) nblocks = nblocks + 1
   ! Local blocks
   nblocks_loc = nblocks / nproc_pool
   rem_loc = MOD(nblocks, nproc_pool)
   IF (me_pool < rem_loc) nblocks_loc = nblocks_loc + 1
   ALLOCATE(psi_gamma(minmn*nblocks_loc,minmn), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psi_gamma', 1)
   ALLOCATE(piv(minmn), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating piv', 1)
#else
   ALLOCATE(piv(npol*nrtot), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating piv', 1)
   ALLOCATE(psi_gamma(npol*nrtot,num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psi_gamma', 1)
   piv(:) = 0
#endif
   !
   minmn2 = MIN(num_bands,n_wannier)
   ALLOCATE(rwork2(5*minmn2), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating rwork2', 1)
   !
   ALLOCATE(piv_pos(n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating piv_pos', 1)
   ALLOCATE(piv_spin(n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating piv_spin', 1)
   ALLOCATE(rpos(3, n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating rpos', 1)
   ALLOCATE(psic_all(nxxs, npol), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psic_all', 1)
   ALLOCATE(psic_1d(npol * nrtot), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psic_1d', 1)
   ALLOCATE(phase(n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phase', 1)
   ALLOCATE(singval(n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating singval', 1)
   ALLOCATE(Umat(num_bands,n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating Umat', 1)
   ALLOCATE(VTmat(n_wannier,n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating VTmat', 1)
   ALLOCATE(Amat(num_bands,n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating Amat', 1)
   ALLOCATE(phase_g(npwx, n_wannier), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating phase_g', 1)
   ALLOCATE(nowfc(n_wannier, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating nowfc', 1)
   !
   IF (wan_mode=='library') THEN
      ALLOCATE(a_mat(num_bands,n_wannier,iknum), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating a_mat', 1)
   ENDIF
   !
   IF (wan_mode == 'standalone') THEN
      ! TODO: append ' with SCDM ' to header
      CALL utility_open_output_file("amn", .TRUE., iun_amn)
      ! CALL date_and_tim( cdate, ctime )
      ! header='Created on '//cdate//' at '//ctime//' with SCDM '
      IF (ionode) THEN
         WRITE (iun_amn, '(3i8,3x,2f10.6)') num_bands, iknum, n_wannier, scdm_mu, scdm_sigma
      ENDIF
   ENDIF
   !
   !vv: Find Gamma-point index in the list of k-vectors
   gamma_idx = -1
   DO ik = 1, nkstot
      IF (ALL(ABS(xk_all(:, ik)) < 1.D-8)) THEN
         gamma_idx = ik
         EXIT
      ENDIF
   ENDDO
   IF (gamma_idx < 0) call errore('compute_amn', 'No Gamma point found.',1)
   !
   ! Find the place of k=Gamma in the pools
   CALL pool_and_local_kpoint_index(nkstot, gamma_idx, ipool_gamma, ik_gamma_loc)
   !
   ! Calculate the pivot points
   !
   CALL start_clock('scdm_QRCP')
   !
   IF (my_pool_id == ipool_gamma) THEN
      !
      ik = ik_gamma_loc
      CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1)
      !
      npw = ngk(ik)
      !
      ! Trim excluded bands from evc and et
      evc_k(:, :) = (0.d0, 0.d0)
      ibnd_m = 0
      DO m = 1, nbnd
         IF (excluded_band(m)) CYCLE
         ibnd_m = ibnd_m + 1
         evc_k(:, ibnd_m) = evc(:, m)
         et_k(ibnd_m) = et(m, ik)
      ENDDO
      !
      DO ibnd = 1, num_bands
         focc = scdm_occupation(et_k(ibnd))
         !
         IF (noncolin) THEN
            psic_nc(:, :) = (0.0_DP, 0.0_DP)
            psic_nc( dffts%nl( igk_k(1:npw, ik) ), 1) = evc_k(1:npw, ibnd)
            psic_nc( dffts%nl( igk_k(1:npw, ik) ), 2) = evc_k(1+npwx:npw+npwx, ibnd)
            CALL invfft('Wave', psic_nc(:,1), dffts)
            CALL invfft('Wave', psic_nc(:,2), dffts)
            !
#if defined(__MPI)
            CALL gather_grid(dffts, psic_nc(:, 1), psic_all(:, 1))
            CALL gather_grid(dffts, psic_nc(:, 2), psic_all(:, 2))
#else
            psic_all = psic_nc
#endif
         ELSE
            ! spin-collinear case
            ! vv: Compute unk's on a real grid (the fft grid)
            psic(:) = (0.0_DP, 0.0_DP)
            psic( dffts%nl( igk_k(1:npw,ik) ) ) = evc_k(1:npw, ibnd)
            CALL invfft ('Wave', psic, dffts)
            !
            psic_all(:, 1) = (0.0_DP, 0.0_DP)
#if defined(__MPI)
            CALL gather_grid(dffts, psic, psic_all(:, 1))
#else
            psic_all(1:nrtot, 1) = psic(1:nrtot)
#endif
         ENDIF ! noncolin
         !
         ! Flatten the real-space and spin indices
         !
         DO ipol = 1, npol
            psic_1d(1+(ipol-1)*nrtot:nrtot+(ipol-1)*nrtot) = psic_all(1:nrtot, ipol)
         ENDDO
         !
         ! vv: Gamma only
         ! vv: Build Psi_k = Unk * focc
         norm_psi = SQRT(SUM( ABS(psic_1d)**2 ))
#if defined(__SCALAPACK_QRCP)
         CALL mp_bcast(norm_psi, root_pool, intra_pool_comm)
         CALL mp_bcast(psic_1d, root_pool, intra_pool_comm)
         DO ibl = 0, nblocks_loc-1
            psi_gamma(minmn*ibl+1:minmn*(ibl+1), ibnd) = &
               psic_1d(minmn*(ibl*nproc_pool+me_pool)+1:minmn*(ibl*nproc_pool+me_pool+1)) * (focc / norm_psi)
         ENDDO
#else
         psi_gamma(1:nrtot, ibnd) = psic_1d(1:nrtot) * (focc / norm_psi)
#endif
         !
      ENDDO
      !
      ! Run QRCP
      !
#if defined(__SCALAPACK_QRCP)
      WRITE(stdout, '(5x,A,I4,A)') "Running QRCP in parallel, using ", nproc_pool, " cores"
      ALLOCATE(piv_p(minmn*nblocks_loc), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating piv_p', 1)
      piv_p(:) = 0
      !
      CALL PZGEQPF( num_bands, npol*nrtot, psi_gamma, 1, 1, descG, piv_p, qr_tau, &
                    tmp_cwork, -1, tmp_rwork, -1, info )
      !
      lcwork = AINT(REAL(tmp_cwork(1)))
      lrwork = AINT(REAL(tmp_rwork(1)))
      ALLOCATE(rwork(lrwork), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating rwork', 1)
      ALLOCATE(cwork(lcwork), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating cwork', 1)
      rwork(:) = 0.0
      cwork(:) = CMPLX(0.0, 0.0)
      !
      CALL PZGEQPF( num_bands, npol*nrtot, TRANSPOSE(CONJG(psi_gamma)), 1, 1, descG, piv_p, qr_tau, &
                    cwork, lcwork, rwork, lrwork, info )
      !
      IF (me_pool == root_pool) piv(1:minmn) = piv_p(1:minmn)
      CALL mp_bcast(info, root_pool, intra_pool_comm)
      CALL mp_bcast(piv, root_pool, intra_pool_comm)
      DEALLOCATE(piv_p)
      DEALLOCATE(cwork)
      DEALLOCATE(rwork)
#else
      WRITE(stdout, '(5x, "Running QRCP in serial")')
#if defined(__SCALAPACK)
      WRITE(stdout, '(10x, A)') "Program compiled with ScaLAPACK but not using it for QRCP."
      WRITE(stdout, '(10x, A)') "To enable ScaLAPACK for QRCP, use valid versions"
      WRITE(stdout, '(10x, A)') "(ScaLAPACK >= 2.1.0 or MKL >= 2020) and set the argument"
      WRITE(stdout, '(10x, A)') "'with-scalapack-qrcp' in configure."
      WRITE(stdout, '(10x, A)') ""
#endif
      ! vv: Perform QR factorization with pivoting on Psi_Gamma
      ! vv: Preliminary call to define optimal values for lwork and cwork size
      ! Perform QR factorization only in a single processer
      IF(me_pool == root_pool) THEN
         ALLOCATE(rwork(2*npol*nrtot), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating rwork', 1)
         rwork(:) = 0.0_DP
         CALL ZGEQP3(num_bands, npol*nrtot, TRANSPOSE(CONJG(psi_gamma)), num_bands, &
            piv, qr_tau, tmp_cwork, -1, rwork, info)
         IF (info/=0) CALL errore('compute_amn', 'Error in priliminary call for the QR factorization', 1)
         lcwork = AINT(REAL(tmp_cwork(1)))
         piv(:) = 0
         ALLOCATE(cwork(lcwork), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating cwork', 1)
         CALL ZGEQP3(num_bands, npol*nrtot, TRANSPOSE(CONJG(psi_gamma)), num_bands, &
            piv, qr_tau, cwork, lcwork, rwork, info)
         DEALLOCATE(cwork)
         DEALLOCATE(rwork)
      ENDIF
      CALL mp_bcast(info, root_pool, intra_pool_comm)
      CALL mp_bcast(piv, root_pool, intra_pool_comm)
#endif
      !
   ENDIF ! ipool_gamma
   !
   CALL mp_bcast(info, ipool_gamma, inter_pool_comm)
   CALL mp_bcast(piv, ipool_gamma, inter_pool_comm)
   IF (info /= 0) CALL errore('compute_amn', 'Error in computing the QR factorization', 1)
   !
   CALL stop_clock('scdm_QRCP')
   !
   ! noncollinear case: calculate position and spin part of piv
   IF (noncolin) THEN
      count_piv_spin_up = 0
      DO iw = 1, n_wannier
         IF (piv(iw) <= nrtot) THEN
            piv_pos(iw) = piv(iw)
            piv_spin(iw) = 1
            count_piv_spin_up = count_piv_spin_up + 1
         ELSE
            piv_pos(iw) = piv(iw) - nrtot
            piv_spin(iw) = -1
         ENDIF
      ENDDO
      WRITE(stdout, '(a,I5)') "  Number of pivot points with spin up  : ", count_piv_spin_up
      WRITE(stdout, '(a,I5)') "  Number of pivot points with spin down: ", n_wannier - count_piv_spin_up
   ELSE
      piv_pos(1:n_wannier) = piv(1:n_wannier)
      piv_spin(:) = 0
      WRITE(stdout, '(a,i8)') '  Number of pivot points: ', n_wannier
   ENDIF
   !
   ! vv: Compute the points in 3d grid
   DO iw = 1, n_wannier
      i = piv_pos(iw) - 1
      k = i / (dffts%nr1 * dffts%nr2)
      i = i - (dffts%nr1 * dffts%nr2) * k
      j = i / dffts%nr1
      i = i - dffts%nr1 * j
      rpos(1, iw) = REAL(i, KIND=DP) / REAL(dffts%nr1, KIND=DP)
      rpos(2, iw) = REAL(j, KIND=DP) / REAL(dffts%nr2, KIND=DP)
      rpos(3, iw) = REAL(k, KIND=DP) / REAL(dffts%nr3, KIND=DP)
      rpos(:, iw) = rpos(:, iw) - ANINT(rpos(:, iw))
   ENDDO
   !
   ! Print pivot positions (and spin indices)
   IF (noncolin) THEN
      WRITE(stdout, '(a)') '  Pivot point positions (alat units) and spin indices:'
   ELSE
      WRITE(stdout, '(a)') '  Pivot point positions (alat units):'
   ENDIF
   DO iw = 1, n_wannier
      rpos_cart(:) = rpos(:, iw)
      CALL cryst_to_cart(1, rpos_cart, at, +1)
      IF (noncolin) THEN
         WRITE(stdout, '(I8,3F12.6,I3)') iw, rpos_cart, piv_spin(iw)
      ELSE
         WRITE(stdout, '(I8,3F12.6)') iw, rpos_cart
      ENDIF
   ENDDO
   WRITE(stdout, *)
   !
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      ik_g_w90 = global_kpoint_index(nkstot, ik) - ikstart + 1
      !
      npw = ngk(ik)
      !
      CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1)
      !
      ! Trim excluded bands from evc and et
      evc_k(:, :) = (0.d0, 0.d0)
      ibnd_m = 0
      DO m = 1, nbnd
         IF (excluded_band(m)) CYCLE
         ibnd_m = ibnd_m + 1
         evc_k(:, ibnd_m) = evc(:, m)
         et_k(ibnd_m) = et(m, ik)
      ENDDO
      !
      ! vv: SCDM method for generating the Amn matrix
      ! jml: calculate of psi_nk at pivot points using slow FT
      !      This is faster than using invfft because the number of pivot
      !      points is much smaller than the number of FFT grid points.
      phase(:) = (0.0_DP,0.0_DP)
      nowfc(:,:) = (0.0_DP,0.0_DP)
      Umat(:,:) = (0.0_DP,0.0_DP)
      VTmat(:,:) = (0.0_DP,0.0_DP)
      Amat(:,:) = (0.0_DP,0.0_DP)
      singval(:) = 0.0_DP
      rwork2(:) = 0.0_DP
      !
      ! jml: calculate phase factors before the loop over bands
      xk_cry = xk(:, ik)
      CALL cryst_to_cart(1, xk_cry, at, -1)
      !
      DO iw = 1, n_wannier
         arg = 2.0_DP * pi * SUM(rpos(:, iw) * xk_cry)
         phase(iw) = CMPLX(COS(arg), SIN(arg), KIND=DP)
         !
         DO ig = 1, npw
            tpi_r_dot_g = 2.0_DP * pi * SUM(rpos(:, iw) * REAL(mill(:, igk_k(ig, ik)), DP))
            phase_g(ig, iw) = CMPLX(COS(tpi_r_dot_g), SIN(tpi_r_dot_g), KIND=DP)
         ENDDO
      ENDDO
      !
      DO ibnd = 1, num_bands
         focc = scdm_occupation(et_k(ibnd))
         !
         norm_psi = SUM( ABS(evc_k(1:npw, ibnd))**2 )
         IF (noncolin) norm_psi = norm_psi + SUM( ABS(evc_k(1+npwx:npw+npwx, ibnd))**2 )
         CALL mp_sum(norm_psi, intra_pool_comm)
         norm_psi = SQRT(norm_psi)
         !
         ! jml: nowfc = sum_G (psi(G) * exp(i*G*r)) * focc  * phase(iw) / norm_psi
         DO iw = 1, n_wannier
            IF (noncolin) THEN
               IF (piv_spin(iw) == 1) THEN
                  ! spin up
                  nowfc_tmp = SUM( evc_k(1:npw, ibnd) * phase_g(1:npw, iw) )
               ELSE
                  ! spin down
                  nowfc_tmp = SUM( evc_k(1+npwx:npw+npwx, ibnd) * phase_g(1:npw, iw) )
               ENDIF
            ELSE
               ! spin collinear
               nowfc_tmp = SUM( evc_k(1:npw, ibnd) * phase_g(1:npw, iw) )
            ENDIF
            nowfc(iw, ibnd) = nowfc_tmp * phase(iw) * focc / norm_psi
         ENDDO
         !
      ENDDO
      CALL mp_sum(nowfc, intra_pool_comm) ! jml
      !
      IF (me_pool == root_pool) THEN
         CALL ZGESVD('S', 'S', num_bands, n_wannier, TRANSPOSE(CONJG(nowfc)), num_bands, &
            singval, Umat, num_bands, VTmat, n_wannier, tmp_cwork, -1, rwork2, info)
         lcwork = AINT(REAL(tmp_cwork(1)))
         ALLOCATE(cwork(lcwork), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating cwork', 1)
         ! vv: SVD to generate orthogonal projections
         CALL ZGESVD('S', 'S', num_bands, n_wannier, TRANSPOSE(CONJG(nowfc)), num_bands, &
            singval, Umat, num_bands, VTmat, n_wannier, cwork, lcwork, rwork2, info)
         DEALLOCATE(cwork)
      ENDIF
      CALL mp_bcast(info, root_pool, intra_pool_comm)
      IF (info /= 0) CALL errore('compute_amn', &
         'Error in computing the SVD of the PSI matrix in the SCDM method', 1)
      !
      IF (me_pool == root_pool) THEN
         Amat = MATMUL(Umat, VTmat)
         DO iw = 1, n_wannier
            DO ibnd = 1, num_bands
               WRITE(iun_amn,'(3i10,2f18.12)') ibnd, iw, ik_g_w90, REAL(Amat(ibnd,iw)), AIMAG(Amat(ibnd,iw))
            ENDDO
         ENDDO
      ENDIF ! root_pool
   ENDDO ! k-points
   !
   IF (me_pool == root_pool) CLOSE(iun_amn, STATUS="KEEP")
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   !
   CALL utility_merge_files("amn", .TRUE.)
   !
   ! vv: Deallocate all the variables for the SCDM method
   DEALLOCATE(psi_gamma)
   DEALLOCATE(nowfc)
   DEALLOCATE(piv)
   DEALLOCATE(piv_pos)
   DEALLOCATE(piv_spin)
   DEALLOCATE(qr_tau)
   DEALLOCATE(rwork2)
   DEALLOCATE(rpos)
   DEALLOCATE(Umat)
   DEALLOCATE(VTmat)
   DEALLOCATE(Amat)
   DEALLOCATE(singval)
   DEALLOCATE(phase_g)
   DEALLOCATE(psic_all)
   !
   WRITE(stdout,*) ' AMN calculated'
   CALL stop_clock('compute_amn')
   !
   CONTAINS
   !
   FUNCTION scdm_occupation(e) RESULT(focc)
      !! vv: Generate the occupation numbers matrix according to scdm_entanglement
      REAL(DP), INTENT(IN) :: e
      !! Energy eigenvalue
      REAL(DP) :: focc
      !! SCDM occupation for constructing quasi-density matrix
      IF(TRIM(scdm_entanglement) == 'isolated') THEN
         focc = 1.0_DP
      ELSEIF (TRIM(scdm_entanglement) == 'erfc') THEN
         focc = 0.5_DP * ERFC((e * rytoev - scdm_mu) / scdm_sigma)
      ELSEIF (TRIM(scdm_entanglement) == 'gaussian') THEN
         focc = EXP(-1.0_DP * ((e * rytoev - scdm_mu)**2) / (scdm_sigma**2))
      ELSE
         CALL errore('compute_amn_with_scdm', 'scdm_entanglement value not recognized.', 1)
      ENDIF
   END FUNCTION scdm_occupation
   !
END SUBROUTINE compute_amn_with_scdm

SUBROUTINE compute_amn_with_atomproj
   ! Use internal UPF atomic projectors or external projectors
   ! to compute amn matrices.
   !
   ! The code is roughly the same as projwfc.x with some entensions:
   ! 1. allow using external projectors (i.e. custom radial functions)
   ! 2. allow skipping orthogonalization of projectors
   ! 3. allow excluding projectors specified by user
   ! 4. allow excluding bands specified by user
   !
   USE kinds, ONLY: DP
   USE mp, ONLY: mp_bcast, mp_barrier
   USE mp_pools, ONLY: me_pool, root_pool, intra_pool_comm
   USE mp_world, ONLY: world_comm
   USE io_global, ONLY: stdout, ionode, ionode_id
   USE ions_base, ONLY: nat, ityp, atm, nsp
   USE basis, ONLY: natomwfc, swfcatom
   USE klist, ONLY: xk, nks, nkstot, nelec, ngk, igk_k
   USE lsda_mod, ONLY: nspin, lsda, isk
   USE noncollin_module, ONLY: noncolin, npol, lspinorb, domag
   USE wvfct, ONLY: npwx, nbnd
   USE uspp, ONLY: nkb, vkb
   USE uspp_init, ONLY : init_us_2
   USE becmod, ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
   USE io_files, ONLY: prefix, restart_dir, tmp_dir, nwordwfc, iunwfc
   USE control_flags, ONLY: gamma_only, use_para_diag
   USE pw_restart_new, ONLY: read_collected_wfc
   USE wavefunctions, ONLY: evc
   USE projections, ONLY: nlmchi, fill_nlmchi, compute_mj, &
                          sym_proj_g, sym_proj_k, sym_proj_nc, sym_proj_so, &
                          compute_zdistmat, compute_ddistmat, &
                          wf_times_overlap, wf_times_roverlap
   USE atproj, ONLY: atom_proj_dir, atom_proj_ext, atom_proj_ortho, &
                     atom_proj_sym, natproj, nexatproj, nexatproj_max, &
                     atproj_excl, atproj_typs, atom_proj_exclude, &
                     allocate_atproj_type, read_atomproj, init_tab_atproj, &
                     deallocate_atproj, atomproj_wfc
   USE wannier
   !
   IMPLICIT NONE
   !
   INCLUDE 'laxlib.fh'
   !
   INTEGER :: npw, npw_, ik, ibnd, nwfc, lmax_wfc
   INTEGER :: i, j, k, it, l, m, ib, ip, ik_g_w90, ibnd1
   REAL(DP), ALLOCATABLE :: e(:)
   COMPLEX(DP), ALLOCATABLE :: wfcatom(:, :), wfcatomall(:, :)
   COMPLEX(DP), ALLOCATABLE :: proj0(:, :), proj0all(:, :), proj(:, :)
   COMPLEX(DP), ALLOCATABLE :: e_work_d(:, :)
   ! Some workspace for gamma-point calculation ...
   REAL(DP), ALLOCATABLE :: rproj0(:, :), rproj0all(:, :)
   COMPLEX(DP), ALLOCATABLE :: overlap_d(:, :), work_d(:, :), diag(:, :), vv(:, :)
   REAL(DP), ALLOCATABLE :: roverlap_d(:, :)
   COMPLEX(DP), ALLOCATABLE :: evc_k(:, :)
   !
   LOGICAL :: freeswfcatom
   !
   INTEGER :: idesc(LAX_DESC_SIZE)
   INTEGER, ALLOCATABLE :: idesc_ip(:, :, :)
   INTEGER, ALLOCATABLE :: rank_ip(:, :)
   ! matrix distribution descriptors
   INTEGER :: nx, nrl, nrlx
   ! maximum local block dimension
   LOGICAL :: la_proc
   ! flag to distinguish procs involved in linear algebra
   INTEGER, ALLOCATABLE :: notcnv_ip(:)
   INTEGER, ALLOCATABLE :: ic_notcnv(:)
   LOGICAL :: do_distr_diag_inside_bgrp
   INTEGER :: nproc_ortho
   ! distinguishes active procs in parallel linear algebra
   CHARACTER(len=256) :: err_str
   LOGICAL :: has_excl_proj
   INTEGER :: ierr
   !
   INTEGER, EXTERNAL :: global_kpoint_index
   !
   CALL start_clock('compute_amn')
   !
   IF (wan_mode == 'library') THEN
     CALL errore('pw2wannier90', 'have not tested with library mode', 1)
   ENDIF
   !
   IF (atom_proj_ext) THEN
      IF (domag) CALL errore('pw2wannier90', &
                   'does not support magnetism with external projectors', 1)
      IF (noncolin) CALL errore('pw2wannier90', &
                   'does not support non-collinear magnetism with external projectors', 1)
      IF (atom_proj_sym) CALL errore('pw2wannier90', &
                   'does not support symmetrization with external projectors', 1)
   ENDIF
   !
   ALLOCATE(evc_k(npol*npwx, num_bands), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_k', 1)
   !
   IF (atom_proj_ext) THEN
      ALLOCATE (atproj_typs(nsp), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating atproj_typs', 1)
   ENDIF
   !
   IF (ionode) THEN
      IF (atom_proj_ext) THEN
         WRITE (stdout, '(a)') '  Using atomic projectors from dir '//TRIM(atom_proj_dir)
         WRITE (stdout, *) ''
         !
         CALL read_atomproj(atproj_typs)
         !
         n_proj = 0
         DO i = 1, nat
            it = ityp(i)
            DO nwfc = 1, atproj_typs(it)%nproj
               l = atproj_typs(it)%l(nwfc)
               DO m = 1, 2*l+1
                  n_proj = n_proj + 1
                  WRITE (stdout, 1000, ADVANCE="no") n_proj, i, atproj_typs(it)%atsym, nwfc, l
                  WRITE (stdout, '(" m=", i2, ")")') m
               ENDDO
            ENDDO
         ENDDO
         WRITE (stdout, *) ''
      ELSE
         WRITE (stdout, '(a)') '  Use atomic projectors from UPF'
         WRITE (stdout, *) ''
         WRITE (stdout, '( 5x,"(read from pseudopotential files):"/)')
         CALL fill_nlmchi(natomwfc, lmax_wfc)
         DO nwfc = 1, natomwfc
            WRITE (stdout, 1000, ADVANCE="no") &
               nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
               nlmchi(nwfc)%n, nlmchi(nwfc)%l
            IF (lspinorb) THEN
               WRITE (stdout, '(" j=", f3.1, " m_j=", f4.1, ")")') &
                  nlmchi(nwfc)%jj, compute_mj(nlmchi(nwfc)%jj, nlmchi(nwfc)%l, nlmchi(nwfc)%m)
            ELSE IF (noncolin) THEN
                WRITE (stdout, '(" m=", i2, " s_z=", f4.1, ")")') &
                  nlmchi(nwfc)%m, 0.5D0 - INT(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l + 2))
            ELSE
               WRITE (stdout, '(" m=", i2, ")")') nlmchi(nwfc)%m
           ENDIF
         ENDDO
         WRITE (stdout, *) ''
         n_proj = natomwfc
      ENDIF ! atom_proj_ext
1000  FORMAT(5X, "state #", i4, ": atom ", i3, " (", a3, "), wfc ", i2, &
             " (l=", i1)
      !
      IF (n_proj <= 0) CALL errore('pw2wannier90', &
            'Cannot project on zero atomic projectors!', 1)
      !
      ! check exclude
      ALLOCATE(atproj_excl(n_proj), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating atproj_excl', 1)
      atproj_excl = .false.
      !
      DO i = 1, nexatproj_max
         IF (atom_proj_exclude(i) > n_proj) THEN
            WRITE (err_str, *) 'atom_proj_exclude(', i, ') = ', atom_proj_exclude(i), &
                               '> total number of projectors (', n_proj, ')'
            CALL errore('pw2wannier90', err_str, i)
         ELSEIF (atom_proj_exclude(i) < 0) THEN
           CYCLE
         ELSE
            atproj_excl(atom_proj_exclude(i)) = .true.
         ENDIF
      ENDDO
      !
      nexatproj = COUNT(atproj_excl)
      has_excl_proj = (nexatproj > 0)
      !
      natproj = n_proj
      !
      IF (has_excl_proj) THEN
         WRITE (stdout, *) '    excluded projectors: '
         j = 0 ! how many elements have been written
         DO i = 1, n_proj
            IF (atproj_excl(i)) THEN
               WRITE (stdout, '(i8)', advance='no') i
               j = j + 1
               IF (MOD(j, 10) == 0) WRITE (stdout, *)
           ENDIF
         ENDDO
         WRITE (stdout, *) ''
         n_proj = n_proj - nexatproj
      ENDIF
      !
      IF (gamma_only) WRITE (stdout, '(5x,"gamma-point specific algorithms are used")')
      !
   ENDIF ! ionode
   !
   ! MPI related calls
   IF (atom_proj_ext) THEN
      DO it = 1, nsp
         i = atproj_typs(it)%ngrid
         j = atproj_typs(it)%nproj
         CALL mp_bcast(i, ionode_id, world_comm)
         CALL mp_bcast(j, ionode_id, world_comm)
         IF (.NOT. ionode) CALL allocate_atproj_type(atproj_typs(it), i, j)
         !
         CALL mp_bcast(atproj_typs(it)%atsym, ionode_id, world_comm)
         CALL mp_bcast(atproj_typs(it)%xgrid, ionode_id, world_comm)
         CALL mp_bcast(atproj_typs(it)%rgrid, ionode_id, world_comm)
         CALL mp_bcast(atproj_typs(it)%l, ionode_id, world_comm)
         CALL mp_bcast(atproj_typs(it)%radial, ionode_id, world_comm)
      ENDDO
      !
      CALL init_tab_atproj(world_comm)
   ELSE
      ! need to access nlmchi, natomwfc, lmax_wfc on each core,
      ! the root node has been filled already
      IF (.NOT. ionode) CALL fill_nlmchi(natomwfc, lmax_wfc)
   ENDIF ! atom_proj_ext
   !
   CALL mp_bcast(natproj, ionode_id, world_comm)
   CALL mp_bcast(n_proj, ionode_id, world_comm)
   CALL mp_bcast(nexatproj, ionode_id, world_comm)
   CALL mp_bcast(has_excl_proj, ionode_id, world_comm)
   IF (.NOT. ionode) THEN
      ALLOCATE(atproj_excl(n_proj+nexatproj), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating atproj_excl', 1)
   ENDIF
   CALL mp_bcast(atproj_excl, ionode_id, world_comm)
   !
   !   Initialize parallelism for linear algebra
   !
   CALL set_para_diag(n_proj, use_para_diag)
   !
   CALL desc_init(n_proj, nx, la_proc, idesc, rank_ip, idesc_ip)
   CALL laxlib_getval(nproc_ortho=nproc_ortho)
   use_para_diag = (nproc_ortho > 1)
   IF (ionode .AND. use_para_diag) THEN
      WRITE (stdout, '(5x,"linear algebra parallelized on ",i3," procs")') nproc_ortho
   ENDIF
   !
   IF (ionode) THEN
      !
      !   nbnd = num_bands + nexband
      ! For UPF projectors:
      !   natomwfc = n_proj + nexatproj
      !
      WRITE (stdout, *)
      WRITE (stdout, *) '    Problem Sizes '
      WRITE (stdout, *) '      n_proj    = ', n_proj
      IF (atom_proj_ext) THEN
        WRITE (stdout, *) '      natproj   = ', natproj
      ELSE
        WRITE (stdout, *) '      natomwfc  = ', natomwfc
      END IF
      WRITE (stdout, *) '      num_bands = ', num_bands
      WRITE (stdout, *) '      nbnd      = ', nbnd
      WRITE (stdout, *) '      nkstot    = ', nkstot
      IF (use_para_diag) WRITE (stdout, *) '      nx        = ', nx
      WRITE (stdout, *) '      npwx      = ', npwx
      WRITE (stdout, *) '      nkb       = ', nkb
      WRITE (stdout, *)
   END IF
   !
   ALLOCATE (proj(num_bands, n_proj), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating proj', 1)
   !
   IF (.NOT. ALLOCATED(swfcatom)) THEN
      ALLOCATE (swfcatom(npwx*npol, n_proj), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating swfcatom', 1)
      freeswfcatom = .TRUE.
   ELSE
      freeswfcatom = .FALSE.
   ENDIF
   !
   ALLOCATE (wfcatom(npwx*npol, n_proj), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wfcatom', 1)
   IF (has_excl_proj) THEN
      ! additional space for excluded projectors, for atomic_wfc(), etc.
      IF (atom_proj_ext) THEN
         ALLOCATE (wfcatomall(npwx*npol, natproj), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wfcatomall', 1)
      ELSE
         ALLOCATE (wfcatomall(npwx*npol, natomwfc), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating wfcatomall', 1)
      ENDIF
   ENDIF
   ALLOCATE (e(n_proj), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating e', 1)
   !
   CALL utility_open_output_file("amn", .TRUE., iun_amn)
   IF (ionode) WRITE(iun_amn, *) num_bands, iknum, n_proj
   !
   ! loop on k points
   !
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      !
      ik_g_w90 = global_kpoint_index(nkstot, ik) - ikstart + 1
      npw = ngk(ik)
      !
      ! Read wavefunctions at k, exclude the excluded bands
      !
      CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1 )
      !
      ibnd1 = 0
      DO ibnd = 1, nbnd
         IF (excluded_band(ibnd)) CYCLE
         ibnd1 = ibnd1 + 1
         evc_k(:, ibnd1) = evc(:, ibnd)
      ENDDO
      !
      wfcatom(:, :) = (0.0_DP, 0.0_DP)
      IF (atom_proj_ext) THEN
         IF (.NOT. has_excl_proj) THEN
            CALL atomproj_wfc(ik, wfcatom)
         ELSE
            wfcatomall(:, :) = (0.0_DP, 0.0_DP)
            CALL atomproj_wfc(ik, wfcatomall)
            ! exclude projectors
            i = 1 ! counter for wfcatom
            DO j = 1, natproj ! counter for wfcatomall
               IF (atproj_excl(j)) CYCLE
               wfcatom(:, i) = wfcatomall(:, j)
               i = i + 1
            END DO
            IF ((i - 1) /= n_proj) THEN
               CALL errore('compute_amn_with_atomproj', &
                          'internal error when excluding projectors', 1)
            ENDIF
         ENDIF
      ELSE
         IF (.NOT. has_excl_proj) THEN
            IF (noncolin) THEN
               CALL atomic_wfc_nc_proj(ik, wfcatom)
            ELSE
               CALL atomic_wfc(ik, wfcatom)
            END IF
         ELSE
            wfcatomall(:, :) = (0.0_DP, 0.0_DP)
            IF (noncolin) THEN
               CALL atomic_wfc_nc_proj(ik, wfcatomall)
            ELSE
               CALL atomic_wfc(ik, wfcatomall)
            ENDIF
            ! exclude projectors
            i = 1 ! counter for wfcatom
            DO j = 1, natomwfc ! counter for wfcatomall
               IF (atproj_excl(j)) CYCLE
               wfcatom(:, i) = wfcatomall(:, j)
               i = i + 1
            END DO
            IF ((i - 1) /= n_proj) CALL errore('compute_amn_with_atomproj', &
                        'internal error when excluding projectors', 1)
            !
         ENDIF
      ENDIF ! atom_proj_ext
      !
      CALL allocate_bec_type(nkb, n_proj, becp)
      !
      CALL init_us_2(npw, igk_k(1, ik), xk(1, ik), vkb)
      CALL calbec(npw, vkb, wfcatom, becp)
      CALL s_psi(npwx, npw, n_proj, wfcatom, swfcatom)
      !
      CALL deallocate_bec_type(becp)
      !
      ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>
      ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j>
      !
      npw_ = npw
      IF (noncolin) npw_ = npol*npwx
      !
      IF (atom_proj_ortho) THEN
         IF (la_proc) THEN
            ALLOCATE (overlap_d(nx, nx), stat=ierr)
            IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating overlap_d', 1)
         ELSE
            ALLOCATE (overlap_d(1, 1), stat=ierr)
            IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating overlap_d', 1)
         ENDIF
         overlap_d = (0.D0, 0.D0)
         IF (gamma_only) THEN
            !
            ! in the Gamma-only case the overlap matrix (real) is copied
            ! to a complex one as for the general case - easy but wasteful
            !
            IF (la_proc) THEN
               ALLOCATE (roverlap_d(nx, nx), stat=ierr)
               IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating roverlap_d', 1)
            ELSE
               ALLOCATE (roverlap_d(1, 1), stat=ierr)
               IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating roverlap_d', 1)
            END IF
            roverlap_d = 0.D0
            CALL compute_ddistmat(npw, n_proj, nx, wfcatom, swfcatom, roverlap_d, &
                                 idesc, rank_ip, idesc_ip)
            overlap_d(:, :) = CMPLX(roverlap_d(:, :), 0.0_DP, kind=dp)
         ELSE
            CALL compute_zdistmat(npw_, n_proj, nx, wfcatom, swfcatom, overlap_d, &
                                 idesc, rank_ip, idesc_ip)
         ENDIF
         !
         ! diagonalize the overlap matrix
         !
         IF (la_proc) THEN
            !
            ALLOCATE (work_d(nx, nx), stat=ierr)
            IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating work_d', 1)
            !
            nrl = idesc(LAX_DESC_NRL)
            nrlx = idesc(LAX_DESC_NRLX)
            !
            ALLOCATE (diag(nrlx, n_proj), stat=ierr)
            IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating diag', 1)
            ALLOCATE (vv(nrlx, n_proj), stat=ierr)
            IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating vv', 1)
            !
            !  re-distribute the overlap matrix for parallel diagonalization
            !
            CALL blk2cyc_redist(n_proj, diag, nrlx, n_proj, overlap_d, nx, nx, idesc)
            !
            ! parallel diagonalization
            !
            CALL zhpev_drv('V', diag, nrlx, e, vv, nrlx, nrl, n_proj, &
                           idesc(LAX_DESC_NPC)*idesc(LAX_DESC_NPR), &
                           idesc(LAX_DESC_MYPE), idesc(LAX_DESC_COMM))
            !
            !  bring distributed eigenvectors back to original distribution
            !
            CALL cyc2blk_redist(n_proj, vv, nrlx, n_proj, work_d, nx, nx, idesc)
            !
            DEALLOCATE (vv)
            DEALLOCATE (diag)
            !
         ELSE
            ALLOCATE (work_d(1, 1), stat=ierr)
            IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating work_d', 1)
         END IF
         !
         CALL mp_bcast(e, root_pool, intra_pool_comm)
         !
         ! calculate O^{-1/2} (actually, its transpose)
         !
         DO i = 1, n_proj
            e(i) = 1.D0/dsqrt(e(i))
         END DO
         !
         IF (la_proc) THEN
            ALLOCATE (e_work_d(nx, nx), stat=ierr)
            IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating e_work_d', 1)
            DO j = 1, idesc(LAX_DESC_NC)
               DO i = 1, idesc(LAX_DESC_NR)
                  e_work_d(i, j) = e(j + idesc(LAX_DESC_IC) - 1)*work_d(i, j)
               END DO
            END DO
            CALL sqr_mm_cannon('N', 'C', n_proj, (1.0_DP, 0.0_DP), e_work_d, &
                               nx, work_d, nx, (0.0_DP, 0.0_DP), overlap_d, nx, idesc)
            CALL laxlib_zsqmher(n_proj, overlap_d, nx, idesc)
            DEALLOCATE (e_work_d)
         END IF
         !
         DEALLOCATE (work_d)
         !
         ! calculate wfcatom = O^{-1/2} \hat S | phi>
         !
         IF (gamma_only) THEN
            roverlap_d(:, :) = REAL(overlap_d(:, :), DP)
            CALL wf_times_roverlap(nx, npw, swfcatom, roverlap_d, wfcatom, &
                                   idesc, rank_ip, idesc_ip, la_proc)
            DEALLOCATE (roverlap_d)
         ELSE
            CALL wf_times_overlap(nx, npw_, swfcatom, overlap_d, wfcatom, &
                                  idesc, rank_ip, idesc_ip, la_proc)
         END IF
         DEALLOCATE (overlap_d)
      ELSE
         wfcatom = swfcatom
      ENDIF ! atom_proj_ortho
      !
      ! make the projection
      !   <psi_i| O^{-1/2} \hat S | phi_j> (if atom_proj_ortho)
      !   <psi_i| \hat S | phi_j> (if not atom_proj_ortho)
      ! symmetrize the projections if required (not implemented yet)
      !
      IF (gamma_only) THEN
         !
         ALLOCATE (rproj0(n_proj, num_bands), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating rproj0', 1)
         CALL calbec(npw_, wfcatom, evc_k, rproj0, nbnd=num_bands)
         ! haven't tested symmetrization with external projectors, so
         ! I disable these for now.
         ! IF ((.NOT. atom_proj_ext) .AND. atom_proj_sym) THEN
         !   IF (has_excl_proj) THEN
         !     ALLOCATE (rproj0all(natomwfc, num_bands))
         !     rproj0all = 0.0_DP
         !     ! expand the size to natomwfc so I can call sym_proj_g
         !     ! the excluded part is just 0.0
         !     i = 1 ! counter for rproj0
         !     DO j = 1, natomwfc ! counter for rproj0all
         !       IF (atproj_excl(j)) CYCLE
         !       rproj0all(j, :) = rproj0(i, :)
         !       i = i + 1
         !     END DO
         !     !
         !     CALL sym_proj_g(rproj0all)
         !     !
         !     ! exclude projectors
         !     i = 1 ! counter for rproj0
         !     DO j = 1, natomwfc ! counter for rproj0all
         !       IF (atproj_excl(j)) CYCLE
         !       rproj0(i, :) = rproj0all(j, :)
         !       i = i + 1
         !     END DO
         !     !
         !     DEALLOCATE (rproj0all)
         !   ELSE
         !     CALL sym_proj_g(rproj0)
         !   END IF
         ! END IF
         !
         ! Note the CONJG, I need <psi|g>, while rpoj0 = <g|psi>
         proj(:, :) = TRANSPOSE(rproj0(:, :))
         DEALLOCATE (rproj0)
         !
      ELSE
         !
         ALLOCATE (proj0(n_proj, num_bands), stat=ierr)
         IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating proj0', 1)
         CALL calbec(npw_, wfcatom, evc_k, proj0, nbnd=num_bands)
         !
         ! IF ((.NOT. atom_proj_ext) .AND. atom_proj_sym) THEN
         !   IF (has_excl_proj) THEN
         !     ALLOCATE (proj0all(natomwfc, num_bands))
         !     proj0all = (0.0_DP, 0.0_DP)
         !     ! expand the size to natomwfc so I can call sym_proj_*
         !     ! the exclude part is just 0.0
         !     i = 1 ! counter for proj0
         !     DO j = 1, natomwfc ! counter for proj0all
         !       IF (atproj_excl(j)) CYCLE
         !       proj0all(j, :) = proj0(i, :)
         !       i = i + 1
         !     END DO
         !     !
         !     IF (lspinorb) THEN
         !       CALL sym_proj_so(domag, proj0all)
         !     ELSE IF (noncolin) THEN
         !       CALL sym_proj_nc(proj0all)
         !     ELSE
         !       CALL sym_proj_k(proj0all)
         !     END IF
         !     !
         !     ! exclude projectors
         !     i = 1 ! counter for proj0
         !     DO j = 1, natomwfc ! counter for proj0all
         !       IF (atproj_excl(j)) CYCLE
         !       proj0(i, :) = proj0all(j, :)
         !       i = i + 1
         !     END DO
         !     !
         !     DEALLOCATE (proj0all)
         !   ELSE
         !     IF (lspinorb) THEN
         !       CALL sym_proj_so(domag, proj0)
         !     ELSE IF (noncolin) THEN
         !       CALL sym_proj_nc(proj0)
         !     ELSE
         !       CALL sym_proj_k(proj0)
         !     END IF
         !   END IF
         ! END IF
         !
         ! Note the CONJG, I need <psi|g>, while proj0 = <g|psi>
         proj(:, :) = TRANSPOSE(CONJG(proj0(:, :)))
         DEALLOCATE (proj0)
         !
      ENDIF ! gamma_only
      !
      ! Write amn to file
      !
      IF (me_pool == root_pool) THEN
         DO ip = 1, n_proj
            DO ib = 1, num_bands
               WRITE (iun_amn, '(3i10,2f18.12)') ib, ip, ik_g_w90, proj(ib, ip)
            ENDDO
         ENDDO
      ENDIF
      !
   ENDDO ! on k-points
   !
   IF (me_pool == root_pool) CLOSE (iun_amn, STATUS="KEEP")
   !
   CALL mp_barrier(world_comm)
   !
   ! If using pool parallelization, concatenate files written by other nodes
   ! to the main output.
   !
   CALL utility_merge_files("amn", .TRUE.)
   !
   CALL deallocate_atproj()
   DEALLOCATE (e)
   DEALLOCATE (wfcatom)
   IF (freeswfcatom) DEALLOCATE (swfcatom)
   IF (has_excl_proj) DEALLOCATE (wfcatomall)
   DEALLOCATE (idesc_ip)
   DEALLOCATE (rank_ip)
   !
   ! write to standard output and to file
   !
   IF (ionode) THEN
      WRITE (stdout, '(/)')
      WRITE (stdout, *) ' AMN calculated'
   END IF
   !
   DEALLOCATE (proj)
   CALL laxlib_end()
   !
   CALL stop_clock('compute_amn')
   !
END SUBROUTINE compute_amn_with_atomproj

subroutine orient_gf_spinor(npw)
   use constants, only: eps6
   use noncollin_module, only: npol
   use wvfct,           ONLY : npwx
   use wannier

   implicit none

   integer :: npw, iw, ipol, istart, iw_spinor
   logical :: spin_z_pos, spin_z_neg
   complex(dp) :: fac(2)


   gf_spinor = (0.0d0, 0.0d0)
  DO iw = 1,n_proj
     spin_z_pos=.false.;spin_z_neg=.false.
     ! detect if spin quantisation axis is along z
     if((abs(spin_qaxis(1,iw)-0.0d0)<eps6).and.(abs(spin_qaxis(2,iw)-0.0d0)<eps6) &
        .and.(abs(spin_qaxis(3,iw)-1.0d0)<eps6)) then
        spin_z_pos=.true.
     elseif(abs(spin_qaxis(1,iw)-0.0d0)<eps6.and.abs(spin_qaxis(2,iw)-0.0d0)<eps6 &
        .and.abs(spin_qaxis(3,iw)+1.0d0)<eps6) then
        spin_z_neg=.true.
     endif
     if(spin_z_pos .or. spin_z_neg) then
        if(spin_z_pos) then
           ipol=(3-spin_eig(iw))/2
        else
           ipol=(3+spin_eig(iw))/2
        endif
        istart = (ipol-1)*npwx + 1
        gf_spinor(istart:istart+npw-1, iw) = gf(1:npw, iw)
     else
       if(spin_eig(iw)==1) then
          fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw)+1)*cmplx(1.0d0,0.0d0,dp)
          fac(2)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
       else
          fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw))*cmplx(1.0d0,0.0d0,dp)
          fac(2)=(1.0_dp/sqrt(1-spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
       endif
       gf_spinor(1:npw, iw) = gf(1:npw, iw) * fac(1)
       gf_spinor(npwx + 1:npwx + npw, iw) = gf(1:npw, iw) * fac(2)
     endif
  enddo
end subroutine orient_gf_spinor
!
SUBROUTINE generate_guiding_functions(ik)
   !! gf should not be normalized at each k point because the atomic orbitals are
   !! not orthonormal so that their Bloch representation is not normalized.
   !
   USE io_global,  ONLY : stdout
   USE constants, ONLY : pi, tpi, fpi, eps8
   USE control_flags, ONLY : gamma_only
   USE gvect, ONLY : g, gstart
   USE cell_base,  ONLY : tpiba
   USE wannier
   USE klist,      ONLY : xk, ngk, igk_k
   USE cell_base, ONLY : bg
   USE mp, ONLY : mp_sum
   USE mp_pools,  ONLY : intra_pool_comm

   IMPLICIT NONE

   INTEGER, INTENT(in) :: ik
   INTEGER, PARAMETER :: lmax=3, lmax2=(lmax+1)**2
   INTEGER :: npw, iw, ig, l, ierr
   INTEGER :: lmax_iw, lm, ipol, n1, n2, n3, nr1, nr2, nr3, iig
   real(DP) :: arg
   COMPLEX(DP) :: lphase
   real(DP), ALLOCATABLE :: gk(:,:), qg(:), ylm(:,:), radial(:,:)
   COMPLEX(DP), ALLOCATABLE :: sk(:)
   !
   npw = ngk(ik)
   ALLOCATE( gk(3,npw), qg(npw), ylm(npw,lmax2), sk(npw), radial(npw,0:lmax), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating gk/qg/ylm/sk/radial', 1)
   !
   DO ig = 1, npw
      gk (1,ig) = xk(1, ik) + g(1, igk_k(ig,ik) )
      gk (2,ig) = xk(2, ik) + g(2, igk_k(ig,ik) )
      gk (3,ig) = xk(3, ik) + g(3, igk_k(ig,ik) )
      qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
   ENDDO

   CALL ylmr2 (lmax2, npw, gk, qg, ylm)
   ! define qg as the norm of (k+g) in a.u.
   qg(:) = sqrt(qg(:)) * tpiba

   DO iw = 1, n_proj
      !
      gf(:,iw) = (0.d0,0.d0)

      CALL radialpart(npw, qg, alpha_w(iw), r_w(iw), lmax, radial)

      DO lm = 1, lmax2
         IF ( abs(csph(lm,iw)) < eps8 ) CYCLE
         l = int (sqrt( lm-1.d0))
         lphase = (0.d0,-1.d0)**l
         !
         DO ig=1,npw
            gf(ig,iw) = gf(ig,iw) + csph(lm,iw) * ylm(ig,lm) * radial(ig,l) * lphase
         ENDDO !ig
      ENDDO ! lm
      DO ig=1,npw
         iig = igk_k(ig,ik)
         arg = ( gk(1,ig)*center_w(1,iw) + gk(2,ig)*center_w(2,iw) + &
                                           gk(3,ig)*center_w(3,iw) ) * tpi
         ! center_w are cartesian coordinates in units of alat
         sk(ig) = cmplx(cos(arg), -sin(arg) ,kind=DP)
         gf(ig,iw) = gf(ig,iw) * sk(ig)
      ENDDO
   ENDDO
   !
   DEALLOCATE ( gk, qg, ylm, sk, radial)
   RETURN
END SUBROUTINE generate_guiding_functions

SUBROUTINE write_band
   USE io_global,  ONLY : stdout, ionode
   USE mp,         ONLY : mp_barrier, mp_sum
   USE mp_world,   ONLY : world_comm
   USE mp_pools,   ONLY : me_pool, root_pool, my_pool_id, inter_pool_comm
   USE constants,  ONLY : rytoev
   USE wvfct,      ONLY : nbnd, et
   USE klist,      ONLY : nkstot, nks
   USE lsda_mod,   ONLY : lsda, isk
   USE wannier
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=256) :: filename
   INTEGER :: ik, ibnd, ibnd1, ikevc, ikevc_g, ierr
   !
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   INTEGER, EXTERNAL :: global_kpoint_index
   !
   IF (wan_mode == 'standalone') THEN
      IF (me_pool == root_pool) THEN
         IF (irr_bz) THEN
            filename = TRIM(seedname) // ".ieig"
         ELSE
            filename = TRIM(seedname) // ".eig"
         ENDIF
         IF (.NOT. ionode) filename = TRIM(filename) // TRIM(int_to_char(my_pool_id+1))
         OPEN(NEWUNIT=iun_band, FILE=TRIM(filename), FORM='FORMATTED', STATUS='REPLACE')
      ENDIF
   ELSEIF (wan_mode == 'library') THEN
      ALLOCATE(eigval(num_bands, iknum), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating eigval', 1)
      eigval = 0.0_DP
   ELSE
      CALL errore('write_band', 'value of wan_mode not recognised', 1)
   ENDIF
   !
   !
   DO ik = 1, nks
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      ikevc = ik - ikstart + 1
      ikevc_g = global_kpoint_index(nkstot, ik) - ikstart + 1
      !
      ibnd1 = 0
      DO ibnd = 1, nbnd
         IF (excluded_band(ibnd)) CYCLE
         ibnd1 = ibnd1 + 1
         IF (wan_mode == 'standalone') THEN
            IF (me_pool == root_pool) WRITE (iun_band,'(2i10,f18.12)') ibnd1, ikevc_g, et(ibnd,ik)*rytoev
         ELSEIF (wan_mode == 'library') THEN
            eigval(ibnd1,ikevc_g) = et(ibnd,ik)*rytoev
         ENDIF
      ENDDO
   ENDDO
   !
   IF (wan_mode == 'standalone') THEN
      IF (me_pool == root_pool) CLOSE(iun_band, STATUS="KEEP")
   ENDIF
   IF (wan_mode == 'library') CALL mp_sum(eigval, inter_pool_comm)
   !
   CALL mp_barrier(world_comm)
   !
   CALL utility_merge_files("eig", .TRUE.)
   !
END SUBROUTINE write_band

SUBROUTINE write_plot
   USE io_global,  ONLY : stdout, ionode
   USE mp_pools,        ONLY : me_pool, root_pool
   USE wvfct, ONLY : nbnd, npwx
   USE control_flags, ONLY : gamma_only
   USE wavefunctions, ONLY : evc
   USE io_files, ONLY : nwordwfc, iunwfc
   USE wannier
   USE klist,           ONLY : nkstot, xk, ngk, igk_k, nks
   USE gvect,           ONLY : g, ngm
   USE fft_base,        ONLY : dffts
   USE scatter_mod,     ONLY : gather_grid
   USE fft_interfaces,  ONLY : invfft
   USE noncollin_module,ONLY : noncolin, npol
   USE lsda_mod,        ONLY : lsda, isk
   !
   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: global_kpoint_index
   !
   CHARACTER(LEN=20) :: wfnname
   INTEGER :: ik, npw, ibnd, ik_g_w90, i1, j, spin, ipol, nxxs
   INTEGER :: nr1, nr2, nr3
   !! Real space grid sizes for the wavefunction data written to file
   INTEGER :: i, k, idx, pos, ierr
   !! aam: 1/5/06: for writing smaller unk files
   COMPLEX(DP), ALLOCATABLE :: evc_r(:, :)
   !! Distributed real-space wavefunction.
   COMPLEX(DP), ALLOCATABLE, TARGET :: psic_all(:, :)
   !! Gathered real-space wavefunction.
   COMPLEX(DP), POINTER :: psic_small(:, :)
   !! Real-space wavefunction written to file.
   !
   CALL start_clock( 'write_unk' )
   !
   IF (reduce_unk .AND. (reduce_unk_factor == 1)) THEN
      WRITE(stdout, *)
      WRITE(stdout, '(5x,a)') "WARNING: reduce_unk is .TRUE. but reduce_unk_factor is 1."
      WRITE(stdout, '(5x,a)') "The UNK file size is not reduced."
      WRITE(stdout, '(5x,a)') "To enable reduction, set reduce_unk_factor to a value greater than 1."
      WRITE(stdout, *)
   ENDIF
   IF ((.NOT. reduce_unk) .AND. (reduce_unk_factor > 1)) THEN
      WRITE(stdout, *)
      WRITE(stdout, '(5x,a)') "WARNING: reduce_unk_factor is not 1 but reduce_unk is .FALSE."
      WRITE(stdout, '(5x,a)') "The UNK file size is not reduced."
      WRITE(stdout, '(5x,a)') "To enable reduction, set reduce_unk to .TRUE."
      WRITE(stdout, *)
   ENDIF
   !
   nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
   ALLOCATE(psic_all(nxxs, npol), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psic_all', 1)
   ALLOCATE(evc_r(dffts%nnr, npol), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_r', 1)
   psic_all = (0.0_DP, 0.0_DP)
   evc_r = (0.0_DP, 0.0_DP)
   !
   IF (reduce_unk) THEN
      WRITE(stdout,'(3(a,i5))') 'nr1s=',dffts%nr1,' nr2s=',dffts%nr2,' nr3s=',dffts%nr3
      nr1 = (dffts%nr1 - 1) / reduce_unk_factor + 1
      nr2 = (dffts%nr2 - 1) / reduce_unk_factor + 1
      nr3 = (dffts%nr3 - 1) / reduce_unk_factor + 1
      WRITE(stdout,'(3(a,i5))') 'n1by2=', nr1, ' n2by2=', nr2, ' n3by2=', nr3
      ALLOCATE(psic_small(nr1*nr2*nr3, npol), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating psic_small', 1)
   ELSE
      psic_small => psic_all
      nr1 = dffts%nr1
      nr2 = dffts%nr2
      nr3 = dffts%nr3
   ENDIF
   !
   WRITE(stdout, '(a,i8)') '  Number of local k points = ', nks
   !
   DO ik = 1, nks
      !
      CALL print_progress(ik, nks)
      !
      IF (lsda .AND. isk(ik) /= ispinw) CYCLE
      ik_g_w90 = global_kpoint_index(nkstot, ik) - ikstart + 1
      !
      npw = ngk(ik)
      !
      IF (.NOT.noncolin) THEN
         spin = ispinw
         IF (ispinw == 0) spin = 1
         WRITE(wfnname, '("UNK", i5.5, ".", i1)') ik_g_w90, spin
      ELSE
         WRITE(wfnname, '("UNK", i5.5, ".", "NC")') ik_g_w90
      ENDIF
      !
      IF (me_pool == root_pool) THEN
         IF(wvfn_formatted) THEN
            OPEN(NEWUNIT=iun_plot, FILE=wfnname, FORM='formatted')
            WRITE(iun_plot, *) nr1, nr2, nr3, ik_g_w90, num_bands
         ELSE
            OPEN(NEWUNIT=iun_plot, FILE=wfnname, FORM='unformatted')
            WRITE(iun_plot) nr1, nr2, nr3, ik_g_w90, num_bands
         ENDIF
      ENDIF
      !
      CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1 )
      !
      DO ibnd = 1, nbnd
         IF (excluded_band(ibnd)) CYCLE
         !
         ! Transform wavefunction to real space
         !
         evc_r(:, :) = (0.d0, 0.d0)
         DO ipol = 1, npol
            evc_r(dffts%nl(igk_k(1:npw,ik)), ipol) = evc(1+npwx*(ipol-1):npw+npwx*(ipol-1), ibnd)
            CALL invfft('Wave', evc_r(:, ipol), dffts)
         ENDDO
         !
#if defined(__MPI)
         DO ipol = 1, npol
            CALL gather_grid(dffts, evc_r(:, ipol), psic_all(:,ipol))
         ENDDO
#else
         psic_all(1:dffts%nnr, :) = evc_r(1:dffts%nnr, :)
#endif
         !
         IF (reduce_unk) THEN
            psic_small = (0.0_DP, 0.0_DP)
            pos = 0
            DO k = 1, dffts%nr3, reduce_unk_factor
               DO j = 1, dffts%nr2, reduce_unk_factor
                  DO i = 1, dffts%nr1, reduce_unk_factor
                     idx = (k-1)*dffts%nr2*dffts%nr1 + (j-1)*dffts%nr1 + i
                     pos = pos + 1
                     DO ipol = 1, npol
                        psic_small(pos, ipol) = psic_all(idx, ipol)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         !
         IF (me_pool == root_pool) THEN
            DO ipol = 1, npol
               IF(wvfn_formatted) THEN
                  WRITE(iun_plot, '(2ES20.10)') (psic_small(j, ipol), j = 1, nr1*nr2*nr3)
               ELSE
                  WRITE(iun_plot) (psic_small(j, ipol), j = 1, nr1*nr2*nr3)
               ENDIF
            ENDDO
         ENDIF
      ENDDO !ibnd
      !
      IF (me_pool == root_pool) CLOSE(iun_plot)
      !
   ENDDO  !ik
   !
   IF (reduce_unk) THEN
      DEALLOCATE(psic_small)
   ENDIF
   DEALLOCATE(evc_r)
   DEALLOCATE(psic_all)
   !
   WRITE(stdout, *) ' UNK written'
   !
   CALL stop_clock( 'write_unk' )
   !
   RETURN
END SUBROUTINE write_plot

SUBROUTINE write_parity
   USE mp_pools,             ONLY : intra_pool_comm, me_pool, root_pool, nproc_pool
   USE mp,                   ONLY : mp_sum
   USE io_global,            ONLY : stdout, ionode
   USE wvfct,                ONLY : nbnd
   USE control_flags,        ONLY : gamma_only
   USE wavefunctions,        ONLY : evc
   USE io_files,             ONLY : nwordwfc, iunwfc
   USE klist,                ONLY : nkstot, xk, igk_k, ngk, nks
   USE gvect,                ONLY : g, ngm, mill
   USE cell_base,            ONLY : at
   USE constants,            ONLY : eps6
   USE lsda_mod,             ONLY : lsda, isk
   USE wannier
   !
   IMPLICIT NONE
   !
   INTEGER :: npw, ibnd, ig, kgamma, ik, i, ig_target, num_G, ierr
   !!
   INTEGER :: g_target(3, 32)
   !! List of G vectors to find
   COMPLEX(KIND=DP), ALLOCATABLE :: evc_target(:, :)
   !! evc values at the target G vectors
   !
   CALL start_clock('write_parity')
   !
   WRITE(stdout, *) "Finding the 32 unkg's per band required for parity signature."
   !
   ! List of target G vectors
   g_target(:,  1) = (/  0,  0,  0 /) ! 1
   g_target(:,  2) = (/  1,  0,  0 /) ! x
   g_target(:,  3) = (/  0,  1,  0 /) ! y
   g_target(:,  4) = (/  0,  0,  1 /) ! z
   g_target(:,  5) = (/  2,  0,  0 /) ! x^2
   g_target(:,  6) = (/  1,  1,  0 /) ! xy
   g_target(:,  7) = (/  1, -1,  0 /) ! xy
   g_target(:,  8) = (/  1,  0,  1 /) ! xz
   g_target(:,  9) = (/  1,  0, -1 /) ! xz
   g_target(:, 10) = (/  0,  2,  0 /) ! y^2
   g_target(:, 11) = (/  0,  1,  1 /) ! yz
   g_target(:, 12) = (/  0,  1, -1 /) ! yz
   g_target(:, 13) = (/  0,  0,  2 /) ! z^2
   g_target(:, 14) = (/  3,  0,  0 /) ! x^3
   g_target(:, 15) = (/  2,  1,  0 /) ! x^2y
   g_target(:, 16) = (/  2, -1,  0 /) ! x^2y
   g_target(:, 17) = (/  2,  0,  1 /) ! x^2z
   g_target(:, 18) = (/  2,  0, -1 /) ! x^2z
   g_target(:, 19) = (/  1,  2,  0 /) ! xy^2
   g_target(:, 20) = (/  1, -2,  0 /) ! xy^2
   g_target(:, 21) = (/  1,  1,  1 /) ! xyz
   g_target(:, 22) = (/  1,  1, -1 /) ! xyz
   g_target(:, 23) = (/  1, -1,  1 /) ! xyz
   g_target(:, 24) = (/  1, -1, -1 /) ! xyz
   g_target(:, 25) = (/  1,  0,  2 /) ! xz^2
   g_target(:, 26) = (/  1,  0, -2 /) ! xz^2
   g_target(:, 27) = (/  0,  3,  0 /) ! y^3
   g_target(:, 28) = (/  0,  2,  1 /) ! y^2z
   g_target(:, 29) = (/  0,  2, -1 /) ! y^2z
   g_target(:, 30) = (/  0,  1,  2 /) ! yz^2
   g_target(:, 31) = (/  0,  1, -2 /) ! yz^2
   g_target(:, 32) = (/  0,  0,  3 /) ! z^3
   !
   ! getting the ik index corresponding to the Gamma point
   ! ... and the spin channel (fix due to N Poilvert, Feb 2011)
   !
   IF (.NOT. gamma_only) THEN
      kgamma = -1
      DO ik = 1, nks
         IF (lsda .AND. isk(ik) /= ispinw) CYCLE
         IF (ALL(ABS(xk(:, ik)) < eps6)) THEN
            kgamma = ik
            EXIT
         ENDIF
         IF (ik == ikstop) CALL errore('write_parity',&
              ' parity calculation may only be performed at the gamma point',1)
      ENDDO
   ELSE
      ! NP: spin unpolarized or "up" component of spin
      IF (ispinw == 0 .OR. ispinw == 1) THEN
         kgamma = 1
      ELSE ! NP: "down" component
         kgamma = 2
      ENDIF
   ENDIF
   !
   ! Run calculation only on the pool with k = Gamma.
   IF (kgamma /= -1) THEN
      ALLOCATE(evc_target(32, nbnd), stat=ierr)
      IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating evc_target', 1)
      evc_target = (0.d0, 0.d0)
      !
      ! building the evc array corresponding to the Gamma point
      !
      CALL davcio(evc, 2*nwordwfc, iunwfc, kgamma, -1)
      npw = ngk(kgamma)
      !
      ! Count and identify the G vectors we will be extracting for each cpu.
      ! Fill evc_target with required fourier component from each cpu dependent evc
      !
      num_G = 0
      DO ig = 1, npw
         DO ig_target = 1, 32
            IF ( ALL(mill(:, igk_k(ig, kgamma)) == g_target(:, ig_target)) ) THEN
               num_G = num_G + 1
               evc_target(ig_target, :) = evc(ig, :)
               EXIT
            ENDIF
         ENDDO
      ENDDO
      CALL mp_sum(evc_target, intra_pool_comm)
      !
      ! Sum laterally across cpus num_G, so it contains
      ! the number of g_vectors on each node, and known to all cpus
      ! Check if all target G vectors are found.
      !
      CALL mp_sum(num_G, intra_pool_comm)
      IF (num_G /= 32) CALL errore('write_parity', 'incorrect number of g-vectors extracted',1)
      !
      ! Write to file
      !
      IF (me_pool == root_pool) THEN
         OPEN(NEWUNIT=iun_parity, FILE=TRIM(seedname)//".unkg", FORM='formatted')
         WRITE(iun_parity, *) num_G ! this value is always 32
         DO ibnd = 1, nbnd
            DO ig_target = 1, 32
               WRITE(iun_parity, '(5i5,2f12.7)') ibnd, ig_target, g_target(:, ig_target), &
                                                 REAL(evc_target(ig_target, ibnd)), &
                                                 AIMAG(evc_target(ig_target, ibnd))
            ENDDO
         ENDDO
         CLOSE(iun_parity, STATUS="KEEP")
      ENDIF
      !
      DEALLOCATE(evc_target)
      !
   ENDIF ! kgamma /= -1
   !
   CALL stop_clock('write_parity')
   !
END SUBROUTINE write_parity


SUBROUTINE wan2sic

  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE io_files, ONLY : iunwfc, nwordwfc, nwordwann
  USE gvect, ONLY : g, ngm
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : nbnd, npwx
  USE gvecw, ONLY : gcutw
  USE klist, ONLY : nkstot, xk, wk, ngk
  USE wannier

  IMPLICIT NONE

  INTEGER :: i, j, nn, ik, ibnd, iw, ikevc, npw, ierr
  COMPLEX(DP), ALLOCATABLE :: orbital(:,:), u_matrix(:,:,:)
  INTEGER :: iunatsicwfc = 31 ! unit for sic wfc

  OPEN (20, file = trim(seedname)//".dat" , form = 'formatted', status = 'unknown')
  WRITE(stdout,*) ' wannier plot '

  ALLOCATE ( u_matrix( n_wannier, n_wannier, nkstot), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating u_matrix', 1)
  ALLOCATE ( orbital( npwx, n_wannier), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating orbital', 1)

  !
  DO i = 1, n_wannier
     DO j = 1, n_wannier
        DO ik = 1, nkstot
           READ (20, * ) u_matrix(i,j,ik)
           !do nn = 1, nnb(ik)
           DO nn = 1, nnb
              READ (20, * ) ! m_matrix (i,j,nkp,nn)
           ENDDO
        ENDDO  !nkp
     ENDDO !j
  ENDDO !i
  !
  DO ik=1,iknum
     ikevc = ik + ikstart - 1
     CALL davcio (evc, 2*nwordwfc, iunwfc, ikevc, -1)
     npw = ngk(ik)
     WRITE(stdout,*) 'npw ',npw
     DO iw = 1, n_wannier
        DO j=1,npw
           orbital(j,iw) = (0.0d0,0.0d0)
           DO ibnd=1,n_wannier
              orbital(j,iw) = orbital(j,iw) + u_matrix(iw,ibnd,ik)*evc(j,ibnd)
              WRITE(stdout,*) j, iw, ibnd, ik, orbital(j,iw), &
                              u_matrix(iw,ibnd,ik), evc(j,ibnd)
           ENDDO !ibnd
        ENDDO  !j
     ENDDO !wannier
     CALL davcio (orbital, 2*nwordwann, iunatsicwfc, ikevc, +1)
  ENDDO ! k-points

  DEALLOCATE ( u_matrix)
  WRITE(stdout,*) ' dealloc u '
  DEALLOCATE (  orbital)
  WRITE(stdout,*) ' dealloc orbital '
  !
END SUBROUTINE wan2sic

SUBROUTINE ylm_expansion
   USE io_global,  ONLY : stdout
   USE kinds, ONLY :  DP
   USE random_numbers,  ONLY : randy
   USE matrix_inversion
   USE wannier
   IMPLICIT NONE
   ! local variables
   INTEGER, PARAMETER :: lmax2=16
   INTEGER ::  lm, i, ir, iw, m, ierr
   real(DP), ALLOCATABLE :: r(:,:), rr(:), rp(:,:), ylm_w(:), ylm(:,:), mly(:,:)
   real(DP) :: u(3,3)

   ALLOCATE (r(3,lmax2), rp(3,lmax2), rr(lmax2), ylm_w(lmax2), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating r/rp/rr/ylm_w', 1)
   ALLOCATE (ylm(lmax2,lmax2), mly(lmax2,lmax2), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating ylm/mly', 1)

   ! generate a set of nr=lmax2 random vectors
   DO ir=1,lmax2
      DO i=1,3
         r(i,ir) = randy() -0.5d0
      ENDDO
   ENDDO
   rr(:) = r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:)
   !- compute ylm(ir,lm)
   CALL ylmr2(lmax2, lmax2, r, rr, ylm)
   !- store the inverse of ylm(ir,lm) in mly(lm,ir)
   CALL invmat(lmax2, ylm, mly)
   !- check that r points are independent
   CALL check_inverse(lmax2, ylm, mly)

   DO iw=1, n_proj

      !- define the u matrix that rotate the reference frame
      CALL set_u_matrix (xaxis(:,iw),zaxis(:,iw),u)
      !- find rotated r-vectors
      rp(:,:) = MATMUL ( u(:,:) , r(:,:) )
      !- set ylm funtion according to wannier90 (l,mr) indexing in the rotaterd points
      CALL ylm_wannier(ylm_w,l_w(iw),mr_w(iw),rp,lmax2)

      csph(:,iw) = MATMUL (mly(:,:), ylm_w(:))

!      write (stdout,*)
!      write (stdout,'(2i4,2(2x,3f6.3))') l_w(iw), mr_w(iw), xaxis(:,iw), zaxis(:,iw)
!      write (stdout,'(16i6)')   (lm, lm=1,lmax2)
!      write (stdout,'(16f6.3)') (csph(lm,iw), lm=1,lmax2)

   ENDDO
   DEALLOCATE (r, rp, rr, ylm_w, ylm, mly )

   RETURN
END SUBROUTINE ylm_expansion

SUBROUTINE check_inverse(lmax2, ylm, mly)
   USE kinds, ONLY :  DP
   USE constants, ONLY :  eps8
   IMPLICIT NONE
   ! I/O variables
   INTEGER :: lmax2
   real(DP) :: ylm(lmax2,lmax2), mly(lmax2,lmax2)
   ! local variables
   real(DP), ALLOCATABLE :: uno(:,:)
   real(DP) :: capel
   INTEGER :: lm, ierr
   !
   ALLOCATE (uno(lmax2,lmax2), stat=ierr)
   IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating uno', 1)
   uno = MATMUL(mly, ylm)
   capel = 0.d0
   DO lm = 1, lmax2
      uno(lm,lm) = uno(lm,lm) - 1.d0
   ENDDO
   capel = capel + sum ( abs(uno(1:lmax2,1:lmax2) ) )
!   write (stdout,*) "capel = ", capel
   IF (capel > eps8) CALL errore('ylm_expansion', &
                    ' inversion failed: r(*,1:nr) are not all independent !!',1)
   DEALLOCATE (uno)
   RETURN
END SUBROUTINE check_inverse

SUBROUTINE set_u_matrix(x,z,u)
   USE kinds, ONLY :  DP
   USE constants, ONLY : eps6
   IMPLICIT NONE
   ! I/O variables
   real(DP) :: x(3),z(3),u(3,3)
   ! local variables
   real(DP) :: xx, zz, y(3), coseno

   xx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
   IF (xx < eps6) CALL errore ('set_u_matrix',' |xaxis| < eps ',1)
!   x(:) = x(:)/xx
   zz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
   IF (zz < eps6) CALL errore ('set_u_matrix',' |zaxis| < eps ',1)
!   z(:) = z(:)/zz

   coseno = (x(1)*z(1) + x(2)*z(2) + x(3)*z(3))/xx/zz
   IF (abs(coseno) > eps6) CALL errore('set_u_matrix',' xaxis and zaxis are not orthogonal !',1)

   y(1) = (z(2)*x(3) - x(2)*z(3))/xx/zz
   y(2) = (z(3)*x(1) - x(3)*z(1))/xx/zz
   y(3) = (z(1)*x(2) - x(1)*z(2))/xx/zz

   u(1,:) = x(:)/xx
   u(2,:) = y(:)
   u(3,:) = z(:)/zz

!   write (stdout,'(3f10.7)') u(:,:)

   RETURN

END SUBROUTINE set_u_matrix

SUBROUTINE ylm_wannier(ylm,l,mr,r,nr)
!
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
! of the spherical harmonic identified  by indices (l,mr)
! in table 3.1 of the wannierf90 specification.
!
! No reference to the particular ylm ordering internal to Quantum ESPRESSO
! is assumed.
!
! If ordering in wannier90 code is changed or extended this should be the
! only place to be modified accordingly
!
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi, fpi, eps8
   IMPLICIT NONE
! I/O variables
!
   INTEGER :: l, mr, nr
   real(DP) :: ylm(nr), r(3,nr)
!
! local variables
!
   real(DP), EXTERNAL :: s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy
   real(DP), EXTERNAL :: fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(DP) :: rr, cost, phi
   INTEGER :: ir
   real(DP) :: bs2, bs3, bs6, bs12
   bs2 = 1.d0/sqrt(2.d0)
   bs3=1.d0/sqrt(3.d0)
   bs6 = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)
!
   IF (l > 3 .or. l < -5 ) CALL errore('ylm_wannier',' l out of range ', 1)
   IF (l>=0) THEN
      IF (mr < 1 .or. mr > 2*l+1) CALL errore('ylm_wannier','mr out of range' ,1)
   ELSE
      IF (mr < 1 .or. mr > abs(l)+1 ) CALL errore('ylm_wannier','mr out of range',1)
   ENDIF

   DO ir=1, nr
      rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      IF (rr < eps8) CALL errore('ylm_wannier',' rr too small ',1)

      cost =  r(3,ir) / rr
      !
      !  beware the arc tan, it is defined modulo pi
      !
      IF (r(1,ir) > eps8) THEN
         phi = atan( r(2,ir)/r(1,ir) )
      ELSEIF (r(1,ir) < -eps8 ) THEN
         phi = atan( r(2,ir)/r(1,ir) ) + pi
      ELSE
         phi = sign( pi/2.d0,r(2,ir) )
      ENDIF


      IF (l==0) THEN   ! s orbital
                    ylm(ir) = s(cost,phi)
      ENDIF
      IF (l==1) THEN   ! p orbitals
         IF (mr==1) ylm(ir) = p_z(cost,phi)
         IF (mr==2) ylm(ir) = px(cost,phi)
         IF (mr==3) ylm(ir) = py(cost,phi)
      ENDIF
      IF (l==2) THEN   ! d orbitals
         IF (mr==1) ylm(ir) = dz2(cost,phi)
         IF (mr==2) ylm(ir) = dxz(cost,phi)
         IF (mr==3) ylm(ir) = dyz(cost,phi)
         IF (mr==4) ylm(ir) = dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = dxy(cost,phi)
      ENDIF
      IF (l==3) THEN   ! f orbitals
         IF (mr==1) ylm(ir) = fz3(cost,phi)
         IF (mr==2) ylm(ir) = fxz2(cost,phi)
         IF (mr==3) ylm(ir) = fyz2(cost,phi)
         IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
         IF (mr==5) ylm(ir) = fxyz(cost,phi)
         IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
         IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
      ENDIF
      IF (l==-1) THEN  !  sp hybrids
         IF (mr==1) ylm(ir) = bs2 * ( s(cost,phi) + px(cost,phi) )
         IF (mr==2) ylm(ir) = bs2 * ( s(cost,phi) - px(cost,phi) )
      ENDIF
      IF (l==-2) THEN  !  sp2 hybrids
         IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
      ENDIF
      IF (l==-3) THEN  !  sp3 hybrids
         IF (mr==1) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+p_z(cost,phi))
         IF (mr==2) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-p_z(cost,phi))
         IF (mr==3) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-p_z(cost,phi))
         IF (mr==4) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+p_z(cost,phi))
      ENDIF
      IF (l==-4) THEN  !  sp3d hybrids
         IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
         IF (mr==4) ylm(ir) = bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
         IF (mr==5) ylm(ir) =-bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
      ENDIF
      IF (l==-5) THEN  ! sp3d2 hybrids
         IF (mr==1) ylm(ir) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
         IF (mr==2) ylm(ir) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
         IF (mr==3) ylm(ir) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
         IF (mr==4) ylm(ir) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = bs6*s(cost,phi)-bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
         IF (mr==6) ylm(ir) = bs6*s(cost,phi)+bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
      ENDIF

   ENDDO

   RETURN

END SUBROUTINE ylm_wannier

!======== l = 0 =====================================================================
FUNCTION s(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) :: s, cost,phi
   s = 1.d0/ sqrt(fpi)
   RETURN
END FUNCTION s
!======== l = 1 =====================================================================
FUNCTION p_z(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::p_z, cost,phi
   p_z =  sqrt(3.d0/fpi) * cost
   RETURN
END FUNCTION p_z
FUNCTION px(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::px, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   px =  sqrt(3.d0/fpi) * sint * cos(phi)
   RETURN
END FUNCTION px
FUNCTION py(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py =  sqrt(3.d0/fpi) * sint * sin(phi)
   RETURN
END FUNCTION py
!======== l = 2 =====================================================================
FUNCTION dz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dz2, cost, phi
   dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)
   RETURN
END FUNCTION dz2
FUNCTION dxz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)
   RETURN
END FUNCTION dxz
FUNCTION dyz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz =  sqrt(15.d0/fpi) * sint*cost * sin(phi)
   RETURN
END FUNCTION dyz
FUNCTION dx2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)
   RETURN
END FUNCTION dx2my2
FUNCTION dxy(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)
   RETURN
END FUNCTION dxy
!======== l = 3 =====================================================================
FUNCTION fz3(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fz3, cost, phi
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
   RETURN
END FUNCTION fz3
FUNCTION fxz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
   RETURN
END FUNCTION fxz2
FUNCTION fyz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
   RETURN
END FUNCTION fyz2
FUNCTION fzx2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
   RETURN
END FUNCTION fzx2my2
FUNCTION fxyz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
   RETURN
END FUNCTION fxyz
FUNCTION fxx2m3y2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
   RETURN
END FUNCTION fxx2m3y2
FUNCTION fy3x2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
   RETURN
END FUNCTION fy3x2my2
!
!
!-----------------------------------------------------------------------
SUBROUTINE radialpart(ng, q, alfa, rvalue, lmax, radial)
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform
  ! of the radial functions.
  !
  USE kinds,      ONLY : dp
  USE constants,  ONLY : fpi
  USE cell_base,  ONLY : omega
  !
  IMPLICIT NONE
  ! I/O
  INTEGER :: ng, rvalue, lmax
  real(DP) :: q(ng), alfa, radial(ng,0:lmax)
  ! local variables
  real(DP), PARAMETER :: xmin=-6.d0, dx=0.025d0, rmax=10.d0

  real(DP) :: rad_int, pref, x
  INTEGER :: l, lp1, ir, ig, mesh_r, ierr
  real(DP), ALLOCATABLE :: bes(:), func_r(:), r(:), rij(:), aux(:)

  mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
  ALLOCATE ( bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating bes/func_r/r/rij', 1)
  ALLOCATE ( aux(mesh_r), stat=ierr)
  IF (ierr /= 0) CALL errore('pw2wannier90', 'Error allocating aux', 1)
  !
  !    compute the radial mesh
  !
  DO ir = 1, mesh_r
     x = xmin  + dble (ir - 1) * dx
     r (ir) = exp (x) / alfa
     rij (ir) = dx  * r (ir)
  ENDDO
  !
  IF (rvalue==1) func_r(:) = 2.d0 * alfa**(3.d0/2.d0) * exp(-alfa*r(:))
  IF (rvalue==2) func_r(:) = 1.d0/sqrt(8.d0) * alfa**(3.d0/2.d0) * &
                     (2.0d0 - alfa*r(:)) * exp(-alfa*r(:)*0.5d0)
  IF (rvalue==3) func_r(:) = sqrt(4.d0/27.d0) * alfa**(3.0d0/2.0d0) * &
                     (1.d0 - 2.0d0/3.0d0*alfa*r(:) + 2.d0*(alfa*r(:))**2/27.d0) * &
                                           exp(-alfa*r(:)/3.0d0)
  pref = fpi/sqrt(omega)
  !
  DO l = 0, lmax
     DO ig=1,ng
       CALL sph_bes (mesh_r, r(1), q(ig), l, bes)
       aux(:) = bes(:) * func_r(:) * r(:) * r(:)
       ! second r factor added upo suggestion by YY Liang
       CALL simpson (mesh_r, aux, rij, rad_int)
       radial(ig,l) = rad_int * pref
     ENDDO
  ENDDO

  DEALLOCATE (bes, func_r, r, rij, aux )
  RETURN
END SUBROUTINE radialpart
