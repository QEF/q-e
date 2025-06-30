  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-----------------------------------------------------------------------
  SUBROUTINE ep_coarse_unfolding(nqc, xqc, w_centers)
  !-----------------------------------------------------------------------
  !!
  !! Electron-phonon calculation with Wannier functions: load all phonon q's
  !!
  !! This routine is the main driver of the electron-phonon
  !! calculation. It first calculates the electron-phonon matrix elements
  !! on the coarse mesh and then passes the data off to [[ephwann_shuffle]]
  !! to perform the interpolation.
  !!
  !! SP - Apr 2021 - Addition of quadrupoles
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE mp_global,     ONLY : my_pool_id, world_comm, npool
  USE mp_images,     ONLY : my_image_id, nimage
  USE mp_world,      ONLY : mpime
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_bcast
  USE io_global,     ONLY : stdout, meta_ionode, meta_ionode_id, ionode_id,     &
                            ionode
  USE qrad_mod,      ONLY : init_tab_qrad, deallocate_tab_qrad
  USE gvect,         ONLY : gcutm
  USE cellmd,        ONLY : cell_factor
  USE io_files,      ONLY : prefix, tmp_dir, diropn
  USE wavefunctions, ONLY : evc
  USE wvfct,         ONLY : npwx
  USE eqv,           ONLY : vlocq, dmuxc
  USE ions_base,     ONLY : nat, tau, ityp, amass
  USE control_flags, ONLY : iverbosity
  USE io_var,        ONLY : iuepb, iuqpeig, crystal, iunpattern, iuquad, iuxqc, &
                            iuqmap
  USE pwcom,         ONLY : nks, nbnd, nkstot, nelec
  USE cell_base,     ONLY : at, bg, alat, omega, tpiba
  USE symm_base,     ONLY : irt, s, nsym, ft, sname, invs, s_axis_to_cart,      &
                            sr, nrot, set_sym_bl, find_sym, inverse_s, t_rev,   &
                            time_reversal
  USE phcom,         ONLY : evq
  USE qpoint,        ONLY : igkq, xq, eigqts
  USE modes,         ONLY : nmodes, u, npert, nirr
  USE lr_symm_base,  ONLY : minus_q, rtau, gi, gimq, irotmq, nsymq, invsymq
  USE input,         ONLY : epbread, epbwrite, epwread, lifc,                   &
                            nbndsub, iswitch, kmaps, eig_read, dvscf_dir,       &
                            nkc1, nkc2, nkc3, nqc1, nqc2, nqc3,                 &
                            fixsym, epw_noinv, system_2d, compute_dmat, lwfpt,  &
                            exciton             
                            ! ZD: last line for ex-plrn
  USE global_var,    ONLY : epmatq, dynq, et_ks, xkq, ifc, veff,&
                            zstar, epsi, cu, cuq, lwin, lwinq, nbndep, ibndkept,&
                            ngxx, exband, wscache, area, ngxxf, shift,  &
                            gmap, qrpl, Qmat, L, do_cutoff_2D_epw, alph, nbndskip,&
                            wf_temp, wf                                         
                            ! ZD: last line for ex-plrn
  USE input,         ONLY : et_loc, et_all, xk_all
  USE ep_constants,  ONLY : ryd2ev, zero, two, czero, eps6, eps8
  USE fft_base,      ONLY : dfftp
  USE control_ph,    ONLY : u_from_file
  USE stop,          ONLY : stop_epw
  USE parallelism,   ONLY : fkbounds
  USE uspp,          ONLY : okvan
  USE lrus,          ONLY : becp1
  USE becmod,        ONLY : deallocate_bec_type
  USE phus,          ONLY : int1, int1_nc, int2, int2_so, alphap
  USE kfold,         ONLY : createkmap_pw2, createkmap
  USE low_lvl,       ONLY : set_ndnmbr, eqvect_strict, copy_sym_epw, fix_sym
  USE ph_restart,    ONLY : read_disp_pattern_only
  USE io,            ONLY : read_ifc_epw, readgmap, dynmat_asr, loadumat
  USE parallelism,   ONLY : poolgather
  USE symmetry,      ONLY : rotate_epmat, rotate_eigenm, star_q2, gmap_sym, &
                            setup_rotate_wavefunction, write_qmap,          &
                            calc_rotation_gauge, rotate_wfn_deallocate
  USE pw2wan,        ONLY : compute_pmn_para
  USE ep_coarse,     ONLY : ep_coarse_compute, dvscf_read
  USE noncollin_module, ONLY : m_loc, npol, noncolin, lspinorb, nspin_mag
  USE exphonon,      ONLY : calc_G_epmat_subbnd, find_iq, prepare_ex_ph_g, release_ex_ph_g
  USE klist,         ONLY : degauss, ngauss
  ! ZD: the last line for ex-ph calculations

#if defined(__NAG)
  USE f90_unix_io,   ONLY : flush
#endif
  !
  ! --------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(inout) :: nqc
  !! number of qpoints in the coarse grid
  REAL(KIND = DP), INTENT(inout) :: xqc(3, nqc1 * nqc2 * nqc3)
  !! qpoint list, coarse mesh
  REAL(KIND = DP), INTENT(inout) :: w_centers(3, nbndsub)
  !! Wannier centers
  !
  CHARACTER(LEN = 256) :: tempfile
  !! Temporary .eig file
  CHARACTER(LEN = 256) :: dirname
  !! Name of the directory
  CHARACTER(LEN = 256) :: filename
  !! Name of the file
  CHARACTER(LEN = 256) :: dummy
  !! Dummy character reading
  CHARACTER(LEN = 4) :: filelab
  !! Append the number of the core that works on that file
  CHARACTER(LEN = 80)   :: line
  !! Use to read external eigenvalues
  CHARACTER(LEN = 6), EXTERNAL :: int_to_char
  !! Transform an INTEGER into a character
  LOGICAL :: timerev
  !!  true if we are using time reversal
  LOGICAL :: sym(48)
  !! Logical vectors that say which crystal symmetries exist in our system
  LOGICAL :: nog
  !! Find if G=0 or not in $$S(q_0)+G=q$$
  LOGICAL :: non_symmorphic
  !! Check whether the symmetry belongs to a non-symmorphic group
  !! non_symmorphic == TRUE if it has fractional translation.
  LOGICAL :: exst
  !! Find if a file exists.
  LOGICAL :: kmesh_symmetry
  !! Whether the k mesh have the symmetry
  LOGICAL, EXTERNAL :: check_q_points_sym
  !! Check symmetry of the mesh (Here, applied to the k point mesh)
  INTEGER :: sym_smallq(48)
  !! Set of all symmetries for the small group of one q.
  !! This is a subset of total crystal symmetries that remains
  !! after the q-point pertubation.
  INTEGER :: imode0
  !! Counter on modes
  INTEGER :: nqc_irr
  !! Number of qpoints in the irreducible wedge
  INTEGER :: irr
  !! Counter on representations
  INTEGER :: npe
  !! Number of perturbations for irr representation
  INTEGER :: ik
  !! Total k-point index
  INTEGER :: ios
  !! Contains the state of the opened file
  INTEGER :: ik_start
  !! Lower bound for the k-point of the coarse grid in parallel
  INTEGER :: ik_stop
  !! Higher bound for the k-point of the coarse grid in parallel
  INTEGER :: nq
  !! Degeneracy of the star of q
  INTEGER :: isq(48)
  !! Index of q in the star of a given sym.op.
  INTEGER :: imq
  !! Index of -q in the star of q (0 if not present)
  INTEGER :: sym_sgq(48)
  !! The symmetries giving the q point iq in the star
  INTEGER :: indsym(48)
  !! The correspondence between the original sym. indices and the reshuffled indices
  INTEGER :: i
  !! Index for the star of q points
  INTEGER :: j
  !! Cartesian index
  INTEGER :: iq
  !! q-point index
  INTEGER :: iq_irr
  !! Counter on irreducible q-points
  INTEGER :: isym
  !! Index of symmetry
  INTEGER :: isym1
  !! Index of symmetry
  INTEGER :: iq_first
  !! First q in the star of q
  INTEGER :: jsym
  !! Symmetry index
  INTEGER :: ism1
  !! Inverse of the symmetry
  INTEGER :: nsq
  !! The number of degeneracy of the small group for this iq in the star
  INTEGER :: ipol
  !! Polarization index
  INTEGER :: jpol
  !! Polarization index
  INTEGER :: ierr
  !! Error index when reading/writing a file
  INTEGER :: na
  !! Atom index
  INTEGER :: idir
  !! Cartesian direction
  !INTEGER :: iunpun
  !! Unit of the file
  INTEGER, ALLOCATABLE :: gmapsym(:, :)
  !! Correspondence G -> S(G)
  REAL(KIND = DP) :: sxq(3, 48)
  !! List of vectors in the star of q
  REAL(KIND = DP) :: et_tmp(nbnd, nkstot)
  !! Temporary array containing the eigenvalues (KS or GW) when read from files
  REAL(KIND = DP) :: xq0(3)
  !! Current coarse q-point coords.
  REAL(KIND = DP) :: aq(3)
  !! Store the current q-point for symmetry multiplication
  REAL(KIND = DP) :: saq(3)
  !! Rotated q-point
  REAL(KIND = DP) :: raq(3)
  !! Rotate q-point in cartesian coordinate
  REAL(KIND = DP) :: ft1
  !! Fractional translation x
  REAL(KIND = DP) :: ft2
  !! Fractional translation y
  REAL(KIND = DP) :: ft3
  !! Fractional translation z
  REAL(KIND = DP) :: qnorm_tmp
  !! Absolute value of xqc_irr
  REAL(KIND = DP) :: qmax
  !! Max value of q for interpoplation table
  REAL(KIND = DP) :: sumr(2, 3, nat, 3)
  !! Sum to impose the ASR
  REAL(KIND = DP) :: Qxx, Qyy, Qzz, Qyz, Qxz, Qxy
  !! Specific quadrupole value read from file.
  REAL(KIND = DP), ALLOCATABLE :: xqc_irr(:, :)
  !! The qpoints in the irr wedge
  REAL(KIND = DP), ALLOCATABLE :: wqlist(:)
  !! The corresponding weigths
  COMPLEX(KIND = DP) :: cz1(nmodes, nmodes)
  !! The eigenvectors for the first q in the star
  COMPLEX(KIND = DP) :: cz2(nmodes, nmodes)
  !! The rotated eigenvectors, for the current q in the star
  COMPLEX(KIND = DP), ALLOCATABLE :: eigv(:, :)
  !! $e^{ iGv}$ for 1...nsym (v the fractional translation)
  COMPLEX(KIND = DP), ALLOCATABLE :: dvscfin(:, :, :)
  !! Change of the scf potential
  COMPLEX(KIND = DP), ALLOCATABLE :: el_ph_mat(:, :, :, :)
  !! ep matrix element for this q-point
  INTEGER :: iq_order
  !! The index of current q in the order of the klist in nscf calculations
  REAL(KIND = DP) :: xq0_cryst(3)
  !! The current coords of iq in the crystal coordinate
  COMPLEX(KIND = DP) :: cz2_out(nmodes, nmodes)
  !! The rotated eigenvectors. 
  !! ZD: This will be set to equal cz2 before rotate_epmat
  !! to avoid possible mess-up from multiply mass factors
  !! Current coarse q-point coords in crystal coordinate.
  !
  CALL start_clock('elphon_wrap')
  !
#if defined(__HDF5)
  WRITE(stdout, '(5x, a)') 'HDF5 is used in the current build. &
                            Exciton-phonon coupling calculations are enabled.'
#else
  exciton = .FALSE.
  WRITE(stdout, '(5x, a)') 'HDF5 is NOT used in the current build. &
                            Exciton-phonon coupling calculations are disabled.'
#endif

  alph = 1.0
  !
  ! S. Tiwari: Initializing the Ewald parameter whose units are (2pi/alat)^{2}
  ! in case we don't use dynamical matrix read from IFC
  !
  IF (epwread .AND. .NOT. epbread) THEN
    CONTINUE
  ELSE
    CALL compute_pmn_para(include_nonlocal = .TRUE.)
    IF (lwfpt) CALL compute_pmn_para(include_nonlocal = .FALSE.)
    !
    ! Regenerate qpoint list
    !
    ALLOCATE(wqlist(nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating wqlist', 1)
    wqlist(:) = zero
    CALL kpoint_grid(nsym, time_reversal, .FALSE., s, t_rev, bg, nqc1 * nqc2 * nqc3, &
                     0, 0, 0, nqc1, nqc2, nqc3, nqc_irr, xqc, wqlist)
    ALLOCATE(xqc_irr(3, nqc_irr), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating xqc_irr', 1)
    xqc_irr(:, :) = zero
    DEALLOCATE(wqlist, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating wqlist', 1)
    xqc_irr(:, :) = xqc(:, 1:nqc_irr)
    xqc(:, :) = zero
    !
    ! ensure that the size of the interpolation table 'tab_qrad' for the Q(r)
    ! Q(r) of ultrasoft pseudopotentials is sufficient; re-initialize otherwise
    !
    IF (.NOT. epbread .AND. .NOT. epwread) THEN
      qnorm_tmp = 0.0_dp
      DO iq_irr = 1, nqc_irr
        qnorm_tmp = MAX( qnorm_tmp, xqc_irr(1, iq_irr)**2 + xqc_irr(2, iq_irr)**2 + xqc_irr(3, iq_irr)**2)
      ENDDO
      qmax = (SQRT(gcutm) + SQRT(qnorm_tmp)) * tpiba * cell_factor
      ! FIXME: I don't think cell_factor should be there
      CALL init_tab_qrad(qmax, omega, intra_bgrp_comm, ierr)

    ENDIF
    !
    IF (nkstot /= nkc1 * nkc2 * nkc3) CALL errore('ep_coarse_unfolding', 'nscf run inconsistent with epw input', 1)
    !
  ENDIF
  !
  ! Read in external electronic eigenvalues. e.g. GW
  !
  IF (.NOT. epwread) THEN
    ALLOCATE(et_ks(nbnd, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating et_ks', 1)
    et_ks(:, :) = zero
    IF (eig_read) THEN
      IF (meta_ionode) THEN
        WRITE(stdout, '(5x, a, i5, a, i5, a)') "Reading external electronic eigenvalues (", &
              nbnd, ",", nkstot,")"
        tempfile = TRIM(prefix) // '.eig'
        OPEN(iuqpeig, FILE = tempfile, FORM = 'formatted', ACTION = 'read', IOSTAT = ios)
        IF (ios /= 0) CALL errore('ep_coarse_unfolding', 'error opening' // tempfile, 1)
        READ(iuqpeig, '(a)') line
        DO ik = 1, nkstot
          ! We do not save the k-point for the moment ==> should be read and
          ! tested against the current one
          READ(iuqpeig, '(a)') line
          READ(iuqpeig, *) et_tmp(:, ik)
        ENDDO
        CLOSE(iuqpeig)
        ! from eV to Ryd
        et_tmp = et_tmp / ryd2ev
      ENDIF
      CALL mp_bcast(et_tmp, meta_ionode_id, world_comm)
      !
      CALL fkbounds(nkstot, ik_start, ik_stop)
      et_ks(:, :)  = et_loc(:, :)
      et_loc(:, :) = et_tmp(:, ik_start:ik_stop)
    ENDIF
  ELSE
    ! if starting from epwread, do not need to get external eigs from file.
    ! allocate zero sized array so no issues with deallocation at end of execution
    ALLOCATE(et_ks(0, 0), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating et_ks', 1)
  ENDIF
  !
  !  gather electronic eigenvalues for subsequent shuffle
  !
  IF (eig_read .AND. .NOT. epwread) THEN
    et_all(:, :) = zero
    CALL poolgather(nbnd, nkstot, nks, et_loc(1:nbnd, 1:nks), et_all)
  ENDIF
  !
  IF (.NOT. epbread .AND. .NOT. epwread) THEN
    IF (.NOT. kmaps) THEN
      CALL start_clock('kmaps')
      CALL createkmap_pw2
      CALL stop_clock('kmaps')
      CALL print_clock('kmaps')
    ELSE
      !
      ! 26/06/2012 RM
      ! if we do not have epmatq already on file then epbread=.FALSE.
      ! .kgmap is used from disk and .kmap is regenerated for each q-point
      !
      WRITE(stdout, '(/5x, a)') 'Using kmap and kgmap from disk'
      !
      ! gmap gets allocated inside readgmap
      !
      CALL readgmap(nkstot)
      !
    ENDIF
    !
    IF (iverbosity == 1) WRITE(stdout, 15) ngxx
15  FORMAT(5x,'Estimated size of gmap: ngxx =', i5)
    !
    ALLOCATE(gmapsym(ngxxf, nsym), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating gmapsym', 1)
    gmapsym = 0
    ALLOCATE(eigv(ngxxf, nsym), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating eigv', 1)
    eigv = czero
    !
  ENDIF
  !
  IF (epwread) THEN
    !
    ! We need some crystal info
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = crystal, FILE = 'crystal.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore('ep_coarse_unfolding', 'error opening crystal.fmt', crystal)
      READ(crystal, *) nat
      READ(crystal, *) nmodes
      READ(crystal, *) nelec, nbndskip
      READ(crystal, *) at
      READ(crystal, *) bg
      READ(crystal, *) omega
      READ(crystal, *) alat
      ALLOCATE(tau(3, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating tau', 1)
      READ(crystal, *) tau
      READ(crystal, *) amass
      ALLOCATE(ityp(nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating ityp', 1)
      READ(crystal, *) ityp
      READ(crystal, *) noncolin
      READ(crystal, *) do_cutoff_2D_epw
      READ(crystal, *) w_centers
      READ(crystal, *) L
      READ(crystal, *) degauss
      READ(crystal, *) ngauss
      !
    ENDIF ! mpime == ionode_id
    CALL mp_bcast(nat      , ionode_id, world_comm)
    IF (mpime /= ionode_id) ALLOCATE(ityp(nat))
    CALL mp_bcast(nmodes   , ionode_id, world_comm)
    CALL mp_bcast(nelec    , ionode_id, world_comm)
    CALL mp_bcast(nbndskip , ionode_id, world_comm)
    CALL mp_bcast(at       , ionode_id, world_comm)
    CALL mp_bcast(bg       , ionode_id, world_comm)
    CALL mp_bcast(omega    , ionode_id, world_comm)
    CALL mp_bcast(alat     , ionode_id, world_comm)
    IF (mpime /= ionode_id) ALLOCATE(tau(3, nat))
    CALL mp_bcast(tau      , ionode_id, world_comm)
    CALL mp_bcast(amass    , ionode_id, world_comm)
    CALL mp_bcast(ityp     , ionode_id, world_comm)
    CALL mp_bcast(noncolin , ionode_id, world_comm)
    CALL mp_bcast(do_cutoff_2D_epw , ionode_id, world_comm)
    CALL mp_bcast(w_centers, ionode_id, world_comm)
    CALL mp_bcast(L        , ionode_id, world_comm)
    CALL mp_bcast(degauss  , ionode_id, world_comm)
    CALL mp_bcast(ngauss   , ionode_id, world_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(crystal)
    ENDIF
  ENDIF ! epwread
  !
  IF (system_2d /= 'no') THEN
    area = omega * bg(3, 3) / alat
    WRITE(stdout, '(5x, a, F14.8, a)') 'Area is', area, ' [Bohr^2]'
  ENDIF
  !
  ALLOCATE(Qmat(nat, 3, 3, 3), STAT = ierr)
  IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating Qmat', 1)
  Qmat(:, :, :, :) = zero
  IF (qrpl) THEN
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iuquad, FILE = 'quadrupole.fmt', STATUS = 'old', IOSTAT = ios)
      READ(iuquad, *) dummy
      DO i = 1, 3 * nat
        READ(iuquad, *) na, idir, Qxx, Qyy, Qzz, Qyz, Qxz, Qxy
        Qmat(na, idir, 1, 1) = Qxx
        Qmat(na, idir, 2, 2) = Qyy
        Qmat(na, idir, 3, 3) = Qzz
        Qmat(na, idir, 2, 3) = Qyz
        Qmat(na, idir, 3, 2) = Qyz
        Qmat(na, idir, 1, 3) = Qxz
        Qmat(na, idir, 3, 1) = Qxz
        Qmat(na, idir, 1, 2) = Qxy
        Qmat(na, idir, 2, 1) = Qxy
      ENDDO
      CLOSE(iuquad)
    ENDIF ! mpime == ionode_id
    CALL mp_bcast(Qmat, ionode_id, world_comm)
    WRITE(stdout, '(a)') '     '
    WRITE(stdout, '(a)') '     ------------------------------------ '
    WRITE(stdout, '(a)') '     Quadrupole tensor is correctly read: '
    WRITE(stdout, '(a)') '     ------------------------------------ '
    WRITE(stdout, '(a)') '     atom   dir        Qxx       Qyy      Qzz        Qyz       Qxz       Qxy'
    DO na = 1, nat
      WRITE(stdout, '(i8, a,6f10.5)' ) na, '        x    ', Qmat(na, 1, 1, 1), Qmat(na, 1, 2, 2), Qmat(na, 1, 3, 3), &
                                                            Qmat(na, 1, 2, 3), Qmat(na, 1, 1, 3), Qmat(na, 1, 1, 2)
      WRITE(stdout, '(i8, a,6f10.5)' ) na, '        y    ', Qmat(na, 2, 1, 1), Qmat(na, 2, 2, 2), Qmat(na, 2, 3, 3), &
                                                            Qmat(na, 2, 2, 3), Qmat(na, 2, 1, 3), Qmat(na, 2, 1, 2)
      WRITE(stdout, '(i8, a,6f10.5)' ) na, '        z    ', Qmat(na, 3, 1, 1), Qmat(na, 3, 2, 2), Qmat(na, 3, 3, 3), &
                                                            Qmat(na, 3, 2, 3), Qmat(na, 3, 1, 3), Qmat(na, 3, 1, 2)
    ENDDO
    WRITE(stdout, '(a)') '     '
  ENDIF ! exst
  !
  IF (lifc) THEN
    ALLOCATE(ifc(nqc1, nqc2, nqc3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating ifc', 1)
    ifc(:, :, :, :, :, :, :) = zero
  ENDIF
  !
  ! SP: Symmetries needs to be consistent with QE so that the order of the q in the star is the
  !     same as in the .dyn files produced by QE.
  !
  ! Initialize symmetries and create the s matrix
  s(:, :, :) = 0 ! Symmetry in crystal axis with dim: 3,3,48
  CALL set_sym_bl()
  !
  ! Setup Bravais lattice symmetry
  WRITE(stdout,'(5x,a,i3)') "Symmetries of Bravais lattice: ", nrot
  !
  ! Setup crystal symmetry
  IF (.NOT. ALLOCATED(m_loc)) ALLOCATE(m_loc(3, nat))
  CALL find_sym(nat, tau, ityp, .FALSE., m_loc)
  IF (fixsym) CALL fix_sym(.FALSE.)
  WRITE(stdout, '(5x, a, i3)') "Symmetries of crystal:         ", nsym
  !
  IF (epwread .AND. .NOT. epbread) THEN
    CONTINUE
  ELSE
    !
    !  allocate dynamical matrix and ep matrix for all q's
    !
    ALLOCATE(dynq(nmodes, nmodes, nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating dynq', 1)
    ALLOCATE(epmatq(nbndep, nbndep, nks, nmodes, nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating epmatq', 1)
    ALLOCATE(epsi(3, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating epsi', 1)
    ALLOCATE(zstar(3, 3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating zstar', 1)
    ALLOCATE(cu(nbndep, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating cu', 1)
    ALLOCATE(cuq(nbndep, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating cuq', 1)
    ALLOCATE(lwin(nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating lwin', 1)
    ALLOCATE(lwinq(nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating lwinq', 1)
    ALLOCATE(exband(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating exband', 1)
    dynq(:, :, :)         = czero
    epmatq(:, :, :, :, :) = czero
    epsi(:, :)            = zero
    zstar(:, :, :)        = zero
    cu(:, :, :)           = czero
    cuq(:, :, :)          = czero
    sxq(:, :)             = zero
    !
    ! read interatomic force constat matrix from q2r
    IF (lifc) THEN
      CALL read_ifc_epw
    ENDIF
    !
    ! The following loop is required to propertly set up the symmetry matrix s.
    ! We here copy the calls made in PHonon/PH/init_representations.f90 to have the same s as in QE 5.
    DO iq_irr = 1, nqc_irr
      xq = xqc_irr(:, iq_irr)
      ! search for the small group of q
      CALL set_small_group_of_q(nsymq, invsymq, minus_q)
      ! calculate rtau with the new symmetry order
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      ! calculate the vectors G associated to the symmetry Sq = q + G
      ! if minus_q is true calculate also irotmq and the G associated to Sq=-g+G
      CALL set_giq(xq, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
    ENDDO
  ENDIF ! epwread .AND. .NOT. epbread
  !
  ! CV: if we read the .fmt files we don't need to read the .epb anymore
  !
  IF (.NOT. epbread .AND. .NOT. epwread) THEN
    !
    ALLOCATE(evq(npwx * npol, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating evq', 1)
    ALLOCATE(xkq(3, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating xkq', 1)
    ALLOCATE(shift(nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating shift', 1)
    IF (lifc) THEN
      ALLOCATE(wscache(-2 * nqc3:2 * nqc3, -2 * nqc2:2 * nqc2, -2 * nqc1:2 * nqc1, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating wscache', 1)
      wscache(:, :, :, :, :) = zero
    ENDIF
    evq(:, :) = zero
    xkq(:, :) = zero
    shift(:)  = 0
    !
    ! In the loop over irr q-point, we need to read the pattern that
    ! corresponds to the dvscf file computed with QE 5.
    !
    sumr(:, :, :, :) = zero
    iq_first = 1
    !
    DO i = 1, 48
      indsym(i) = i
    ENDDO
    !
    !  determine the G vector map S(G) -> G
    !
    CALL gmap_sym(nsym, s, ft, gmapsym, eigv)
    !
    ! ZD: initialize arrays for ex-ph matrix
    IF (exciton) THEN
      CALL prepare_ex_ph_g()
    ENDIF
    !
    ! Setup for compute_dmat
    !
    IF (compute_dmat) THEN
      !
      ! Check k mesh have symmetry. If not, raise error.
      kmesh_symmetry = check_q_points_sym(nkstot, xk_all, at, bg, nsym, &
                                          s, invs, nkc1, nkc2, nkc3)
      !
      IF (.NOT. kmesh_symmetry) THEN
        CALL errore('ep_coarse_unfolding', &
            'k mesh breaks the symmetry. Cannot compute dmat.', 1)
      ENDIF
      !
      ! Load lwin from ukk file
      !
      xq(:) = zero
      CALL loadumat(nbndep, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq, exband, w_centers)
      !
      ! Group et_all into degenerate groups
      !
      CALL setup_rotate_wavefunction(nsym, s, ft, invs)
      !
      ! Compute gauge matrix for symmetry operation
      !
      CALL calc_rotation_gauge(nsym, s, invs, gmapsym, eigv, lwin)
      !
      IF (meta_ionode) OPEN(iuqmap, FILE = TRIM(prefix) // '.qmap', FORM = 'formatted')
      !
    ENDIF
    !
    DO iq_irr = 1, nqc_irr
      u_from_file = .TRUE.
      !
      !  read the displacement patterns
      !
      IF (u_from_file) THEN
         ierr = 0
         dirname = TRIM(dvscf_dir) // TRIM(prefix) // '.phsave'
         filename = TRIM(dirname) // '/patterns.' // TRIM(int_to_char(iq_irr)) // '.xml'
         INQUIRE(FILE = TRIM(filename), EXIST = exst)
         IF (.NOT. exst) CALL errore('ep_coarse_unfolding', &
                   'cannot open file for reading or writing', 1)
         CALL read_disp_pattern_only(iunpattern, filename, iq_irr, ierr)
         IF (ierr /= 0) CALL errore('ep_coarse_unfolding', ' Problem with modes file', 1)
      ENDIF
      !
      WRITE(stdout, '(//5x, a)') REPEAT('=', 67)
      WRITE(stdout, '(5x, "irreducible q point # ", i4)') iq_irr
      WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
      FLUSH(stdout)
      !
      xq = xqc_irr(:, iq_irr)
      !
      ! SP : The following is largely inspiered by PHonon/PH/q2qstar.f90
      !
      ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~
      !
      minus_q = .TRUE.
      sym = .FALSE.
      sym(1:nsym) = .TRUE.
      CALL smallg_q(xq, 0, at, bg, nsym, s, sym, minus_q) ! s is intent(in)
      !
      nsymq = copy_sym_epw(nsym, sym, indsym)
      !
      ! Recompute the inverses as the order of sym.ops. has changed
      CALL inverse_s()
      CALL s_axis_to_cart()
      !
      ! This computes gi, gimq
      CALL set_giq(xq, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
      WRITE(stdout, '(5x, a, i3)') "Symmetries of small group of q:", nsymq
      IF(minus_q) WRITE(stdout, '(10x, a)') "in addition sym. q -> -q+G:"
      !
      ! Finally this does some of the above again and also computes rtau...
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      ! ######################### star of q #########################
      !
      sym_smallq(:) = 0
      CALL star_q2(xq, at, bg, nsym, s, invs, t_rev, nq, sxq, isq, imq, .TRUE., sym_smallq)
      IF (fixsym) THEN
        IF (epw_noinv) imq = 1 ! Any non-zero integer is ok.
      ENDIF
      !
      ! The reason for xq instead of xq0 in the above is because xq is passed to QE through module
      xq0 = xq
      !
      !  Re-set the variables needed for the pattern representation
      !  and the symmetries of the small group of irr-q
      !  (from phq_setup.f90)
      !
      DO isym = 1, nsym
        sym(isym) = .TRUE.
      ENDDO
      !
      ! ZD: allocate array for temporary store of phonon frequencies
      IF(exciton) THEN 
        ALLOCATE(wf_temp(nmodes,nq),STAT=ierr)
        IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error allocating wf_temp', 1) 
      ENDIF     
      !
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      CALL dynmat_asr(iq_irr, nqc_irr, nq, iq_first, sxq, imq, isq, invs, s, irt, rtau, sumr)
      !
      ! now dynq is the cartesian dyn mat (not divided by the masses)
      !
      minus_q = (iswitch > -3)
      !
      !  loop over the q points of the star
      !
      DO iq = 1, nq
        ! SP: First the vlocq needs to be initialized properly with the first
        !     q in the star
        xq = xq0
        CALL init(.FALSE.)
        !
        ! retrieve the q in the star
        xq = sxq(:, iq)
        !
        ! and populate the uniform grid
        nqc = nqc + 1
        xqc(:, nqc) = xq
        !
        IF (iq == 1) WRITE(stdout, *)
        WRITE(stdout, 5) nqc, xq
        ! ZD: determine the nscf-iq index of the current (iq_irr, iq) pair
        IF(exciton) THEN
          xq0_cryst = xq
          CALL cryst_to_cart(1, xq0_cryst, at, -1)
          CALL find_iq(xq0_cryst, nqc1, nqc2, nqc3, iq_order)
          ! WRITE(stdout, '(5x, a, I5, 3F10.5)') 'Current q in cryst. coord.: ', iq_order, xq0_cryst
          wf(:,iq_order)=wf_temp(:,iq)
        ENDIF
        !
        ! Prepare the kmap for the refolding
        !
        CALL createkmap(xq)
        !
        IF (iverbosity == 1) THEN
          !
          ! Description of symmetries
          !
          WRITE(stdout, '(36x, "s", 24x, "frac. trans.")')
          CALL s_axis_to_cart() ! give sr(:,:, isym)
          DO isym = 1, nsym
            WRITE(stdout, '(/6x, "isym = ", i2, 5x, a45/)') isym, sname(isym)
            IF (ft(1, isym)**two + ft(2, isym)**two + ft(3, isym)**two > eps8) THEN
              ft1 = at(1, 1) * ft(1, isym) + at(1, 2) * ft(2, isym) + at(1, 3) * ft(3, isym)
              ft2 = at(2, 1) * ft(1, isym) + at(2, 2) * ft(2, isym) + at(2, 3) * ft(3, isym)
              ft3 = at(3, 1) * ft(1, isym) + at(3, 2) * ft(2, isym) + at(3, 3) * ft(3, isym)
              WRITE(stdout, '(1x, "cryst.", 3x, "s(", i2, ") = (", 3(i6, 5x), &
                            " )    f =( ", f10.7, " )")') &
                            isym, (s(1, ipol, isym), ipol = 1, 3), ft(1, isym)
              WRITE(stdout, '(17x, " (", 3(i6, 5x), " )       ( ", f10.7, " )")') &
                                  (s(2, ipol, isym), ipol = 1, 3), ft(2, isym)
              WRITE(stdout, '(17x, " (", 3(i6, 5x), " )       ( ", f10.7, " )"/)') &
                                  (s(3, ipol, isym), ipol = 1, 3), ft(3, isym)
              WRITE(stdout, '(1x, "cart. ", 3x, "s(", i2, ") = (", 3f11.7, &
                            " )    f =( ", f10.7, " )")') &
                            isym, (sr(1, ipol, isym), ipol = 1, 3), ft1
              WRITE(stdout, '(17x, " (", 3f11.7, " )       ( ", f10.7, " )")') &
                                  (sr(2, ipol, isym), ipol = 1, 3), ft2
              WRITE(stdout, '(17x, " (", 3f11.7, " )       ( ", f10.7, " )"/)') &
                                  (sr(3, ipol, isym), ipol = 1, 3), ft3
            ELSE
              WRITE(stdout, '(1x, "cryst.", 3x, "s(", i2, ") = (", 3(i6, 5x), " )")') &
                                                      isym,  (s (1, ipol, isym) , ipol = 1,3)
              WRITE(stdout, '(17x, " (", 3(i6, 5x), " )")')  (s(2, ipol, isym), ipol = 1, 3)
              WRITE(stdout, '(17x, " (", 3(i6, 5x), " )"/)') (s(3, ipol, isym), ipol = 1, 3)
              WRITE(stdout, '(1x, "cart. ", 3x, "s(", i2,") = (", 3f11.7, " )")') &
                                                    isym, (sr(1, ipol, isym), ipol = 1, 3)
              WRITE(stdout, '(17x, " (", 3f11.7, " )")')  (sr(2, ipol, isym), ipol = 1, 3)
              WRITE(stdout, '(17x, " (", 3f11.7, " )"/)') (sr(3, ipol, isym), ipol = 1, 3)
            ENDIF
            !
          ENDDO
          !
        ENDIF
        !
        ! isq(isym)=iq means: when we apply symmetry isym to the originating q
        ! of the star, we get the iq-th member of the star. There are as many
        ! matches as the degeneracy of the star.
        !
        ! We now need to pick up the q in the small group of q* so that Sxq0+G=iq with G=0.
        ! If we choose another element in the small group
        ! the actual q-point may be Sq+G and we screw up the q-vector below to generate
        ! k+q from k and for the KB projectors
        !
        nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
        !
        DO jsym = 1, nsym
          IF (isq(jsym) == iq) THEN
            nsq = nsq + 1
            sym_sgq(nsq) = jsym
          ENDIF
        ENDDO
        IF (nsq * nq /= nsym ) CALL errore('ep_coarse_unfolding', 'wrong degeneracy', iq)
        !
        IF (iverbosity == 1) THEN
          !
          WRITE(stdout,*) 'iq, i, isym, nog, non_symmorphic'
          DO i = 1, nsq
            !
            isym = sym_sgq(i)
            ism1 = invs (isym)
            !
            !  check for G such that Sq = q* + G
            !
            aq  = xq0
            saq = xq
            CALL cryst_to_cart(1, aq, at, -1)
            DO j = 1, 3
              raq(j) = s(j, 1, ism1) * aq(1) &
                     + s(j, 2, ism1) * aq(2) &
                     + s(j, 3, ism1) * aq(3)
            ENDDO
            CALL cryst_to_cart(1, saq, at, -1)
            nog = eqvect_strict(raq, saq, eps6)
            !
            !  check whether the symmetry belongs to a symmorphic group
            !
            !symmo = (ft(1, isym)**2 + ft(2, isym)**2 + ft(3, isym)**2 > 1.0d-8)
            non_symmorphic = (ft(1, isym) /= 0.0d0 .OR. ft(2, isym) /= 0.0d0 .OR. ft(3, isym) /= 0.0d0)
            !
            WRITE(stdout, '(3i5, L3, L3)') iq, i, isym, nog, non_symmorphic
            !
          ENDDO
          !
        ENDIF
        !
        ! SP: We now need to select one symmetry among the small group of q that has G=0
        !     (i.e. that respects Sq0+G=q ). There should always be such symmetry.
        !     We enforce this for later easiness.
        !
        aq = xq0
        saq = xq
        CALL cryst_to_cart(1, aq,  at, -1)
        CALL cryst_to_cart(1, saq, at, -1)
        DO jsym = 1, nsq
          ism1 = invs(sym_sgq(jsym))
          raq = zero
          DO ipol = 1, 3
            DO jpol = 1, 3
              raq(ipol) = raq(ipol) + s(ipol, jpol, ism1) * aq(jpol)
            ENDDO
          ENDDO
          nog = eqvect_strict(raq, saq, eps6)
          IF (nog) THEN ! This is the symmetry such that Sq=q
            isym = sym_sgq(jsym)
            isym1 = indsym(isym)
            EXIT
          ENDIF
          ! If we enter into that loop it means that we have not found
          ! such symmetry within the small group of q.
          IF (jsym == nsq) THEN
            CALL errore('ep_coarse_unfolding ', 'No sym. such that Sxq0=iq was found in the sgq !', 1)
          ENDIF
        ENDDO ! jsym
        !
        ! Write symmetry data to .qmap file.
        !
        IF (compute_dmat) CALL write_qmap(nqc, iq_irr, iq_first, isym, isym1, s, ft, timerev)
        !
        CALL loadumat(nbndep, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq, exband, w_centers)
        !
        !   calculate the sandwiches
        !
        ! A more accurate way of doing this is to symmetrize the matrix element w.r.t.
        ! the small group of the given q in the star. I'm not doing this here.
        ! (but I checked that even without symm the result of full zone and irr zone
        ! are equal to 5+ digits).
        ! For any volunteers, please write to the EPW team.
        !
        timerev = .FALSE.
        imode0 = 0
        ALLOCATE(el_ph_mat(nbndep, nbndep, nks, 3 * nat), STAT = ierr)
        IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating el_ph_mat', 1)
        el_ph_mat(:, :, :, :) = czero
        DO irr = 1, nirr
          npe = npert(irr)
          ALLOCATE(dvscfin(dfftp%nnr, nspin_mag, npe), STAT = ierr)
          IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating dvscfin', 1)
          !
          ! Read the dvscf from file and add the Ewald part.
          CALL dvscf_read(imode0, iq_irr, nqc_irr, npe, xq0, dvscfin, timerev)
          ! Reconstruct the EPC on the full BZ from the IBZ read from file.
          CALL ep_coarse_compute(irr, npe, imode0, dvscfin, gmapsym(:, isym1), eigv(:, isym1), &
                                 isym, xq0, nqc, el_ph_mat, epmatq, timerev)
          !
          DEALLOCATE(dvscfin, STAT = ierr)
          IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating dvscfin', 1)
          imode0 = imode0 + npe
        ENDDO ! irr
        DEALLOCATE(el_ph_mat, STAT = ierr)
        IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating el_ph_mat', 1)
        !
        !  bring epmatq in the mode representation of iq_first,
        !  and then in the cartesian representation of iq
        !
        CALL rotate_eigenm(iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2, timerev)
        !
        ! ZD: Store cz2 for ex-ph calculation when nedded
        cz2_out(:,:) = cz2(:,:)
        !
        CALL rotate_epmat(cz1, cz2, xq, nqc, lwin, lwinq, exband, timerev)
        !DBSP
        !write(*,*)'epmatq(:,:,2,:,nqc)',SUM(epmatq(:,:,2,:,nqc))
        !write(*,*)'epmatq(:,:,2,:,nqc)**2',SUM((REAL(REAL(epmatq(:,:,2,:,nqc))))**2)+&
        !  SUM((REAL(AIMAG(epmatq(:,:,2,:,nqc))))**2)
        !print*,'dynq ', SUM(dynq(:,:,nqc))
        !print*,'et ', et_loc(:,2)
        !END
        ! ZD: Calculate the ex-ph matrix
        IF (exciton) THEN
          CALL calc_G_epmat_subbnd(iq_order, cz2_out, nqc)
        ENDIF
        ! 
        ! SP: Now we treat separately the case imq == 0
        IF (imq == 0) THEN
          timerev = .TRUE.
          !
          ! SP: First the vlocq need to be initialized propertly with the first
          !     q in the star
          xq = -xq0
          CALL init(.FALSE.)
          !
          ! retrieve the q in the star
          xq = -sxq(:, iq)
          !
          ! and populate the uniform grid
          nqc = nqc + 1
          xqc(:,nqc) = xq
          !
          ! ZD: Determine the nscf-iq index of the current (iq_irr, iq) pair
          IF(exciton) THEN
            xq0_cryst = xq
            CALL cryst_to_cart(1, xq0_cryst, at, -1)
            CALL find_iq(xq0_cryst, nqc1, nqc2, nqc3, iq_order)
            ! WRITE(stdout, '(5x, a, I5, 3F10.5)') 'Current q in cryst. coord.: ', iq_order, xq0_cryst
            wf(:,iq_order)=wf_temp(:,iq)
          ENDIF
          !
          IF (iq == 1) WRITE(stdout, *)
          WRITE(stdout,5) nqc, xq
          !
          ! Write symmetry data to .qmap file.
          !
          IF (compute_dmat) CALL write_qmap(nqc, iq_irr, iq_first, isym, isym1, s, ft, timerev)
          !
          !  prepare the kmap for the refolding
          !
          CALL createkmap(xq)
          !
          CALL loadumat(nbndep, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq, exband, w_centers)
          !
          xq0 = -xq0
          !
          imode0 = 0
          ALLOCATE(el_ph_mat(nbndep, nbndep, nks, 3 * nat), STAT = ierr)
          IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating el_ph_mat', 1)
          el_ph_mat(:, :, :, :) = czero
          DO irr = 1, nirr
            npe = npert(irr)
            ALLOCATE(dvscfin(dfftp%nnr, nspin_mag, npe), STAT = ierr)
            IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating dvscfin', 1)
            !
            ! Read the dvscf from file and add the Ewald part.
            CALL dvscf_read(imode0, iq_irr, nqc_irr, npe, xq0, dvscfin, timerev)
            ! Reconstruct the EPC on the full BZ from the IBZ read from file.
            CALL ep_coarse_compute(irr, npe, imode0, dvscfin, gmapsym(:, isym1), eigv(:, isym1), &
                                   isym, xq0, nqc, el_ph_mat, epmatq, timerev)
            !
            DEALLOCATE(dvscfin, STAT = ierr)
            IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating dvscfin', 1)
            imode0 = imode0 + npe
          ENDDO ! irr
          DEALLOCATE(el_ph_mat, STAT = ierr)
          IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating el_ph_mat', 1)
          !
          !  bring epmatq in the mode representation of iq_first,
          !  and then in the cartesian representation of iq
          !
          CALL rotate_eigenm(iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2, timerev)
          !
          ! ZD: Store cz2 for ex-ph calculation when nedded
          cz2_out(:,:) = cz2(:,:)
          !
          CALL rotate_epmat(cz1, cz2, xq, nqc, lwin, lwinq, exband, timerev)
          !
          ! ZD: Calculate the ex-ph matrix
          IF (exciton) THEN
            CALL calc_G_epmat_subbnd(iq_order, cz2_out, nqc)
          ENDIF
          ! 
          xq0 = -xq0
        ENDIF ! end imq == 0
        !
      ENDDO
      ! ZD: Deallocate the temporary array for phonon frequencies
      IF(exciton) THEN 
        DEALLOCATE(wf_temp,STAT=ierr)
        IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error deallocating wf_temp', 1)   
      ENDIF  
      !
      iq_first = iq_first + nq
      if (imq == 0) iq_first = iq_first + nq
      !
    ENDDO ! irr-q loop
    !
    ! ZD: Deallocate all the arrays related to ex-ph coupling
    IF (exciton) THEN 
      CALL release_ex_ph_g()
    ENDIF
    !
    IF (nqc /= nqc1 * nqc2 * nqc3) CALL errore('ep_coarse_unfolding', 'nqc /= nq1*nq2*nq3', nqc)
    !
    ! Write xqc to file. This is needed for wfpt
    !
    IF (lwfpt .AND. ionode) THEN
      CALL diropn(iuxqc, 'xqc', 3 * nqc, exst)
      CALL davcio(xqc, 3 * nqc, iuxqc, 1, +1)
      CLOSE(iuxqc)
    ENDIF
    !
    IF (lifc) THEN
      DEALLOCATE(wscache, STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating wscache', 1)
    ENDIF
    DEALLOCATE(evc, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating evc', 1)
    DEALLOCATE(evq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating evq', 1)
    DEALLOCATE(xkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating xkq', 1)
    DEALLOCATE(shift, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating shift', 1)
    DEALLOCATE(gmap, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating gmap', 1)
    DEALLOCATE(gmapsym, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating gmapsym', 1)
    DEALLOCATE(eigv, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating eigv', 1)
    DEALLOCATE(vlocq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating vlocq', 1)
    DEALLOCATE(dmuxc, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating dmuxc', 1)
    DEALLOCATE(eigqts, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating eigqts', 1)
    DEALLOCATE(rtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating rtau', 1)
    DEALLOCATE(u, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating u', 1)
    DEALLOCATE(npert, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating npert', 1)
    IF (okvan) THEN
      DEALLOCATE(veff, STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating veff', 1)
      DEALLOCATE(int1, STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating int1', 1)
      DEALLOCATE(int2, STAT = ierr)
      IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating int2', 1)
      IF (noncolin) THEN
        DEALLOCATE(int1_nc, STAT = ierr)
        IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating int1_nc', 1)
        IF (lspinorb) THEN
          DEALLOCATE(int2_so, STAT = ierr)
          IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating int2_so', 1)
        ENDIF
      ENDIF
    ENDIF
    DO ik = 1, nks
      DO ipol = 1, 3
        CALL deallocate_bec_type(alphap(ipol, ik))
      ENDDO
    ENDDO
    DEALLOCATE(alphap, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating alphap', 1)
    DO ik = 1, SIZE(becp1)
      CALL deallocate_bec_type(becp1(ik))
    ENDDO
    DEALLOCATE(becp1, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating becp1', 1)
    IF (compute_dmat) THEN
      CALL rotate_wfn_deallocate()
      IF (meta_ionode) CLOSE(iuqmap)
    ENDIF
  ENDIF ! IF (.NOT. epbread .AND. .NOT. epwread) THEN
  !
  IF (my_image_id == 0) THEN
    IF (epbread .OR. epbwrite) THEN
      !
      ! read/write the e-ph matrix elements and other info in the Bloch representation
      ! (coarse mesh) from/to .epb files (one for each pool)
      !
      tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.epb'
      CALL set_ndnmbr(0, my_pool_id + 1, 1, npool, filelab)
      tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.epb' // filelab
      !
      IF (epbread) THEN
        INQUIRE(FILE = tempfile, EXIST = exst)
        IF (.NOT.  exst) CALL errore('ep_coarse_unfolding', 'epb files not found ', 1)
        OPEN(iuepb, FILE = tempfile, FORM = 'unformatted')
        WRITE(stdout, '(/5x, "Reading epmatq from .epb files"/)')
        READ(iuepb) nqc, xqc, et_loc, dynq, epmatq, zstar, epsi
        CLOSE(iuepb)
        WRITE(stdout, '(/5x, "The .epb files have been correctly read"/)')
      ENDIF
      !
      IF (epbwrite) THEN
        OPEN(iuepb, FILE = tempfile, FORM = 'unformatted')
        WRITE(stdout, '(/5x, "Writing epmatq on .epb files"/)')
        WRITE(iuepb) nqc, xqc, et_loc, dynq, epmatq, zstar, epsi
        CLOSE(iuepb)
        WRITE(stdout, '(/5x, "The .epb files have been correctly written"/)')
      ENDIF
    ENDIF
  ENDIF
  !
  ! In case of image parallelization we want to stop after writing the .epb file
  ! S. Tiwari: Since the image parallelization is enabled, we allow the check to
  ! fail
  !  IF (nimage > 1) THEN
  !    WRITE(stdout, '(/5x, "Image parallelization. The code will stop now. "/)')
  !    WRITE(stdout, '(/5x, "You need to restart a calculation by reading the .epb "/)')
  !    WRITE(stdout, '(/5x, "                       with pool parallelization only. "/)')
  !    CALL stop_epw()
  !  ENDIF
  !
  IF (.NOT. epbread .AND. epwread) THEN
    ! CV: need dummy nqc, xqc for the ephwann_shuffle call
    nqc = 1
    xqc = zero
    WRITE(stdout, '(/5x, "Do not need to read .epb files; read .fmt files"/)')
    !
  ENDIF
  !
  ! now dynq is the cartesian dyn mat ( not divided by the masses)
  ! and epmatq is the epmat in cartesian representation (rotation in elphon_shuffle)
  !
  ! free up some memory
  !
  NULLIFY(igkq)
  IF (.NOT. (epwread .AND. .NOT. epbread)) THEN
    DEALLOCATE(xqc_irr, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating xqc_irr', 1)
  ENDIF
  !
  IF ( .NOT. epbread .AND. .NOT. epwread ) THEN
     ! FIXME: the original instruction was "deallocate qrad if re-allocated above"
     ! FIXME: I don't see any reason not to deallocate here if qrad no longer needed
     CALL deallocate_tab_qrad ( )
  ENDIF
  !
  IF (.NOT. (epwread .AND. .NOT. epbread)) THEN
    DEALLOCATE(cu, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating cu', 1)
    DEALLOCATE(cuq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating cuq', 1)
    DEALLOCATE(lwin, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating lwin', 1)
    DEALLOCATE(lwinq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating lwinq', 1)
    DEALLOCATE(exband, STAT = ierr)
    IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error deallocating exband', 1)
  ENDIF
  !
  DEALLOCATE(ibndkept, STAT = ierr)
  IF (ierr /= 0) CALL errore('ep_coarse_unfolding', 'Error allocating ibndkept', 1)
  !
  CALL stop_clock('elphon_wrap')
  !DBSP
  !  DO iq = 1, nqc
  !    write(*,*) iq, xqc(:,iq)
  !    write(*,*)'epmatq(:,:,2,:,iq)',SUM(epmatq(:,:,2,:,iq))
  !    write(*,*)'epmatq(:,:,2,:,iq)**2',SUM((REAL(REAL(epmatq(:,:,2,:,iq))))**2)+&
  !               SUM((REAL(AIMAG(epmatq(:,:,2,:,iq))))**2)
  !  ENDDO
  !END
  !
5 FORMAT (8x, "q(", i5, " ) = (", 3f12.7, " )")
  !
  RETURN
  !---------------------------------------------------------------------------
  END SUBROUTINE ep_coarse_unfolding
  !---------------------------------------------------------------------------
