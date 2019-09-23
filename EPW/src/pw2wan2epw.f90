  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  MODULE pw2wan2epw
  !----------------------------------------------------------------------
  !! 
  !! This module contains routine to go from PW results to Wannier90 to EPW
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw2wan90epw 
    !------------------------------------------------------------------------
    !! 
    !! This is the interface to the Wannier90 code: see http://www.wannier.org
    !! Adapted from the code PP/pw2wannier - Quantum-ESPRESSO group               
    !!
    !! 10/2008  Parellel computation of Amn and Mmn 
    !! 12/2008  Added phase setting of overlap matrix elements
    !! 02/2009  works with standard nkc1*nkc2*nkc3 grids
    !! 12/2009  works with USPP 
    !! 12/2014  RM: Imported the noncolinear case implemented by xlzhang
    !! 06/2016  SP: Debug of SOC + print/reading of nnkp file
    !!
    !------------------------------------------------------------------------
    USE io_global,        ONLY : stdout
    USE klist,            ONLY : nkstot
    USE io_files,         ONLY : prefix
    USE epwcom,           ONLY : write_wfn, scdm_proj, scdm_entanglement, & 
                                 scdm_mu, scdm_sigma
    USE wannierEPW,       ONLY : seedname2, wvfn_formatted, reduce_unk, ispinw, &
                                 ikstart, ikstop, iknum
    USE constants_epw,    ONLY : zero, one
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 4) :: spin_component
    !! Determine the spin case
    CHARACTER(LEN = 256) :: outdir
    !! Name of the output directory
    !
    outdir = './'
    seedname2 = prefix
    spin_component = 'none'
    wvfn_formatted = .FALSE.
    reduce_unk = .FALSE.
    !
    IF (scdm_proj) THEN
      IF ((TRIM(scdm_entanglement) /= 'isolated') .AND. &
          (TRIM(scdm_entanglement) /= 'erfc')     .AND. &
          (TRIM(scdm_entanglement) /= 'gaussian')) THEN
          CALL errore('pw2wan90epw', &
               'Can not recognize the choice for scdm_entanglement. ' &
                //'Valid options are: isolated, erfc and gaussian')
      ENDIF
    ENDIF
    IF (scdm_sigma <= zero) &
      CALL errore('pw2wan90epw', 'Sigma in the SCDM method must be positive.')
    !
    WRITE(stdout, *)
    SELECT CASE(TRIM(spin_component))
    CASE ('up')
      WRITE(stdout, *) '    Spin CASE ( up )'
      ispinw  = 1
      ikstart = 1
      ikstop  = nkstot / 2
      iknum   = nkstot / 2
    CASE ('down')
      WRITE(stdout, *) '    Spin CASE ( down )'
      ispinw  = 2
      ikstart = nkstot / 2 + 1
      ikstop  = nkstot
      iknum   = nkstot / 2
    CASE DEFAULT
      IF (noncolin) THEN
        WRITE(stdout, *) '    Spin CASE ( non-collinear )'
      ELSE
        WRITE(stdout, *) '    Spin CASE ( default = unpolarized )'
      ENDIF
      ispinw  = 0
      ikstart = 1
      ikstop  = nkstot
      iknum   = nkstot
    END SELECT
    !
    WRITE(stdout, *) 
    WRITE(stdout, *) '    Initializing Wannier90'
    WRITE(stdout, *) 
    !
    CALL setup_nnkp()
    CALL ylm_expansion()
    IF (scdm_proj) THEN
      CALL compute_amn_with_scdm()
    ELSE
      CALL compute_amn_para()
    ENDIF
    CALL compute_mmn_para()
    CALL phases_a_m()
    CALL write_band()
    !
    IF (write_wfn) CALL write_plot()
    !
    WRITE(stdout, *)
    WRITE(stdout, *) '    Running Wannier90'
    !
    CALL run_wannier()
    !
    CALL lib_dealloc()
    !
    !-------------------------------------------------------------------------
    END SUBROUTINE pw2wan90epw
    !-------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------
    SUBROUTINE lib_dealloc()
    !-----------------------------------------------------------------------
    !! 
    !! Routine to de-allocate Wannier related matrices. 
    !!
    USE wannierEPW, ONLY : atcart, atsym, kpb, g_kpb, center_w, alpha_w, &
                           l_w, mr_w, r_w, zaxis, xaxis, excluded_band,  &
                           m_mat, u_mat, u_mat_opt, a_mat, eigval,       &
                           lwindow, gf, ig_, zerophase, wann_centers,    &
                           wann_spreads
    !
    IMPLICIT NONE
    ! 
    ! Local variables
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(atcart, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating atcart', 1)
    DEALLOCATE(atsym, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating atsym', 1)
    DEALLOCATE(kpb, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating kpb', 1)
    DEALLOCATE(g_kpb, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating g_kpb', 1)
    DEALLOCATE(center_w, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating center_w', 1)
    DEALLOCATE(alpha_w, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating alpha_w', 1)
    DEALLOCATE(l_w, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating l_w', 1)
    DEALLOCATE(mr_w, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating mr_w', 1)
    DEALLOCATE(r_w, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating r_w', 1)
    DEALLOCATE(zaxis, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating zaxis', 1)
    DEALLOCATE(xaxis, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating xaxis', 1)
    DEALLOCATE(excluded_band, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating excluded_band', 1)
    DEALLOCATE(m_mat, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating m_mat', 1)
    DEALLOCATE(u_mat, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating u_mat', 1)
    DEALLOCATE(u_mat_opt, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating u_mat_opt', 1)
    DEALLOCATE(a_mat, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating a_mat', 1)
    DEALLOCATE(eigval, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating eigval', 1)
    DEALLOCATE(lwindow, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating lwindow', 1)
    DEALLOCATE(gf, STAT = ierr) 
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating gf', 1)
    DEALLOCATE(ig_, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating ig_', 1)
    DEALLOCATE(zerophase, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating zerophase', 1)
    DEALLOCATE(wann_centers, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating wann_centers', 1)
    DEALLOCATE(wann_spreads, STAT = ierr)
    IF (ierr /= 0) CALL errore('lib_dealloc', 'Error deallocating wann_spreads', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lib_dealloc
    !-----------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------
    SUBROUTINE setup_nnkp()
    !-----------------------------------------------------------------------
    !! 
    !! This routine writes and reads the .nnkp file. 
    !! The file specifies
    !! 1) The initial projections functions in the format
    !! num_proj
    !! proj_site(1,i), proj_site(2,i), proj_site(3,i), proj_l(i), proj_m(i), proj_radial(i)
    !! proj_z(1,i), proj_z(2,i), proj_z(3,i),proj_x(1,i), proj_x(2,i), proj_x(3,i), proj_zona(i)
    !! proj_s(i), proj_s_qaxis(1,i), proj_s_qaxis(2,i), proj_s_qaxis(3,i) 
    !!
    !! 2) begin nnkpts: the nearest neighbours of each k-point,         
    !! and therefore provides the information required to      
    !! calculate the M_mn(k,b) matrix elements 
    !! 
    ! ---------------------------------------------------------------------- 
    USE kinds,     ONLY : DP
    USE io_global, ONLY : meta_ionode, stdout, meta_ionode_id
    USE mp_world,  ONLY : world_comm
    USE cell_base, ONLY : at, bg, alat
    USE gvect,     ONLY : g, gg
    USE ions_base, ONLY : nat, tau, ityp, atm
    USE mp,        ONLY : mp_bcast, mp_sum
    USE wvfct,     ONLY : nbnd, npwx
    USE wannierEPW, ONLY : num_nnmax, mp_grid, atcart, atsym, kpb, g_kpb, &
                           center_w, alpha_w, l_w, mr_w, r_w, zaxis,      &
                           xaxis, excluded_band, rlatt, glatt, gf,        &
                           csph, ig_, iknum, seedname2, kpt_latt, nnb,    &
                           num_bands, n_wannier, nexband, nnbx, n_proj,   &
                           spin_eig, spin_qaxis, zerophase
    USE noncollin_module, ONLY : noncolin
    USE constants_epw,    ONLY : bohr, eps6, twopi
    USE mp_pools,         ONLY : intra_pool_comm
    USE epwcom,           ONLY : nbndskip, scdm_proj
    USE w90_io,           ONLY : post_proc_flag
    USE io_var,           ONLY : iunnkp
    ! 
    IMPLICIT NONE
    !
    LOGICAL  :: have_nnkp
    !! Check if the .nnkp file exists.
    LOGICAL  :: found
    !! Check if the section in the .nnkp file is found.
    ! 
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ib
    !! Counter on b-vectors
    INTEGER :: ig
    !! Index on G_k+b vectors
    INTEGER :: iw
    !! Counter on number of projections 
    INTEGER :: ia
    !! Counter on number of atoms
    INTEGER  :: type
    !! Type of atom
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER  :: indexb
    !! Index of exclude_bands 
    INTEGER  :: ipol
    !! Counter on polarizations
    INTEGER  :: idum
    !! Dummy index for reading nnkp file
    INTEGER :: tmp_auto
    !! Needed for the selection of projections with SCDM
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: ig_check(:, :)
    !! Temporary index on G_k+b vectors
    INTEGER  :: exclude_bands(nbnd)
    !! Bands excluded from the calcultion of WFs
    REAL(KIND = DP) :: xnorm
    !! Norm of xaxis
    REAL(KIND = DP) :: znorm
    !! Norm of zaxis
    REAL(KIND = DP) :: coseno
    !! Cosine between xaxis and zaxis
    REAL(KIND = DP) :: gg_
    !! Square of g_(3)
    REAL(KIND = DP) :: g_(3)
    !! Temporary vector G_k+b, g_(:) = g_kpb(:,ik,ib)
    !
    INTERFACE
      !!
      !! SP: An interface is required because the Wannier routine has optional arguments
      !! 
      !-----------------------------------------------------------------------
      SUBROUTINE wannier_setup(seed__name, mp_grid_loc, num_kpts_loc, &
        real_lattice_loc, recip_lattice_loc, kpt_latt_loc, num_bands_tot, &
        num_atoms_loc, atom_symbols_loc, atoms_cart_loc, gamma_only_loc, spinors_loc, &
        nntot_loc, nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc, &
        proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc, proj_z_loc, &
        proj_x_loc, proj_zona_loc, exclude_bands_loc, proj_s_loc, proj_s_qaxis_loc)
      !-----------------------------------------------------------------------
      !
      USE kinds,       ONLY : DP
      USE wannierEPW,  ONLY : num_nnmax
      ! 
      IMPLICIT NONE
      !
      CHARACTER(LEN = *), INTENT(in) :: seed__name
      !! Name
      CHARACTER(LEN = *), INTENT(in) :: atom_symbols_loc(num_atoms_loc)
      !! Symbol for atoms
      LOGICAL, INTENT(in) :: gamma_only_loc
      !! Gamma point only
      LOGICAL, INTENT(in) :: spinors_loc
      !! Local spinor
      INTEGER, INTENT(in) :: mp_grid_loc(3)
      !! MP grid
      INTEGER, INTENT(in) :: num_kpts_loc
      !! Local number of k-points
      INTEGER, INTENT(in) :: num_bands_tot
      !! Number of bands
      INTEGER, INTENT(in) :: num_atoms_loc
      !! Number of atoms
      INTEGER, INTENT(out) :: nntot_loc
      !! Local nntot
      INTEGER, INTENT(out) :: nnlist_loc(num_kpts_loc, num_nnmax)
      !! Local nnlist
      INTEGER, INTENT(out) :: nncell_loc(3, num_kpts_loc, num_nnmax)
      !! Local ncell
      INTEGER, INTENT(out) :: num_bands_loc
      !! Number of bands
      INTEGER, INTENT(out) :: num_wann_loc
      !! Number of Wannier functions
      INTEGER, INTENT(out) :: proj_l_loc(num_bands_tot)
      !! Projection of l local momentum 
      INTEGER, INTENT(out) :: proj_m_loc(num_bands_tot) 
      !! Local projection
      INTEGER, INTENT(out) :: proj_radial_loc(num_bands_tot)
      !! Radial projection
      INTEGER, INTENT(out) :: exclude_bands_loc(num_bands_tot)
      !! Exclude bands
      INTEGER, OPTIONAL, INTENT(out) :: proj_s_loc(num_bands_tot)
      !! Projection on s
      REAL(KIND = DP), INTENT(in) :: REAL_lattice_loc(3, 3)
      !! Lattice 
      REAL(KIND = DP), INTENT(in) :: recip_lattice_loc(3, 3)
      !! Reciprocal lattice
      REAL(KIND = DP), INTENT(in) :: kpt_latt_loc(3, num_kpts_loc)
      !! kpt grid
      REAL(KIND = DP), INTENT(in) :: atoms_cart_loc(3, num_atoms_loc)
      !! Cartesian location of atoms
      REAL(KIND = DP), INTENT(out) :: proj_site_loc(3, num_bands_tot)
      !! Site projection
      REAL(KIND = DP), INTENT(out) :: proj_z_loc(3, num_bands_tot)
      !! Projection on z
      REAL(KIND = DP), INTENT(out) :: proj_x_loc(3, num_bands_tot)
      !! Projection on x
      REAL(KIND = DP), INTENT(out) :: proj_zona_loc(num_bands_tot)
      !! Projection on z 
      REAL(KIND = DP), OPTIONAL, INTENT(out) :: proj_s_qaxis_loc(3, num_bands_tot)
      !! Projection s q-axis
      !
      !-----------------------------------------------------------------------
      END SUBROUTINE wannier_setup
      !-----------------------------------------------------------------------
    END INTERFACE
    ! 
    num_nnmax = 32
    ! 
    ! aam: translations between PW2Wannier90 and Wannier90
    ! pw2wannier90   <==>   Wannier90
    !    nbnd                num_bands_tot
    !    n_wannier           num_wann
    !    num_bands           num_bands
    !    nat                 num_atoms
    !    iknum               num_kpts
    !    rlatt               TRANSPOSE(real_lattice)
    !    glatt               TRANSPOSE(recip_lattice)
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
    ! 
    ALLOCATE(atcart(3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating atcart', 1)
    ALLOCATE(atsym(nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating atsym', 1)
    ALLOCATE(kpb(iknum, num_nnmax), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating kpb', 1)
    ALLOCATE(g_kpb(3, iknum, num_nnmax), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating b_kpb', 1)
    ALLOCATE(center_w(3, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating center_w', 1)
    ALLOCATE(alpha_w(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating alpha_w', 1)
    ALLOCATE(l_w(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating l_w', 1)
    ALLOCATE(mr_w(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating mr_w', 1)
    ALLOCATE(r_w(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating r_w', 1)
    ALLOCATE(zaxis(3, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating zaxis', 1)
    ALLOCATE(xaxis(3, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating xaxis', 1)
    ALLOCATE(excluded_band(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating excluded_band', 1)
    ! 
    ! real lattice (Cartesians, Angstrom)
    rlatt(:, :) = TRANSPOSE(at(:, :)) * alat * bohr
    ! reciprocal lattice (Cartesians, Angstrom)
    glatt(:, :) = TRANSPOSE(bg(:, :)) * twopi / ( alat * bohr )
    ! atom coordinates in Cartesian coords and Angstrom units
    atcart(:, :) = tau(:, :) * bohr * alat
    ! atom symbols
    DO ia = 1, nat
      type = ityp(ia)
      atsym(ia) = atm(type)
    ENDDO
    ! 
    IF (meta_ionode) THEN
      !postproc_setup = .TRUE.
      post_proc_flag = .TRUE.
      CALL wannier_setup(seedname2, mp_grid, iknum,        &  ! input
             rlatt, glatt, kpt_latt, nbnd,                 &  ! input
             nat, atsym, atcart, .FALSE., noncolin,        &  ! input
             nnb, kpb, g_kpb, num_bands, n_wannier,        &  ! output
             center_w, l_w, mr_w, r_w, zaxis,              &  ! output
             xaxis, alpha_w, exclude_bands)                   ! output
     ! SP: In wannier_setup, the .nnkp file is produced.
    ENDIF
    ! 
    CALL mp_bcast(nnb,           meta_ionode_id, world_comm)
    CALL mp_bcast(kpb,           meta_ionode_id, world_comm)
    CALL mp_bcast(g_kpb,         meta_ionode_id, world_comm)
    CALL mp_bcast(num_bands,     meta_ionode_id, world_comm)
    CALL mp_bcast(n_wannier,     meta_ionode_id, world_comm)
    CALL mp_bcast(center_w,      meta_ionode_id, world_comm)
    CALL mp_bcast(l_w,           meta_ionode_id, world_comm)
    CALL mp_bcast(mr_w,          meta_ionode_id, world_comm)
    CALL mp_bcast(r_w,           meta_ionode_id, world_comm)
    CALL mp_bcast(zaxis,         meta_ionode_id, world_comm)
    CALL mp_bcast(xaxis,         meta_ionode_id, world_comm)
    CALL mp_bcast(alpha_w,       meta_ionode_id, world_comm)
    CALL mp_bcast(exclude_bands, meta_ionode_id, world_comm)
    CALL mp_bcast(noncolin,      meta_ionode_id, world_comm)
    !
    ! SP: Commented because we now write on file the .nnkp file and read from it.
    ! 
    ! n_proj = nr. of projections (=#WF unless spinors then =#WF/2) 
    !IF (noncolin) THEN
    !   n_proj = n_wannier/2
    !ELSE
       n_proj = n_wannier
    !ENDIF
    !
    WRITE(stdout, *)
    WRITE(stdout, *) '    Initial Wannier projections'
    WRITE(stdout, *)
    !
    DO iw = 1, n_proj
      WRITE(stdout, '(5x, "(", 3f10.5, ") :  l = ", i3, " mr = ", i3)') & 
                       center_w(:, iw), l_w(iw), mr_w(iw)
    ENDDO
    !
    excluded_band(1:nbnd) = .FALSE.
    nexband = 0
    band_loop: DO ibnd = 1, nbnd
      indexb = exclude_bands(ibnd)
      IF (indexb > nbnd .OR. indexb < 0) THEN
        CALL errore('setup_nnkp', ' wrong excluded band index ', 1)
      ELSEIF (indexb == 0) THEN
        EXIT band_loop
      ELSE
        nexband = nexband + 1
        excluded_band(indexb) = .TRUE.
      ENDIF
    ENDDO band_loop
    !  
    WRITE(stdout, '(/, "      - Number of bands is (", i3, ")")') num_bands
    WRITE(stdout, '("      - Number of total bands is (", i3, ")")') nbnd
    WRITE(stdout, '("      - Number of excluded bands is (", i3, ")")') nexband
    WRITE(stdout, '("      - Number of wannier functions is (", i3, ")")') n_wannier
    IF ((nexband > 0) .AND. (nbndskip /= nexband)) THEN
      WRITE(stdout, '(/5x, "Warning: check if nbndskip = ", i3 " makes sense since ", i3, &
                    &" bands are excluded from wannier projection")') nbndskip, nexband
    ENDIF
    !
    IF ((nbnd-nexband) /= num_bands ) &
        CALL errore('setup_nnkp', ' something wrong with num_bands', 1)
    ! 
    ! Now we read the .nnkp file 
    !
    IF (meta_ionode) THEN  ! Read nnkp file on ionode only
      INQUIRE(FILE = TRIM(seedname2)//".nnkp", EXIST = have_nnkp)
      IF (.NOT. have_nnkp) THEN
         CALL errore('setup_nnkp', 'Could not find the file '&
            &//TRIM(seedname2)//'.nnkp', 1 )
      ENDIF
      OPEN(UNIT = iunnkp, FILE = TRIM(seedname2)//".nnkp", FORM = 'formatted')
    ENDIF
    !
    IF (meta_ionode) THEN   ! read from ionode only
      IF (noncolin) THEN
        CALL scan_file_to(iunnkp, 'spinor_projections', found)
        IF (.NOT. found) THEN
          CALL errore('setup_nnkp', 'Could not find projections block in '&
             &//TRIM(seedname2)//'.nnkp', 1)
        ENDIF
      ELSE
        CALL scan_file_to(iunnkp, 'projections', found)
        IF (.NOT. found) THEN
          CALL errore('setup_nnkp', 'Could not find projections block in '&
             &//TRIM(seedname2)//'.nnkp', 1)
        ENDIF
      ENDIF
      READ(iunnkp, *) n_proj
    ENDIF
    CALL mp_bcast(n_proj, meta_ionode_id, world_comm)
    ! 
    ALLOCATE(gf(npwx, n_proj), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating gf', 1)
    ALLOCATE(csph(16, n_proj), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating csph', 1)
    IF (noncolin) THEN 
      ALLOCATE(spin_eig(n_proj), STAT = ierr)
      IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating spin_eig', 1)
      ALLOCATE(spin_qaxis(3, n_proj), STAT = ierr)
      IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating spin_qaxis', 1)
    ENDIF
    ! 
    IF (meta_ionode) THEN   ! read from ionode only
      DO iw = 1, n_proj
        READ(iunnkp, *) (center_w(ipol, iw), ipol = 1, 3), l_w(iw), mr_w(iw), r_w(iw)
        READ(iunnkp, *) (zaxis(ipol, iw), ipol = 1, 3), (xaxis(ipol, iw), ipol = 1, 3), alpha_w(iw)
        xnorm = SQRT(SUM(xaxis(:, iw) * xaxis(:, iw)))
        IF (xnorm < eps6) CALL errore('setup_nnkp', '|xaxis| < eps', 1)
        znorm = SQRT(SUM(zaxis(:, iw) * zaxis(:, iw)))
        IF (znorm < eps6) CALL errore('setup_nnkp', '|zaxis| < eps', 1)
        coseno = SUM(xaxis(:, iw) * zaxis(:, iw)) / xnorm / znorm
        IF (ABS(coseno) > eps6) CALL errore('setup_nnkp', 'xaxis and zaxis are not orthogonal!', 1)
        IF (alpha_w(iw) < eps6) CALL errore('setup_nnkp', 'zona value must be positive', 1)
        ! convert wannier center in cartesian coordinates (in unit of alat)
        CALL cryst_to_cart(1, center_w(:, iw), at, 1)
        IF (noncolin) THEN
          READ(iunnkp, *) spin_eig(iw), (spin_qaxis(ipol, iw), ipol = 1, 3)
          xnorm = SQRT(SUM(spin_qaxis(:, iw) * spin_qaxis(:, iw)))
          IF (xnorm < eps6) CALL errore('setup_nnkp', '|xaxis| < eps', 1)
          spin_qaxis(:, iw) = spin_qaxis(:, iw) / xnorm
        ENDIF
      ENDDO
    ENDIF
    ! 
    ! automatic projections
    IF (meta_ionode) THEN
      CALL scan_file_to(iunnkp, 'auto_projections', found)
      IF (found) THEN
        READ(iunnkp, *) n_wannier
        READ(iunnkp, *) tmp_auto
        !
        IF (scdm_proj) THEN
          IF (n_proj > 0) THEN
            WRITE(stdout, '(//, " ****** begin Error message ******",/)')
            WRITE(stdout, '(/, " Found a projection block, an auto_projections block", /)')
            WRITE(stdout, '(/, " and scdm_proj = T in the input file. These three options are inconsistent.", /)')
            WRITE(stdout, '(/, " Please refer to the Wannier90 User guide for correct use of these flags.",/)')
            WRITE(stdout, '(/, " ****** end Error message ******",//)')
            CALL errore('setup_nnkp', 'Inconsistent options for projections.', 1)
          ELSE
            IF (tmp_auto /= 0) CALL errore('setup_nnkp', 'Second entry in auto_projections block is not 0. ' // &
               'See Wannier90 User Guide in the auto_projections section for clarifications.', 1)
          ENDIF
        ELSE
          ! Fire an error whether or not a projections block is found
          CALL errore('setup_nnkp', 'scdm_proj = F but found an auto_projections block in '&
                 &//trim(seedname2)//'.nnkp', 1)
        ENDIF
      ELSE
        IF (scdm_proj) THEN
          ! Fire an error whether or not a projections block is found
          CALL errore('setup_nnkp', 'scdm_proj = T but cannot find an auto_projections block in '&
                 &//trim(seedname2)//'.nnkp', 1)
        ENDIF
      ENDIF
    ENDIF
    !
    IF (.NOT. scdm_proj) WRITE(stdout,*) '     - All guiding functions are given '
    ! 
    ! Broadcast
    CALL mp_bcast(center_w, meta_ionode_id, world_comm)
    CALL mp_bcast(l_w,      meta_ionode_id, world_comm)
    CALL mp_bcast(mr_w,     meta_ionode_id, world_comm)
    CALL mp_bcast(r_w,      meta_ionode_id, world_comm)
    CALL mp_bcast(zaxis,    meta_ionode_id, world_comm)
    CALL mp_bcast(xaxis,    meta_ionode_id, world_comm)
    CALL mp_bcast(alpha_w,  meta_ionode_id, world_comm)
    IF (noncolin) THEN
      CALL mp_bcast(spin_eig,   meta_ionode_id, world_comm)
      CALL mp_bcast(spin_qaxis, meta_ionode_id, world_comm)
    ENDIF
    !
    IF (meta_ionode) THEN   ! read from ionode only
      CALL scan_file_to(iunnkp, 'nnkpts', found)
      IF (.NOT. found) THEN
         CALL errore('setup_nnkp', 'Could not find nnkpts block in '&
            &//TRIM(seedname2)//'.nnkp', 1)
      ENDIF
      READ(iunnkp, *) nnb
    ENDIF
    !
    ! Broadcast
    CALL mp_bcast(nnb, meta_ionode_id, world_comm)
    !
    nnbx = 0
    nnbx = MAX(nnbx, nnb)
    ALLOCATE(ig_(iknum, nnbx), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating ig_', 1)
    ALLOCATE(ig_check(iknum, nnbx), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating ig_check', 1)
    ALLOCATE(zerophase(iknum, nnb), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error allocating zerophase', 1)
    zerophase = .FALSE.
    !
    ! Read data about neighbours
    WRITE(stdout, *)
    WRITE(stdout, *) ' Reading data about k-point neighbours '
    WRITE(stdout, *)
    IF (meta_ionode) THEN
      DO ik = 1, iknum
        DO ib = 1, nnb
          READ(iunnkp, *) idum, kpb(ik, ib), (g_kpb(ipol, ik, ib), ipol = 1, 3)
        ENDDO
      ENDDO
    ENDIF
    !
    ! Broadcast
    CALL mp_bcast(kpb,   meta_ionode_id, world_comm)
    CALL mp_bcast(g_kpb, meta_ionode_id, world_comm)
    ! 
    DO ik  = 1, iknum
      DO ib = 1, nnb
        IF ((g_kpb(1, ik, ib) == 0) .AND.  &
            (g_kpb(2, ik, ib) == 0) .AND.  &
            (g_kpb(3, ik, ib) == 0) ) zerophase(ik, ib) = .TRUE.
        g_(:) = REAL(g_kpb(:, ik, ib))
        CALL cryst_to_cart(1, g_, bg, 1)
        gg_ = g_(1) * g_(1) + g_(2) * g_(2) + g_(3) * g_(3)
        ig_(ik, ib) = 0
        ig = 1
        DO WHILE (gg(ig) <= gg_ + eps6)
          IF ((ABS(g(1, ig) - g_(1)) < eps6) .AND.  &
              (ABS(g(2, ig) - g_(2)) < eps6) .AND.  &
              (ABS(g(3, ig) - g_(3)) < eps6) ) ig_(ik, ib) = ig
          ig = ig + 1
        ENDDO
      ENDDO
    ENDDO
    !
    ig_check(:, :) = ig_(:, :)
    CALL mp_sum(ig_check, intra_pool_comm)
    DO ik = 1, iknum
      DO ib = 1, nnb
        IF (ig_check(ik, ib) == 0) &
          CALL errore('setup_nnkp', 'g_kpb vector is not in the list of Gs', 100 * ik + ib)
      ENDDO
    ENDDO
    DEALLOCATE(ig_check, STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_nnkp', 'Error deallocating ig_check', 1)
    !
    WRITE(stdout, *) '     - All neighbours are found '
    WRITE(stdout, *)
    !
    IF (meta_ionode) THEN 
      CALL scan_file_to(iunnkp, 'exclude_bands', found)
      IF (.NOT. found) THEN
        CALL errore('setup_nnkp', 'Could not find exclude_bands block in '&
           &//TRIM(seedname2)//'.nnkp', 1)
      ENDIF
      READ(iunnkp, *) nexband
      excluded_band(1:nbnd) = .FALSE.
      DO ibnd = 1, nexband
        READ(iunnkp, *) indexb
        IF (indexb < 1 .OR. indexb > nbnd) &
          CALL errore('setup_nnkp', ' wrong excluded band index ', 1)
        excluded_band(indexb) = .TRUE.
      ENDDO
    ENDIF
    num_bands = nbnd - nexband
    !
    ! Broadcast
    CALL mp_bcast(nexband,       meta_ionode_id, world_comm)
    CALL mp_bcast(excluded_band, meta_ionode_id, world_comm)
    CALL mp_bcast(num_bands,     meta_ionode_id, world_comm)
    ! 
    IF (meta_ionode) THEN
      CLOSE(iunnkp)
    ENDIF
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE setup_nnkp
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE scan_file_to(iunr, keyword, found)
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: keyword
    !! Keyword searched for
    LOGICAL, INTENT(out)         :: found
    !! Check if the section in the .nnkp file is found.
    INTEGER, INTENT(in)          :: iunr 
    !! Unit number for file
    !  
    ! Local variables
    CHARACTER(LEN = 80) :: line1, line2
    !!
    INTEGER :: ierr
    !! Error status
    !
    ! Uncommenting the following line the file scan restarts every time
    ! from the beginning thus making the reading independent on the order
    ! of data-blocks
    ! REWIND(iunr)
    !
    10 CONTINUE
    READ(iunr, *, END = 20) line1, line2
    IF (line1 /= 'begin') GOTO 10
    IF (line2 /= keyword) GOTO 10
    found = .TRUE.
    RETURN
    20 found = .FALSE.
    REWIND(iunr)
    !
    !------------------------------------------------------------------------
    END SUBROUTINE scan_file_to
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    SUBROUTINE run_wannier()
    !-----------------------------------------------------------------------
    !
    USE io_global, ONLY : stdout, meta_ionode_id, meta_ionode
    USE ions_base, ONLY : nat
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : world_comm
    USE cell_base, ONLY : alat
    USE io_files,  ONLY : prefix
    USE io_var,    ONLY : iuqpeig
    USE pwcom,     ONLY : nkstot
    USE wannierEPW,ONLY : u_mat, lwindow, wann_centers, wann_spreads, eigval,  &
                          n_wannier, spreads, nnb, rlatt, glatt, kpt_latt,     &
                          iknum, seedname2, num_bands, u_mat_opt, atsym, a_mat,&
                          atcart, m_mat, mp_grid
    USE epwcom,    ONLY : eig_read
    USE wvfct,     ONLY : nbnd
    USE constants_epw, ONLY : zero, czero, bohr
    !
    IMPLICIT NONE
    ! 
    CHARACTER(LEN = 256) :: tempfile
    !! Temporary file
    CHARACTER (LEN = 80) :: line
    !! Temporary character
    INTEGER :: iw 
    !! Counter on wannier functions
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: ios
    !! Integer variable for I/O control
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(u_mat(n_wannier, n_wannier, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('run_wannier', 'Error allocating u_mat', 1)
    ALLOCATE(u_mat_opt(num_bands, n_wannier, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('run_wannier', 'Error allocating u_mat_opt', 1)
    ALLOCATE(lwindow(num_bands, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('run_wannier', 'Error allocating lwindow', 1)
    ALLOCATE(wann_centers(3, n_wannier), STAT = ierr)
    IF (ierr /= 0) CALL errore('run_wannier', 'Error allocating wann_centers', 1)
    ALLOCATE(wann_spreads(n_wannier), STAT = ierr)
    IF (ierr /= 0) CALL errore('run_wannier', 'Error allocating wann_spreads', 1)
    !
    u_mat(:, :, :) = czero
    u_mat_opt(:, :, :) = czero
    wann_centers(:, :) = zero
    wann_spreads(:) = zero
    !
    IF (meta_ionode) THEN
      ! read in external eigenvalues, e.g.  GW
      IF (eig_read) THEN
        WRITE (stdout, '(5x, a, i5, a, i5, a)') "Reading external electronic eigenvalues (", &
             nbnd, ",", nkstot,")"
        tempfile = TRIM(prefix)//'.eig'
        OPEN(iuqpeig, FILE = tempfile, FORM = 'formatted', ACTION = 'read', IOSTAT = ios)
        IF (ios /= 0) CALL errore('run_wannier', 'error opening'//tempfile, 1)
        READ(iuqpeig, '(a)') line
        DO ik = 1, nkstot
          ! We do not save the k-point for the moment ==> should be read and
          ! tested against the current one  
          READ(iuqpeig, '(a)') line
          READ(iuqpeig, *) eigval(:, ik)
        ENDDO
        CLOSE(iuqpeig)
      ENDIF
  
  ! SP : This file is not used for now. Only required to build the UNK file
  !      tempfile = TRIM(prefix)//'.mmn'
  !      OPEN(iummn, FILE = tempfile, IOSTAT = ios, FORM = 'unformatted')
  !      WRITE(iummn) m_mat
  !      CLOSE(iummn)
  
      CALL wannier_run(seedname2, mp_grid, iknum,   &              ! input
                       rlatt, glatt, kpt_latt, num_bands,       &  ! input
                       n_wannier, nnb, nat, atsym,              &  ! input
                       atcart, .FALSE., m_mat, a_mat, eigval,   &  ! input
                       u_mat, u_mat_opt, lwindow, wann_centers, &  ! output
                       wann_spreads, spreads)                      ! output
    ENDIF
    !
    CALL mp_bcast(u_mat,        meta_ionode_id, world_comm)
    CALL mp_bcast(u_mat_opt,    meta_ionode_id, world_comm)
    CALL mp_bcast(lwindow,      meta_ionode_id, world_comm)
    CALL mp_bcast(wann_centers, meta_ionode_id, world_comm)
    CALL mp_bcast(wann_spreads, meta_ionode_id, world_comm)
    CALL mp_bcast(spreads,      meta_ionode_id, world_comm)
    !
    !
    ! output the results of the wannierization
    !
    WRITE(stdout, *)
    WRITE(stdout, *) '    Wannier Function centers (cartesian, alat) and spreads (ang):'
    WRITE(stdout, *)
    DO iw = 1, n_wannier
      WRITE(stdout, '(5x, "(", 3f10.5, ") :  ",f8.5)') & 
           wann_centers(:, iw) / alat / bohr, wann_spreads(iw)
    ENDDO
    WRITE(stdout, *)
    !
    ! store the final minimisation matrix on disk for later use
    !
    CALL write_filukk
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE run_wannier
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_amn_para()
    !-----------------------------------------------------------------------
    !!  adapted from compute_amn in pw2wannier90.f90
    !!  parallelization on k-points has been added
    !!  10/2008 Jesse Noffsinger UC Berkeley
    !!  01/2018 Roxana Margine updated 
    !!
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout 
    USE klist,           ONLY : nks, igk_k
    USE klist_epw,       ONLY : xk_loc
    USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
    USE wavefunctions,   ONLY : evc
    USE gvect,           ONLY : g, ngm
    USE gvecw,           ONLY : gcutw
    USE uspp,            ONLY : nkb, vkb
    USE becmod,          ONLY : becp, calbec, deallocate_bec_type, & 
                                allocate_bec_type
    USE wannierEPW,      ONLY : csph, excluded_band, gf, num_bands, &
                                n_wannier, iknum, n_proj, a_mat, &  
                                spin_qaxis, spin_eig, gf_spinor, sgf_spinor
    USE uspp_param,      ONLY : upf
    USE noncollin_module,ONLY : noncolin, npol
    USE constants_epw,   ONLY : czero, cone, zero, eps6
    USE mp_global,       ONLY : my_pool_id, npool, intra_pool_comm, inter_pool_comm
    USE mp,              ONLY : mp_sum
    USE kfold,           ONLY : ktokpmq
    USE mp_world,        ONLY : mpime
    USE io_epw,          ONLY : readwfc
    ! 
    IMPLICIT NONE
    !
    LOGICAL :: any_uspp
    !! Check if uspp are used
    LOGICAL :: spin_z_pos
    !! Detect if spin quantisation axis is along +z
    LOGICAL :: spin_z_neg
    !! Detect if spin quantisation axis is along -z
    INTEGER :: iw
    !! Counter on number of projections
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: ibnd1
    !! Band index
    INTEGER :: ipool
    !! Index of current pool
    INTEGER ::  nkq
    !! Index of k-point in the pool 
    INTEGER ::  nkq_abs
    !! Absolute index of k-point 
    INTEGER ::  ik_g
    !! Temporary index of k-point, ik_g = nkq_abs
    INTEGER  :: ipol
    !! Index of spin-up/spin-down polarizations
    INTEGER ::  istart
    !! Starting index on plane waves
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: zero_vect(3)
    !! Temporary zero vector
    COMPLEX(KIND = DP) :: amn
    !! Element of A_mn matrix
    COMPLEX(KIND = DP), EXTERNAL :: ZDOTC
    !! <psi_mk|g_n>
    COMPLEX(KIND = DP) :: amn_tmp
    !! Element of A_mn matrix
    COMPLEX(KIND = DP) :: fac(2)
    !! Factor for spin quantization axis
    COMPLEX(KIND = DP), ALLOCATABLE :: sgf(:, :)
    !! Guiding functions
    !
    !nocolin: we have half as many projections g(r) defined as wannier
    !         functions. We project onto (1,0) (ie up spin) and then onto
    !         (0,1) to obtain num_wann projections. jry
    !
    any_uspp = ANY(upf(:)%tvanp)
    !
    ALLOCATE(a_mat(num_bands, n_wannier, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_para', 'Error allocating a_mat', 1)
    ALLOCATE(sgf(npwx, n_proj), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_para', 'Error allocating sgf', 1)
    a_mat(:, :, :) = czero
    sgf(:, :)      = czero
    zero_vect(:)   = zero
    !
    IF (noncolin) THEN
      ALLOCATE( gf_spinor(2 * npwx, n_proj), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_para', 'Error allocating gf_spinor', 1)
      ALLOCATE(sgf_spinor(2 * npwx, n_proj), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_para', 'Error allocating sgf_spinor', 1)
      gf_spinor(:, :)  = czero
      sgf_spinor(:, :) = czero
    ENDIF
    !
    WRITE(stdout, '(5x, a)') 'AMN'
    !
    IF (any_uspp) THEN
      CALL deallocate_bec_type(becp)
      CALL allocate_bec_type(nkb, n_wannier, becp)
    ENDIF
    !
#if defined(__MPI)
    WRITE(stdout, '(6x, a, i5, a, i4, a)') 'k points = ', iknum, ' in ', npool, ' pools'
#endif
    ! 
    DO ik = 1, nks
      !
      ! returns in-pool index nkq and absolute index nkq_abs of xk
      CALL ktokpmq(xk_loc(:, ik), zero_vect, +1, ipool, nkq, nkq_abs)
      ik_g = nkq_abs
      !
      WRITE(stdout, '(5x, i8, " of ", i4, a)') ik , nks, ' on ionode'
      FLUSH(stdout)
      !
      ! read wfc at k
      CALL readwfc(my_pool_id + 1, ik, evc) 
      !
      ! sorts k+G vectors in order of increasing magnitude, up to ecut
      CALL gk_sort(xk_loc(1, ik), ngm, g, gcutw, npw, igk_k(1, ik), g2kin)
      !
      CALL generate_guiding_functions(ik)   ! they are called gf(npw, n_proj)
      !
      IF (noncolin) THEN
        CALL orient_gf_spinor(npw)
      ENDIF
      !
      !  USPP
      !
      IF (any_uspp) THEN
        CALL init_us_2(npw, igk_k(1,ik), xk_loc(1,ik), vkb)
        ! below we compute the product of beta functions with trial func.
        IF (noncolin) THEN
          CALL calbec(npw, vkb, gf_spinor, becp, n_proj)
        ELSE
          CALL calbec(npw, vkb, gf, becp, n_proj)
        ENDIF
        ! and we use it for the product S|trial_func>
        IF (noncolin) THEN
          CALL s_psi(npwx, npw, n_proj, gf_spinor, sgf_spinor)
        ELSE
          CALL s_psi(npwx, npw, n_proj, gf, sgf)
        ENDIF
      ELSE
        sgf(:, :) = gf(:, :)
      ENDIF
      !
      IF (noncolin) THEN
        ! we do the projection as g(r)*a(r) and g(r)*b(r)
        DO iw = 1, n_proj
          !
          spin_z_pos = .FALSE.
          spin_z_neg = .FALSE.
          ! detect if spin quantisation axis is along z
          IF ((ABS(spin_qaxis(1, iw) - 0.0d0) < eps6) .AND. & 
              (ABS(spin_qaxis(2, iw) - 0.0d0) < eps6) .AND. &
              (ABS(spin_qaxis(3, iw) - 1.0d0) < eps6)) THEN
            spin_z_pos = .TRUE.
          ELSEIF ((ABS(spin_qaxis(1, iw) - 0.0d0) < eps6) .AND. & 
                  (ABS(spin_qaxis(2, iw) - 0.0d0) < eps6) .AND. &
                  (ABS(spin_qaxis(3, iw) + 1.0d0) < eps6)) THEN
            spin_z_neg = .TRUE.
          ENDIF
          !
          IF (spin_z_pos .OR. spin_z_neg) THEN
            !
            ibnd1 = 0
            DO ibnd = 1, nbnd
              IF (excluded_band(ibnd)) CYCLE
              ibnd1 = ibnd1 + 1
              !
              IF (spin_z_pos) THEN
                ipol = (3 - spin_eig(iw)) / 2
              ELSE
                ipol = (3 + spin_eig(iw)) / 2
              ENDIF
              istart = (ipol - 1) * npwx + 1
              ! 
              amn = czero
              IF (any_uspp) THEN
                amn = ZDOTC(npw, evc(1, ibnd), 1, sgf_spinor(1, iw), 1)
                amn = amn + ZDOTC(npw, evc(npwx + 1, ibnd), 1, sgf_spinor(npwx + 1, iw), 1)
              ELSE
                amn = ZDOTC(npw, evc(istart, ibnd), 1, sgf(1, iw), 1)
              ENDIF
              CALL mp_sum(amn, intra_pool_comm)
              !
              ! changed since n_proj is now read from .nnkp file (n_proj=n_wannier) 
              ! a_mat(ibnd1, iw + n_proj * (ipol - 1), ik_g) = amn
              a_mat(ibnd1, iw, ik_g) = amn
            ENDDO
          ELSE
            ! general routine for quantisation axis (a,b,c) 
            ! 'up'    eigenvector is 1/SQRT(1+c) [c+1,a+ib]
            ! 'down'  eigenvector is 1/SQRT(1-c) [c-1,a+ib]
            IF (spin_eig(iw) == 1) THEN
              fac(1) = (1.0d0 / SQRT(1 + spin_qaxis(3, iw))) & 
                     * (spin_qaxis(3, iw) + 1) * cone
              fac(2) = (1.0d0 / SQRT(1 + spin_qaxis(3, iw))) & 
                     * CMPLX(spin_qaxis(1, iw), spin_qaxis(2, iw), KIND = DP)
            ELSE
              fac(1) = (1.0d0 / SQRT(1 + spin_qaxis(3, iw))) & 
                     * (spin_qaxis(3, iw)) * cone
              fac(2) = (1.0d0 / SQRT(1 - spin_qaxis(3, iw))) & 
                     * CMPLX(spin_qaxis(1, iw), spin_qaxis(2, iw), KIND = DP)
            ENDIF
            !
            ibnd1 = 0
            DO ibnd = 1, nbnd
              IF (excluded_band(ibnd)) CYCLE
              ibnd1 = ibnd1 + 1
              !
              amn = czero
              DO ipol = 1, npol
                istart = (ipol - 1) * npwx + 1
                amn_tmp = czero
                IF (any_uspp) THEN 
                  amn_tmp = ZDOTC(npw, evc(istart, ibnd), 1, sgf_spinor(istart, iw), 1)
                  CALL mp_sum(amn_tmp, intra_pool_comm)
                  amn = amn + amn_tmp
                ELSE
                  amn_tmp = ZDOTC(npw, evc(istart, ibnd), 1, sgf(1, iw), 1)
                  CALL mp_sum(amn_tmp, intra_pool_comm)
                  amn = amn + fac(ipol) * amn_tmp
                ENDIF
              ENDDO
              !
              ! changed since n_proj is now read from .nnkp file (n_proj=n_wannier)
              !a_mat(ibnd1, iw + n_proj * (ipol - 1), ik_g) = amn
              a_mat(ibnd1, iw, ik_g) = amn
            ENDDO
          ENDIF ! spin_z_pos
        ENDDO
      ELSE ! scalar wavefunction
        DO iw = 1, n_proj
          !
          ibnd1 = 0
          DO ibnd = 1, nbnd
            IF (excluded_band(ibnd)) CYCLE
            ibnd1 = ibnd1 + 1
            !
            amn = czero
            amn = ZDOTC(npw, evc(1, ibnd), 1, sgf(1, iw), 1)
            CALL mp_sum(amn, intra_pool_comm)
            ! 
            a_mat(ibnd1, iw, ik_g) = amn
          ENDDO !bands
        ENDDO !wannier fns
      ENDIF
    ENDDO  ! k-points
    !
    DEALLOCATE(csph, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_para', 'Error deallocating csph', 1)
    DEALLOCATE(sgf, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_para', 'Error deallocating sgf', 1)
    IF (noncolin) THEN  
      DEALLOCATE(gf_spinor, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_para', 'Error deallocating gf_spinor', 1)
      DEALLOCATE(sgf_spinor, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_para', 'Error deallocating sgf_spinor', 1)
    ENDIF
    !
    IF (any_uspp) CALL deallocate_bec_type(becp)
    !
    CALL mp_sum(a_mat, inter_pool_comm)
    !
    !IF (mpime == 0) WRITE(900,*) a_mat
    !IF (mpime == 0) FLUSH(900)
    !IF (mpime == 1) WRITE(901,*) a_mat
    !IF (mpime == 1) FLUSH(901)
  
    WRITE(stdout,*)
    WRITE(stdout,'(5x,a)') 'AMN calculated'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_amn_para
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_amn_with_scdm()
    !-----------------------------------------------------------------------
    !
    !! This subroutine computes automatic Wannier functions using 
    !! the selected columns of density matrix (SCDM) method. 
    !! The subroutine is adapted from compute_amn_with_scdm and 
    !! compute_amn_with_scdm_spinor subroutines in pw2wannier90.f90. 
    !
    USE kinds,           ONLY : DP
    USE constants,       ONLY : rytoev
    USE io_global,       ONLY : stdout, meta_ionode, meta_ionode_id
    USE wvfct,           ONLY : nbnd, npw, npwx, et, g2kin
    USE gvecw,           ONLY : gcutw
    USE wavefunctions,   ONLY : evc, psic, psic_nc
    USE noncollin_module, ONLY : noncolin
    USE wannierEPW,      ONLY : n_wannier, iknum, num_bands, nexband, &
                                excluded_band, kpt_latt, a_mat
    USE klist,           ONLY : nks, igk_k
    USE klist_epw,       ONLY : xk_all, xk_loc
    USE gvect,           ONLY : g, ngm, mill
    USE fft_base,        ONLY : dffts
    USE scatter_mod,     ONLY : gather_grid
    USE fft_interfaces,  ONLY : invfft
    USE mp,              ONLY : mp_bcast, mp_sum
    USE mp_world,        ONLY : world_comm
    USE mp_pools,        ONLY : intra_pool_comm
    USE mp_global,       ONLY : my_pool_id, npool, inter_pool_comm
    USE ions_base,       ONLY : ntyp => nsp
    USE constants_epw,   ONLY : zero, czero, one, twopi
    USE kfold,           ONLY : ktokpmq
    USE epwcom,          ONLY : scdm_entanglement, scdm_mu, scdm_sigma
    USE io_epw,          ONLY : readwfc
    !
    IMPLICIT NONE
    !
    ! Local variables
    LOGICAL :: found_gamma
    !! find G-point index
    INTEGER :: ik
    !! Counter over k-points
    INTEGER :: ibnd
    !! Counter over bands
    INTEGER :: ibnd1
    !! Band index
    INTEGER :: iw
    !! Counter over the nr. of projections
    INTEGER :: nrtot
    !! Size of FFT grid
    INTEGER :: ierr
    !! Error status
    INTEGER :: info
    !! if info=0 succesful exit of QR factorization
    INTEGER :: lcwork
    !!
    INTEGER :: ipt, jpt, kpt, lpt
    !! Counter over FFT grid
    INTEGER :: ispin
    !! Spin index
    INTEGER :: count_piv_spin
    !!
    INTEGER :: minmn
    !!
    INTEGER :: minmn2
    !!
    INTEGER :: maxmn2
    !!
    INTEGER :: numbands
    !! Nr. of bands
    INTEGER :: nbtot
    !! Total nr. of bands
    INTEGER :: ig
    !! Counter on G vectors
    INTEGER :: ig_local
    !! Counter on G vectors
    INTEGER :: ipool
    !! Index of current pool
    INTEGER ::  nkq
    !! Index of k-point in the pool
    INTEGER ::  nkq_abs
    !! Absolute index of k-point
    INTEGER ::  ik_g
    !! Temporary index of k-point, ik_g = nkq_abs
    INTEGER, ALLOCATABLE :: piv(:)
    !! vv: Pivot array in the QR factorization
    INTEGER, ALLOCATABLE :: piv_pos(:)
    !! jml: Position of piv
    INTEGER, ALLOCATABLE :: piv_spin(:)
    !! jml: Spin index of piv
    !
    REAL(KIND = DP):: pnorm
    !! 
    REAL(KIND = DP):: norm_psi
    !! norm of wave function |Psi|
    REAL(KIND = DP):: f_gamma
    !! generalization of Fermi-Dirac probability for occupied states at G-point
    REAL(KIND = DP):: tpi_r_dot_k
    !! scalar product i*2*pi*k*r
    REAL(KIND = DP):: tpi_r_dot_g
    !! scalar prodoct i*2*pi*G*r
    REAL(KIND = DP), EXTERNAL :: DDOT
    !! Scalar product of two vectors
    REAL(KIND = DP) :: zero_vect(3)
    !! Temporary zero vector
    REAL(KIND = DP) :: g_(3)
    !! Temporary vector G
    REAL(KIND = DP), ALLOCATABLE :: focc(:)
    !! generalization of Fermi-Dirac probability for occupied states
    !
    ! vv: Real array for the QR factorization and SVD
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: rwork2(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: singval(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: rpos(:, :)
    !! real space coords. of FFT grid
    REAL(KIND = DP), ALLOCATABLE :: cpos(:, :)
    !! real space coords. of FFT grid
    !
    COMPLEX(KIND = DP) :: nowfc_tmp
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: nowfc(:, :)
    !! sum_G (psi(G) * e^(i*G*r)) * focc  * phase / norm_psi
    COMPLEX(KIND = DP), ALLOCATABLE :: psi_gamma(:, :)
    !! psi * f_gamma / norm_psi
    COMPLEX(KIND = DP) :: tmp_cwork(2)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: phase(:)
    !! exp(i*k*r)
    COMPLEX(KIND = DP), ALLOCATABLE :: phase_g(:, :)
    !! exp(i*G*r)
    !
    ! complex arrays in QR factorization
    COMPLEX(KIND = DP), ALLOCATABLE :: qr_tau(:)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: cwork(:)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: cwork2(:)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: Umat(:, :)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: VTmat(:, :)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: Amat(:, :)
    !! Elements of A_mn matrix
    !
#if defined(__MPI)
    INTEGER nxxs
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: psic_all(:)
    !! 
    COMPLEX(KIND = DP), ALLOCATABLE :: psic_nc_all(:, :)
    !!
    nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
    IF (noncolin) THEN
      ALLOCATE(psic_nc_all(nxxs, 2), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating psic_nc_all', 1)
      psic_nc_all(:, :) = czero
    ELSE
      ALLOCATE(psic_all(nxxs), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating psic_all', 1)
      psic_all(:) = czero
    ENDIF
#endif
    !
    ! vv: Write info about SCDM in output
    IF (TRIM(scdm_entanglement) == 'isolated') THEN
      WRITE(stdout, '(1x, a, a/)') 'Case  : ', TRIM(scdm_entanglement)
    ELSEIF ((TRIM(scdm_entanglement) == 'erfc') .OR. &
            (TRIM(scdm_entanglement) == 'gaussian')) THEN
      WRITE(stdout, '(1x, a, a)') 'Case  : ', TRIM(scdm_entanglement)
      WRITE(stdout, '(1x, a, f10.3, a/, 1x, a, f10.3, a/)') 'mu    = ', scdm_mu, ' eV', 'sigma =', scdm_sigma, ' eV'
    ENDIF
    !
    ! vv: Allocate all the variables for the SCDM method:
    !     1)For the QR decomposition
    !     2)For the unk's on the real grid
    !     3)For the SVD
    !
    IF (TRIM(scdm_entanglement) == 'isolated') THEN
      numbands = n_wannier
      nbtot = n_wannier + nexband
    ELSE
      numbands = nbnd - nexband
      nbtot = nbnd
    ENDIF
    !
    nrtot = dffts%nr1 * dffts%nr2 * dffts%nr3
    IF (noncolin) THEN
      minmn = MIN(numbands, 2 * nrtot)
      ALLOCATE(qr_tau(2 * minmn), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating qr_tau', 1)
      ALLOCATE(piv(2 * nrtot), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating piv', 1)
      ALLOCATE(piv_pos(n_wannier), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating piv_pos', 1)
      ALLOCATE(piv_spin(n_wannier), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating piv_spin', 1)
      ALLOCATE(rwork(2 * 2 * nrtot), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating rwork', 1)
      ALLOCATE(psi_gamma(2 * nrtot, numbands), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating psi_gamma', 1)
    ELSE
      minmn = MIN(numbands, nrtot)
      ALLOCATE(qr_tau(2 * minmn), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating qr_tau', 1)
      ALLOCATE(piv(nrtot), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating piv', 1)
      ALLOCATE(rwork(2 * nrtot), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating rwork', 1)
      ALLOCATE(psi_gamma(nrtot, numbands), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating psi_gamma', 1)
    ENDIF
    !
    ALLOCATE(nowfc(n_wannier, numbands), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating nowfc', 1)
    ALLOCATE(focc(numbands), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating focc', 1)
    minmn2 = MIN(numbands, n_wannier)
    maxmn2 = MAX(numbands, n_wannier)
    ALLOCATE(rwork2(5 * minmn2), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating rwork2', 1)
    !
    ALLOCATE(rpos(nrtot, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating rpos', 1)
    ALLOCATE(cpos(n_wannier, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating cpos', 1)
    ALLOCATE(phase(n_wannier), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating phase', 1)
    ALLOCATE(singval(n_wannier), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating singval', 1)
    ALLOCATE(Umat(numbands, n_wannier), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating Umat', 1)
    ALLOCATE(VTmat(n_wannier, n_wannier), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating VTmat', 1)
    ALLOCATE(Amat(numbands, n_wannier), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating Amat', 1)
    !
    ALLOCATE(a_mat(num_bands, n_wannier, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating a_mat', 1)
    a_mat(:, :, :) = czero
    zero_vect(:) = zero
    !
    WRITE(stdout, '(5x, a)') 'AMN'
    !
    ! Check that Gamma-point is first in the list of k-vectors
    !
    found_gamma = .FALSE.
    found_gamma = kpt_latt(1, 1) == 0.0d0 .AND. &
                  kpt_latt(2, 1) == 0.0d0 .AND. &
                  kpt_latt(3, 1) == 0.0d0
    IF (.NOT. found_gamma) CALL errore('compute_amn_with_scdm', 'No Gamma point found.', 1)
    !
    ! read wfc at G-point
    ik = 1
    CALL readwfc(my_pool_id + 1, ik, evc)
    !
    ! sorts k+G vectors in order of increasing magnitude, up to ecut
    CALL gk_sort(xk_all(1, ik), ngm, g, gcutw, npw, igk_k(1, ik), g2kin)
    !
    ibnd1 = 0
    f_gamma = zero
    psi_gamma(:, :) = czero
    IF (noncolin) THEN
      DO ibnd = 1, nbtot
        IF (excluded_band(ibnd)) CYCLE
        ibnd1 = ibnd1 + 1
        !
        ! check ibnd1 <= numbands
        IF (ibnd1 > numbands) CALL errore('compute_amn_with_scdm', &
                                         &'Something wrong with the number of bands. Check exclude_bands.', 1)
        IF (TRIM(scdm_entanglement) == 'isolated') THEN
          f_gamma = 1.0d0
        ELSEIF (TRIM(scdm_entanglement) == 'erfc') THEN
          f_gamma = 0.5d0 * ERFC((et(ibnd, ik) * rytoev - scdm_mu) / scdm_sigma)
        ELSEIF (TRIM(scdm_entanglement) == 'gaussian') THEN
          f_gamma = EXP(- 1.0d0 * ((et(ibnd, ik) * rytoev - scdm_mu)**2) / (scdm_sigma**2))
        ELSE
          CALL errore('compute_amn_with_scdm', 'scdm_entanglement value not recognized.', 1)
        ENDIF
        !
        ! vv: Compute unk's on a real grid (the fft grid)
        !
        psic_nc(:, :) = czero
        psic_nc(dffts%nl(igk_k(1:npw, ik)), 1) = evc(       1:npw,        ibnd)
        psic_nc(dffts%nl(igk_k(1:npw, ik)), 2) = evc(1 + npwx:npw + npwx, ibnd)
        CALL invfft('Wave', psic_nc(:, 1), dffts)
        CALL invfft('Wave', psic_nc(:, 2), dffts)
        !
        ! vv: Build Psi_k = Unk * focc at G-point only
#if defined(__MPI)
        CALL gather_grid(dffts, psic_nc(:, 1), psic_nc_all(:, 1))
        CALL gather_grid(dffts, psic_nc(:, 2), psic_nc_all(:, 2))
        pnorm = REAL(SUM(psic_nc_all(1:nrtot, 1) * CONJG(psic_nc_all(1:nrtot, 1))), KIND = DP) + &
                REAL(SUM(psic_nc_all(1:nrtot, 2) * CONJG(psic_nc_all(1:nrtot, 2))), KIND = DP)
        norm_psi = SQRT(pnorm)
        psi_gamma(        1:nrtot,     ibnd1) = psic_nc_all(1:nrtot, 1) * f_gamma / norm_psi
        psi_gamma(1 + nrtot:2 * nrtot, ibnd1) = psic_nc_all(1:nrtot, 2) * f_gamma / norm_psi
#else
        pnorm = REAL(SUM(psic_nc(1:nrtot, 1) * CONJG(psic_nc(1:nrtot, 1))), KIND = DP) + &
                REAL(SUM(psic_nc(1:nrtot, 2) * CONJG(psic_nc(1:nrtot, 2))), KIND = DP)
        norm_psi = SQRT(pnorm) 
        psi_gamma(        1:nrtot,     ibnd1) = psic_nc(1:nrtot, 1) * f_gamma / norm_psi
        psi_gamma(1 + nrtot:2 * nrtot, ibnd1) = psic_nc(1:nrtot, 2) * f_gamma / norm_psi
#endif
      ENDDO ! ibnd
    ELSE
      DO ibnd = 1, nbtot
        IF (excluded_band(ibnd)) CYCLE
        ibnd1 = ibnd1 + 1
        !
        ! check ibnd1 <= numbands
        IF (ibnd1 > numbands) CALL errore('compute_amn_with_scdm', &
                                          &'Something wrong with the number of bands. Check exclude_bands.', 1)
        IF (TRIM(scdm_entanglement) == 'isolated') THEN
          f_gamma = 1.0d0
        ELSEIF (TRIM(scdm_entanglement) == 'erfc') THEN
          f_gamma = 0.5d0 * ERFC((et(ibnd, ik) * rytoev - scdm_mu) / scdm_sigma)
        ELSEIF (TRIM(scdm_entanglement) == 'gaussian') THEN
          f_gamma = EXP(- 1.0d0 * ((et(ibnd, ik) * rytoev - scdm_mu)**2) / (scdm_sigma**2))
        ELSE
          CALL errore('compute_amn_with_scdm', 'scdm_entanglement value not recognized.', 1)
        ENDIF
        !
        ! vv: Compute unk's on a real grid (the fft grid)
        !
        psic(:) = czero
        psic(dffts%nl(igk_k(1:npw, ik))) = evc(1:npw, ibnd)
        CALL invfft('Wave', psic, dffts)
        !
        ! vv: Build Psi_k = Unk * focc at G-point only
        psi_gamma(:, :) = czero
#if defined(__MPI)
        CALL gather_grid(dffts, psic, psic_all)
        pnorm = REAL(SUM(psic_all(1:nrtot) * CONJG(psic_all(1:nrtot))), KIND = DP)
        norm_psi = SQRT(pnorm)
        psi_gamma(1:nrtot, ibnd1) = psic_all(1:nrtot) * f_gamma / norm_psi
#else
        pnorm = REAL(SUM(psic(1:nrtot) * CONJG(psic(1:nrtot))), KIND = DP)
        norm_psi = SQRT(pnorm)
        psi_gamma(1:nrtot, ibnd1) = psic(1:nrtot) * f_gamma / norm_psi
#endif
      ENDDO ! ibnd
    ENDIF !noncolin
    !
    ! vv: Perform QR factorization with pivoting on Psi_Gamma
    ! vv: Preliminary call to define optimal values for lwork and cwork size
    !
    piv(:)       = 0
    qr_tau(:)    = czero
    tmp_cwork(:) = czero
    rwork(:)     = zero
    IF (noncolin) THEN
      CALL ZGEQP3(numbands, 2 * nrtot, TRANSPOSE(CONJG(psi_gamma)), numbands, &
                  piv, qr_tau, tmp_cwork, -1, rwork, info)
    ELSE
      CALL ZGEQP3(numbands, nrtot, TRANSPOSE(CONJG(psi_gamma)), numbands, &
                  piv, qr_tau, tmp_cwork, -1, rwork, info)
    ENDIF
    IF (info /= 0) CALL errore('compute_amn_with_scdm', 'Error in computing the QR factorization', 1)
    !
    lcwork = AINT(REAL(tmp_cwork(1)))
    ALLOCATE(cwork(lcwork), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating cwork', 1)
    piv(:)    = 0
    qr_tau(:) = czero
    cwork(:)  = czero
    rwork(:)  = zero
    !
#if defined(__MPI)
    IF (meta_ionode) THEN
      IF (noncolin) THEN 
        CALL ZGEQP3(numbands, 2 * nrtot, TRANSPOSE(CONJG(psi_gamma)), numbands, &
                    piv, qr_tau, cwork, lcwork, rwork, info)
      ELSE
        CALL ZGEQP3(numbands, nrtot, TRANSPOSE(CONJG(psi_gamma)), numbands, &
                    piv, qr_tau, cwork, lcwork, rwork, info)
      ENDIF
      IF (info /= 0) CALL errore('compute_amn_with_scdm', 'Error in computing the QR factorization', 1)
    ENDIF
    CALL mp_bcast(piv, meta_ionode_id, world_comm)
#else
    IF (noncolin) THEN 
      CALL ZGEQP3(numbands, 2 * nrtot, TRANSPOSE(CONJG(psi_gamma)), numbands, &
                  piv, qr_tau, cwork, lcwork, rwork, info)
    ELSE
      CALL ZGEQP3(numbands, nrtot, TRANSPOSE(CONJG(psi_gamma)), numbands, &
                  piv, qr_tau, cwork, lcwork, rwork, info)
    ENDIF
    IF (info /= 0) CALL errore('compute_amn_with_scdm', 'Error in computing the QR factorization', 1)
#endif
    !
    DEALLOCATE(cwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating cwork', 1)
    !
    ! jml: calculate position and spin part of piv
    IF (noncolin) THEN
      count_piv_spin = 0
      piv_pos(:) = 0
      piv_spin(:) = 0
      DO iw = 1, n_wannier
        IF (piv(iw) <= nrtot) THEN
          piv_pos(iw) = piv(iw)
          piv_spin(iw) = 1
          count_piv_spin = count_piv_spin + 1
        ELSE
          piv_pos(iw) = piv(iw) - nrtot
          piv_spin(iw) = 2
        ENDIF
      ENDDO
      WRITE(stdout, '(a, i5)') " Number of pivot points with spin up  : ", count_piv_spin
      WRITE(stdout, '(a, i5)') " Number of pivot points with spin down: ", n_wannier - count_piv_spin
    ENDIF
    !
    ! vv: Compute the points
    lpt = 0
    rpos(:, :) = zero
    DO kpt = 0, dffts%nr3 - 1
      DO jpt = 0, dffts%nr2 - 1
        DO ipt = 0, dffts%nr1 - 1
          lpt = lpt + 1
          rpos(lpt, 1) = REAL(ipt, KIND = DP) / REAL(dffts%nr1, KIND = DP)
          rpos(lpt, 2) = REAL(jpt, KIND = DP) / REAL(dffts%nr2, KIND = DP)
          rpos(lpt, 3) = REAL(kpt, KIND = DP) / REAL(dffts%nr3, KIND = DP)
        ENDDO
      ENDDO
    ENDDO
    !
    cpos(:, :) = zero
    IF (noncolin) THEN 
      DO iw = 1, n_wannier
        cpos(iw, :) = rpos(piv_pos(iw), :)
        cpos(iw, :) = cpos(iw, :) - ANINT(cpos(iw, :))
      ENDDO
    ELSE
      DO iw = 1, n_wannier
        cpos(iw, :) = rpos(piv(iw), :)
        cpos(iw, :) = cpos(iw, :) - ANINT(cpos(iw, :))
      ENDDO
    ENDIF
    !
#if defined(__MPI)
    WRITE(stdout,'(6x,a,i5,a,i4,a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
    !
    DO ik = 1, nks
      !
      ! returns in-pool index nkq and absolute index nkq_abs of xk
      CALL ktokpmq(xk_loc(:,ik), zero_vect, +1, ipool, nkq, nkq_abs)
      ik_g = nkq_abs
      !
      WRITE(stdout,'(5x,i8, " of ", i4,a)') ik , nks, ' on ionode'
      FLUSH(stdout)
      !
      ! read wfc at k
      CALL readwfc(my_pool_id + 1, ik, evc)
      !
      ! sorts k+G vectors in order of increasing magnitude, up to ecut
      CALL gk_sort(xk_loc(1, ik), ngm, g, gcutw, npw, igk_k(1, ik), g2kin)
      !
      ! vv: SCDM method for generating the Amn matrix
      ! jml: calculate of psi_nk at pivot points using slow FT
      !      This is faster than using invfft because the number of pivot
      !      points is much smaller than the number of FFT grid points.
      !
      ! jml: calculate phase factors before the loop over bands
      ALLOCATE(phase_g(npw, n_wannier), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating phase_g', 1)
      phase_g(:, :) = czero
      phase(:)      = czero
      !
      DO iw = 1, n_wannier
        !tpi_r_dot_k = twopi * DDOT(3, cpos(iw, :), 1, kpt_latt(:, ik_g), 1)
        tpi_r_dot_k = twopi * (cpos(iw, 1) * kpt_latt(1, ik_g) + &
                               cpos(iw, 2) * kpt_latt(2, ik_g) + &
                               cpos(iw, 3) * kpt_latt(3, ik_g))
        phase(iw) = CMPLX(COS(tpi_r_dot_k), SIN(tpi_r_dot_k), KIND = DP)
        DO ig_local = 1, npw
          ig = igk_k(ig_local, ik)
          g_(:) = REAL(mill(:, ig), KIND = DP)
          !tpi_r_dot_g = twopi * DDOT(3, cpos(iw, :), 1, g_(:), 1)
          tpi_r_dot_g = twopi * (cpos(iw, 1) * g_(1) + &
                                 cpos(iw, 2) * g_(2) + &
                                 cpos(iw, 3) * g_(3))
          phase_g(ig_local, iw) = CMPLX(COS(tpi_r_dot_g), SIN(tpi_r_dot_g), KIND = DP)
        ENDDO
      ENDDO
      !
      ! vv: Generate the occupation numbers matrix according to scdm_entanglement
      focc(:)     = zero
      nowfc(:, :) = czero
      ibnd1 = 0
      IF (noncolin) THEN
        DO ibnd = 1, nbtot
          IF (excluded_band(ibnd)) CYCLE
          ibnd1 = ibnd1 + 1
          !
          IF (TRIM(scdm_entanglement) == 'isolated') THEN
            focc(ibnd1) = 1.0d0
          ELSEIF (TRIM(scdm_entanglement) == 'erfc') THEN
            focc(ibnd1) = 0.5d0 * ERFC((et(ibnd, ik_g) * rytoev - scdm_mu) / scdm_sigma)
          ELSEIF (TRIM(scdm_entanglement) == 'gaussian') THEN
            focc(ibnd1) = EXP( - 1.0d0 * ((et(ibnd, ik_g) * rytoev - scdm_mu)**2) / (scdm_sigma**2))
          ELSE
            CALL errore('compute_amn_with_scdm', 'scdm_entanglement value not recognized.', 1)
          ENDIF
          !
          norm_psi = REAL(SUM(evc(       1:npw,        ibnd) * CONJG(evc(       1:npw,        ibnd)))) + &
                     REAL(SUM(evc(1 + npwx:npw + npwx, ibnd) * CONJG(evc(1 + npwx:npw + npwx, ibnd)) ))
          CALL mp_sum(norm_psi, intra_pool_comm)
          norm_psi = SQRT(norm_psi)
          !
          ! jml: nowfc = sum_G (psi(G) * exp(i*G*r)) * focc  * phase(iw) / norm_psi
          !
          DO iw = 1, n_wannier
            IF (piv_spin(iw) == 1) THEN ! spin up
              nowfc_tmp = SUM(evc(       1:npw,        ibnd) * phase_g(1:npw, iw))
            ELSE
              nowfc_tmp = SUM(evc(1 + npwx:npw + npwx, ibnd) * phase_g(1:npw, iw) )
            ENDIF
            nowfc(iw, ibnd1) = nowfc_tmp * phase(iw) * focc(ibnd1) / norm_psi
          ENDDO
        ENDDO ! ibnd
      ELSE
        DO ibnd = 1, nbtot
          IF (excluded_band(ibnd)) CYCLE
          ibnd1 = ibnd1 + 1
          !
          IF (TRIM(scdm_entanglement) == 'isolated') THEN
            focc(ibnd1) = 1.0d0
          ELSEIF (TRIM(scdm_entanglement) == 'erfc') THEN
            focc(ibnd1) = 0.5d0 * ERFC((et(ibnd, ik_g) * rytoev - scdm_mu) / scdm_sigma)
          ELSEIF (TRIM(scdm_entanglement) == 'gaussian') THEN
            focc(ibnd1) = EXP( - 1.0d0 * ((et(ibnd, ik_g) * rytoev - scdm_mu)**2) / (scdm_sigma**2))
          ELSE
            CALL errore('compute_amn_with_scdm', 'scdm_entanglement value not recognized.', 1)
          ENDIF
          !
          norm_psi = REAL(SUM(evc(1:npw, ibnd) * CONJG(evc(1:npw, ibnd))))
          CALL mp_sum(norm_psi, intra_pool_comm)
          norm_psi = SQRT(norm_psi)
          !
          ! jml: nowfc = sum_G (psi(G) * exp(i*G*r)) * focc  * phase(iw) / norm_psi
          !
          DO iw = 1, n_wannier
            nowfc_tmp = SUM(evc(1:npw, ibnd) * phase_g(1:npw, iw))
            nowfc(iw, ibnd1) = nowfc_tmp * phase(iw) * focc(ibnd1) / norm_psi
          ENDDO 
        ENDDO ! ibnd
      ENDIF
      CALL mp_sum(nowfc, intra_pool_comm) 
      !
      DEALLOCATE(phase_g, STAT = ierr) 
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating phase_g', 1) 
      !
      singval(:)   = zero
      Umat(:, :)   = czero
      VTmat(:, :)  = czero
      tmp_cwork(:) = czero
      rwork2(:)    = zero
      !
      CALL ZGESVD('s', 's', numbands, n_wannier, TRANSPOSE(CONJG(nowfc)), numbands, &
                  singval, Umat, numbands, VTmat, n_wannier, tmp_cwork, -1, rwork2, info)
      !
      lcwork = AINT(REAL(tmp_cwork(1)))
      ALLOCATE(cwork(lcwork), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error allocating cwork', 1)
      singval(:)  = zero
      Umat(:, :)  = czero
      VTmat(:, :) = czero
      cwork(:)    = czero
      rwork2(:)   = zero
      !
      ! vv: SVD to generate orthogonal projections
#if defined(__MPI)
      IF (meta_ionode) THEN
        CALL ZGESVD('s', 's', numbands, n_wannier, TRANSPOSE(CONJG(nowfc)), numbands, &
                    singval, Umat, numbands, VTmat, n_wannier, cwork, lcwork, rwork2, info)
        IF(info /= 0) CALL errore('compute_amn_with_scdm', 'Error in computing the SVD of the PSI matrix in the SCDM method', 1)
      ENDIF
      CALL mp_bcast(Umat,  meta_ionode_id, world_comm)
      CALL mp_bcast(VTmat, meta_ionode_id, world_comm)
#else
      CALL ZGESVD('s', 's', numbands, n_wannier, TRANSPOSE(CONJG(nowfc)), numbands, &
                  singval, Umat, numbands, VTmat, n_wannier, cwork, lcwork, rwork2, info)
      IF(info /= 0) CALL errore('compute_amn_with_scdm', 'Error in computing the SVD of the PSI matrix in the SCDM method', 1)
#endif
      DEALLOCATE(cwork, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating cwork', 1)
      !
      Amat(:, :) = czero
      Amat = MATMUL(Umat, VTmat)
      DO iw = 1, n_wannier
        ibnd1 = 0
        DO ibnd = 1, nbtot
          IF (excluded_band(ibnd)) CYCLE
          ibnd1 = ibnd1 + 1
          !
          a_mat(ibnd1, iw, ik_g) = Amat(ibnd1, iw)
        ENDDO ! bands
      ENDDO ! wannier fns
      !
    ENDDO ! k-points
    !
    CALL mp_sum(a_mat, inter_pool_comm)
    !
    DEALLOCATE(piv, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating piv', 1)
    DEALLOCATE(qr_tau, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating qr_tau', 1)
    DEALLOCATE(rwork, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating rwork', 1)
    DEALLOCATE(nowfc, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating nowfc', 1)
    DEALLOCATE(psi_gamma, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating psi_gamma', 1)
    DEALLOCATE(focc, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating focc', 1)
    DEALLOCATE(rwork2, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating rwork2', 1)
    DEALLOCATE(rpos, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating rpos', 1)
    DEALLOCATE(cpos, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating cpos', 1)
    DEALLOCATE(phase, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating phase', 1)
    DEALLOCATE(singval, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating singval', 1)
    DEALLOCATE(Umat, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating Umat', 1)
    DEALLOCATE(VTmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating VTmat', 1)
    DEALLOCATE(Amat, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating Amat', 1)
    !
    IF (noncolin) THEN 
      DEALLOCATE(piv_pos, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating piv_pos', 1)
      DEALLOCATE(piv_spin, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating piv_spin', 1)
    ENDIF
    !
#if defined(__MPI)
    IF (noncolin) THEN
      DEALLOCATE(psic_nc_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating psic_nc_all', 1)
    ELSE
      DEALLOCATE(psic_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_amn_with_scdm', 'Error deallocating psic_all', 1)
    ENDIF
#endif
    !
    WRITE(stdout, *)
    WRITE(stdout, '(5x, a)') 'AMN calculated'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_amn_with_scdm
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE orient_gf_spinor(npw)
    !-----------------------------------------------------------------------
    !! 
    !! FIXME
    !! 
    USE kinds,            ONLY : DP
    USE constants_epw,    ONLY : czero, cone, eps6
    USE noncollin_module, ONLY : npol
    USE wvfct,            ONLY : npwx
    USE wannierEPW,       ONLY : gf, n_proj, spin_qaxis, spin_eig, gf_spinor
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: npw
    !! Number of plane waves 
    !
    LOGICAL :: spin_z_pos
    !! Detect if spin quantisation axis is along +z
    LOGICAL :: spin_z_neg
    !! Detect if spin quantisation axis is along -z
    !
    INTEGER :: iw
    !! Counter on number of projections
    INTEGER  :: ipol
    !! Index of spin-up/spin-down polarizations
    INTEGER ::  istart, iend
    !! Starting/ending index on plane waves
    INTEGER :: ierr
    !! Error status
    !
    COMPLEX(KIND = DP) :: fac(2)
    !!
    !
    gf_spinor(:, :) = czero
    DO iw = 1, n_proj
      spin_z_pos = .FALSE.
      spin_z_neg = .FALSE.
      ! detect if spin quantisation axis is along z
      IF ((ABS(spin_qaxis(1, iw) - 0.0d0) < eps6) .AND. & 
          (ABS(spin_qaxis(2, iw) - 0.0d0) < eps6) .AND. &
          (ABS(spin_qaxis(3, iw) - 1.0d0) < eps6)) THEN 
        spin_z_pos = .TRUE.
      ELSEIF ((ABS(spin_qaxis(1, iw) - 0.0d0) < eps6) .AND. &
              (ABS(spin_qaxis(2, iw) - 0.0d0) < eps6) .AND. &
              (ABS(spin_qaxis(3, iw) + 1.0d0) < eps6)) THEN
        spin_z_neg = .TRUE.
      ENDIF
      IF (spin_z_pos .OR. spin_z_neg) THEN
        IF (spin_z_pos) THEN
          ipol = (3 - spin_eig(iw)) / 2
        ELSE
          ipol = (3 + spin_eig(iw)) / 2
        ENDIF
        istart = (ipol - 1) * npwx + 1
        iend   = (ipol - 1) * npwx + npw
        gf_spinor(istart:iend, iw) = gf(1:npw, iw)
      ELSE
        IF (spin_eig(iw) == 1) THEN
          fac(1) = (1.0d0 / SQRT(1 + spin_qaxis(3, iw))) &
                 * (spin_qaxis(3, iw) + 1 ) * cone
          fac(2) = (1.0d0 / SQRT(1 + spin_qaxis(3, iw))) &
                 * CMPLX(spin_qaxis(1, iw), spin_qaxis(2, iw), KIND = DP)
        ELSE
          fac(1) = (1.0d0 / SQRT(1 + spin_qaxis(3, iw))) &
                 * ( spin_qaxis(3, iw) ) * cone
          fac(2) = (1.0d0 / SQRT(1 - spin_qaxis(3, iw))) &
                 * CMPLX(spin_qaxis(1, iw), spin_qaxis(2,iw), KIND = DP)
        ENDIF
        gf_spinor(1:npw, iw) = gf(1:npw, iw) * fac(1)
        gf_spinor(1 + npwx:npw + npwx, iw) = gf(1:npw, iw) * fac(2)
      ENDIF
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE orient_gf_spinor
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_mmn_para()
    !-----------------------------------------------------------------------
    !!
    !!  adapted from compute_mmn in pw2wannier90.f90
    !!  parallelization on k-points has been added
    !!  10/2008 Jesse Noffsinger UC Berkeley
    !!
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout, meta_ionode
    USE io_files,        ONLY : diropn, prefix
    USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
    USE wavefunctions,   ONLY : evc, psic, psic_nc
    USE units_lr,        ONLY : lrwfc, iuwfc
    USE fft_base,        ONLY : dffts
    USE fft_interfaces,  ONLY : fwfft, invfft
    USE klist,           ONLY : nkstot, nks, igk_k
    USE klist_epw,       ONLY : xk_all, xk_loc
    USE gvect,           ONLY : g, ngm, gstart
    USE gvecw,           ONLY : gcutw
    USE cell_base,       ONLY : omega, tpiba, bg
    USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
    USE uspp,            ONLY : nkb, vkb
    USE uspp_param,      ONLY : upf, lmaxq, nh, nhm
    USE becmod,          ONLY : becp, calbec, allocate_bec_type, & 
                                deallocate_bec_type
    USE noncollin_module,ONLY : noncolin, npol
    USE spin_orb,        ONLY : lspinorb
    USE wannierEPW,      ONLY : m_mat, num_bands, nnb, iknum, g_kpb, kpb, ig_, &
                                excluded_band, write_mmn, zerophase
    USE constants_epw,   ONLY : czero, cone, twopi, zero
    USE io_var,          ONLY : iummn
#if defined(__NAG)
    USE f90_unix_io,     ONLY : flush
#endif
    USE mp,              ONLY : mp_sum
    USE mp_global,       ONLY : my_pool_id, npool, intra_pool_comm, inter_pool_comm
    USE mp_world,        ONLY : mpime
    USE kfold,           ONLY : ktokpmq
    USE io_epw,          ONLY : readwfc
    ! 
    IMPLICIT NONE
    !
    CHARACTER (LEN=80) :: filmmn
    !! file containing M_mn(k,b) matrix
    LOGICAL :: any_uspp
    !! Check if uspp are used
    LOGICAL :: exst
    !! Does the file exist
    !
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ikp
    !! Index of k-point nearest neighbour
    INTEGER :: ib
    !! Counter on b-vectors
    INTEGER :: npwq
    !! Number of planes waves at k+b
    INTEGER :: m
    !! Counter over bands
    INTEGER :: n
    !! Counter over bands
    INTEGER :: ibnd_m
    !! Band index
    INTEGER :: ibnd_n
    !! Band index
    INTEGER :: ih, jh
    !! Counters over the beta functions per atomic type
    INTEGER :: ikb, jkb, ijkb0
    !! Indexes over beta functions
    INTEGER :: na
    !! Counter over atoms
    INTEGER :: nt
    !! Number of type of atoms
    INTEGER :: ind
    !! Counter over nearest neighbours of all k-points
    INTEGER :: nbt
    !!
    INTEGER :: ipool
    !! Index of current pool
    INTEGER ::  nkq
    !! Index of k-point in the pool
    INTEGER ::  nkq_abs
    !! Absolute index of k-point
    INTEGER  :: ipol
    !! Index of spin-up/spin-down polarizations
    INTEGER ::  istart, iend
    !! Starting/ending index on plane waves 
    INTEGER ::  ik_g
    !! Temporary index of k-point, ik_g = nkq_abs
    INTEGER :: ind0 
    !! Starting index for k-point nearest neighbours in each pool
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE  :: igkq(:)
    !!
    REAL(KIND = DP) :: arg
    !! 2*pi*(k+b - k)*R
    COMPLEX(KIND = DP) :: mmn
    !! Element of M_mn matrix
    REAL(KIND = DP), ALLOCATABLE :: dxk(:, :)
    !! Temporary vector k+b - k 
    REAL(KIND = DP), ALLOCATABLE :: qg(:)
    !!  Square of dxk
    REAL(KIND = DP), ALLOCATABLE :: ylm(:, :)
    !! Spherical harmonics
    REAL(KIND = DP) :: zero_vect(3)
    !! Temporary zero vector
    REAL(KIND = DP) :: g_(3)
    !! Temporary vector G_k+b, g_(:) = g_kpb(:,ik,ib)
    COMPLEX(KIND = DP), EXTERNAL :: ZDOTC
    !! Scalar product of complex vectors
    COMPLEX(KIND = DP) :: phase1
    !! e^{i*2*pi*(k+b - k)*R}
    COMPLEX(KIND = DP), ALLOCATABLE :: phase(:)
    !! Phase
    COMPLEX(KIND = DP), ALLOCATABLE :: aux(:)
    !! Auxillary variable
    COMPLEX(KIND = DP), ALLOCATABLE :: aux_nc(:, :)
    !! NC auxillary variable
    COMPLEX(KIND = DP), ALLOCATABLE :: evcq(:, :)
    !! Wave functions psi_k+b
    COMPLEX(KIND = DP), ALLOCATABLE :: Mkb(:, :)
    !! Element of M_mn matrix
    COMPLEX(KIND = DP), ALLOCATABLE :: becp2(:, :)
    !! Beta functions
    COMPLEX(KIND = DP), ALLOCATABLE :: becp2_nc(:, :, :)
    !! Beta functions, noncolin
    COMPLEX(KIND = DP), ALLOCATABLE :: qb(:, :, :, :)
    !! Local variables for uspp 
    COMPLEX(KIND = DP), ALLOCATABLE :: qgm(:)
    !! Local variables for uspp
    COMPLEX(KIND = DP), ALLOCATABLE :: qq_so(:, :, :, :)
    !! Local variables for uspp
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mat_tmp(:, :, :, :)
    !! M_mn matrix
    !
    any_uspp = ANY(upf(:)%tvanp)
    !
    ALLOCATE(phase(dffts%nnr), STAT = ierr) 
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating phase', 1)
    ALLOCATE(igkq(npwx), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating igkq', 1)
    ALLOCATE(evcq(npol*npwx, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating evcq', 1)
    !
    IF (noncolin) THEN
      ALLOCATE(aux_nc(npwx, npol), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating aux_nc', 1)
    ELSE
      ALLOCATE(aux(npwx), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating aux', 1)
    ENDIF
    !
    ALLOCATE(m_mat(num_bands, num_bands, nnb, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating m_mat', 1)
    !
    ! close all the wfc files to allow access for each pool to all wfs
    CLOSE(UNIT = iuwfc, STATUS = 'keep')
    !
    WRITE(stdout, *)
    WRITE(stdout, '(5x, a)') 'MMN'
    !
    zero_vect(:) = zero
    m_mat(:, :, :, :) = czero
    !
    !   USPP
    !
    IF (any_uspp) THEN
      CALL allocate_bec_type(nkb, nbnd, becp)
      IF (noncolin) THEN
        ALLOCATE(becp2_nc(nkb, 2, nbnd), STAT = ierr)
        IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating becp2_nc', 1)
      ELSE
        ALLOCATE(becp2(nkb, nbnd), STAT = ierr)
        IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating becp2', 1)
      ENDIF
      !
      !     qb is  FT of Q(r)
      !
      nbt = nnb * iknum
      !
      ALLOCATE(qg(nbt), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating qg', 1)
      ALLOCATE(dxk(3, nbt), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating dxk', 1)
      !
      ind = 0
      DO ik = 1, iknum ! loop over k-points
        DO ib = 1, nnb ! loop over nearest neighbours for each k-point
          ind = ind + 1
          ikp = kpb(ik, ib) 
          !
          g_(:) = REAL(g_kpb(:, ik, ib))
          ! bring g_ from crystal to cartesian
          CALL cryst_to_cart(1, g_, bg, 1)
          dxk(:, ind) = xk_all(:, ikp) + g_(:) - xk_all(:, ik) 
          qg(ind) = SUM(dxk(:, ind) * dxk(:, ind)) 
        ENDDO
      ENDDO
      !
      ALLOCATE(ylm(nbt, lmaxq * lmaxq), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating ylm', 1)
      ALLOCATE(qgm(nbt), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating qgm', 1)
      ALLOCATE(qb(nhm, nhm, ntyp, nbt), STAT = ierr) 
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating qb', 1)
      ALLOCATE(qq_so(nhm, nhm, 4, ntyp), STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating qq_so', 1)
      !
      ! get spherical harmonics ylm 
      CALL ylmr2(lmaxq * lmaxq, nbt, dxk, qg, ylm)
      qg(:) = SQRT(qg(:)) * tpiba
      !
      DO nt = 1, ntyp
        IF (upf(nt)%tvanp) THEN
          DO ih = 1, nh(nt)
            DO jh = 1, nh(nt)
              CALL qvan2(nbt, ih, jh, nt, qg, qgm, ylm)
              qb(ih, jh, nt, 1:nbt) = omega * qgm(1:nbt)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      !
      DEALLOCATE(qg, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating qg', 1)
      DEALLOCATE(qgm, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating qgm', 1)
      DEALLOCATE(ylm, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating ylm', 1)
      !
    ENDIF
    !
    ALLOCATE(Mkb(nbnd, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating Mkb', 1)
    !
#if defined(__MPI)
    WRITE(stdout, '(6x, a, i5, a, i4, a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
    !
    ! returns in-pool index nkq and absolute index nkq_abs of first k-point in this pool 
    CALL ktokpmq( xk_loc(:,1), zero_vect, +1, ipool, nkq, nkq_abs )
    ind0 = (nkq_abs - 1) * nnb
    !
    ind = ind0
    DO ik = 1, nks 
      !
      ! returns in-pool index nkq and absolute index nkq_abs of xk
      CALL ktokpmq(xk_loc(:, ik), zero_vect, +1, ipool, nkq, nkq_abs)
      ik_g = nkq_abs
      !
      WRITE(stdout, '(5x, i8, " of ", i4, a)') ik , nks, ' on ionode'
      FLUSH(stdout)
      !
      ! read wfc at k
      CALL readwfc(my_pool_id + 1, ik, evc)
      !
      ! sorts k+G vectors in order of increasing magnitude, up to ecut
      CALL gk_sort(xk_loc(1, ik), ngm, g, gcutw, npw, igk_k(1, ik), g2kin)
      !
      !  USPP
      !
      IF (any_uspp) THEN
        CALL init_us_2(npw, igk_k(1, ik), xk_loc(1, ik), vkb)
        ! below we compute the product of beta functions with |psi>
        CALL calbec(npw, vkb, evc, becp)
      ENDIF
      !
      !
      DO ib = 1, nnb  ! loop on finite diff b-vectors
        ind = ind + 1
        !
        ikp = kpb(ik_g, ib)
        !
        CALL ktokpmq(xk_all(:, ikp), zero_vect, +1, ipool, nkq, nkq_abs)
        !
        ! read wfc at k+b
        CALL readwfc(ipool, nkq, evcq)
        !
        ! sorts k+G vectors in order of increasing magnitude, up to ecut
        CALL gk_sort(xk_all(1,ikp), ngm, g, gcutw, npwq, igkq, g2kin)
        !
        ! compute the phase
        IF (.NOT. zerophase(ik_g, ib)) THEN
          phase(:) = czero
          IF (ig_(ik_g, ib) > 0) phase(dffts%nl(ig_(ik_g, ib))) = cone
          CALL invfft('Wave', phase, dffts)
        ENDIF
        !
        !  USPP
        !
        IF (any_uspp) THEN
          CALL init_us_2(npwq, igkq, xk_all(1, ikp), vkb)
          ! below we compute the product of beta functions with |psi>
          IF (noncolin) THEN
            CALL calbec(npwq, vkb, evcq, becp2_nc)
            !
            IF (lspinorb) THEN
              qq_so(:, :, :, :) = czero
              CALL transform_qq_so(qb(:, :, :, ind), qq_so)
            ENDIF
            !
          ELSE
            CALL calbec(npwq, vkb, evcq, becp2)
          ENDIF
        ENDIF
        !
        Mkb(:, :) = czero
        !
        IF (any_uspp) THEN
          ijkb0 = 0
          DO nt = 1, ntyp
            IF (upf(nt)%tvanp) THEN
              DO na = 1, nat
                !
                arg = twopi * DOT_PRODUCT(dxk(:, ind), tau(:, na))
                phase1 = CMPLX(COS(arg), -SIN(arg), KIND = DP)
                !
                IF (ityp(na) == nt) THEN
                  DO jh = 1, nh(nt)
                    jkb = ijkb0 + jh
                    DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      !
                      DO m = 1, nbnd
                        IF (excluded_band(m)) CYCLE
                        IF (noncolin) THEN
                          DO n = 1, nbnd
                            IF (excluded_band(n)) CYCLE
                            IF (lspinorb) THEN
                              Mkb(m, n) = Mkb(m, n) + phase1 * &
                                 (qq_so(ih, jh, 1, nt) * CONJG(becp%nc(ikb, 1, m)) * becp2_nc(jkb, 1, n) + &
                                  qq_so(ih, jh, 2, nt) * CONJG(becp%nc(ikb, 1, m)) * becp2_nc(jkb, 2, n) + &
                                  qq_so(ih, jh, 3, nt) * CONJG(becp%nc(ikb, 2, m)) * becp2_nc(jkb, 1, n) + &
                                  qq_so(ih, jh, 4, nt) * CONJG(becp%nc(ikb, 2, m)) * becp2_nc(jkb, 2, n))
                            ELSE
                              Mkb(m, n) = Mkb(m, n) + phase1 * qb(ih, jh, nt, ind) * &
                                 (CONJG(becp%nc(ikb, 1, m)) * becp2_nc(jkb, 1, n) + &
                                  CONJG(becp%nc(ikb, 2, m)) * becp2_nc(jkb, 2, n))
                            ENDIF
                          ENDDO
                        ELSE
                          DO n = 1, nbnd
                            IF (excluded_band(n)) CYCLE
                            Mkb(m, n) = Mkb(m, n) + phase1 * qb(ih, jh, nt, ind) * &
                                        CONJG(becp%k(ikb, m)) * becp2(jkb, n)
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
                IF (ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
              ENDDO
            ENDIF !tvanp
          ENDDO !ntyp
        ENDIF ! any_uspp
        !
        DO m = 1, nbnd  ! loop on band m
          IF (excluded_band(m)) CYCLE
          !
          IF (noncolin) THEN
            psic_nc(:, :) = czero
            DO ipol = 1, 2 !npol
              istart = (ipol - 1) * npwx + 1
              iend   = (ipol - 1) * npwx + npw
              psic_nc(dffts%nl(igk_k(1:npw, ik)), ipol) = evc(istart:iend, m)
              IF (.NOT. zerophase(ik_g, ib)) THEN
                CALL invfft('Wave', psic_nc(:, ipol), dffts)
                psic_nc(1:dffts%nnr, ipol) = psic_nc(1:dffts%nnr, ipol) * &
                                             phase(1:dffts%nnr)
                CALL fwfft('Wave', psic_nc(:, ipol), dffts)
              ENDIF
              aux_nc(1:npwq, ipol) = psic_nc(dffts%nl(igkq(1:npwq)), ipol)
            ENDDO
          ELSE
            psic(:) = czero
            psic(dffts%nl(igk_k(1:npw, ik))) = evc(1:npw, m)
            IF (.NOT. zerophase(ik_g,ib)) THEN
              CALL invfft('Wave', psic, dffts)
              psic(1:dffts%nnr) = psic(1:dffts%nnr) * phase(1:dffts%nnr)
              CALL fwfft('Wave', psic, dffts)
            ENDIF
            aux(1:npwq) = psic(dffts%nl(igkq(1:npwq)))
          ENDIF
          !  aa = 0.d0
          !
          !  Mkb(m,n) = Mkb(m,n) + \sum_{ijI} qb_{ij}^I * e^-i(b*tau_I)
          !             <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 > 
          !
          IF (noncolin) THEN
            DO n = 1, nbnd
              IF (excluded_band(n)) CYCLE
              mmn = czero
              mmn = mmn + ZDOTC(npwq, aux_nc(1, 1), 1, evcq(       1, n), 1) &
                        + ZDOTC(npwq, aux_nc(1, 2), 1, evcq(1 + npwx, n), 1)
              CALL mp_sum(mmn, intra_pool_comm)
              Mkb(m, n) = mmn + Mkb(m, n)
              !  aa = aa + ABS(mmn)**2
            ENDDO
          ELSE
            DO n = 1, nbnd
              IF (excluded_band(n)) CYCLE
              mmn = ZDOTC(npwq, aux, 1, evcq(1, n), 1)
              CALL mp_sum(mmn, intra_pool_comm)
              Mkb(m, n) = mmn + Mkb(m, n)
              !  aa = aa + ABS(mmn)**2
            ENDDO
          ENDIF
        ENDDO   ! m
        !
        ibnd_n = 0
        DO n = 1, nbnd
          IF (excluded_band(n)) CYCLE
          ibnd_n = ibnd_n + 1
          !
          ibnd_m = 0
          DO m = 1, nbnd
            IF (excluded_band(m)) CYCLE
            ibnd_m = ibnd_m + 1
            !
            m_mat(ibnd_m, ibnd_n, ib, ik_g) = Mkb(m, n)
          ENDDO
        ENDDO
        !
      ENDDO !ib
    ENDDO !ik
    !
    CALL mp_sum(m_mat, inter_pool_comm)
    !
    ! RM - write mmn to file (file needed with vme = true)
    IF (meta_ionode) THEN
      write_mmn = .TRUE.
      IF (write_mmn) THEN
        ALLOCATE(m_mat_tmp(nbnd, nbnd, nnb, nkstot), STAT = ierr)
        IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error allocating m_mat_tmp', 1)
        m_mat_tmp(:, :, :, :) = czero
        !
        filmmn = TRIM(prefix)//'.mmn'
        OPEN(UNIT = iummn, FILE = filmmn, FORM = 'formatted')
        DO ik = 1, nkstot
          DO ib = 1, nnb
            !
            ibnd_n = 0
            DO n = 1, nbnd
              IF (excluded_band(n)) CYCLE
              ibnd_n = ibnd_n + 1
              !
              ibnd_m = 0
              DO m = 1, nbnd
                IF (excluded_band(m)) CYCLE
                ibnd_m = ibnd_m + 1
                !   
                m_mat_tmp(m, n, ib, ik) =  m_mat(ibnd_m, ibnd_n, ib, ik)
              ENDDO
            ENDDO
            !
            DO n = 1, nbnd
              DO m = 1, nbnd
                WRITE(iummn,*) m_mat_tmp(m, n, ib, ik)
              ENDDO
            ENDDO
            !
          ENDDO
        ENDDO
        !
        CLOSE(iummn)
        DEALLOCATE(m_mat_tmp, STAT = ierr)
        IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating m_mat_tmp', 1)
      ENDIF
    ENDIF
    !
    DEALLOCATE(Mkb, STAT = ierr) 
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating Mkb', 1)
    DEALLOCATE(phase, STAT = ierr) 
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating phase', 1)
    DEALLOCATE(evcq, STAT = ierr) 
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating evcq', 1)
    DEALLOCATE(igkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating igkq', 1)
    IF (noncolin) THEN
      DEALLOCATE(aux_nc, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating aux_nc', 1)
    ELSE
      DEALLOCATE(aux, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating aux', 1)
    ENDIF
    !
    IF (any_uspp) THEN
      DEALLOCATE(dxk, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating dxk', 1) 
      DEALLOCATE(qb, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating qb', 1)
      DEALLOCATE(qq_so, STAT = ierr)
      IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating qq_so', 1)
      CALL deallocate_bec_type(becp)
      IF (noncolin) THEN
        DEALLOCATE(becp2_nc, STAT = ierr)
        IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating becp2_nc', 1)
      ELSE
        DEALLOCATE(becp2, STAT = ierr)
        IF (ierr /= 0) CALL errore('compute_mmn_para', 'Error deallocating becp2', 1)
      ENDIF
    ENDIF
    !
    WRITE(stdout, '(5x, a)') 'MMN calculated'
    !
    ! reopen wfc here, leaving UNIT = 20 in the same state
    iuwfc = 20
    CALL diropn(iuwfc, 'wfc', lrwfc, exst)  
    !
    RETURN
    !
    !------------------------------------------------------------------------
    END SUBROUTINE compute_mmn_para
    !------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_pmn_para()
    !-----------------------------------------------------------------------
    !
    !!
    !!  Computes dipole matrix elements.
    !!  This can be used to compute the velocity in the local approximation.
    !!  The commutator with the non-local pseudopotetial is neglected.
    !!
    !!  06/2010  Jesse Noffsinger: adapted from compute_amn_para
    !!  08/2016  Samuel Ponce: adapted to work with SOC
    !! 
    !
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout
    USE mp,              ONLY : mp_sum
    USE mp_global,       ONLY : my_pool_id
    USE klist,           ONLY : nks, igk_k
    USE klist_epw,       ONLY : xk_loc
    USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
    USE gvecw,           ONLY : gcutw
    USE wavefunctions,   ONLY : evc
    USE units_lr,        ONLY : lrwfc, iuwfc
    USE gvect,           ONLY : g, ngm
    USE cell_base,       ONLY : tpiba
    USE noncollin_module,ONLY : noncolin
    USE elph2,           ONLY : dmec
    USE constants_epw,   ONLY : czero
    USE uspp_param,      ONLY : upf
    USE becmod,          ONLY : becp, deallocate_bec_type, allocate_bec_type
    USE uspp,            ONLY : nkb
    USE wannierEPW,      ONLY : n_wannier
    USE kfold,           ONLY : ktokpmq
    USE io_epw,          ONLY : readwfc
    !
    IMPLICIT NONE
    !  
    LOGICAL :: any_uspp
    !! Check if USPP is present
    !
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd, jbnd
    !! Counter on bands
    INTEGER :: ig
    !! Counter on G vector
    INTEGER :: ierr
    !! Error status
    ! 
    COMPLEX(KIND = DP) :: dipole_aux(3, nbnd, nbnd)
    !! Auxilary dipole
    COMPLEX(KIND = DP) :: caux 
    !! Wavefunction squared
    !
    any_uspp = ANY(upf(:)%tvanp)
    !
    IF (any_uspp) CALL errore('pw2wan90epw', &
      'dipole matrix calculation not implimented with USPP - set vme=.TRUE.',1)
    !
    ALLOCATE(dmec(3, nbnd, nbnd, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_pmn_para', 'Error allocating dmec', 1)
    !
    ! initialize
    dmec(:, :, :, :) = czero
    dipole_aux(:, :, :) = czero
    !
    IF (any_uspp) THEN
      CALL deallocate_bec_type(becp)
      CALL allocate_bec_type(nkb, n_wannier, becp)
    ENDIF
    !
    ! computes velocity dmec = v_mn(\alpha,k) in the local approximation 
    ! acording to [Eqn. 60 of Comp. Phys. Commun. 209, 116 (2016)] 
    ! v_mn(\alpha,k) = k \delta_mn + \sum_G * CONJG(c_mk(G)) * c_nk(G) * G 
    !
    DO ik = 1, nks
      !
      ! read wfc for the given kpt
      CALL readwfc(my_pool_id + 1, ik, evc)
      !
      ! sorts k+G vectors in order of increasing magnitude, up to ecut
      CALL gk_sort(xk_loc(:, ik), ngm, g, gcutw, npw, igk_k(:, ik), g2kin)
      !
      dipole_aux = czero
      DO jbnd = 1, nbnd
        DO ibnd = 1, nbnd
          !
          IF (ibnd == jbnd) CYCLE
          !
          ! taken from PP/epsilon.f90 SUBROUTINE dipole_calc
          DO ig = 1, npw
            IF (igk_k(ig, ik) > SIZE(g, 2) .OR. igk_k(ig, ik) < 1) CYCLE
            !
            caux = CONJG(evc(ig, ibnd)) * evc(ig, jbnd) 
            !
            IF (noncolin) THEN
              !
              caux = caux + CONJG(evc(ig + npwx, ibnd)) * evc(ig + npwx, jbnd)
              !
            ENDIF
            !
            dipole_aux(:, ibnd, jbnd) = dipole_aux(:, ibnd, jbnd) + &
                                        (g(:, igk_k(ig, ik))) * caux
            !
          ENDDO
          !
        ENDDO !bands i
      ENDDO ! bands j
      !
      ! metal diagonal part
      DO ibnd = 1, nbnd
        DO ig = 1, npw
          IF (igk_k(ig, ik) > SIZE(g, 2) .OR. igk_k(ig, ik) < 1) CYCLE
          !
          caux = CONJG(evc(ig, ibnd)) * evc(ig, ibnd) 
          !
          IF (noncolin) THEN
            !
            caux = caux + CONJG(evc(ig + npwx, ibnd)) * evc(ig + npwx, ibnd)
            !
          ENDIF
          !
          dipole_aux(:, ibnd, ibnd) = dipole_aux(:, ibnd, ibnd) + &
                                      (g(:, igk_k(ig, ik)) + xk_loc(:, ik)) * caux
          !
        ENDDO
      ENDDO
      ! need to divide by 2pi/a to fix the units
      dmec(:, :, :, ik) = dipole_aux(:, :, :) * tpiba
      !
    ENDDO  ! k-points
    !
    WRITE(stdout, '(/5x, a)') 'Dipole matrix elements calculated'
    WRITE(stdout, *)
    !DBSP 
    !WRITE(stdout, *) 'dmec ', SUM(dmec)
    !IF (meta_ionode) THEN
    !   WRITE(stdout, *) 'dmec(:, :, :, 1) ',sum(dmec(:, :, :, 1))
    !   WRITE(stdout, *) 'dmec(:, :, :, 2) ',sum(dmec(:, :, :, 2))
    !ENDIF
    !
    RETURN
    ! 
    !------------------------------------------------------------------------
    END SUBROUTINE compute_pmn_para
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_filukk
    !-----------------------------------------------------------------------
    !
    ! Here we compute and write out the final ukk matrix which is used by
    ! epw.x to localize the electron wavefuctions (and therefore the ep-matrix 
    ! elements)
    ! 10/2008 Jesse Noffsinger UC Berkeley
    ! 07/2010 Fixed the rotation for ndimwin when lower bands are not included
    !
    USE kinds,        ONLY : DP
    USE io_var,       ONLY : iunukk
    USE wvfct,        ONLY : nbnd
    USE wannierEPW,   ONLY : n_wannier, iknum, u_mat, u_mat_opt, lwindow, &
                             excluded_band, num_bands, wann_centers
    USE epwcom,       ONLY : filukk
    USE constants_epw,ONLY : czero, bohr
    USE io_global,    ONLY : meta_ionode
    USE cell_base,    ONLY : alat
    !
    IMPLICIT NONE
    !
    LOGICAL, ALLOCATABLE :: lwindow_tmp(:, :)
    !! .TRUE. if the band ibnd lies within the outer window at k-point ik
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: ibnd1
    !! Band index
    INTEGER :: iw
    !! Counter on number of Wannier functions
    INTEGER :: ierr
    !! Error status
    INTEGER :: ndimwin(iknum)
    !! Number of bands within outer window at each k-point
    !
    COMPLEX(KIND = DP), ALLOCATABLE :: u_kc_tmp(:, :, :)
    !! Temporaty rotation matrix
    COMPLEX(KIND = DP), ALLOCATABLE :: u_kc(:, :, :)
    !! Rotation matrix
    !
    ! RM: Band-dimension of u_mat_opt and lwindow is num_bands while
    !     that of u_kc is nbnd to be compatible when reading umat in loadumat. 
    !     Care needs to be taken if exclude_bands is used because
    !     num_bands = nbnd - nexband
    !
    !     u_mat_opt(num_bands, n_wannier, nkstot) in wannier_run
    !     u_mat(n_wannier, n_wannier, nkstot) in wannier_run
    !     lwindow(num_bands, nkstot) in wannier_run
    !     u_kc_tmp(num_bands, n_wannier, nkstot) in write_filukk
    !     u_kc(nbnd, n_wannier, nkstot) in write_filukk
    !
    IF (meta_ionode) THEN
      !
      ndimwin(:) = 0
      DO ik = 1, iknum
        DO ibnd = 1, num_bands
          IF (lwindow(ibnd, ik)) ndimwin(ik) = ndimwin(ik) + 1
        ENDDO
      ENDDO
      !
      ! get the final rotation matrix, which is the product of the optimal
      ! subspace and the rotation among the n_wannier wavefunctions
      !
      ALLOCATE(u_kc_tmp(num_bands, n_wannier, iknum), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_filukk', 'Error allocating u_kc_tmp', 1)
      u_kc_tmp(:, :, :) = czero
      !
      DO ik = 1, iknum
        !
        u_kc_tmp(1:ndimwin(ik), 1:n_wannier, ik) = &
             MATMUL(u_mat_opt(1:ndimwin(ik), :, ik), u_mat(:, 1:n_wannier, ik))
        !
      ENDDO
      !
      ALLOCATE(u_kc(nbnd, n_wannier, iknum), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_filukk', 'Error allocating u_kc', 1)
      u_kc(:, :, :) = czero
      !
      OPEN(UNIT = iunukk, FILE = filukk, FORM = 'formatted')
      ! u_kc(1:num_bands, :, :) = u_kc_tmp(1:num_bands, :, :)
      ! u_kc(num_bands+1:nbnd, :, :) = czero
      DO ik = 1, iknum
        DO iw = 1, n_wannier
          ibnd1 = 0
          DO ibnd = 1, nbnd
            IF (excluded_band(ibnd)) CYCLE
            ibnd1 = ibnd1 + 1
            u_kc(ibnd, iw, ik) = u_kc_tmp(ibnd1, iw, ik)
          ENDDO
        ENDDO
        !
        DO ibnd = 1, nbnd
          DO iw = 1, n_wannier
            WRITE(iunukk, *) u_kc(ibnd, iw, ik)
          ENDDO
        ENDDO
      ENDDO
      !
      ! needs also lwindow when disentanglement is used
      ALLOCATE(lwindow_tmp(nbnd, iknum), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_filukk', 'Error allocating lwindow_tmp', 1)
      lwindow_tmp(:, :) = .FALSE.
      !
      DO ik = 1, iknum
        ibnd1 = 0
        DO ibnd = 1, nbnd
          IF (excluded_band(ibnd)) CYCLE
          ibnd1 = ibnd1 + 1
          lwindow_tmp(ibnd, ik) = lwindow(ibnd1, ik) 
        ENDDO
        !
        DO ibnd = 1, nbnd
          WRITE(iunukk, *) lwindow_tmp(ibnd, ik)
        ENDDO
      ENDDO
      !
      DO ibnd = 1, nbnd
        WRITE(iunukk, *) excluded_band(ibnd)
      ENDDO
      ! 
      ! Now write the Wannier centers to files
      DO iw = 1, n_wannier
        ! SP : Need more precision other WS are not determined properly. 
        !WRITE (iuukk, '(3f12.8)') wann_centers(:, iw) / alat / bohr
        WRITE (iunukk, '(3E22.12)') wann_centers(:, iw) / alat / bohr
      ENDDO
      !
      CLOSE(iunukk)
      !
      DEALLOCATE(u_kc_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_filukk', 'Error deallocating u_kc_tmp', 1)
      DEALLOCATE(u_kc, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_filukk', 'Error deallocating u_kc', 1)
      DEALLOCATE(lwindow_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_filukk', 'Error deallocating lwindow_tmp', 1) 
      !
    ENDIF
    !
    RETURN
    !
    !------------------------------------------------------------------------
    END SUBROUTINE write_filukk
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE phases_a_m
    !-----------------------------------------------------------------------
    !!
    !! We will set phases here on the matrices. It should not affect 
    !! the spreads and centers found in w90, but it will leave 
    !! u_mat_opt and u_mat to reflect the known phases. 
    !!
    !
    USE kinds,           ONLY : DP
    USE mp_global,       ONLY : inter_pool_comm
    USE mp,              ONLY : mp_sum
    USE io_global,       ONLY : meta_ionode
    USE klist,           ONLY : nkstot, nks
    USE klist_epw,       ONLY : xk_loc
    USE wvfct,           ONLY : nbnd
    USE wannierEPW,      ONLY : a_mat, m_mat, num_bands, n_wannier, n_proj, & 
                                nnb, kpb, iknum, excluded_band
    USE elph2,           ONLY : umat, umat_all
    USE constants_epw,   ONLY : czero, cone, zero
    USE kfold,           ONLY : ktokpmq
    !
    IMPLICIT NONE
    !
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ikb
    !! Index of k-point nearest neighbour
    INTEGER :: ib
    !! Counter on b-vectors
    INTEGER :: ipool
    !! Index of current pool
    INTEGER ::  nkq
    !! Index of k-point in the pool
    INTEGER ::  nkq_abs
    !! Absolute index of k-point
    INTEGER ::  ik_g
    !! Temporary index of k-point, ik_g = nkq_abs
    INTEGER :: m, n
    !! Counter over bands
    INTEGER :: ibnd_m, ibnd_n
    !! Band index
    INTEGER :: iw
    !! Counter on number of projections
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: zero_vect(3)
    !! Temporary zero vector
    !
    COMPLEX(KIND = DP), ALLOCATABLE :: a_mat_tmp(:, :, :)
    !! Temporary a_mat matrices
    COMPLEX(KIND = DP), ALLOCATABLE :: a_mat_tmp1(:, :, :)
    !! Temporary a_mat matrices
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mn_tmp1(:, :)
    !! Temporary m_mat matrices
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mn_tmp2(:, :)
    !! Temporary m_mat matrices
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mn_tmp3(:, :, :, :)
    !! Temporary m_mat matrices
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mat_tmp(:, :, :, :)
    !! Temporary m_mat matrices
    !
    ! RM: Band-dimension of a_mat and m_mat is num_bands while that of
    !     umat and umat_all is nbnd. This causes a problem if exclude_bands 
    !     is used because num_bands = nbnd - nexband 
    !     a_mat(num_bands,n_wannier,nkstot) in compute_amn_para    
    !     m_mat(num_bands,n_wannier,nnb,nkstot) in compute_mmn_para
    !     umat(nbnd,nbnd,nks) and umat_mat(nbnd,nbnd,nkstot) in setphases_wrap
    !
    ALLOCATE(a_mat_tmp(nbnd, n_wannier, iknum), STAT = ierr) 
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error allocating a_mat_tmp', 1)
    ALLOCATE(a_mat_tmp1(nbnd, n_wannier, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error allocating a_mat_tmp1', 1)
    ALLOCATE(m_mat_tmp(nbnd, nbnd, nnb, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error allocating m_mat_tmp', 1)
    ALLOCATE(m_mn_tmp1(nbnd, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error allocating m_mn_tmp1', 1)
    ALLOCATE(m_mn_tmp2(nbnd, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error allocating m_mn_tmp2', 1)
    ALLOCATE(m_mn_tmp3(nbnd, nbnd, nnb, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error allocating m_mn_tmp3', 1)
    !
    ! zero all temporary/work quantities
    !
    zero_vect(:) = zero
    a_mat_tmp(:, :, :) = czero
    a_mat_tmp1(:, :, :) = czero
    m_mn_tmp1(:, :) = czero
    m_mn_tmp2(:, :) = czero
    m_mn_tmp3(:, :, :, :) = czero
    m_mat_tmp(:, :, :, :) = czero
    !
    ! full size a_mat_tmp1 and m_mat_tmp matrices to nbnd bands
    !
    DO ik = 1, nkstot
      DO iw = 1, n_proj
        ibnd_m = 0
        DO m = 1, nbnd
          IF (excluded_band(m)) CYCLE
          ibnd_m = ibnd_m + 1
          a_mat_tmp1(m, iw, ik) = a_mat(ibnd_m, iw, ik)
        ENDDO
      ENDDO
      !
      DO ib = 1, nnb
        ibnd_n = 0
        DO n = 1, nbnd
          IF (excluded_band(n)) CYCLE
          ibnd_n = ibnd_n + 1
          ibnd_m = 0
          DO m = 1, nbnd
            IF (excluded_band(m)) CYCLE
            ibnd_m = ibnd_m + 1
            m_mat_tmp(m, n, ib, ik) = m_mat(ibnd_m, ibnd_n, ib, ik)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    DO ik = 1, nks
      !
      ! returns in-pool index nkq and absolute index nkq_abs of xk
      CALL ktokpmq(xk_loc(:, ik), zero_vect, +1, ipool, nkq, nkq_abs)
      ik_g = nkq_abs
      !
      !  GF_n are the guiding functions which are our initial guesses 
      !  Amn(k) = <psi_k,m|GF_n>.  
      !  We want U(k)^\dagger<psi_k,m|GF_m>
      !
      ! CALL ZGEMM('c', 'n', nbnd, n_wannier, nbnd, cone, umat(:, :, ik), & 
      !      nbnd, a_mat(:, :, ik_g), nbnd, czero, a_mat_tmp(:, :, ik_g), nbnd)
      CALL ZGEMM('c', 'n', nbnd, n_wannier, nbnd, cone, umat(:, :, ik), &
           nbnd, a_mat_tmp1(:, :, ik_g), nbnd, czero, a_mat_tmp(:, :, ik_g), nbnd)
      !
      DO ib = 1,nnb
        ikb = kpb(ik_g, ib)
        !
        ! Mmn(k,k+b)  = <psi_k_m| psi_(k+b)_n> so we need
        !  (U(k)^\dagger <psi_k_m| ) * (|psi_k+b_n> U(k+b)
        ! = U(k)^\dagger (M_mn) = m_mat_tmp, 
        ! Mmn(k,k+b)' = m_mat_tmp*U(k+b) 
        !
        ! CALL ZGEMM('c', 'n', nbnd, nbnd, nbnd, cone, umat(:, :, ik), & 
        !      nbnd, m_mat(:, :, ib, ik_g), nbnd, czero, m_mn_tmp1(:, :), nbnd)
        CALL ZGEMM('c', 'n', nbnd, nbnd, nbnd, cone, umat(:, :, ik), &
             nbnd, m_mat_tmp(:, :, ib, ik_g), nbnd, czero, m_mn_tmp1(:, :), nbnd)
        CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, m_mn_tmp1(:, :), & 
             nbnd, umat_all(:, :, ikb), nbnd, czero, m_mn_tmp2(:, :), nbnd)
        ! 
        ! m_mn_tmp1 = MATMUL(CONJG(TRANSPOSE(umat(:, :, ik))), m_mat(:, :, ib, ik_g))
        ! m_mn_tmp2 = MATMUL(m_mn_tmp1, umat_g(:, :, ikb))
        !
        m_mn_tmp3(:, :, ib, ik_g) = m_mn_tmp2(:, :)
      ENDDO
    ENDDO
    CALL mp_sum(a_mat_tmp, inter_pool_comm)
    CALL mp_sum(m_mn_tmp3, inter_pool_comm)
    !
    ! a_mat(:, :, :) = a_mat_tmp(:, :, :)
    ! m_mat(:, :, :, :) = m_mn_tmp3(:, :, :, :)
    !
    ! slim down a_mat and m_mat matrices to num_bands=nbnd-nexband bands
    ! 
    DO ik = 1, nkstot
      DO iw = 1, n_proj
        ibnd_m = 0
        DO m = 1, nbnd
          IF (excluded_band(m)) CYCLE
          ibnd_m = ibnd_m + 1
          a_mat(ibnd_m, iw, ik) = a_mat_tmp(m, iw, ik)
        ENDDO
      ENDDO
      !
      DO ib = 1, nnb
        ibnd_n = 0
        DO n = 1, nbnd
          IF (excluded_band(n)) CYCLE
          ibnd_n = ibnd_n + 1
          ibnd_m = 0
          DO m = 1, nbnd
            IF (excluded_band(m)) CYCLE
            ibnd_m = ibnd_m + 1
            m_mat(ibnd_m, ibnd_n, ib, ik) = m_mn_tmp3(m, n, ib, ik) 
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    DEALLOCATE(a_mat_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error deallocating a_mat_tmp', 1)
    DEALLOCATE(a_mat_tmp1, STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error deallocating a_mat_tmp1', 1)
    DEALLOCATE(m_mat_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error deallocating m_mat_tmp', 1)
    DEALLOCATE(m_mn_tmp1, STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error deallocating m_mn_tmp1', 1)
    DEALLOCATE(m_mn_tmp2, STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error deallocating m_mn_tmp2', 1)
    DEALLOCATE(m_mn_tmp3, STAT = ierr)
    IF (ierr /= 0) CALL errore('phases_a_m', 'Error deallocating m_mn_tmp3', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE phases_a_m
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE generate_guiding_functions(ik)
    !-----------------------------------------------------------------------
    !! 
    !! Guiding functions
    !!  
    !
    USE kinds,          ONLY : DP
    USE mp_global,      ONLY : intra_pool_comm
    USE mp,             ONLY : mp_sum
    USE wvfct,          ONLY : npw
    USE gvect,          ONLY : g
    USE cell_base,      ONLY : tpiba
    USE wannierEPW,     ONLY : n_proj, gf, center_w, csph, alpha_w, r_w
    USE klist,          ONLY : igk_k 
    USE klist_epw,      ONLY : xk_loc
    USE constants_epw,  ONLY : zero, czero, ci, twopi, eps8
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Index of k-point
    !
    ! Local variables
    INTEGER, PARAMETER :: lmax = 3
    !! Nr of spherical harmonics
    INTEGER, PARAMETER :: lmax2 = (lmax + 1)**2
    !! Total nr. of spherical harmonics
    INTEGER :: iw
    !! Counter on wannier projections
    INTEGER :: ig
    !! Counter on G-vectors
    INTEGER :: iig
    !! Index of G vector
    INTEGER :: lm
    !! Counter on spherical harmonics
    INTEGER :: l
    !! Momentum l
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: arg
    !! 2*pi*(k+G)*r 
    REAL(KIND = DP) :: anorm
    !! Anormal
    REAL(KIND = DP), EXTERNAL :: DDOT
    !! Scalar product of two vectors
    REAL(KIND = DP), ALLOCATABLE :: gk(:, :)
    !! k+G vectors
    REAL(KIND = DP), ALLOCATABLE :: qg(:)
    !! Norm of k+G
    REAL(KIND = DP), ALLOCATABLE :: ylm(:, :)
    !! Spherical harmonics
    REAL(KIND = DP), ALLOCATABLE :: radial(:, :)
    !! Radiam
    COMPLEX(KIND = DP) :: lphase
    !! (-i)^l
    COMPLEX(KIND = DP), EXTERNAL :: ZDOTC
    !! Scalar product of two complex vectors
    COMPLEX(KIND = DP), ALLOCATABLE :: sk(:) 
    !! e^{-i*2*pi*(k+G)*r}
    !
    ALLOCATE(gk(3, npw), STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error allocating gk', 1)
    ALLOCATE(qg(npw), STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error allocating qg', 1)
    ALLOCATE(ylm(npw, lmax2), STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error allocating ylm', 1)
    ALLOCATE(sk(npw), STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error allocating sk', 1)
    ALLOCATE(radial(npw, 0:lmax), STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error allocating radial', 1)
    gk(:, :) = zero
    qg(:) = zero
    ylm(:, :) = zero
    sk(:) = czero
    radial(:, 0:lmax) = zero
    !
    DO ig = 1, npw
      gk(:, ig) = xk_loc(:, ik) + g(:, igk_k(ig, ik))
      qg(ig) = SUM(gk(:, ig) * gk(:, ig))
    ENDDO
    !
    ! get spherical harmonics ylm up to lmax
    CALL ylmr2(lmax2, npw, gk, qg, ylm)
    ! define qg as the norm of (k+g) in a.u.
    qg(:) = SQRT(qg(:)) * tpiba
    !
    ! RM changed according to QE4.0.3/PP/pw2wannier90 
    DO iw = 1, n_proj
      !
      gf(:, iw) = czero
      !
      CALL radialpart(npw, qg, alpha_w(iw), r_w(iw), lmax, radial) 
      !
      DO lm = 1, lmax2
        IF (ABS(csph(lm, iw)) < eps8) CYCLE
        l = INT(SQRT(lm - 1.d0))
        lphase = (- ci)**l
        !
        DO ig = 1, npw
          gf(ig, iw) = gf(ig, iw) + csph(lm, iw) * ylm(ig, lm) * radial(ig, l) * lphase
        ENDDO !ig
      ENDDO ! lm
      !
      DO ig = 1, npw
        iig = igk_k(ig, ik)
        arg = twopi * DDOT(3, gk(:, ig), 1, center_w(:, iw), 1)
        ! center_w are cartesian coordinates in units of alat 
        sk(ig) = CMPLX(COS(arg), - SIN(arg), KIND = DP)
        gf(ig, iw) = gf(ig, iw) * sk(ig) 
      ENDDO
      anorm = REAL(ZDOTC(npw, gf(1, iw), 1, gf(1, iw), 1))
      CALL mp_sum(anorm, intra_pool_comm)
      gf(:, iw) = gf(:, iw) / SQRT(anorm)
    ENDDO
    !
    DEALLOCATE(gk, STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error deallocating gk', 1)
    DEALLOCATE(qg, STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error deallocating qg', 1)
    DEALLOCATE(ylm, STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error deallocating ylm', 1)
    DEALLOCATE(sk, STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error deallocating sk', 1)
    DEALLOCATE(radial, STAT = ierr)
    IF (ierr /= 0) CALL errore('generate_guiding_functions', 'Error deallocating radial', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE generate_guiding_functions
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_band()
    !-----------------------------------------------------------------------
    !! 
    !! Write bands
    !! 
    USE wvfct,         ONLY : nbnd, et
    USE constants,     ONLY : rytoev
    USE wannierEPW,    ONLY : ikstart, ikstop, iknum, num_bands, eigval, &
                              excluded_band
    !
    IMPLICIT NONE
    ! 
    ! Local variables
    INTEGER :: ik
    !! Counter on k-point
    INTEGER ikevc
    !! k-point index
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ibnd1
    !! Band index
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(eigval(num_bands, iknum), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_band', 'Error deallocating eigval', 1)
    !
    DO ik = ikstart, ikstop
      ikevc = ik - ikstart + 1
      ibnd1 = 0
      DO ibnd = 1, nbnd
        IF (excluded_band(ibnd)) CYCLE
        ibnd1 = ibnd1 + 1
        ! 
        ! RM - same value for rytoev as in wannier90
        ! eigval(ibnd1, ikevc) = et(ibnd, ik) * ryd2ev
        eigval(ibnd1, ikevc) = et(ibnd, ik) * rytoev
      ENDDO
    ENDDO
    !
    RETURN
    !
    !------------------------------------------------------------------------
    END SUBROUTINE write_band
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_plot()
    !-----------------------------------------------------------------------
    !!
    !! JN 06/2009: 
    !! added a couple of calls -- now works with multiple
    !! pools/procs (but one proc per pool)
    !!
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout, meta_ionode
    USE io_var,          ONLY : iun_plot
    USE wvfct,           ONLY : nbnd, npw, npwx, g2kin
    USE gvecw,           ONLY : gcutw
    USE wavefunctions,   ONLY : evc, psic, psic_nc
    USE wannierEPW,      ONLY : reduce_unk, wvfn_formatted, ispinw, nexband, &
                                excluded_band 
    USE klist,           ONLY : nks, igk_k
    USE klist_epw,       ONLY : xk_loc
    USE gvect,           ONLY : g, ngm 
    USE fft_base,        ONLY : dffts
    USE fft_interfaces,  ONLY : invfft
    USE noncollin_module,ONLY : noncolin
    USE scatter_mod,     ONLY : gather_grid
    USE constants_epw,   ONLY : czero, zero
    USE mp_global,       ONLY : my_pool_id
    USE kfold,           ONLY : ktokpmq
    USE io_epw,          ONLY : readwfc
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 20) wfnname
    !! Name of WF
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ibnd1
    !! Band index
    INTEGER :: ipool
    !! Index of current pool
    INTEGER ::  nkq
    !! Index of k-point in the pool
    INTEGER ::  nkq_abs
    !! Absolute index of k-point
    INTEGER ::  ik_g
    !! Temporary index of k-point, ik_g = nkq_abs
    INTEGER :: spin
    !! Spin case
    INTEGER :: ipol
    !! Counter on polarizations
    INTEGER :: istart, iend
    !! Starting/ ending index of plane waves
    INTEGER :: i, j, k
    !!
    INTEGER :: n1by2, n2by2, n3by2
    !!
    INTEGER :: idx
    !!
    INTEGER :: pos
    !!
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: zero_vect(3)
    !! Temporary zero vector
    !
    COMPLEX(KIND = DP), ALLOCATABLE :: psic_small(:)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: psic_nc_small(:, :)
    !!
#if defined(__MPI)
    INTEGER nxxs
    !! 
    COMPLEX(KIND = DP), ALLOCATABLE :: psic_all(:) 
    !! 
    COMPLEX(KIND = DP), ALLOCATABLE :: psic_nc_all(:, :) 
    !! 
    !
    nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
    IF (noncolin) THEN
      ALLOCATE(psic_nc_all(nxxs, 2), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_band', 'Error allocating psic_nc_all', 1)
      psic_nc_all(:, :) = czero
    ELSE
      ALLOCATE(psic_all(nxxs), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_band', 'Error allocating psic_all', 1)
      psic_all(:) = czero
    ENDIF
#endif
    !
    CALL start_clock('write_unk')
    ! 
    zero_vect(:) = zero
    WRITE(stdout,*) '    Writing out UNK plot files'
    !
    IF (reduce_unk) THEN
      WRITE(stdout, '(3(a, i5))') 'nr1s =', dffts%nr1, 'nr2s =', dffts%nr2, 'nr3s =', dffts%nr3
      n1by2 = (dffts%nr1 + 1) / 2
      n2by2 = (dffts%nr2 + 1) / 2
      n3by2 = (dffts%nr3 + 1) / 2
      WRITE(stdout, '(3(a, i5))') 'n1by2 =', n1by2, 'n2by2 =', n2by2, 'n3by2 =', n3by2
      IF (noncolin) THEN                                                               
        ALLOCATE(psic_nc_small(n1by2 * n2by2 * n3by2, 2), STAT = ierr)                                            
        IF (ierr /= 0) CALL errore('write_band', 'Error allocating psic_nc_small', 1)
        psic_nc_small(:, :) = czero                                                      
      ELSE                                                                                  
        ALLOCATE(psic_small(n1by2 * n2by2 * n3by2), STAT = ierr)      
        IF (ierr /= 0) CALL errore('write_band', 'Error allocating psic_small', 1)                              
        psic_small(:) = czero                                                   
      ENDIF 
    ENDIF
    !
    DO ik = 1, nks
      !
      ! returns in-pool index nkq and absolute index nkq_abs of xk
      CALL ktokpmq(xk_loc(:, ik), zero_vect, +1, ipool, nkq, nkq_abs)
      ik_g = nkq_abs
      !
      spin = ispinw
      IF (ispinw == 0) spin = 1
      IF (noncolin) THEN                                                               
        WRITE(wfnname, 201) ik_g                            
      ELSE                                                                                  
        WRITE(wfnname, 200) ik_g, spin                                                           
      ENDIF       
200 FORMAT('UNK', i5.5, '.', i1)
201 FORMAT('UNK', i5.5, '.', 'NC')                                                          
      !
      IF (meta_ionode) THEN 
        IF (wvfn_formatted) THEN 
          OPEN(UNIT = iun_plot, FILE = wfnname, FORM = 'formatted')
           IF (reduce_unk) THEN
             WRITE(iun_plot, *) n1by2, n2by2, n3by2, ik_g, nbnd - nexband
           ELSE
             WRITE(iun_plot, *) dffts%nr1, dffts%nr2, dffts%nr3, ik_g, nbnd - nexband
           ENDIF
        ELSE
           OPEN(UNIT = iun_plot, FILE = wfnname, FORM = 'unformatted')
           IF (reduce_unk) THEN
             WRITE(iun_plot) n1by2, n2by2, n3by2, ik_g, nbnd - nexband
           ELSE
             WRITE(iun_plot) dffts%nr1, dffts%nr2, dffts%nr3, ik_g, nbnd - nexband
           ENDIF
        ENDIF
      ENDIF
      !
      ! read wfc at ik
      CALL readwfc(my_pool_id + 1, ik, evc)
      !
      ! sorts k+G vectors in order of increasing magnitude, up to ecut
      CALL gk_sort(xk_loc(1, ik), ngm, g, gcutw, npw, igk_k(1, ik), g2kin)
      !
      ibnd1 = 0
      DO ibnd = 1, nbnd
        IF (excluded_band(ibnd)) CYCLE
        ibnd1 = ibnd1 + 1
        !
        IF (noncolin) THEN
          psic_nc(:, :) = czero
          DO ipol = 1, 2                                                              
            istart = 1 + npwx * (ipol - 1)
            iend   = npw + npwx * (ipol - 1) 
            psic_nc(dffts%nl(igk_k(1:npw, ik)), ipol) = evc(istart:iend, ibnd) 
            CALL invfft('Wave', psic_nc(:,ipol), dffts)                                 
          ENDDO 
        ELSE
          psic(:) = czero
          psic(dffts%nl(igk_k(1:npw, ik))) = evc(1:npw, ibnd)
          CALL invfft('Wave', psic, dffts)
        ENDIF
        IF (reduce_unk) pos = 0
#if defined(__MPI)
        IF (noncolin) THEN
          DO ipol = 1, 2
            CALL gather_grid(dffts, psic_nc(:, ipol), psic_nc_all(:, ipol))
          ENDDO
        ELSE
          CALL gather_grid(dffts, psic, psic_all)
        ENDIF
        !
        IF (reduce_unk) THEN
          DO k = 1, dffts%nr3, 2
            DO j = 1, dffts%nr2, 2
              DO i = 1, dffts%nr1, 2
                idx = (k - 1) * dffts%nr3 * dffts%nr2 + (j - 1) * dffts%nr2 + i
                pos = pos + 1
                IF (noncolin) THEN
                  DO ipol = 1, 2
                    psic_nc_small(pos, ipol) = psic_nc_all(idx, ipol)
                  ENDDO
                ELSE
                  psic_small(pos) = psic_all(idx)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        !
        IF (meta_ionode) THEN
          IF (wvfn_formatted) THEN
            IF (reduce_unk) THEN
              IF (noncolin) THEN
                DO ipol = 1, 2
                  WRITE(iun_plot, '(2ES20.10)') (psic_nc_small(j, ipol), j = 1, n1by2 * n2by2 * n3by2)
                ENDDO
              ELSE
                WRITE(iun_plot, '(2ES20.10)') (psic_small(j), j = 1, n1by2 * n2by2 * n3by2)
              ENDIF
            ELSE
              IF (noncolin) THEN
                DO ipol = 1, 2
                  WRITE(iun_plot, '(2ES20.10)') (psic_nc_all(j, ipol), j = 1, n1by2 * n2by2 * n3by2)
                ENDDO
              ELSE
                WRITE(iun_plot, '(2ES20.10)') (psic_all(j), j = 1, n1by2 * n2by2 * n3by2)
              ENDIF
            ENDIF
          ELSE
            IF (reduce_unk) THEN
              IF (noncolin) THEN            
                DO ipol = 1, 2
                  WRITE(iun_plot) (psic_nc_small(j, ipol), j = 1, n1by2 *n2by2 * n3by2)
                ENDDO
              ELSE
                WRITE(iun_plot) (psic_small(j), j = 1, n1by2 * n2by2 * n3by2)
              ENDIF
            ELSE
              IF (noncolin) THEN
                DO ipol = 1, 2
                  WRITE(iun_plot) (psic_nc_all(j, ipol), j = 1, n1by2 *n2by2 * n3by2)
                ENDDO
              ELSE
                WRITE(iun_plot) (psic_all(j), j = 1, n1by2 * n2by2 * n3by2)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
#else
        IF (reduce_unk) THEN
          DO k = 1, dffts%nr3, 2
            DO j = 1, dffts%nr2, 2
              DO i = 1, dffts%nr1, 2
                idx = (k-1) * dffts%nr3 * dffts%nr2 + (j-1) * dffts%nr2 + i
                pos = pos + 1
                IF (noncolin) THEN
                  DO ipol = 1, 2
                    psic_nc_small(pos, ipol) = psic_nc(idx, ipol)
                  ENDDO
                ELSE
                  psic_small(pos) = psic(idx)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        !
        IF (meta_ionode) THEN
          IF (wvfn_formatted) THEN
            IF (noncolin) THEN
              DO ipol = 1, 2
                IF (reduce_unk) THEN
                  WRITE(iun_plot, '(2ES20.10)') (psic_nc_small(j, ipol), j = 1, n1by2 * n2by2 * n3by2)
                ELSE
                  WRITE(iun_plot, '(2ES20.10)') (psic_nc(j, ipol), j = 1, dffts%nr1 * dffts%nr2 * dffts%nr3)
                ENDIF
              ENDDO
            ELSE
              IF (reduce_unk) THEN
                WRITE(iun_plot, '(2ES20.10)') (psic_small(j), j = 1, n1by2 * n2by2 * n3by2)
              ELSE
                WRITE(iun_plot, '(2ES20.10)') (psic(j), j = 1, dffts%nr1 * dffts%nr2 * dffts%nr3)
              ENDIF
            ENDIF
          ELSE
            IF (noncolin) THEN
              DO ipol = 1, 2
                IF (reduce_unk) THEN
                  WRITE(iun_plot) (psic_nc_small(j, ipol), j = 1, n1by2 * n2by2 * n3by2)
                ELSE
                  WRITE(iun_plot) (psic_nc(j, ipol), j = 1, dffts%nr1 * dffts%nr2 * dffts%nr3)
                ENDIF
              ENDDO
            ELSE
              IF (reduce_unk) THEN
                WRITE(iun_plot) (psic_small(j), j = 1, n1by2 * n2by2 * n3by2)
              ELSE
                WRITE(iun_plot) (psic(j), j = 1, dffts%nr1 * dffts%nr2 * dffts%nr3)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
#endif
      ENDDO ! ibnd
      !
      IF (meta_ionode) CLOSE(iun_plot)
      !
    ENDDO ! ik
    !
    IF (reduce_unk) THEN
      IF (noncolin) THEN
        DEALLOCATE(psic_nc_small, STAT = ierr)
        IF (ierr /= 0) CALL errore('write_band', 'Error deallocating psic_nc_small', 1)
      ELSE
        DEALLOCATE(psic_small, STAT = ierr)
        IF (ierr /= 0) CALL errore('write_band', 'Error deallocating psic_small', 1)
      ENDIF
    ENDIF
    !
#if defined(__MPI)
    IF (noncolin) THEN
      DEALLOCATE(psic_nc_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_band', 'Error deallocating psic_nc_all', 1)
    ELSE
      DEALLOCATE(psic_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_band', 'Error deallocating psic_all', 1)
    ENDIF
#endif
    !
    WRITE(stdout, '(/)')
    WRITE(stdout,*) ' UNK written'
    !
    CALL stop_clock('write_unk')
    !
    RETURN
    !
    !------------------------------------------------------------------------
    END SUBROUTINE write_plot
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE ylm_expansion() 
    !-----------------------------------------------------------------------
    !!
    !! Expansion 
    !!  
    USE kinds,         ONLY : DP
    USE random_numbers,ONLY : randy
    USE wannierEPW,       ONLY : n_proj, xaxis, zaxis, csph, l_w, mr_w
    USE matrix_inversion, ONLY: invmat
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    ! local variables
    INTEGER, PARAMETER :: lmax2 = 16
    !!
    INTEGER ::  i
    !! Couter on polarizations
    INTEGER ::  ir
    !! Counter on spherical harmonics
    INTEGER :: iw
    !! Counter on wannier projections
    INTEGER :: ierr
    !! Error status
    ! 
    REAL(KIND = DP) :: u(3, 3)
    !!
    REAL(KIND = DP), ALLOCATABLE :: r(:, :)
    !!
    REAL(KIND = DP), ALLOCATABLE :: rr(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: rp(:, :)
    !!
    REAL(KIND = DP), ALLOCATABLE :: ylm_w(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: ylm(:, :)
    !!
    REAL(KIND = DP), ALLOCATABLE :: mly(:, :)
    !!
    !
    ALLOCATE(r(3, lmax2), STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error allocating r', 1)
    ALLOCATE(rp(3, lmax2), STAT = ierr)  
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error allocating rp', 1)
    ALLOCATE(rr(lmax2), STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error allocating rr', 1)
    ALLOCATE(ylm_w(lmax2), STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error allocating ylm_w', 1)
    ALLOCATE(ylm(lmax2, lmax2), STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error allocating ylm', 1)
    ALLOCATE(mly(lmax2, lmax2), STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error allocating mly', 1)
    r(:, :) = zero
    rr(:) = zero
    rp(:, :) = zero
    ylm_w(:) = zero
    ylm(:, :) = zero
    mly(:, :) = zero
    !
    ! generate a set of nr=lmax2 random vectors
    DO ir = 1, lmax2
      DO i = 1, 3
        r(i, ir) = randy() - 0.5d0
      ENDDO
    ENDDO
    rr(:) = r(1, :) * r(1, :) + r(2, :) * r(2, :) + r(3, :) * r(3, :)
    !- compute ylm(ir,lm)
    CALL ylmr2(lmax2, lmax2, r, rr, ylm)
    !- store the inverse of ylm(ir,lm) in mly(lm,ir)
    CALL invmat(lmax2, ylm, mly)
    !- check that r points are independent
    CALL check_inverse(lmax2, ylm, mly)
    !
    DO iw = 1, n_proj 
      !
      !- define the u matrix that rotate the reference frame
      CALL set_u_matrix(xaxis(:, iw), zaxis(:, iw), u)
      !- find rotated r-vectors 
      rp(:, :) = MATMUL(u(:, :), r(:, :))
      !- set ylm funtion according to wannier90 (l,mr) indexing in the rotaterd points
      CALL ylm_wannier(ylm_w, l_w(iw), mr_w(iw), rp, lmax2) 
      !
      csph(:, iw) = MATMUL(mly(:, :), ylm_w(:))
      !
      !
    ENDDO
    !
    DEALLOCATE(r, STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error deallocating r', 1)
    DEALLOCATE(rp, STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error deallocating rp', 1)
    DEALLOCATE(rr, STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error deallocating rr', 1)
    DEALLOCATE(ylm_w, STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error deallocating ylm_w', 1)
    DEALLOCATE(ylm, STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error deallocating ylm', 1)
    DEALLOCATE(mly, STAT = ierr)
    IF (ierr /= 0) CALL errore('ylm_expansion', 'Error deallocating mly', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE ylm_expansion
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE check_inverse(lmax2, ylm, mly)
    !-----------------------------------------------------------------------
    !! 
    !! Check the inverse 
    !!
    USE kinds,     ONLY : DP
    USE constants_epw, ONLY : zero, eps8
    !
    IMPLICIT NONE
    !
    ! I/O variables
    INTEGER, INTENT(in) :: lmax2
    !!
    REAL(KIND = DP), INTENT(in) :: ylm(lmax2, lmax2)
    !!
    REAL(KIND = DP), INTENT(in) :: mly(lmax2, lmax2)
    !!
    !
    ! local variables
    INTEGER :: lm
    !! Counter on spherical harmonics
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: capel
    !!
    REAL(KIND = DP), ALLOCATABLE :: uno(:, :)
    !!
    !
    ALLOCATE(uno(lmax2, lmax2), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_inverse', 'Error allocating uno', 1)
    uno(:, :) = zero
    !
    uno = MATMUL(mly, ylm)
    capel = zero
    DO lm = 1, lmax2
      uno(lm, lm) = uno(lm, lm) - 1.d0
    ENDDO
    capel = capel + SUM(ABS(uno(1:lmax2, 1:lmax2)))
    IF (capel > eps8) CALL errore('ylm_expansion', &
                     'Inversion failed: r(*, 1:nr) are not all independent!!', 1)
    DEALLOCATE(uno, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_inverse', 'Error deallocating uno', 1)
    !
    RETURN
    !
    !------------------------------------------------------------------------
    END SUBROUTINE check_inverse
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE set_u_matrix(x, z, u)
    !-----------------------------------------------------------------------
    !!
    !! Set the U matrix
    !!  
    USE kinds,     ONLY : DP
    USE constants_epw, ONLY : eps8
    !
    IMPLICIT NONE
    !
    ! I/O variables
    REAL(KIND = DP), INTENT(in) :: x(3)
    !!
    REAL(KIND = DP), INTENT(in) :: z(3)
    !!
    REAL(KIND = DP), INTENT(out) :: u(3, 3)
    !!
    ! local variables
    REAL(KIND = DP) :: coseno
    !! cosine between x(3) and z(3)
    REAL(KIND = DP) :: xx
    !! Norm of x(3)
    REAL(KIND = DP) :: zz
    !! Norm of z(3)
    REAL(KIND = DP) :: y(3)
    !!
    !
    xx = SQRT(SUM(x(:) * x(:)))
    IF (xx < eps8) CALL errore('set_u_matrix', '|xaxis| < eps8', 1)
    !
    zz = SQRT(SUM(z(:) * z(:)))
    IF (zz < eps8) CALL errore ('set_u_matrix', '|zaxis| < eps8', 1)
    !
    coseno = SUM(x(:) * z(:)) / xx / zz
    IF (ABS(coseno) > eps8) CALL errore('set_u_matrix', 'xaxis and zaxis are not orthogonal!', 1)
    !
    y(1) = (z(2) * x(3) - x(2) * z(3)) / xx / zz
    y(2) = (z(3) * x(1) - x(3) * z(1)) / xx / zz
    y(3) = (z(1) * x(2) - x(1) * z(2)) / xx / zz
    !
    u(1, :) = x(:) / xx
    u(2, :) = y(:)
    u(3, :) = z(:) / zz
    !
    !
    RETURN
    !
    !------------------------------------------------------------------------
    END SUBROUTINE set_u_matrix
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE ylm_wannier(ylm, l, mr, r, nr) 
    !-----------------------------------------------------------------------
    !!
    !! This routine returns in ylm(r) the values at the nr points r(1:3,1:nr) 
    !! of the spherical harmonic identified  by indices (l,mr) 
    !! in table 3.1 of the wannierf90 specification.
    !! 
    !! No reference to the particular ylm ordering internal to quantum-espresso
    !! is assumed. 
    !!
    !! If ordering in wannier90 code is changed or extended this should be the 
    !! only place to be modified accordingly
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi, pibytwo, eps8
    USE low_lvl,       ONLY : fy3x2my2, fxx2m3y2, fxx2m3y2, fxyz, fzx2my2, &
                              fyz2, fxz2, fz3, dxy, dx2my2, dyz, dxz, dz2, &
                              p_z, py, px, s 
    !
    IMPLICIT NONE
    !
    ! I/O variables
    INTEGER, INTENT(in) :: l, mr
    !! quantum numbers
    INTEGER, INTENT(in) :: nr
    !!
    !
    REAL(KIND = DP), INTENT(in) :: r(3, nr)
    !!
    REAL(KIND = DP), INTENT(out) :: ylm(nr)
    !! spherical harmonics
    !
    ! local variables
    INTEGER :: ir
    !!
    REAL(KIND = DP) :: rr
    !!
    REAL(KIND = DP) :: cost
    !!
    REAL(KIND = DP) :: phi
    !!
    REAL(KIND = DP) :: bs2, bs3, bs6, bs12
    !! Inverse SQRT of 2, 3, 6, and 12
    !
    bs2  = 1.d0 / SQRT(2.d0)
    bs3  = 1.d0 / SQRT(3.d0)
    bs6  = 1.d0 / SQRT(6.d0)
    bs12 = 1.d0 / SQRT(12.d0)
    !
    IF (l > 3 .OR. l < -5) CALL errore('ylm_wannier', 'l out of range ', 1)
    IF (l >= 0) THEN
      IF (mr < 1 .OR. mr > 2 * l + 1) CALL errore('ylm_wannier', 'mr out of range' ,1)
    ELSE
      IF (mr < 1 .OR. mr > ABS(l) + 1 ) CALL errore('ylm_wannier', 'mr out of range', 1)
    ENDIF
    !
    DO ir  = 1, nr
      rr = SQRT(SUM(r(:, ir) * r(:, ir)))
      IF (rr < eps8) CALL errore('ylm_wannier', 'rr too small', 1)
      !
      cost =  r(3, ir) / rr
      !
      !  beware the arc tan, it is defined modulo pi
      !
      IF (r(1, ir) > eps8) THEN
        phi = ATAN(r(2, ir) / r(1, ir))
      ELSEIF (r(1, ir) < -eps8) THEN
        phi = ATAN(r(2, ir) / r(1, ir)) + pi
      ELSE
        phi = SIGN(pibytwo, r(2, ir))
      ENDIF
      !
      IF (l == 0) THEN ! s orbital
        ylm(ir) = s()  
      ENDIF
      IF (l == 1) THEN   ! p orbitals
        IF (mr == 1) ylm(ir) = p_z(cost) 
        IF (mr == 2) ylm(ir) = px(cost, phi)
        IF (mr == 3) ylm(ir) = py(cost, phi)
      ENDIF
      IF (l == 2) THEN   ! d orbitals
        IF (mr == 1) ylm(ir) = dz2(cost)
        IF (mr == 2) ylm(ir) = dxz(cost, phi)
        IF (mr == 3) ylm(ir) = dyz(cost, phi)
        IF (mr == 4) ylm(ir) = dx2my2(cost, phi)
        IF (mr == 5) ylm(ir) = dxy(cost, phi)
      ENDIF
      IF (l == 3) THEN   ! f orbitals
        IF (mr == 1) ylm(ir) = fz3(cost)
        IF (mr == 2) ylm(ir) = fxz2(cost, phi)
        IF (mr == 3) ylm(ir) = fyz2(cost, phi)
        IF (mr == 4) ylm(ir) = fzx2my2(cost, phi)
        IF (mr == 5) ylm(ir) = fxyz(cost, phi)
        IF (mr == 6) ylm(ir) = fxx2m3y2(cost, phi)
        IF (mr == 7) ylm(ir) = fy3x2my2(cost, phi)
      ENDIF
      IF (l == -1) THEN  !  sp hybrids
        IF (mr == 1) ylm(ir) = bs2 * (s() + px(cost, phi)) 
        IF (mr == 2) ylm(ir) = bs2 * (s() - px(cost, phi)) 
      ENDIF
      IF (l == -2) THEN  !  sp2 hybrids 
        IF (mr == 1) ylm(ir) = bs3 * s() - bs6 * px(cost, phi) + bs2 * py(cost, phi)
        IF (mr == 2) ylm(ir) = bs3 * s() - bs6 * px(cost, phi) - bs2 * py(cost, phi)
        IF (mr == 3) ylm(ir) = bs3 * s() + 2.d0 *bs6 * px(cost, phi) 
      ENDIF
      IF (l == -3) THEN  !  sp3 hybrids
        IF (mr == 1) ylm(ir) = 0.5d0 *(s() + px(cost, phi) + py(cost, phi) + p_z(cost))
        IF (mr == 2) ylm(ir) = 0.5d0 *(s() + px(cost, phi) - py(cost, phi) - p_z(cost))
        IF (mr == 3) ylm(ir) = 0.5d0 *(s() - px(cost, phi) + py(cost, phi) - p_z(cost))
        IF (mr == 4) ylm(ir) = 0.5d0 *(s() - px(cost, phi) - py(cost, phi) + p_z(cost))
      ENDIF
      IF (l == -4) THEN  !  sp3d hybrids
        IF (mr == 1) ylm(ir) = bs3 * s() -bs6 * px(cost, phi) + bs2 * py(cost, phi)
        IF (mr == 2) ylm(ir) = bs3 * s() -bs6 * px(cost, phi) - bs2 * py(cost, phi)
        IF (mr == 3) ylm(ir) = bs3 * s() +2.d0 * bs6 * px(cost, phi) 
        IF (mr == 4) ylm(ir) = bs2 * p_z(cost) + bs2 * dz2(cost)
        IF (mr == 5) ylm(ir) =-bs2 * p_z(cost) + bs2 * dz2(cost)
      ENDIF
      IF (l == -5) THEN  ! sp3d2 hybrids
        IF (mr == 1) ylm(ir) = bs6 *s() - bs2 * px(cost, phi) - bs12 * dz2(cost) + .5 * dx2my2(cost, phi)
        IF (mr == 2) ylm(ir) = bs6 *s() + bs2 * px(cost, phi) - bs12 * dz2(cost) + .5 * dx2my2(cost, phi)
        IF (mr == 3) ylm(ir) = bs6 *s() - bs2 * py(cost, phi) - bs12 * dz2(cost) - .5 * dx2my2(cost, phi)
        IF (mr == 4) ylm(ir) = bs6 *s() + bs2 * py(cost, phi) - bs12 * dz2(cost) - .5 * dx2my2(cost, phi)
        IF (mr == 5) ylm(ir) = bs6 *s() - bs2 * p_z(cost) + bs3 * dz2(cost)
        IF (mr == 6) ylm(ir) = bs6 *s() + bs2 * p_z(cost) + bs3 * dz2(cost)
      ENDIF
      !
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE ylm_wannier
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE radialpart(ng, q, alfa, rvalue, lmax, radial)
    !-----------------------------------------------------------------------
    !!
    !! This routine computes a table with the radial Fourier transform 
    !! of the radial functions.
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : fpi
    USE cell_base, ONLY : omega
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    ! I/O
    INTEGER, INTENT(in) :: ng
    !! Number of plane waves
    INTEGER, INTENT(in) :: rvalue
    !!
    INTEGER, INTENT(in) :: lmax
    !! Maxium angular momentum
    !
    REAL(KIND = DP), INTENT(in) :: alfa
    !!
    REAL(KIND = DP), INTENT(in) :: q(ng)
    !!
    REAL(KIND = DP), INTENT(out) :: radial(ng, 0:lmax)
    !!
    !
    ! local variables
    INTEGER :: l
    !! Counter on angular momentum
    INTEGER :: ig
    !! Counter on plane waves
    INTEGER :: ir
    !! Counter on mesh_r
    INTEGER :: mesh_r
    !! Size of the radial FFT mesh
    INTEGER :: ierr
    !! Error status
    ! 
    REAL(KIND = DP), PARAMETER :: xmin = -6.d0
    !!
    REAL(KIND = DP), PARAMETER :: dx = 0.025d0
    !!
    REAL(KIND = DP), PARAMETER :: rmax = 10.d0
    !!
    REAL(KIND = DP) :: rad_int
    !!
    REAL(KIND = DP) :: pref
    !!
    REAL(KIND = DP) :: x
    !!
    REAL(KIND = DP), ALLOCATABLE :: bes(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: func_r(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: r(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: rij(:)
    !!
    REAL(KIND = DP), ALLOCATABLE :: aux(:)
    !! 
    ! 
    mesh_r = NINT((log(rmax) - xmin) / dx + 1)
    !
    ALLOCATE(bes(mesh_r), STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error allocating bes', 1)
    ALLOCATE(func_r(mesh_r), STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error allocating func_r', 1)
    ALLOCATE(r(mesh_r), STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error allocating r', 1)
    ALLOCATE(rij(mesh_r), STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error allocating rij', 1)
    ALLOCATE(aux(mesh_r), STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error allocating aux', 1)
    bes(:) = zero
    func_r(:) = zero
    r(:) = zero
    rij(:) = zero
    aux(:) = zero
    !
    ! compute the radial mesh
    !
    DO ir = 1, mesh_r
      x = xmin  + DBLE(ir - 1) * dx 
      r(ir) = EXP(x) / alfa
      rij(ir) = dx  * r(ir)
    ENDDO
    !
    IF (rvalue == 1) func_r(:) = 2.d0 * alfa**(3.d0 / 2.d0) * EXP(-alfa * r(:))
    IF (rvalue == 2) func_r(:) = 1.d0 / SQRT(8.d0) * alfa**(3.d0 / 2.d0) * & 
                                (2.0d0 - alfa * r(:)) * EXP(-alfa * r(:) * 0.5d0)
    IF (rvalue == 3) func_r(:) = SQRT(4.d0 / 27.d0) * alfa**(2.0d0 / 3.0d0) * &
                                (1.d0 - 1.5d0 * alfa * r(:) + & 
                                2.d0 * (alfa * r(:))**2 / 27.d0) * & 
                                EXP(-alfa * r(:) / 3.0d0)
    pref = fpi / SQRT(omega)
    !
    DO l = 0, lmax
      DO ig = 1, ng
        CALL sph_bes(mesh_r, r(1), q(ig), l, bes)
        aux(:) = bes(:) * func_r(:) * r(:) * r(:)
        CALL simpson(mesh_r, aux, rij, rad_int)
        radial(ig, l) = rad_int * pref
      ENDDO
    ENDDO
    !
    DEALLOCATE(bes, STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error deallocating bes', 1)
    DEALLOCATE(func_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error deallocating func_r', 1)
    DEALLOCATE(r, STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error deallocating r', 1)
    DEALLOCATE(rij, STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error deallocating rij', 1)
    DEALLOCATE(aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('radialpart', 'Error deallocating aux', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE radialpart
    !-----------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE pw2wan2epw
  !------------------------------------------------------------------------------
