  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  SUBROUTINE ephwann_shuffle_mem( nqc, xqc )
  !---------------------------------------------------------------------
  !!
  !!  Wannier interpolation of electron-phonon vertex
  !!  Here we do an additional mode parallelism to save memory
  !!  (only fine interpolation part)
  !!
  !!  Scalar implementation   Feb 2006
  !!  Parallel version        May 2006
  !!  Disentenglement         Oct 2006
  !!  Compact formalism       Dec 2006
  !!  Phonon irreducible zone Mar 2007
  !!
  !!  No ultrasoft now
  !!  No spin polarization
  !!
  !!  RM - add noncolin case
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : nbnd, nks, nkstot, isk, et, xk, ef, nelec
  USE cell_base,     ONLY : at, bg, omega, alat
  USE start_k,       ONLY : nk1, nk2, nk3
  USE ions_base,     ONLY : nat, amass, ityp, tau
  USE phcom,         ONLY : nq1, nq2, nq3, nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange,               &
                            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen, &
                            nbndskip, scr_typ,                                  &
                            elecselfen, phonselfen, nest_fn, a2f, specfun_ph,   &
                            vme, eig_read, ephwrite, nkf1, nkf2, nkf3,          & 
                            efermi_read, fermi_energy, specfun_el, band_plot,   &
                            scattering, nstemp, int_mob, scissor, carrier,      &
                            iterative_bte, longrange, scatread, nqf1, prtgkk,   &
                            nqf2, nqf3, mp_mesh_k, restart, ncarrier, plselfen, &
                            specfun_pl, use_ws, epmatkqread, selecqread
  USE noncollin_module, ONLY : noncolin
  USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, zero, czero, cone,       & 
                            twopi, ci, kelvin2eV, eps8
  USE io_files,      ONLY : prefix, diropn, tmp_dir
  USE io_global,     ONLY : stdout, ionode
  USE io_epw,        ONLY : lambda_phself, linewidth_phself, iunepmatwe,        &
                            iunepmatwp, crystal, iunepmatwp2, iunrestart
  USE elph2,         ONLY : cu, cuq, lwin, lwinq,                               &
                            chw, chw_ks, cvmew, cdmew, rdw,                     &
                            epmatq, wf, etf, etf_k, etf_ks, xqf, xkf,           &
                            wkf, dynq, nkqf, epf17, nkf, nqf, et_ks,    &
                            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef,     &
                            sigmai_all, sigmai_mode, gamma_all, epsi, zstar,    &
                            efnew, sigmar_all, zi_all, nkqtotf, eps_rpa,        &
                            nkqtotf, sigmar_all, zi_allvb, inv_tau_all,         &
                            inv_tau_allcb, zi_allcb, exband
  USE wigner,        ONLY : wigner_seitz_wrap
  USE transportcom,  ONLY : transp_temp, mobilityh_save, mobilityel_save, lower_bnd, &
                            upper_bnd 
  USE wan2bloch,     ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch,   &
                            vmewan2bloch, dynifc2blochf, dynifc2blochc, & 
                            ephwan2blochp_mem, ephwan2bloch_mem
  USE bloch2wan,     ONLY : hambloch2wan, dmebloch2wan, dynbloch2wan,   &
                            vmebloch2wan, ephbloch2wane, ephbloch2wanp_mem
  USE io_eliashberg, ONLY : write_ephmat, count_kpoints, kmesh_fine,    &
                            kqmap_fine
  USE transport,     ONLY : transport_coeffs, scattering_rate_q, qwindow
  USE transport_iter,ONLY : iter_restart
  USE printing,      ONLY : print_gkk
  USE io_scattering, ONLY : electron_read, tau_read, iter_open
  USE close_epw,     ONLY : iter_close
  USE division,      ONLY : fkbounds
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool
  USE mp_world,      ONLY : mpime, world_comm
#if defined(__MPI)
  USE parallel_include, ONLY : MPI_MODE_RDONLY, MPI_INFO_NULL, MPI_OFFSET_KIND, &
                               MPI_OFFSET
#endif
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nqc
  !! number of qpoints in the coarse grid
  !
  REAL(kind=DP), INTENT (in) :: xqc(3,nqc)
  !! qpoint list, coarse mesh
  ! 
  ! Local  variables
  LOGICAL :: already_skipped
  !! Skipping band during the Wannierization
  LOGICAL :: exst
  !! If the file exist
  LOGICAL :: first_cycle
  !! Check wheter this is the first cycle after a restart. 
  LOGICAL :: first_time
  !! Check wheter this is the first timeafter a restart. 
  LOGICAL :: homogeneous
  !! Check if the k and q grids are homogenous and commensurate.
  !
  CHARACTER (len=256) :: filint
  !! Name of the file to write/read 
  CHARACTER (len=30)  :: myfmt
  !! Variable used for formatting output
  ! 
  INTEGER :: ios
  !! integer variable for I/O control
  INTEGER :: iq 
  !! Counter on coarse q-point grid
  INTEGER :: iqq 
  !! Counter on coarse q-point grid
  INTEGER :: iq_restart
  !! Counter on coarse q-point grid
  INTEGER :: ik
  !! Counter on coarse k-point grid
  INTEGER :: ikk
  !! Counter on k-point when you have paired k and q
  INTEGER :: ikq
  !! Paired counter so that q is adjacent to its k
  INTEGER :: ibnd
  !! Counter on band
  INTEGER :: jbnd
  !! Counter on band
  INTEGER :: imode
  !! Counter on mode
  INTEGER :: na
  !! Counter on atom
  INTEGER :: mu
  !! counter on mode
  INTEGER :: nu
  !! counter on mode
  INTEGER :: fermicount
  !! Number of states at the Fermi level
  INTEGER :: lrepmatw
  !! record length while reading file
  INTEGER :: ikx
  !! Counter on the coase k-grid
  INTEGER :: ikfx 
  !! Counter on the fine k-grid. 
  INTEGER :: xkk1, xkq1
  !! Integer of xkk when multiplied by nkf/nk
  INTEGER :: xkk2, xkq2
  !! Integer of xkk when multiplied by nkf/nk
  INTEGER :: xkk3, xkq3
  !! Integer of xkk when multiplied by nkf/nk
  INTEGER :: ir
  !! Counter for WS loop
  INTEGER :: nrws
  !! Number of real-space Wigner-Seitz
  INTEGER :: valueRSS(2)
  !! Return virtual and resisdent memory from system
  INTEGER :: ierr
  !! Error status
  INTEGER :: nrr_k
  !! Number of WS points for electrons
  INTEGER :: nrr_q
  !! Number of WS points for phonons
  INTEGER :: nrr_g
  !! Number of WS points for electron-phonons
  INTEGER :: dims
  !! Dims is either nbndsub if use_ws or 1 if not
  INTEGER :: dims2
  !! Dims is either nat if use_ws or 1 if not
  INTEGER :: iw
  !! Counter on bands when use_ws == .true.
  INTEGER :: iw2
  !! Counter on bands when use_ws == .true.
  INTEGER :: iter
  !! Current iteration number
  INTEGER :: itemp
  !! Temperature index
  INTEGER :: icbm
  !! Index of the CBM
  INTEGER :: totq
  !! Total number of q-points within the fsthick window. 
  INTEGER, ALLOCATABLE :: irvec_k(:,:)
  !! integer components of the ir-th Wigner-Seitz grid point in the basis
  !! of the lattice vectors for electrons
  INTEGER, ALLOCATABLE :: irvec_q(:,:)
  !! integer components of the ir-th Wigner-Seitz grid point for phonons
  INTEGER, ALLOCATABLE :: irvec_g(:,:)
  !! integer components of the ir-th Wigner-Seitz grid point for electron-phonon
  INTEGER, ALLOCATABLE :: ndegen_k (:,:,:)
  !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
  INTEGER, ALLOCATABLE :: ndegen_q (:,:,:)
  !! Wigner-Seitz weights for the phonon grid that depend on 
  !! atomic positions $R + \tau(nb) - \tau(na)$
  INTEGER, ALLOCATABLE :: ndegen_g (:,:,:,:)
  !! Wigner-Seitz weights for the electron-phonon grid that depend on 
  !! atomic positions $R - \tau(na)$
  INTEGER, ALLOCATABLE :: selecq(:)
  !! Selected q-points within the fsthick window
  INTEGER, PARAMETER :: nrwsx=200
  !! Maximum number of real-space Wigner-Seitz
#if defined(__MPI)
  INTEGER (kind=MPI_OFFSET_KIND) :: ind_tot
  INTEGER (kind=MPI_OFFSET_KIND) :: ind_totcb
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw2
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw4
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw5
  INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw6
  !! Offset to tell where to start reading the file
  INTEGER (kind=MPI_OFFSET_KIND) :: lsize
  !! Offset to tell where to start reading the file
#else
  INTEGER :: ind_tot
  INTEGER :: ind_totcb
  INTEGER :: lrepmatw2
  INTEGER :: lrepmatw4
  INTEGER :: lrepmatw5
  INTEGER :: lrepmatw6
  !! Offset to tell where to start reading the file
  INTEGER :: lsize
#endif
  !  
  REAL(kind=DP) :: rdotk_scal
  !! Real (instead of array) for $r\cdot k$
  REAL(kind=DP) :: xxq(3)
  !! Current q-point 
  REAL(kind=DP) :: xxk(3)
  !! Current k-point on the fine grid
  REAL(kind=DP) :: xkk(3)
  !! Current k-point on the fine grid
  REAL(kind=DP) :: xkq(3)
  !! Current k+q point on the fine grid
  REAL(kind=DP) :: rws(0:3,nrwsx)
  !! Real-space wigner-Seitz vectors
  REAL(kind=DP) :: atws(3,3)
  !! Maximum vector: at*nq
  REAL(kind=DP) :: w_centers(3,nbndsub)
  !! Wannier centers
  REAL(KIND=DP) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND=DP) :: ef0(nstemp)
  !! Fermi level for the temperature itemp  
  REAL(KIND=DP) :: efcb(nstemp)
  !! Second Fermi level for the temperature itemp  
  REAL(KIND=DP) :: dummy(3)
  !! Dummy variable
  REAL(KIND=DP), EXTERNAL :: fermicarrier
  !! Function that returns the Fermi level so that n=p (if int_mob = .true.)
  REAL(kind=DP), EXTERNAL :: efermig
  !! External function to calculate the fermi energy
  REAL(kind=DP), EXTERNAL :: efermig_seq
  !! Same but in sequential
  REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
  REAL(kind=DP), ALLOCATABLE :: w2 (:)
  !! Interpolated phonon frequency
  REAL(kind=DP), ALLOCATABLE :: irvec_r (:,:)
  !! Wigner-Size supercell vectors, store in real instead of integer
  REAL(kind=DP), ALLOCATABLE :: rdotk(:)
  !! $r\cdot k$
  REAL(kind=DP), ALLOCATABLE :: rdotk2(:)
  !! $r\cdot k$
  REAL(kind=DP), ALLOCATABLE :: wslen_k(:)
  !! real-space length for electrons, in units of alat
  REAL(kind=DP), ALLOCATABLE :: wslen_q(:)
  !! real-space length for phonons, in units of alat
  REAL(kind=DP), ALLOCATABLE :: wslen_g(:)
  !! real-space length for electron-phonons, in units of alat
  REAL(kind=DP), ALLOCATABLE :: vkk_all(:,:,:)
  !! velocity from all the k-point
  REAL(kind=DP), ALLOCATABLE :: wkf_all(:)
  !! k-point weights for all the k-points
  !
  COMPLEX(kind=DP), ALLOCATABLE :: epmatwe_mem  (:,:,:,:)
  !! e-p matrix  in wannier basis - electrons (written on disk)
  COMPLEX(kind=DP), ALLOCATABLE :: epmatwef (:,:,:)
  !! e-p matrix  in el wannier - fine Bloch phonon grid
  COMPLEX(kind=DP), ALLOCATABLE :: epmatf( :, :)
  !! e-p matrix  in smooth Bloch basis, fine mesh
  COMPLEX(kind=DP), ALLOCATABLE :: cufkk ( :, :)
  !! Rotation matrix, fine mesh, points k
  COMPLEX(kind=DP), ALLOCATABLE :: cufkq ( :, :)
  !! the same, for points k+q
  COMPLEX(kind=DP), ALLOCATABLE :: uf( :, :)
  !! Rotation matrix for phonons
  COMPLEX(kind=DP), ALLOCATABLE :: bmatf ( :, :)
  !! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
  COMPLEX(kind=DP), ALLOCATABLE :: cfac(:,:,:)
  !! Used to store $e^{2\pi r \cdot k}$ exponential 
  COMPLEX(kind=DP), ALLOCATABLE :: cfacq(:,:,:)
  !! Used to store $e^{2\pi r \cdot k+q}$ exponential
  COMPLEX(kind=DP), ALLOCATABLE :: eptmp(:,:,:,:)
  !! Temporary el-ph matrices. 
  COMPLEX(kind=DP), ALLOCATABLE :: epmatlrT(:,:,:,:)
  !! Long-range temp. save
  ! 
  IF (nbndsub.ne.nbnd) &
       WRITE(stdout, '(/,5x,a,i4)' ) 'Band disentanglement is used:  nbndsub = ', nbndsub
  !
  ALLOCATE ( cu ( nbnd, nbndsub, nks), & 
             cuq ( nbnd, nbndsub, nks), & 
             lwin ( nbnd, nks ), &
             lwinq ( nbnd, nks ), &
             exband ( nbnd ) )
  !
  CALL start_clock ( 'ephwann' )
  !
  IF ( epwread ) THEN
    !
    ! Might have been pre-allocate depending of the restart configuration 
    IF(ALLOCATED(tau))  DEALLOCATE( tau )
    IF(ALLOCATED(ityp)) DEALLOCATE( ityp )
    IF(ALLOCATED(w2))   DEALLOCATE( w2 )
    ! 
    ! We need some crystal info
    IF (mpime.eq.ionode_id) THEN
      !
      OPEN(unit=crystal,file='crystal.fmt',status='old',iostat=ios)
      READ (crystal,*) nat
      READ (crystal,*) nmodes
      READ (crystal,*) nelec
      READ (crystal,*) at
      READ (crystal,*) bg
      READ (crystal,*) omega
      READ (crystal,*) alat
      ALLOCATE( tau( 3, nat ) )
      READ (crystal,*) tau
      READ (crystal,*) amass
      ALLOCATE( ityp( nat ) )
      READ (crystal,*) ityp
      READ (crystal,*) isk
      READ (crystal,*) noncolin
      READ (crystal,*) w_centers
      ! 
    ENDIF
    CALL mp_bcast (nat     , ionode_id, inter_pool_comm)
    CALL mp_bcast (nat     , root_pool, intra_pool_comm)  
    IF (mpime /= ionode_id) ALLOCATE( ityp( nat ) )
    CALL mp_bcast (nmodes  , ionode_id, inter_pool_comm)
    CALL mp_bcast (nmodes  , root_pool, intra_pool_comm)  
    CALL mp_bcast (nelec   , ionode_id, inter_pool_comm)
    CALL mp_bcast (nelec   , root_pool, intra_pool_comm)  
    CALL mp_bcast (at      , ionode_id, inter_pool_comm)
    CALL mp_bcast (at      , root_pool, intra_pool_comm)  
    CALL mp_bcast (bg      , ionode_id, inter_pool_comm)
    CALL mp_bcast (bg      , root_pool, intra_pool_comm)  
    CALL mp_bcast (omega   , ionode_id, inter_pool_comm)
    CALL mp_bcast (omega   , root_pool, intra_pool_comm)  
    CALL mp_bcast (alat    , ionode_id, inter_pool_comm)
    CALL mp_bcast (alat    , root_pool, intra_pool_comm)  
    IF (mpime /= ionode_id) ALLOCATE( tau( 3, nat ) )
    CALL mp_bcast (tau     , ionode_id, inter_pool_comm)
    CALL mp_bcast (tau     , root_pool, intra_pool_comm)  
    CALL mp_bcast (amass   , ionode_id, inter_pool_comm)
    CALL mp_bcast (amass   , root_pool, intra_pool_comm)  
    CALL mp_bcast (ityp    , ionode_id, inter_pool_comm)
    CALL mp_bcast (ityp    , root_pool, intra_pool_comm)  
    CALL mp_bcast (isk     , ionode_id, inter_pool_comm)
    CALL mp_bcast (isk     , root_pool, intra_pool_comm)  
    CALL mp_bcast (noncolin, ionode_id, inter_pool_comm)
    CALL mp_bcast (noncolin, root_pool, intra_pool_comm)  
    CALL mp_bcast (w_centers, ionode_id, inter_pool_comm)
    CALL mp_bcast (w_centers, root_pool, intra_pool_comm)
    IF (mpime.eq.ionode_id) THEN
      CLOSE(crystal)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    ! 
  ELSE
    continue
  ENDIF
  !
  ALLOCATE( w2( 3*nat) )
  !
  ! Determine Wigner-Seitz points
  ! For this we need the Wannier centers
  ! w_centers is allocated inside loadumat
  IF (.not. epwread) THEN
    xxq = 0.d0
    CALL loadumat( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband, w_centers )
  ENDIF
  !
  IF (use_ws) THEN
    ! Use Wannier-centers to contstruct the WS for electonic part and el-ph part
    ! Use atomic position to contstruct the WS for the phonon part and el-ph part
    dims  = nbndsub
    dims2 = nat
    CALL wigner_seitz_wrap ( nk1, nk2, nk3, nq1, nq2, nq3, irvec_k, irvec_q, irvec_g, &
                             ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g,   &
                             w_centers, dims, tau, dims2 )
  ELSE
    ! Center the WS for for electonic part and el-ph part
    ! Use atomic position to contstruct the WS for the phonon part and el-ph part
    dims  = 1
    dims2 = 1
    dummy(:) = (/0.0,0.0,0.0/)
    CALL wigner_seitz_wrap ( nk1, nk2, nk3, nq1, nq2, nq3, irvec_k, irvec_q, irvec_g, &
                             ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, & 
                             dummy, dims, dummy, dims2 )
  ENDIF
  ! 
  ! Determine the size of the respective WS sets based on the length of the matrices
  nrr_k = SIZE(irvec_k(1,:))
  nrr_q = SIZE(irvec_q(1,:))
  nrr_g = SIZE(irvec_g(1,:))
  IF (use_ws) THEN
    WRITE(stdout, '(5x,a)' )   '      '
    WRITE(stdout, '(5x,a)' )   'Construct the Wigner-Seitz cell using Wannier centers and atomic positions '
    WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ',nrr_k
    WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ',nrr_q
    WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electron-phonon ',nrr_g
    WRITE(stdout, '(5x,a,i8)' ) 'Maximum number of cores for efficient parallelization ',nrr_g * nat 
  ELSE
    WRITE(stdout, '(5x,a)' )      '      '
    WRITE(stdout, '(5x,a)' )      'Use zone-centred Wigner-Seitz cells '
    WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ',nrr_k
    WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ',nrr_q
    WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electron-phonon ',nrr_g
    WRITE(stdout, '(5x,a,i8)' ) 'Maximum number of cores for efficient parallelization ',nrr_g * nmodes
    WRITE(stdout, '(5x,a)' )      'Results may improve by using use_ws == .true. '
  ENDIF
  !
#ifndef __MPI  
  ! Open like this only in sequential. Otherwize open with MPI-open
  IF (ionode) THEN
    ! open the .epmatwe file with the proper record length
    lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
    filint    = trim(prefix)//'.epmatwp'
    CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
  ENDIF
#endif
  ! 
  ! At this point, we will interpolate the Wannier rep to the Bloch rep 
  !
  IF ( epwread ) THEN
     !
     !  read all quantities in Wannier representation from file
     !  in parallel case all pools read the same file
     !
     CALL epw_read(nrr_k, nrr_q, nrr_g)
     !
  ELSE !if not epwread (i.e. need to calculate fmt file)
     ! 
     IF (ionode) THEN
       lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
       filint    = trim(prefix)//'.epmatwe'
       CALL diropn (iunepmatwe, 'epmatwe', lrepmatw, exst)
       filint    = trim(prefix)//'.epmatwp'
       CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
     ENDIF
     !
     !xxq = 0.d0 
     !CALL loadumat &
     !     ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband )  
     !
     ! ------------------------------------------------------
     !   Bloch to Wannier transform
     ! ------------------------------------------------------
     !
     ALLOCATE( chw    ( nbndsub, nbndsub, nrr_k ),        &
               chw_ks ( nbndsub, nbndsub, nrr_k ),        &
               rdw    ( nmodes,  nmodes,  nrr_q ) )

     IF (vme) THEN
       ALLOCATE( cvmew ( 3, nbndsub, nbndsub, nrr_k ) )
     ELSE
       ALLOCATE( cdmew ( 3, nbndsub, nbndsub, nrr_k ) )
     ENDIF
     ! 
     ! SP : Let the user chose. If false use files on disk
     ALLOCATE(epmatwe_mem ( nbndsub, nbndsub, nrr_k, nmodes))
     epmatwe_mem(:,:,:,:) = czero
     !
     ! Hamiltonian
     !
     CALL hambloch2wan &
          ( nbnd, nbndsub, nks, nkstot, et, xk, cu, lwin, exband, nrr_k, irvec_k, wslen_k, chw )
     !
     ! Kohn-Sham eigenvalues
     !
     IF (eig_read) THEN
       WRITE (stdout,'(5x,a)') "Interpolating MB and KS eigenvalues"
       CALL hambloch2wan &
            ( nbnd, nbndsub, nks, nkstot, et_ks, xk, cu, lwin, exband, nrr_k, irvec_k, wslen_k, chw_ks )
     ENDIF
     !
     IF (vme) THEN
       ! Transform of position matrix elements
       ! PRB 74 195118  (2006)
       CALL vmebloch2wan &
            ( nbnd, nbndsub, nks, nkstot, xk, cu, nrr_k, irvec_k, wslen_k, lwin, exband )
     ELSE
       ! Dipole
       CALL dmebloch2wan &
            ( nbnd, nbndsub, nks, nkstot, dmec, xk, cu, nrr_k, irvec_k, wslen_k, lwin, exband )
     ENDIF
     !
     ! Dynamical Matrix 
     !
     IF (.not. lifc) CALL dynbloch2wan &
                          ( nmodes, nqc, xqc, dynq, nrr_q, irvec_q, wslen_q )
     !
     !
     ! Electron-Phonon vertex (Bloch el and Bloch ph -> Wannier el and Bloch ph)
     !
     DO iq = 1, nqc
       !
       xxq = xqc (:, iq)
       !
       ! we need the cu again for the k+q points, we generate the map here
       !
       CALL loadumat ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband, w_centers )
       !
       DO imode = 1, nmodes
         !
         CALL ephbloch2wane &
           ( nbnd, nbndsub, nks, nkstot, xk, cu, cuq, &
           epmatq (:,:,:,imode,iq), nrr_k, irvec_k, wslen_k, epmatwe_mem(:,:,:,imode) )
         !
       ENDDO
       ! Only the master node writes 
       IF (ionode) THEN
         ! direct write of epmatwe for this iq 
         CALL rwepmatw ( epmatwe_mem, nbndsub, nrr_k, nmodes, iq, iunepmatwe, +1)       
         !   
       ENDIF   
       !
     ENDDO
     !
     ! Electron-Phonon vertex (Wannier el and Bloch ph -> Wannier el and Wannier ph)
     !
     ! Only master perform this task. Need to be parallelize in the future (SP)
     IF (ionode) THEN
       CALL ephbloch2wanp_mem &
        ( nbndsub, nmodes, xqc, nqc, irvec_k, irvec_g, nrr_k, nrr_g, epmatwe_mem )
     ENDIF
     !
     CALL mp_barrier(inter_pool_comm)
     !
     IF ( epwwrite ) THEN
        CALL epw_write(nrr_k, nrr_q, nrr_g, w_centers)
        CALL epw_read(nrr_k, nrr_q, nrr_g) 
     ENDIF
     !
  ENDIF
  !
  IF ( ALLOCATED (epmatwe_mem) ) DEALLOCATE (epmatwe_mem)
  IF ( ALLOCATED (epmatq) )  DEALLOCATE (epmatq)
  IF ( ALLOCATED (cu) )      DEALLOCATE (cu)
  IF ( ALLOCATED (cuq) )     DEALLOCATE (cuq)
  IF ( ALLOCATED (lwin) )    DEALLOCATE (lwin)
  IF ( ALLOCATED (lwinq) )   DEALLOCATE (lwinq)
  IF ( ALLOCATED (exband) )  DEALLOCATE (exband)
  CLOSE(iunepmatwe, status = 'delete')
  CLOSE(iunepmatwp)
  ! 
  ! Check Memory usage
  CALL system_mem_usage(valueRSS)
  ! 
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a,i10,a)' ) '     Memory usage:  VmHWM =',valueRSS(2)/1024,'Mb'
  WRITE(stdout, '(a,i10,a)' ) '                   VmPeak =',valueRSS(1)/1024,'Mb'
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a)' )             '     '
  
  !
  !  At this point, we will interpolate the Wannier rep to the Bloch rep 
  !  for electrons, phonons and the ep-matrix
  !
  !  need to add some sort of parallelization (on g-vectors?)  what
  !  else can be done when we don't ever see the wfcs??
  !
  ! SP: Only k-para
  CALL loadqmesh_serial
  CALL loadkmesh_para
  !
  ALLOCATE ( epmatwef( nbndsub, nbndsub, nrr_k),         &
             wf( nmodes,  nqf ),                         &
             etf( nbndsub, nkqf),                        &
             etf_ks( nbndsub, nkqf),                     &
             epmatf( nbndsub, nbndsub),                  &
             cufkk( nbndsub, nbndsub),                   &
             cufkq( nbndsub, nbndsub),                   & 
             uf( nmodes, nmodes),                        &
             bmatf( nbndsub, nbndsub),                   & 
             eps_rpa( nmodes) )
  !
  ! Need to be initialized
  etf_ks(:,:) = zero
  epmatf(:,:) = czero
  ! allocate velocity and dipole matrix elements after getting grid size
  !
  IF (vme) THEN 
     ALLOCATE ( vmef(3, nbndsub, nbndsub, 2 * nkf) )
  ELSE
     ALLOCATE ( dmef(3, nbndsub, nbndsub, 2 * nkf) )
  ENDIF
  !
  ALLOCATE(cfac(nrr_k,dims,dims))
  ALLOCATE(cfacq(nrr_k,dims,dims))
  ALLOCATE(rdotk(nrr_k))
  ALLOCATE(rdotk2(nrr_k))
  ! This is simply because dgemv take only real number (not integer)
  ALLOCATE(irvec_r(3,nrr_k))
  irvec_r = REAL(irvec_k,KIND=dp)
  ! 
  ! ------------------------------------------------------
  ! Hamiltonian : Wannier -> Bloch (preliminary)
  ! ------------------------------------------------------
  !
  ! We here perform a preliminary interpolation of the hamiltonian
  ! in order to determine the fermi window ibndmin:ibndmax for later use.
  ! We will interpolate again afterwards, for each k and k+q separately
  !
  xxq = 0.d0
  !
  ! nkqf is the number of kpoints in the pool
  DO ik = 1, nkqf
     !
     xxk = xkf (:, ik)
     !
     IF ( 2*(ik/2).eq.ik ) THEN
        !
        !  this is a k+q point : redefine as xkf (:, ik-1) + xxq
        !
        CALL cryst_to_cart ( 1, xxq, at,-1 )
        xxk = xkf (:, ik-1) + xxq
        CALL cryst_to_cart ( 1, xxq, bg, 1 )
        !
     ENDIF
     !
     ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
     ! + optimize the 2\pi r\cdot k with Blas
     CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xxk, 1, 0.0_DP, rdotk, 1 )
     !  
     DO iw=1, dims
       DO iw2=1, dims
         DO ir=1, nrr_k
           IF (ndegen_k(ir,iw2,iw) > 0 ) &
             cfac(ir,iw2,iw) = exp( ci*rdotk(ir) ) / ndegen_k(ir,iw2,iw)
         ENDDO
       ENDDO
     ENDDO
     ! 
     CALL hamwan2bloch &
          ( nbndsub, nrr_k, cufkk, etf (:, ik), chw, cfac, dims)
     !
     !
  ENDDO
  !
  WRITE(stdout,'(/5x,a,f10.6,a)') 'Fermi energy coarse grid = ', ef * ryd2ev, ' eV'
  !
  IF( efermi_read ) THEN
     !
     ef = fermi_energy
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout, '(/5x,a,f10.6,a)') &
         'Fermi energy is read from the input file: Ef = ', ef * ryd2ev, ' eV'
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     !
     ! SP: even when reading from input the number of electron needs to be correct
     already_skipped = .false.
     IF ( nbndskip .gt. 0 ) THEN
        IF ( .not. already_skipped ) THEN
           IF ( noncolin ) THEN
              nelec = nelec - one * nbndskip
           ELSE
              nelec = nelec - two * nbndskip
           ENDIF
           already_skipped = .true.
           WRITE(stdout,'(/5x,"Skipping the first ",i4," bands:")') nbndskip
           WRITE(stdout,'(/5x,"The Fermi level will be determined with ",f9.5," electrons")') nelec
        ENDIF
     ENDIF
     !      
  ELSEIF( band_plot ) THEN 
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout, '(/5x,"Fermi energy corresponds to the coarse k-mesh")')
     WRITE(stdout,'(/5x,a)') repeat('=',67) 
     !
  ELSE 
     ! here we take into account that we may skip bands when we wannierize
     ! (spin-unpolarized)
     ! RM - add the noncolin case
     already_skipped = .false.
     IF ( nbndskip .gt. 0 ) THEN
        IF ( .not. already_skipped ) THEN
           IF ( noncolin ) THEN 
              nelec = nelec - one * nbndskip
           ELSE
              nelec = nelec - two * nbndskip
           ENDIF
           already_skipped = .true.
           WRITE(stdout,'(/5x,"Skipping the first ",i4," bands:")') nbndskip
           WRITE(stdout,'(/5x,"The Fermi level will be determined with ",f9.5," electrons")') nelec
        ENDIF
     ENDIF
     !
     ! Fermi energy
     !  
     ! since wkf(:,ikq) = 0 these bands do not bring any contribution to Fermi level
     !  
     efnew = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
     !
     WRITE(stdout, '(/5x,a,f10.6,a)') &
         'Fermi energy is calculated from the fine k-mesh: Ef = ', efnew * ryd2ev, ' eV'
     !
     ! if 'fine' Fermi level differs by more than 250 meV, there is probably something wrong
     ! with the wannier functions, or 'coarse' Fermi level is inaccurate
     IF (abs(efnew - ef) * ryd2eV .gt. 0.250d0 .and. (.not.eig_read) ) &
        WRITE(stdout,'(/5x,a)') 'Warning: check if difference with Fermi level fine grid makes sense'
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     !
     ef = efnew
     !
  ENDIF
  !
  ! identify the bands within fsthick from the Fermi level
  ! (in shuffle mode this actually does not depend on q)
  !
  !
  CALL fermiwindow
  ! 
  ! Define it only once for the full run. 
  CALL fkbounds( nkqtotf/2, lower_bnd, upper_bnd )
  !
  ! Re-order the k-point according to weather they are in or out of the fshick
  ! windows
  IF (iterative_bte .and. mp_mesh_k) THEN
    CALL load_rebal()
  ENDIF
  !  xqf must be in crystal coordinates
  !
  ! this loops over the fine mesh of q points.
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  IF (lifc) THEN
    !
    ! build the WS cell corresponding to the force constant grid
    !
    atws(:,1) = at(:,1)*DBLE(nq1)
    atws(:,2) = at(:,2)*DBLE(nq2)
    atws(:,3) = at(:,3)*DBLE(nq3)
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,atws)
  ENDIF
  !
  ! Open the ephmatwp file here
#if defined(__MPI)
  ! Check for directory given by "outdir"
  !      
  filint = trim(tmp_dir)//trim(prefix)//'.epmatwp1'
  CALL MPI_FILE_OPEN(world_comm,filint,MPI_MODE_RDONLY,MPI_INFO_NULL,iunepmatwp2,ierr)
  IF( ierr /= 0 ) CALL errore( 'ephwann_shuffle', 'error in MPI_FILE_OPEN',1 )
#endif
  !
  ! get the size of the matrix elements stored in each pool
  ! for informational purposes.  Not necessary
  !
  CALL mem_size(ibndmin, ibndmax, nmodes, nkf)
  !
  ! Fine mesh set of g-matrices.  It is large for memory storage
  ALLOCATE ( epf17 (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nkf) )
  ALLOCATE ( etf_all (ibndmax-ibndmin+1, nkqtotf/2 ) )
  ALLOCATE ( inv_tau_all (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
  ALLOCATE ( zi_allvb (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
  ALLOCATE ( epmatlrT (nbndsub, nbndsub, nmodes, nkf) )
  ALLOCATE ( eptmp (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nkf) )
  epf17(:,:,:,:)     = czero
  etf_all(:,:)       = zero
  inv_tau_all(:,:,:) = zero
  zi_allvb(:,:,:)    = zero
  ! 
  IF (int_mob .AND. carrier) THEN
    ALLOCATE ( inv_tau_allcb (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
    ALLOCATE ( zi_allcb (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
    inv_tau_allcb(:,:,:) = zero
    zi_allcb(:,:,:) = zero
  ENDIF
  !
  !
  ! ------------------------------------------------
  ! The IBTE implement works in two steps
  ! 1) compute the dominant scattering rates and store them to file
  ! 2) read them from file and solve the IBTE where all important element are in memory
  ! ------------------------------------------------
  !  
  ! Initialization and restart when doing IBTE
  IF (iterative_bte) THEN
998 continue ! Continue after all scattering rates have been computed in print_ibte
    IF (epmatkqread) THEN
      ! 
      ALLOCATE( vkk_all( 3, ibndmax-ibndmin+1, nkqtotf/2 ) )
      ALLOCATE( wkf_all( nkqtotf/2 ) )
      ! 
      CALL iter_restart(etf_all, wkf_all, vkk_all, ind_tot, ind_totcb, ef0, efcb)
      ! 
      DEALLOCATE(vkk_all)
      DEALLOCATE(wkf_all)
      DEALLOCATE(etf_all)
      GOTO 999
      ! 
    ELSE ! epmatkqread
      !  
      ! Open the required files
      CALL iter_open(ind_tot, ind_totcb, lrepmatw2, lrepmatw4, lrepmatw5, lrepmatw6)
      ! 
    ENDIF
  ENDIF
  ! 
  ! -----------------------------------------------------------------------
  ! Determines which q-points falls within the fsthick windows
  ! Store the result in the selecq.fmt file 
  ! If the file exists, automatically restart from the file
  ! -----------------------------------------------------------------------
  ! 
  ! Check if the grids are homogeneous and commensurate
  homogeneous = .FALSE.
  IF ( (nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0) .AND. &
       (nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0) .AND. &
       (MOD(nkf1,nqf1) == 0) .AND. (MOD(nkf2,nqf2) == 0) .AND. (MOD(nkf3,nqf3) == 0) ) THEN
    homogeneous = .TRUE.
  ELSE
    homogeneous = .FALSE.
  ENDIF
  ! 
  totq = 0
  ! 
  ! Check if the file has been pre-computed
  IF (mpime == ionode_id) THEN
    INQUIRE(FILE='selecq.fmt',EXIST=exst)
  ENDIF
  CALL mp_bcast(exst, ionode_id, world_comm)
  ! 
  IF (exst) THEN
    IF (selecqread) THEN
      WRITE(stdout,'(5x,a)')' '
      WRITE(stdout,'(5x,a)')'Reading selecq.fmt file. '
      CALL qwindow(exst, nrr_k, dims, totq, selecq, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
    ELSE
      WRITE(stdout,'(5x,a)')' '
      WRITE(stdout,'(5x,a)')'A selecq.fmt file was found but re-created because selecqread == .false. '
      CALL qwindow(.FALSE., nrr_k, dims, totq, selecq, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
    ENDIF
  ELSE ! exst
    IF (selecqread) THEN
      CALL errore( 'ephwann_shuffle', 'Variable selecqread == .true. but file selecq.fmt not found.',1 )
    ELSE
      CALL qwindow(exst, nrr_k, dims, totq, selecq, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
    ENDIF
  ENDIF
  ! 
  WRITE(stdout,'(5x,a,i8,a)')'We only need to compute ',totq, ' q-points'
  WRITE(stdout,'(5x,a)')' '
  ! 
  ! -----------------------------------------------------------------------
  ! Possible restart during step 1) 
  ! -----------------------------------------------------------------------
  iq_restart = 1
  first_cycle = .FALSE.
  first_time = .TRUE.
  ! 
  ! Restart in SERTA case or self-energy case
  IF (restart) THEN
    IF ( elecselfen ) THEN
      ALLOCATE( sigmar_all(ibndmax-ibndmin+1, nkqtotf/2) )
      ALLOCATE( sigmai_all(ibndmax-ibndmin+1, nkqtotf/2) )
      ALLOCATE( zi_all(ibndmax-ibndmin+1, nkqtotf/2) )
      sigmar_all(:,:) = zero
      sigmai_all(:,:) = zero
      zi_all(:,:)     = zero
      CALL electron_read(iq_restart, totq, nkqtotf/2, sigmar_all, sigmai_all, zi_all)
    ENDIF
    IF ( scattering ) THEN
      IF (int_mob .AND. carrier) THEN
        ! Here inv_tau_all and inv_tau_allcb gets updated
        CALL tau_read(iq_restart, totq, nkqtotf/2, .TRUE.)
      ELSE
        ! Here inv_tau_all gets updated
        CALL tau_read(iq_restart, totq, nkqtotf/2, .FALSE.)
      ENDIF
    ENDIF
    !
    ! If you restart from reading a file. This prevent 
    ! the case were you restart but the file does not exist
    IF (iq_restart > 1) first_cycle = .TRUE.
    ! 
  ENDIF ! restart
  ! 
  ! Scatread assumes that you alread have done the full q-integration
  ! We just do one loop to get interpolated eigenenergies.  
  IF(scatread) iq_restart = totq -1
  ! 
  ! Restart in IBTE case
  IF (iterative_bte) THEN
    IF (mpime == ionode_id) THEN
      INQUIRE(FILE='restart_ibte.fmt',EXIST=exst)
    ENDIF
    CALL mp_bcast(exst, ionode_id, world_comm)
    ! 
    IF (exst) THEN
      IF (mpime.eq.ionode_id) THEN
        OPEN(unit=iunrestart,file='restart_ibte.fmt',status='old',iostat=ios)
        READ (iunrestart,*) iq_restart
        READ (iunrestart,*) ind_tot
        READ (iunrestart,*) ind_totcb
        READ (iunrestart,*) lrepmatw2
        READ (iunrestart,*) lrepmatw4
        READ (iunrestart,*) lrepmatw5
        READ (iunrestart,*) lrepmatw6
        CLOSE(iunrestart)
      ENDIF
      CALL mp_bcast(iq_restart, ionode_id, world_comm )
#if defined(__MPI)
      CALL MPI_BCAST( ind_tot, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
      CALL MPI_BCAST( ind_totcb, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
      CALL MPI_BCAST( lrepmatw2, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
      CALL MPI_BCAST( lrepmatw4, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
      CALL MPI_BCAST( lrepmatw5, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
      CALL MPI_BCAST( lrepmatw6, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
#else
      CALL mp_bcast( ind_tot,   ionode_id, world_comm )
      CALL mp_bcast( ind_totcb, ionode_id, world_comm )
      CALL mp_bcast( lrepmatw2, ionode_id, world_comm )
      CALL mp_bcast( lrepmatw4, ionode_id, world_comm )
      CALL mp_bcast( lrepmatw5, ionode_id, world_comm )
      CALL mp_bcast( lrepmatw6, ionode_id, world_comm )
#endif
      IF( ierr /= 0 ) CALL errore( 'ephwann_shuffle', 'error in MPI_BCAST',1 )
      ! 
      ! Now, the iq_restart point has been done, so we need to do the next one 
      iq_restart = iq_restart + 1
      WRITE(stdout,'(5x,a,i8,a)')'We restart from ',iq_restart, ' q-points'
    ENDIF ! exst
  ENDIF
  ! -----------------------------------------------------------------------------
  !
  DO iqq = iq_restart, totq
    ! This needs to be uncommented. 
    epf17(:,:,:,:) = czero
    ! 
    iq = selecq(iqq)
    !   
    CALL start_clock ( 'ep-interp' )
    !
    ! In case of big calculation, show progression of iq (especially usefull when
    ! elecselfen = true as nothing happen during the calculation otherwise. 
    !
    IF ( .not. phonselfen) THEN 
      IF (MOD(iqq,100) == 0) THEN
        WRITE(stdout, '(a,i10,a,i10)' ) '     Progression iq (fine) = ',iqq,'/',totq
      ENDIF
    ENDIF
    !
    xxq = xqf (:, iq)
    !
    ! ------------------------------------------------------
    ! dynamical matrix : Wannier -> Bloch
    ! ------------------------------------------------------
    !
    IF (.not. lifc) THEN
      CALL dynwan2bloch &
          ( nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf, w2 )
    ELSE
      CALL dynifc2blochf ( nmodes, rws, nrws, xxq, uf, w2 )
    ENDIF
    !
    ! ...then take into account the mass factors and square-root the frequencies...
    !
    DO nu = 1, nmodes
      !
      ! wf are the interpolated eigenfrequencies
      ! (omega on fine grid)
      !
      IF ( w2 (nu) .gt. 0.d0 ) THEN
        wf(nu,iq) =  sqrt(abs( w2 (nu) ))
      ELSE
        wf(nu,iq) = -sqrt(abs( w2 (nu) ))
      ENDIF
      !
      DO mu = 1, nmodes
        na = (mu - 1) / 3 + 1
        uf (mu, nu) = uf (mu, nu) / sqrt(amass(ityp(na)))
      ENDDO
    ENDDO
    !
    ! --------------------------------------------------------------
    ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
    ! --------------------------------------------------------------
    !
    epf17(:,:,:,:) = czero
    eptmp(:,:,:,:) = czero
    epmatlrT(:,:,:,:) = czero
    cufkk(:,:) = czero
    cufkq(:,:) = czero
    DO imode = 1, nmodes 
      epmatwef(:,:,:) = czero
      !DBSP              
      !CALL start_clock ( 'cl2' )
      !write(stdout,*) 'imode, nmodes, xxq, SUM(irvec_g), SUM(ndegen_g), nrr_g, nbndsub, nrr_k, dims ',&
      !       imode, nmodes, xxq, SUM(irvec_g), SUM(ndegen_g), nrr_g, nbndsub, nrr_k, dims
      IF (.NOT. longrange) THEN
        CALL ephwan2blochp_mem &
            (imode, nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatwef, nbndsub, nrr_k, dims, dims2 )
      ENDIF
      !write(stdout,*)'epmatwef ',sum(epmatwef)
      !CALL stop_clock ( 'cl2' )
      !
      !
      !  number of k points with a band on the Fermi surface
      fermicount = 0
      !
      IF (lscreen) THEN
        IF (scr_typ == 0) CALL rpa_epsilon (xxq, wf(:,iq), nmodes, epsi, eps_rpa)
        IF (scr_typ == 1) CALL tf_epsilon (xxq, nmodes, epsi, eps_rpa)
      ENDIF
      ! this is a loop over k blocks in the pool
      ! (size of the local k-set)
      DO ik = 1, nkf
        !print*,'ik in ephwann ', ik
        !
        ! xkf is assumed to be in crys coord
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        xkk = xkf(:, ikk)
        xkq = xkk + xxq
        !
        CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
        CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk2, 1 )
        !
        DO iw=1, dims
          DO iw2=1, dims
            DO ir = 1, nrr_k
              IF (ndegen_k(ir,iw2,iw) > 0) THEN
                cfac(ir,iw2,iw)  = exp( ci*rdotk(ir) ) / ndegen_k(ir,iw2,iw)
                cfacq(ir,iw2,iw) = exp( ci*rdotk2(ir) ) / ndegen_k(ir,iw2,iw)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        !
        ! ------------------------------------------------------        
        ! hamiltonian : Wannier -> Bloch 
        ! ------------------------------------------------------
        !
        ! Kohn-Sham first, then get the rotation matricies for following interp.
        IF (eig_read) THEN
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkk, etf_ks(:, ikk), chw_ks, cfac, dims)
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkq, etf_ks(:, ikq), chw_ks, cfacq, dims)
        ENDIF
        !
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkk, etf(:, ikk), chw, cfac, dims)
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkq, etf(:, ikq), chw, cfacq, dims)
        !
        IF (vme) THEN
           !
           ! ------------------------------------------------------
           !  velocity: Wannier -> Bloch
           ! ------------------------------------------------------
           !
           IF (eig_read) THEN
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_k, cufkk, vmef(:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw_ks, cfac, dims )
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_k, cufkq, vmef(:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks, cfacq, dims )
           ELSE
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_k, cufkk, vmef(:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw, cfac, dims )
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_k, cufkq, vmef(:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw, cfacq, dims )
           ENDIF
        ELSE
           !
           ! ------------------------------------------------------
           !  dipole: Wannier -> Bloch
           ! ------------------------------------------------------
           !
           CALL dmewan2bloch &
                ( nbndsub, nrr_k, cufkk, dmef(:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), cfac, dims )
           CALL dmewan2bloch &
                ( nbndsub, nrr_k, cufkq, dmef(:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), cfacq, dims )
           !
        ENDIF
        !
        IF (.NOT. scatread) THEN
          ! interpolate only when (k,k+q) both have at least one band 
          ! within a Fermi shell of size fsthick 
          !
          !IF (ik==2 .and. imode == 1) print*,iq, etf(:, ikk), etf(:, ikq), ef 
          IF ( (( minval ( abs(etf(:, ikk) - ef) ) < fsthick ) .and. & 
                ( minval ( abs(etf(:, ikq) - ef) ) < fsthick )) ) THEN
            !
            !  fermicount = fermicount + 1
            !
            ! --------------------------------------------------------------
            ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
            ! --------------------------------------------------------------
            !
            !
            ! SP: Note: In case of polar materials, computing the long-range and short-range term 
            !     separately might help speed up the convergence. Indeed the long-range term should be 
            !     much faster to compute. Note however that the short-range term still contains a linear
            !     long-range part and therefore could still be a bit more difficult to converge than 
            !     non-polar materials. 
            ! 
            IF (longrange) THEN
              !      
              epmatf(:,:) = czero
              !
            ELSE
              !
              epmatf(:,:) = czero
              CALL ephwan2bloch_mem &
                ( nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, cfac, dims )
              !
            ENDIF
            !
            IF (lpolar) THEN
              !
              CALL compute_umn_f( nbndsub, cufkk, cufkq, bmatf )
              !
              IF ( (abs(xxq(1)) > eps8) .or. (abs(xxq(2)) > eps8) .or. (abs(xxq(3)) > eps8) ) THEN
                !      
                CALL cryst_to_cart (1, xxq, bg, 1)
                CALL rgd_blk_epw_fine_mem(imode, nq1, nq2, nq3, xxq, uf, epmatlrT(:,:,imode,ik), &
                                      nmodes, epsi, zstar, bmatf, one)
                CALL cryst_to_cart (1, xxq, at, -1)
                !
              ENDIF
              !
            ENDIF
            ! 
            ! Store epmatf in memory
            !
            DO jbnd = ibndmin, ibndmax
              DO ibnd = ibndmin, ibndmax
                ! 
                IF (lscreen) THEN
                  eptmp(ibnd-ibndmin+1,jbnd-ibndmin+1,imode,ik) = epmatf(ibnd,jbnd) / eps_rpa(imode)
                ELSE 
                  eptmp(ibnd-ibndmin+1,jbnd-ibndmin+1,imode,ik) = epmatf(ibnd,jbnd)
                ENDIF
                !
              ENDDO
            ENDDO
            ! 
          ENDIF
        ENDIF ! scatread 
        
      ENDDO  ! end loop over k points
    ENDDO ! modes 
    !
    ! Now do the eigenvector rotation:
    ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
    !
    DO ik=1, nkf
      CALL zgemm( 'n', 'n', (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1), nmodes, nmodes, cone, eptmp(:,:,:,ik),&
            (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1), uf, nmodes, czero, &
            epf17(:,:,:,ik), (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) )
      ! 
    ENDDO
    ! 
    ! After the rotation, add the long-range that is already rotated
    DO jbnd = ibndmin, ibndmax
      DO ibnd = ibndmin, ibndmax
        epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,:) = epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,:) + epmatlrT(ibnd,jbnd,:,:)
      ENDDO
    ENDDO
    !
    !
    !
    IF (prtgkk     ) CALL print_gkk( iq )
    IF (phonselfen ) CALL selfen_phon_q( iqq, iq, totq )
    IF (elecselfen ) CALL selfen_elec_q( iqq, iq, totq, first_cycle )
    IF (plselfen .and. .not. vme ) CALL selfen_pl_q( iqq, iq, totq )
    IF (nest_fn    ) CALL nesting_fn_q( iqq, iq )
    IF (specfun_el ) CALL spectral_func_q( iqq, iq, totq )
    IF (specfun_ph ) CALL spectral_func_ph( iqq, iq, totq )
    IF (specfun_pl .and. .not. vme ) CALL spectral_func_pl_q( iqq, iq, totq )
    IF (ephwrite) THEN
      IF ( iq == 1 ) THEN 
        CALL kmesh_fine
        CALL kqmap_fine
      ENDIF
      CALL write_ephmat( iq ) 
      CALL count_kpoints( iq )
    ENDIF
    ! 
    ! Conductivity ---------------------------------------------------------
    !IF (scattering) CALL scattering_rate_q( iq )
    IF (.NOT. scatread) THEN
      IF (scattering) THEN
        ! Apply a scissor shift to CBM if required by user
        ! The shift is apply to k and k+q
        !IF (scissor > 0.000001) THEN
        IF (ABS(scissor) > 0.000001) THEN
          IF ( noncolin ) THEN
            icbm = FLOOR(nelec/1.0d0) +1
          ELSE
            icbm = FLOOR(nelec/2.0d0) +1
          ENDIF
          !
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            DO ibnd = icbm, nbndsub
              etf (ibnd, ikk) = etf (ibnd, ikk) + scissor
              etf (ibnd, ikq) = etf (ibnd, ikq) + scissor
            ENDDO
          ENDDO
          IF ( iqq == 1 ) THEN
            WRITE(stdout, '(5x,"Applying a scissor shift of ",f9.5," eV to the conduction states")' ) scissor * ryd2ev
          ENDIF
        ENDIF
        !   
        ! If we want to compute intrinsic mobilities, call fermicarrier to 
        ! correctly positionned the ef0 level.
        ! This is only done once for iq = 0 
        IF ( iqq == iq_restart ) THEN
          ! 
          DO itemp = 1, nstemp
            ! 
            etemp = transp_temp(itemp) 
            WRITE(stdout, '(/5x,"Temperature ",f8.3," K")' ) etemp * ryd2ev / kelvin2eV
            ! 
            ! Small gap semiconductor. Computes intrinsic mobility by placing 
            ! the Fermi level such that carrier density is equal for electron and holes
            IF (int_mob .AND. .NOT. carrier) THEN               
              !
              ef0(itemp) = fermicarrier( etemp )
              WRITE(stdout, '(5x,"Mobility Fermi level ",f10.6," eV")' )  ef0(itemp) * ryd2ev  
              ! We only compute 1 Fermi level so we do not need the other
              efcb(itemp) = 0
              !   
            ENDIF
            ! 
            ! Large bandgap semiconductor. Place the gap at the value ncarrier.
            ! The user want both VB and CB mobilities. 
            IF (int_mob .AND. carrier) THEN
              ! 
              ncarrier = - ABS(ncarrier) 
              ef0(itemp) = fermicarrier( etemp )
              WRITE(stdout, '(5x,"Mobility VB Fermi level ",f10.6," eV")' )  ef0(itemp) * ryd2ev 
              ! 
              ncarrier = ABS(ncarrier) 
              efcb(itemp) = fermicarrier( etemp )
              WRITE(stdout, '(5x,"Mobility CB Fermi level ",f10.6," eV")' )  efcb(itemp) * ryd2ev
              !  
            ENDIF   
            ! 
            ! User decide the carrier concentration and choose to only look at VB or CB  
            IF (.NOT. int_mob .AND. carrier) THEN
              ! SP: Determination of the Fermi level for intrinsic or doped carrier 
              !     One also need to apply scissor before calling it.
              ! 
              ef0(itemp) = fermicarrier( etemp )               
              WRITE(stdout, '(5x,"Mobility Fermi level ",f10.6," eV")' )  ef0(itemp) * ryd2ev
              ! We only compute 1 Fermi level so we do not need the other
              efcb(itemp) = 0
              ! 
            ENDIF
            ! 
            IF (.NOT. int_mob .AND. .NOT. carrier ) THEN
              IF ( efermi_read ) THEN
                !
                ef0(itemp) = fermi_energy
                !
              ELSE !SP: This is added for efficiency reason because the efermig routine is slow
                ef0(itemp) = efnew
              ENDIF
              ! We only compute 1 Fermi level so we do not need the other
              efcb(itemp) = 0
              !  
            ENDIF
            ! 
          ENDDO
          !
          ! 
        ENDIF ! iqq=0
        !   
        IF ( .NOT. iterative_bte ) THEN
          CALL scattering_rate_q( iqq, iq, totq, ef0, efcb, first_cycle )
          ! Computes the SERTA mobility
          IF (iqq == totq) CALL transport_coeffs (ef0,efcb)
        ENDIF
        ! 
        IF (iterative_bte) THEN
          CALL print_ibte( iqq, iq, totq, ef0, efcb, first_cycle, ind_tot, ind_totcb, lrepmatw2,&
                           lrepmatw4, lrepmatw5, lrepmatw6 )
          !  
          ! Finished, now compute SERTA and IBTE mobilities
          IF (iqq == totq) THEN
            WRITE(stdout, '(5x,a)')' '
            WRITE(stdout, '(5x,"epmatkqread automatically changed to .true. as all scattering have been computed.")')
            WRITE(stdout, '(5x,a)')' '
            epmatkqread = .true.
            GOTO 998
          ENDIF
        ENDIF
        ! 
      ENDIF ! scattering
      ! --------------------------------------    
      !
      CALL stop_clock ( 'ep-interp' )
      !
    ENDIF ! scatread
  ENDDO  ! end loop over q points
  !
  ! close files
  CALL iter_close()
  ! 
  ! Check Memory usage
  CALL system_mem_usage(valueRSS)
  ! 
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a,i10,a)' ) '     Memory usage:  VmHWM =',valueRSS(2)/1024,'Mb'
  WRITE(stdout, '(a,i10,a)' ) '                   VmPeak =',valueRSS(1)/1024,'Mb'
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a)' )
  ! 
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------  
  !
  ! SP: Added lambda and phonon lifetime writing to file.
  ! 
  CALL mp_barrier(inter_pool_comm)
  IF (mpime.eq.ionode_id) THEN
    !
    IF (phonselfen) THEN
      OPEN(unit=lambda_phself,file='lambda.phself')
      WRITE(lambda_phself, '(/2x,a/)') '#Lambda phonon self-energy'
      WRITE(lambda_phself, *) '#Modes     ',(imode, imode=1,nmodes)
      DO iqq = 1, totq
          !
          !myfmt = "(*(3x,E15.5))"  This does not work with PGI
        myfmt = "(1000(3x,E15.5))"
        WRITE(lambda_phself,'(i9,4x)',advance='no') iqq
        WRITE(lambda_phself, fmt=myfmt) (REAL(lambda_all(imode,iqq,1)),imode=1,nmodes)
          !
      ENDDO
      CLOSE(lambda_phself)
      OPEN(unit=linewidth_phself,file='linewidth.phself')
      WRITE(linewidth_phself, '(a)') '# Phonon frequency and phonon lifetime in meV '
      WRITE(linewidth_phself,'(a)') '# Q-point  Mode   Phonon freq (meV)   Phonon linewidth (meV)'
      DO iqq = 1, totq
        !
        DO imode=1, nmodes
          WRITE(linewidth_phself,'(i9,i6,E20.8,E22.10)') iqq,imode,&
                                 ryd2mev*wf(imode,iqq),ryd2mev*REAL(gamma_all(imode,iqq,1))
        ENDDO
        !
      ENDDO
      CLOSE(linewidth_phself)
    ENDIF
  ENDIF
  IF (band_plot) CALL plot_band
  !
  IF (a2f) CALL eliashberg_a2f
  ! 
  ! if scattering is read then Fermi level and scissor have not been computed.
  IF (scatread) THEN
    IF (ABS(scissor) > 0.000001) THEN
      icbm = FLOOR(nelec/2.0d0) + nbndskip + 1
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        DO ibnd = icbm, nbndsub
          etf (ibnd, ikk) = etf (ibnd, ikk) + scissor
          etf (ibnd, ikq) = etf (ibnd, ikq) + scissor
        ENDDO
      ENDDO
      WRITE( stdout, '(5x,"Applying a scissor shift of ",f9.5," eV to the conduction states")' ) scissor * ryd2ev
    ENDIF          
    DO itemp = 1, nstemp
      etemp = transp_temp(itemp)      
      IF (int_mob .OR. carrier) THEN
        ! SP: Determination of the Fermi level for intrinsic or doped carrier 
        !     One also need to apply scissor before calling it.
        ef0(itemp) = fermicarrier( etemp )
      ELSE
        IF ( efermi_read ) THEN
          ef0(itemp) = fermi_energy
        ELSE !SP: This is added for efficiency reason because the efermig routine is slow
          ef0(itemp) = efnew
        ENDIF
      ENDIF
    ENDDO ! itemp
  ENDIF ! if scattering 
  ! 
  IF ( ALLOCATED(lambda_all) )   DEALLOCATE( lambda_all )
  IF ( ALLOCATED(gamma_all) )    DEALLOCATE( gamma_all )
  IF ( ALLOCATED(sigmai_all) )   DEALLOCATE( sigmai_all )
  IF ( ALLOCATED(sigmai_mode) )  DEALLOCATE( sigmai_mode )
  IF ( ALLOCATED(w2) )           DEALLOCATE( w2 )
  DEALLOCATE(cfac)
  DEALLOCATE(cfacq)
  DEALLOCATE(rdotk)
  DEALLOCATE(irvec_r)
  DEALLOCATE(irvec_k)
  DEALLOCATE(irvec_q)
  DEALLOCATE(ndegen_k)
  DEALLOCATE(ndegen_q)
  DEALLOCATE(ndegen_g)
  DEALLOCATE(wslen_k)
  DEALLOCATE(wslen_q)
  DEALLOCATE(wslen_g)
  ! 
999 continue
  ! 
  CALL stop_clock ( 'ephwann' )
  !
  END SUBROUTINE ephwann_shuffle_mem
  !
