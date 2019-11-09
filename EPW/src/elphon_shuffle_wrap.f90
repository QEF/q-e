  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle_wrap()
  !-----------------------------------------------------------------------
  !!
  !! Electron-phonon calculation with Wannier functions: load all phonon q's
  !!
  !! This routine is the main driver of the electron-phonon 
  !! calculation. It first calculates the electron-phonon matrix elements
  !! on the coarse mesh and then passes the data off to [[ephwann_shuffle]]
  !! to perform the interpolation.
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE mp_global,     ONLY : my_pool_id, world_comm, npool  
  USE mp_images,     ONLY : my_image_id, nimage
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_barrier, mp_bcast
  USE io_global,     ONLY : stdout, meta_ionode, meta_ionode_id, ionode_id
  USE us,            ONLY : nqxq, dq, qrad
  USE gvect,         ONLY : gcutm, ngm
  USE cellmd,        ONLY : cell_factor
  USE uspp_param,    ONLY : lmaxq, nbetam
  USE io_files,      ONLY : prefix, tmp_dir
  USE wavefunctions, ONLY : evc
  USE wvfct,         ONLY : npwx
  USE eqv,           ONLY : vlocq, dmuxc
  USE ions_base,     ONLY : nat, nsp, tau, ityp, amass
  USE control_flags, ONLY : iverbosity
  USE io_var,        ONLY : iuepb, iuqpeig, crystal
  USE pwcom,         ONLY : nks, nbnd, nkstot, nelec
  USE cell_base,     ONLY : at, bg, alat, omega
  USE symm_base,     ONLY : irt, s, nsym, ft, sname, invs, s_axis_to_cart,      &
                            sr, nrot, copy_sym, set_sym_bl, find_sym, inverse_s,& 
                            remove_sym, allfrac
  USE phcom,         ONLY : evq
  USE qpoint,        ONLY : igkq, xq, eigqts
  USE modes,         ONLY : nmodes, u, npert
  USE lr_symm_base,  ONLY : minus_q, rtau, gi, gimq, irotmq, nsymq, invsymq
  USE epwcom,        ONLY : epbread, epbwrite, epwread, lifc, etf_mem, vme,     &
                            nbndsub, iswitch, kmaps, eig_read, dvscf_dir,       & 
                            nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, lpolar
  USE elph2,         ONLY : epmatq, dynq, et_ks, xkq, ifc, umat, umat_all,      &
                            zstar, epsi, cu, cuq, lwin, lwinq, bmat,            &
                            exband, wscache
  USE klist_epw,     ONLY : et_loc, et_all
  USE constants_epw, ONLY : ryd2ev, zero, czero, eps6
  USE fft_base,      ONLY : dfftp
  USE control_ph,    ONLY : u_from_file
  USE noncollin_module, ONLY : m_loc, npol, noncolin
  USE iotk_module,   ONLY : iotk_open_read, iotk_scan_dat, iotk_free_unit,      &
                            iotk_close_read
  USE division,      ONLY : fkbounds
  USE uspp,          ONLY : okvan
  USE spin_orb,      ONLY : lspinorb 
  USE lrus,          ONLY : becp1
  USE becmod,        ONLY : becp, deallocate_bec_type
  USE phus,          ONLY : int1, int1_nc, int2, int2_so, int4, int4_nc, int5,  &
                            int5_so, alphap
  USE kfold,         ONLY : shift, createkmap_pw2, createkmap
  USE low_lvl,       ONLY : set_ndnmbr, eqvect_strict, read_modes
  USE io_epw,        ONLY : read_ifc, readdvscf
  USE poolgathering, ONLY : poolgather
  USE rigid_epw,     ONLY : compute_umn_c
  USE rotate,        ONLY : rotate_epmat, rotate_eigenm, star_q2, gmap_sym
  USE pw2wan2epw,    ONLY : compute_pmn_para
#if defined(__NAG)
  USE f90_unix_io,   ONLY : flush
#endif
  !
  ! --------------------------------------------------------------
  !
  IMPLICIT NONE
  ! 
  CHARACTER(LEN = 256) :: tempfile
  !! Temporary .eig file
  CHARACTER(LEN = 256) :: dirname
  !! Name of the directory
  CHARACTER(LEN = 256) :: filename
  !! Name of the file
  CHARACTER(LEN = 3) :: filelab
  !! Append the number of the core that works on that file
  CHARACTER(LEN = 80)   :: line
  !! Use to read external eigenvalues
  CHARACTER(LEN = 6), EXTERNAL :: int_to_char
  !! Transform an INTEGER into a character
  LOGICAL :: sym(48)
  !! Logical vectors that say which crystal symmetries exist in our system
  LOGICAL :: nog
  !! Find if G=0 or not in $$S(q_0)+G=q$$
  LOGICAL :: non_symmorphic
  !! Check whether the symmetry belongs to a non-symmorphic group
  !! non_symmorphic == TRUE if it has fractional translation. 
  LOGICAL :: exst
  !! Find if a file exists.
  INTEGER :: sym_smallq(48) 
  !! Set of all symmetries for the small group of one q.
  !! This is a subset of total crystal symmetries that remains
  !! after the q-point pertubation.
  INTEGER :: nqc_irr
  !! Number of qpoints in the irreducible wedge
  INTEGER :: nqc
  !! Number of qpoints on the uniform grid
  INTEGER :: maxvalue
  !! Temporary INTEGER for max value
  INTEGER :: nqxq_tmp
  !! Maximum G+q length  
  INTEGER :: ik
  !! Total k-point index
  INTEGER :: ios
  !! Contains the state of the opened file 
  INTEGER :: ik_start
  !! Lower bound for the k-point of the coarse grid in parallel 
  INTEGER :: ik_stop
  !! Higher bound for the k-point of the coarse grid in parallel 
  INTEGER :: gmapsym(ngm, 48)
  !! Correspondence G -> S(G)
  INTEGER :: nq
  !! Degeneracy of the star of q
  INTEGER :: isq(48)
  !! Index of q in the star of a given sym.op.
  INTEGER :: imq              
  !! Index of -q in the star of q (0 if not present)
  INTEGER :: sym_sgq(48)
  !! The symmetries giving the q point iq in the star
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
  INTEGER :: iunpun
  !! Unit of the file
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
  REAL(KIND = DP) :: w_centers(3, nbndsub)
  !! Wannier centers
  REAL(KIND = DP) :: qnorm_tmp
  !! Absolute value of xqc_irr
  REAL(KIND = DP) :: sumr(2, 3, nat, 3)
  !! Sum to impose the ASR
  REAL(KIND = DP), ALLOCATABLE :: xqc_irr(:, :)
  !! The qpoints in the irr wedge
  REAL(KIND = DP), ALLOCATABLE :: xqc(:, :)
  !! The qpoints in the uniform mesh
  REAL(KIND = DP), ALLOCATABLE :: wqlist(:)
  !! The corresponding weigths
  COMPLEX(KIND = DP) :: eigv(ngm, 48)
  !! $e^{ iGv}$ for 1...nsym (v the fractional translation)
  COMPLEX(KIND = DP) :: cz1(nmodes, nmodes)
  !! The eigenvectors for the first q in the star
  COMPLEX(KIND = DP) :: cz2(nmodes, nmodes)
  !! The rotated eigenvectors, for the current q in the star
  !
  CALL start_clock('elphon_wrap')
  !
  ! Read qpoint list from stdin
  !
  IF (meta_ionode) READ(5, *) nqc_irr
  CALL mp_bcast(nqc_irr, meta_ionode_id, world_comm)
  ALLOCATE(xqc_irr(3, nqc_irr), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating xqc_irr', 1)
  ALLOCATE(xqc(3, nqc1 * nqc2 * nqc3), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating xqc', 1)
  ALLOCATE(wqlist(nqc1 * nqc2 * nqc3), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating wqlist', 1)
  xqc_irr(:, :) = zero
  xqc(:, :)     = zero
  wqlist(:)     = zero
  !  
  IF (meta_ionode) THEN
    DO iq_irr = 1, nqc_irr
      READ(5,*) xqc_irr(:, iq_irr)
    ENDDO
  ENDIF
  CALL mp_bcast(xqc_irr, meta_ionode_id, world_comm)
  !
  ! fix for uspp
  ! this is needed to get the correct size of the interpolation table 'qrad' 
  ! for the non-local part of the pseudopotential in PW/src/allocate_nlpot.f90
  !
  maxvalue = nqxq
  DO iq_irr = 1, nqc_irr
    qnorm_tmp = DSQRT(xqc_irr(1, iq_irr)**2 + xqc_irr(2, iq_irr)**2 + xqc_irr(3, iq_irr)**2)
    nqxq_tmp = INT(((DSQRT(gcutm) + qnorm_tmp) / dq + 4) * cell_factor)
    IF (nqxq_tmp > maxvalue)  maxvalue = nqxq_tmp
  ENDDO
  !
  IF (maxvalue > nqxq) THEN
    IF (.NOT. epwread) THEN
      DEALLOCATE(qrad, STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating qrad', 1) 
    ENDIF
    ALLOCATE(qrad(maxvalue, nbetam * (nbetam + 1) / 2, lmaxq, nsp), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating qrad ', 1)
    ! 
    qrad(:, :, :, :) = zero
    ! RM - need to call init_us_1 to re-calculate qrad 
    CALL init_us_1()
  ENDIF
  ! 
  ! do not perform the check if restart
  IF (epwread .AND. .NOT. epbread) THEN
    CONTINUE
  ELSE
    IF (nkstot /= nkc1 * nkc2 * nkc3) CALL errore('elphon_shuffle_wrap', 'nscf run inconsistent with epw input', 1)  
  ENDIF
  !
  ! Read in external electronic eigenvalues. e.g. GW 
  !
  ALLOCATE(et_ks(nbnd, nks), STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating et_ks', 1)
  et_ks(:, :) = zero
  IF (eig_read) THEN
    IF (meta_ionode) THEN
      WRITE (stdout,'(5x,a,i5,a,i5,a)') "Reading external electronic eigenvalues (", &
           nbnd, ",", nkstot,")"
      tempfile = TRIM(prefix)//'.eig'
      OPEN(iuqpeig, FILE = tempfile, FORM = 'formatted', ACTION = 'read', IOSTAT = ios)
      IF (ios /= 0) CALL errore('elphon_shuffle_wrap','error opening' // tempfile, 1)
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
  !
  ! Do not recompute dipole matrix elements
  IF (epwread .AND. .NOT. epbread) THEN 
    CONTINUE
  ELSE
    ! compute coarse grid dipole matrix elements.  Very fast 
    IF (.NOT. vme) CALL compute_pmn_para
  ENDIF
  !
  !  gather electronic eigenvalues for subsequent shuffle
  !  
  IF (eig_read) THEN
    et_all(:, :) = zero
    CALL poolgather(nbnd, nkstot, nks, et_loc(1:nbnd, 1:nks), et_all)
  ENDIF
  !
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
    WRITE(stdout,'(/5x,a)') 'Using kmap and kgmap from disk'
  ENDIF
  ! 
  IF (epwread) THEN
    !
    ! We need some crystal info
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = crystal, FILE = 'crystal.fmt', STATUS = 'old', IOSTAT = ios)
      READ(crystal,*) nat
      READ(crystal,*) nmodes
      READ(crystal,*) nelec
      READ(crystal,*) at
      READ(crystal,*) bg
      READ(crystal,*) omega
      READ(crystal,*) alat
      ALLOCATE(tau(3, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating tau', 1)
      READ(crystal,*) tau
      READ(crystal,*) amass
      ALLOCATE(ityp(nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating ityp', 1)
      READ(crystal,*) ityp
      READ(crystal,*) noncolin
      READ(crystal,*) w_centers
      ! 
    ENDIF ! mpime == ionode_id
    CALL mp_bcast(nat      , ionode_id, world_comm)
    IF (mpime /= ionode_id) ALLOCATE(ityp(nat))
    CALL mp_bcast(nmodes   , ionode_id, world_comm)
    CALL mp_bcast(nelec    , ionode_id, world_comm)
    CALL mp_bcast(at       , ionode_id, world_comm)
    CALL mp_bcast(bg       , ionode_id, world_comm)
    CALL mp_bcast(omega    , ionode_id, world_comm)
    CALL mp_bcast(alat     , ionode_id, world_comm)
    IF (mpime /= ionode_id) ALLOCATE(tau(3, nat) )
    CALL mp_bcast(tau      , ionode_id, world_comm)
    CALL mp_bcast(amass    , ionode_id, world_comm)
    CALL mp_bcast(ityp     , ionode_id, world_comm)
    CALL mp_bcast(noncolin , ionode_id, world_comm)
    CALL mp_bcast(w_centers, ionode_id, world_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(crystal)
    ENDIF
  ENDIF ! epwread
  ! 
  IF (lifc) THEN
    ALLOCATE(ifc(nqc1, nqc2, nqc3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating ifc', 1)
    ifc(:, :, :, :, :, :, :) = zero
  ENDIF
  !
  ! Do not do symmetry stuff 
  IF (epwread .AND. .NOT. epbread) THEN
    CONTINUE
  ELSE
    !
    !  allocate dynamical matrix and ep matrix for all q's
    !
    ALLOCATE(dynq(nmodes, nmodes, nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating dynq', 1)
    ALLOCATE(epmatq(nbnd, nbnd, nks, nmodes, nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating epmatq', 1)
    ALLOCATE(epsi(3, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating epsi', 1)
    ALLOCATE(zstar(3, 3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating zstar', 1)
    ALLOCATE(bmat(nbnd, nbnd, nks, nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating bmat', 1)
    ALLOCATE(cu(nbnd, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating cu', 1)
    ALLOCATE(cuq(nbnd, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating cuq', 1)
    ALLOCATE(lwin(nbnd, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating lwin', 1)
    ALLOCATE(lwinq(nbnd, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating lwinq', 1)
    ALLOCATE(exband(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating exband', 1)
    dynq(:, :, :)         = czero
    epmatq(:, :, :, :, :) = czero
    epsi(:, :)            = zero
    zstar(:, :, :)        = zero
    bmat(:, :, :, :)      = czero
    cu(:, :, :)           = czero
    cuq(:, :, :)          = czero
    !
    ! read interatomic force constat matrix from q2r
    IF (lifc) THEN
      CALL read_ifc
    ENDIF
    !
    ! SP: The symmetries are now consistent with QE 5. This means that the order of the q in the star
    !     should be the same as in the .dyn files produced by QE 5.
    ! 
    !     First we start by setting up the lattice & crystal symm. as done in PHonon/PH/q2qstar.f90
    ! 
    ! ~~~~~~~~ setup Bravais lattice symmetry ~~~~~~~~
    CALL set_sym_bl() ! This should define the s matrix
    WRITE(stdout,'(5x,a,i3)') "Symmetries of Bravais lattice: ", nrot
    !
    ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
    CALL find_sym(nat, tau, ityp, .FALSE., m_loc)
    IF (.NOT. allfrac) CALL remove_sym(dfftp%nr1, dfftp%nr2, dfftp%nr3)
    WRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
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
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating evq', 1)
    ALLOCATE(xkq(3, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating xkq', 1)
    ALLOCATE(shift(nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating shift', 1)
    IF (lifc) THEN
      ALLOCATE(wscache(-2 * nqc3:2 * nqc3, -2 * nqc2:2 * nqc2, -2 * nqc1:2 * nqc1, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error allocating wscache', 1)
      wscache(:, :, :, :, :) = zero      
    ENDIF
    evq(:,:)   = zero
    xkq(:, :)  = zero
    shift(:)   = 0
    ! 
    ! In the loop over irr q-point, we need to read the pattern that
    ! corresponds to the dvscf file computed with QE 5.
    !
    sumr(:, :, :, :) = zero
    nqc = 0
    iq_first = 1
    DO iq_irr = 1, nqc_irr
      u_from_file = .TRUE.
      !tmp_dir_ph = './_ph0/'
      !
      !  read the displacement patterns
      !
      IF (u_from_file) THEN
         ierr = 0
         ! ... look for an empty unit (only ionode needs it)
         IF (meta_ionode) CALL iotk_free_unit(iunpun, ierr)
         dirname = TRIM(dvscf_dir) // TRIM(prefix) // '.phsave'
         filename = TRIM(dirname) // '/patterns.' // TRIM(int_to_char(iq_irr)) // '.xml'
         INQUIRE(FILE = TRIM(filename), EXIST = exst )
         IF (.NOT. exst) CALL errore('elphon_shuffle_wrap', &
                   'cannot open file for reading or writing', ierr)
         CALL iotk_open_read(iunpun, FILE = TRIM(filename), binary = .FALSE., ierr = ierr)
         CALL read_modes(iunpun, iq_irr, ierr)
         IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', ' Problem with modes file', 1)
         IF (meta_ionode) CALL iotk_close_read(iunpun)
      ENDIF
      !  
      WRITE(stdout, '(//5x,a)') REPEAT('=',67) 
      WRITE(stdout, '(5x,"irreducible q point # ",i4)') iq_irr
      WRITE(stdout, '(5x,a/)') REPEAT('=',67) 
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
      ! SP: Notice that the function copy_sym reshuffles the s matrix for each irr_q.  
      !     This is why we then need to call gmap_sym for each irr_q [see below]. 
      nsymq = copy_sym(nsym, sym)    
      !
      ! Recompute the inverses as the order of sym.ops. has changed
      CALL inverse_s()
      CALL s_axis_to_cart()
      !
      ! This computes gi, gimq
      CALL set_giq(xq, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
      WRITE(stdout,'(5x,a,i3)') "Symmetries of small group of q:", nsymq
      IF(minus_q) WRITE(stdout,'(10x,a)') "in addition sym. q -> -q+G:"
      ! 
      ! Finally this does some of the above again and also computes rtau...
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      ! ######################### star of q #########################
      ! 
      sym_smallq(:) = 0
      CALL star_q2(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .TRUE., sym_smallq)
      !
      ! The reason for xq instead of xq0 in the above is because xq is passed to QE through module  
      xq0 = xq
      !
      !  determine the G vector map S(G) -> G 
      !  SP: The mapping needs to be done for each irr_q because the QE 5 symmetry routine
      !      reshuffles the s matrix for each irr_q [putting the sym of the small group of q first].
      !
      !  [I checked that gmapsym(gmapsym(ig,isym),invs(isym)) = ig]
      CALL gmap_sym(nsym, s, ft, gmapsym, eigv, invs)
      !
      !  Re-set the variables needed for the pattern representation
      !  and the symmetries of the small group of irr-q
      !  (from phq_setup.f90)
      !
      DO isym = 1, nsym
        sym(isym) = .TRUE.
      ENDDO
      !
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      IF (meta_ionode) THEN
        CALL dynmat_asr(iq_irr, nqc_irr, nq, iq_first, sxq, imq, isq, invs, s, irt, rtau, sumr)
      ENDIF
      CALL mp_bcast(zstar, meta_ionode_id, world_comm)
      CALL mp_bcast(epsi , meta_ionode_id, world_comm)
      CALL mp_bcast(dynq , meta_ionode_id, world_comm)
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
        CALL epw_init(.FALSE.)
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
        !
        ! Prepare the gmap for the refolding
        !
        CALL createkmap(xq)              
        !
        IF (iverbosity == 1) THEN
          !
          ! Description of symmetries
          !
          WRITE(stdout, '(36x,"s",24x,"frac. trans.")')
          CALL s_axis_to_cart() ! give sr(:,:, isym)
          DO isym = 1, nsym
            WRITE(stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
            IF (ft(1,isym)**2 + ft(2,isym)**2 + ft(3,isym)**2 > 1.0d-8) THEN
                ft1 = at(1,1)*ft(1,isym) + at(1,2)*ft(2,isym) + at(1,3)*ft(3,isym)
                ft2 = at(2,1)*ft(1,isym) + at(2,2)*ft(2,isym) + at(2,3)*ft(3,isym)
                ft3 = at(3,1)*ft(1,isym) + at(3,2)*ft(2,isym) + at(3,3)*ft(3,isym)
                WRITE(stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (s(1,ipol,isym),ipol = 1,3), ft(1,isym)
                WRITE(stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                            (s(2,ipol,isym),ipol = 1,3), ft(2,isym)
                WRITE(stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                            (s(3,ipol,isym),ipol = 1,3), ft(3,isym)
                WRITE(stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (sr(1,ipol,isym),ipol = 1,3), ft1
                WRITE(stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                            (sr(2,ipol,isym),ipol = 1,3), ft2
                WRITE(stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                            (sr(3,ipol,isym),ipol = 1,3), ft3
            ELSE
                WRITE(stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                       isym,  (s (1, ipol, isym) , ipol = 1,3)
                WRITE(stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol = 1,3)
                WRITE(stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol = 1,3)
                WRITE(stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                                   isym, (sr(1,ipol,isym), ipol = 1, 3)
                WRITE(stdout, '(17x," (",3f11.7," )")')  (sr(2,ipol,isym), ipol = 1, 3)
                WRITE(stdout, '(17x," (",3f11.7," )"/)') (sr(3,ipol,isym), ipol = 1, 3)
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
        IF (nsq * nq /= nsym ) CALL errore('elphon_shuffle_wrap', 'wrong degeneracy', iq)
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
            WRITE(stdout,'(3i5,L3,L3)') iq, i, isym, nog, non_symmorphic
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
        CALL cryst_to_cart(1, aq, at, - 1)
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
            EXIT
          ENDIF
          ! If we enter into that loop it means that we have not found 
          ! such symmetry within the small group of q. 
          IF (jsym == nsq) THEN
            CALL errore('elphon_shuffle_wrap ', 'No sym. such that Sxq0=iq was found in the sgq !', 1)
          ENDIF
        ENDDO
        !
        !
        CALL loadumat(nbnd, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq, exband, w_centers)
        !
        ! Calculate overlap U_k+q U_k^\dagger
        IF (lpolar) CALL compute_umn_c(nbnd, nbndsub, nks, cu, cuq, bmat(:, :, :, nqc))
        !
        !   calculate the sandwiches
        !
        ! A more accurate way of doing this is to symmetrize the matrix element w.r.t.
        ! the small group of the given q in the star. I'm not doing this here.
        ! (but I checked that even without symm the result of full zone and irr zone
        ! are equal to 5+ digits).
        ! For any volunteers, please write to giustino@civet.berkeley.edu
        !
        CALL elphon_shuffle(iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, xq0, .FALSE.)
        !
        !  bring epmatq in the mode representation of iq_first, 
        !  and then in the cartesian representation of iq
        !
        CALL rotate_eigenm(iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2)
        !
        CALL rotate_epmat(cz1, cz2, xq, nqc, lwin, lwinq, exband)
        !DBSP
        !write(*,*)'epmatq(:,:,2,:,nqc)',SUM(epmatq(:,:,2,:,nqc))
        !write(*,*)'epmatq(:,:,2,:,nqc)**2',SUM((REAL(REAL(epmatq(:,:,2,:,nqc))))**2)+&
        !  SUM((REAL(AIMAG(epmatq(:,:,2,:,nqc))))**2)
        !print*,'dynq ', SUM(dynq(:,:,nqc))
        !print*,'et ', et_loc(:,2)
        !END
        ! SP: Now we treat separately the case imq == 0
        IF (imq == 0) THEN
          !
          ! SP: First the vlocq need to be initialized propertly with the first
          !     q in the star
          xq = -xq0
          CALL epw_init(.FALSE.)
          !
          ! retrieve the q in the star
          xq = -sxq(:, iq)
          !
          ! and populate the uniform grid
          nqc = nqc + 1
          xqc(:,nqc) = xq
          !
          IF (iq == 1) WRITE(stdout, *)
          WRITE(stdout,5) nqc, xq
          !
          !  prepare the gmap for the refolding
          !
          CALL createkmap(xq)
          !
          CALL loadumat(nbnd, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq, exband, w_centers)
          !
          ! Calculate overlap U_k+q U_k^\dagger
          IF (lpolar) CALL compute_umn_c(nbnd, nbndsub, nks, cu, cuq, bmat(:, :, :, nqc))
          !
          xq0 = -xq0
          !
          CALL elphon_shuffle(iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, xq0, .TRUE.)
          !  bring epmatq in the mode representation of iq_first, 
          !  and then in the cartesian representation of iq
          !
          CALL rotate_eigenm(iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2)
          !
          CALL rotate_epmat(cz1, cz2, xq, nqc, lwin, lwinq, exband)
          !
          xq0 = -xq0
        ENDIF ! end imq == 0  
        !
      ENDDO
      !
      iq_first = iq_first + nq
      if (imq == 0) iq_first = iq_first + nq
      !
    ENDDO ! irr-q loop
    ! 
    IF (nqc /= nqc1 * nqc2 * nqc3) CALL errore('elphon_shuffle_wrap', 'nqc /= nq1*nq2*nq3', nqc)
    wqlist = DBLE(1) / DBLE(nqc)
    !
    IF (lifc) THEN
      DEALLOCATE(wscache, STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating wscache', 1)
    ENDIF
    DEALLOCATE(evc, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating evc', 1)
    DEALLOCATE(evq, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating evq', 1)
    DEALLOCATE(xkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating xkq', 1)
    DEALLOCATE(shift, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating shift', 1)
    DEALLOCATE(vlocq, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating vlocq', 1)
    DEALLOCATE(dmuxc, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating dmuxc', 1)
    DEALLOCATE(eigqts, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating eigqts', 1)
    DEALLOCATE(rtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating rtau', 1)
    DEALLOCATE(u, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating u', 1)
    DEALLOCATE(npert, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating npert', 1)
    IF (okvan) THEN
      DEALLOCATE(int1, STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int1', 1)
      DEALLOCATE(int2, STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int2', 1)
      DEALLOCATE(int4, STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int4', 1)
      DEALLOCATE(int5, STAT = ierr)
      IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int5', 1)
      IF (noncolin) THEN 
        DEALLOCATE(int1_nc, STAT = ierr)
        IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int1_nc', 1)
        DEALLOCATE(int4_nc, STAT = ierr)
        IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int4_nc', 1)
        IF (lspinorb) THEN
          DEALLOCATE(int2_so, STAT = ierr)
          IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int2_so', 1)
          DEALLOCATE(int5_so, STAT = ierr)
          IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating int5_so', 1)
        ENDIF
      ENDIF
    ENDIF
    DO ik = 1, nks
      DO ipol = 1, 3
        CALL deallocate_bec_type(alphap(ipol, ik))
      ENDDO
    ENDDO
    DEALLOCATE(alphap, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating alphap', 1)
    DO ik = 1, SIZE(becp1)
      CALL deallocate_bec_type(becp1(ik))
    ENDDO
    DEALLOCATE(becp1, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating becp1', 1)
    CALL deallocate_bec_type(becp)
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
        IF (.NOT.  exst) CALL errore('elphon_shuffle_wrap', 'epb files not found ', 1)
        OPEN(iuepb, FILE = tempfile, FORM = 'unformatted')
        WRITE(stdout,'(/5x,"Reading epmatq from .epb files"/)') 
        READ(iuepb) nqc, xqc, et_loc, dynq, epmatq, zstar, epsi
        CLOSE(iuepb)
        WRITE(stdout,'(/5x,"The .epb files have been correctly read"/)')
      ENDIF
      !
      IF (epbwrite) THEN
        OPEN(iuepb, FILE = tempfile, FORM = 'unformatted')
        WRITE(stdout, '(/5x,"Writing epmatq on .epb files"/)') 
        WRITE(iuepb) nqc, xqc, et_loc, dynq, epmatq, zstar, epsi
        CLOSE(iuepb)
        WRITE(stdout, '(/5x,"The .epb files have been correctly written"/)')
      ENDIF
    ENDIF
  ENDIF
  !
  ! In case of image parallelization we want to stop after writing the .epb file
  IF (nimage > 1) THEN
    WRITE(stdout, '(/5x,"Image parallelization. The code will stop now. "/)')
    WRITE(stdout, '(/5x,"You need to restart a calculation by reading the .epb "/)')
    WRITE(stdout, '(/5x,"                       with pool parallelization only. "/)')
    CALL stop_epw()
  ENDIF
  !
  IF (.NOT. epbread .AND. epwread) THEN
    ! CV: need dummy nqc, xqc for the ephwann_shuffle call
    nqc = 1
    xqc = zero
    WRITE(stdout, '(/5x,"Do not need to read .epb files; read .fmt files"/)')
    !
  ENDIF
  !
  ! now dynq is the cartesian dyn mat ( not divided by the masses)
  ! and epmatq is the epmat in cartesian representation (rotation in elphon_shuffle)
  !
  ! free up some memory
  !
  NULLIFY(igkq)
  DEALLOCATE(umat_all, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating umat_all', 1)
  DEALLOCATE(umat, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating umat', 1)
  DEALLOCATE(xqc_irr, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating xqc_irr', 1)
  DEALLOCATE(wqlist, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating wqlist', 1)
  ! 
  IF (maxvalue > nqxq) THEN
    DEALLOCATE(qrad, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating qrad', 1)
  ENDIF
  ! 
  IF (.NOT. (epwread .AND. .NOT. epbread)) THEN
    DEALLOCATE(cu, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating cu', 1)
    DEALLOCATE(cuq, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating cuq', 1)
    DEALLOCATE(lwin, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating lwin', 1)
    DEALLOCATE(lwinq, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating lwinq', 1)
    DEALLOCATE(exband, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating exband', 1)
    DEALLOCATE(bmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating bmat', 1)
  ENDIF
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
  ! The electron-phonon wannier interpolation
  IF(etf_mem == 2) THEN
#if defined(__MPI)         
    CALL ephwann_shuffle_mem(nqc, xqc, w_centers)
#else
    WRITE(stdout,'(/5x,a)') 'WARNING: etf_mem==2 only works with MPI'
    WRITE(stdout,'(5x,a)')  '         Changing to etf_mem == 1 and continue ...'
    etf_mem = 1
    CALL ephwann_shuffle(nqc, xqc, w_centers)
#endif
  ELSE ! etf_mem == 0, 1 or 4
    CALL ephwann_shuffle(nqc, xqc, w_centers)
  ENDIF        
  ! 
  DEALLOCATE(xqc, STAT = ierr)
  IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating xqc', 1)
  IF (lifc) THEN
    DEALLOCATE(ifc, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating ifc', 1)    
  ENDIF
  ! 
  IF (epwread) THEN
    DEALLOCATE(tau, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating tau', 1)
    DEALLOCATE(ityp, STAT = ierr)
    IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', 'Error deallocating ityp', 1)
  ENDIF ! epwread
  !
5 FORMAT (8x,"q(",i5," ) = (",3f12.7," )") 
  !
  RETURN
  !---------------------------------------------------------------------------
  END SUBROUTINE elphon_shuffle_wrap
  !---------------------------------------------------------------------------
