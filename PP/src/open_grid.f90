!------------------------------------------------------------------------
PROGRAM open_grid
  !------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup
  USE mp_images,  ONLY : intra_image_comm
  USE mp_pools,   ONLY : npool
  USE mp,         ONLY : mp_bcast
  USE cell_base,  ONLY : at, bg, tpiba2, alat
  USE klist,      ONLY : nks, nkstot, xk, wk, igk_k, ngk, qnorm
  USE io_files,   ONLY : prefix, tmp_dir, nwordwfc, iunwfc, diropn
  USE noncollin_module,   ONLY : noncolin, domag, m_loc, angle1, angle2, nspin_lsda
  USE control_flags,      ONLY : gamma_only
  USE environment,        ONLY : environment_start, environment_end
  USE ions_base,          ONLY : nat, tau, ityp
  USE symm_base,          ONLY : nrot, nsym, s, t_rev, fft_fact, find_sym
  USE parameters,         ONLY : npk
  USE exx_base,           ONLY : nq1,nq2,nq3, xkq_collect, &
                                 nkqs, exx_mp_init, index_xk, exx_grid_init 
  USE exx,                ONLY : exxbuff, exxinit, use_ace, ecutfock
  USE gvecw,              ONLY : ecutwfc, gcutw
  USE gvect,              ONLY : g, ngm
  USE xc_lib,             ONLY : dft_force_hybrid
  USE wvfct,              ONLY : nbnd, npwx, g2kin, et, wg
  USE wavefunctions, ONLY : evc
  USE buffers,            ONLY : save_buffer, open_buffer, close_buffer
  USE scf,                ONLY : rho
  USE lsda_mod,           ONLY : nspin, isk, lsda, starting_magnetization
  USE io_rho_xml,         ONLY : write_scf
  USE noncollin_module,   ONLY : nspin_mag, npol
  USE fft_interfaces,     ONLY : fwfft
  !
  USE fft_base,           ONLY : dffts
  USE control_flags,      ONLY : gamma_only, io_level
  USE start_k, ONLY : init_start_k
  USE extfield,           ONLY : gate
  USE esm,                ONLY : esm_z_inv
  USE rism_module,        ONLY : lrism
  USE command_line_options, ONLY : nband_, ntg_ 
  USE mp_pools,             ONLY : intra_pool_comm 
  USE mp_exx,               ONLY : mp_start_exx
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, ik, ibnd, ik_idx, ik_idx_kpt, ik_idx_exx, is, na
  CHARACTER(len=4) :: spin_component
  CHARACTER(len=256) :: outdir
  INTEGER :: k1, k2, k3
  LOGICAL :: exst, opnd, exst_mem, magnetic_sym
  REAL(DP),ALLOCATABLE :: et0(:,:), wg0(:,:), yk(:,:), wk0(:)
  INTEGER, EXTERNAL  :: n_plane_waves
  COMPLEX(DP),ALLOCATABLE :: psic(:), evx(:,:)
  !
  LOGICAL           :: use_ace_back, exx_status_back
  REAL(DP)          :: ecutfock_back
  INTEGER           :: nq_back(3)
  ! if true do not append '_open' to the prefix
  LOGICAL  :: overwrite_prefix = .false.

  ! these are in wannier module.....-> integer :: ispinw, ikstart, ikstop, iknum
  NAMELIST / inputpp / outdir, prefix, overwrite_prefix !, nq
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  !! not sure if this should be called also in 'library' mode or not !!
  CALL environment_start ( 'OPEN_GRID' )
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
     !nq = 0
     !
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck(outdir)
  ENDIF
  !
  !
  CALL mp_bcast(outdir,ionode_id, intra_image_comm)
  CALL mp_bcast(tmp_dir,ionode_id, intra_image_comm)
  CALL mp_bcast(prefix,ionode_id, intra_image_comm)
  CALL mp_bcast(overwrite_prefix,ionode_id, intra_image_comm)
  !
  WRITE(stdout,*)
  WRITE(stdout,*) ' Reading nscf_save data'
  CALL read_file
  !print*, "initial", nks, nkstot
  CALL open_buffer(iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst)
  !
  WRITE(stdout,*)
  IF ( npool > 1 .and. nspin_mag>1) CALL errore( 'open_grid', &
      'pools+spin not implemented, run this tool with only one pool', npool )
  !
  IF(gamma_only) CALL errore("open_grid", &
      "not implemented, and pointless, for gamma-only",1)
  !
  ! Store some variables related to exx to be put back before writing
  use_ace_back = use_ace
  ecutfock_back = ecutfock
  nq_back = (/ nq1, nq2, nq3 /)
  exx_status_back = .true.
  CALL dft_force_hybrid(exx_status_back)

  magnetic_sym = noncolin .AND. domag 
  ALLOCATE(m_loc(3,nat))
  IF (noncolin.and.domag) THEN
     DO na = 1, nat
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     ENDDO
  ENDIF
  CALL find_sym ( nat, tau, ityp, magnetic_sym, m_loc, &
                & gate .OR. (.NOT. esm_z_inv(lrism)) )

  nq1 = -1
  nq2 = -1
  nq3 = -1
  ecutfock = 4*ecutwfc
  use_ace = .false.
 
  CALL mp_start_exx (nband_, ntg_, intra_pool_comm)
  CALL exx_grid_init()
  CALL exx_mp_init()
  !
  CALL exxinit(.false.)
  qnorm = 0._dp
  !
  CALL close_buffer(iunwfc,'keep')
  !
  ALLOCATE(et0(nbnd,nks), wg0(nbnd,nks), wk0(nks))
  et0 = et
  wg0 = wg
  wk0 = wk(1:nks)
  !
  nks = nkqs
  nkstot = nks
  ! Temporary order of xk points, they will be resorted later (if nspin>1) 
  xk(:,1:nks) = xkq_collect(:,1:nks)
  wk(1:nks) = 1._dp/DBLE(nks) !/DBLE(nspin_mag)
  npwx = n_plane_waves(gcutw, nks, xk, g, ngm)
  IF (nspin==2) THEN
    isk(1:nks/2) = 1
    isk(nks/2+1:nks) = 2
  ENDIF
  !
  DEALLOCATE(igk_k, ngk, et, wg, g2kin)
  ALLOCATE(igk_k(npwx,nks), ngk(nks), g2kin(npwx))
  ALLOCATE(et(nbnd,nks), wg(nbnd,nks))
  !
  DEALLOCATE(evc)
  ALLOCATE(evc(npwx*npol,nbnd))
  !
  if (.not. overwrite_prefix) then
    prefix = TRIM(prefix)//"_open"
  end if
  nwordwfc = nbnd * npwx * npol
  CALL open_buffer(iunwfc, 'wfc', nwordwfc, +1, exst_mem, exst)
  !
  ! Write everything again with the new prefix
  CALL write_scf(rho, nspin)
  !
  ALLOCATE(psic(dffts%nnr), evx(npol*npwx, nbnd))
 
  DO ik = 1, nks
    ik_idx_kpt = ik !+ (is-1)*(nks/nspin_mag) !(ik-1)*nspin_mag + is
    ik_idx_exx = ik !+ (is-1)*(nks/nspin_mag)
    xk(:,ik_idx_kpt) = xkq_collect(:,ik_idx_exx)
!    wk(ik_idx_kpt) = 1._dp/nkqs*DBLE(nspin_mag)
    !
    CALL gk_sort (xk(:,ik_idx_kpt), ngm, g, ecutwfc / tpiba2, &
                  ngk(ik_idx_kpt), igk_k(:,ik_idx_kpt), g2kin)
    evx(:,:) = 0._dp
    DO ibnd = 1, nbnd
      psic(1:dffts%nnr) = exxbuff(1:dffts%nnr,ibnd,ik_idx_exx)
      CALL fwfft('Wave', psic, dffts)
      evx(1:ngk(ik_idx_kpt),ibnd) = psic(dffts%nl(igk_k(1:ngk(ik_idx_kpt),ik_idx_kpt)))
      !
      IF(noncolin)THEN
        psic(1:dffts%nnr) = exxbuff(dffts%nnr+1:2*dffts%nnr,ibnd,ik_idx_exx)
        CALL fwfft('Wave', psic, dffts)
        evx(npwx+1:npwx+ngk(ik_idx_kpt),ibnd) = &
                                  psic(dffts%nl(igk_k(1:ngk(ik_idx_kpt),ik_idx_kpt)))
      ENDIF
      !
    ENDDO
    !
    CALL save_buffer( evx, nwordwfc, iunwfc, ik_idx_kpt )
    et(:,ik_idx_kpt) = et0(:,index_xk(ik_idx_kpt))
    wg(:,ik_idx_kpt) = wg0(:,index_xk(ik_idx_kpt))/wk0(index_xk(ik_idx_kpt))*wk(ik_idx_kpt)
  !ENDDO
  ENDDO
  DEALLOCATE(psic, et0, wg0)
  !
  ! reconstruct input variables
  k1 = 0
  k2 = 0
  k3 = 0
  CALL init_start_k(nq1,nq2,nq3, k1, k2, k3, "automatic",nks/nspin_lsda, xk, wk)
  !
  ! Restore EXX variables
  use_ace = use_ace_back
  ecutfock = ecutfock_back
  nq1 = nq_back(1)
  nq2 = nq_back(2)
  nq3 = nq_back(3)
  CALL dft_force_hybrid(exx_status_back)
  !
  CALL punch('all')
  !
  ALLOCATE(yk(3,nks))
  yk = xk
  CALL cryst_to_cart(nks, yk, at, -1)
  WRITE(stdout,'(5x,a)') "Grid of q-points"
  WRITE(stdout,'(5x,a,3i4)') "Dimensions:", nq1, nq2, nq3
  WRITE(stdout,'(5x,a,3i4)') "Shift:     ", k1,k2,k3
  WRITE(stdout,'(5x,a)') "List to be put in the .win file of wannier90: &
                          &(already in crystal/fractionary coordinates):"

  DO ik = 1, nks/nspin_lsda
    WRITE(stdout,'(3f21.15,3x,f13.10)') yk(:,ik), wk(ik)
  ENDDO
  DEALLOCATE(yk)
  !
  CALL environment_end ( 'OPEN_GRID' )
  WRITE( stdout, *  )
  CALL stop_pp
  !
END PROGRAM open_grid
!
