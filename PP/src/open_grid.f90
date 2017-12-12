!------------------------------------------------------------------------
PROGRAM open_grid
  !------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup, npool, nproc_pool, nproc_pool_file
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE cell_base,  ONLY : at, bg, tpiba2, alat
  USE lsda_mod,   ONLY : nspin, isk
  USE klist,      ONLY : nks, nkstot, xk, wk, igk_k, ngk, qnorm
  USE io_files,   ONLY : prefix, tmp_dir, nwordwfc, iunwfc, diropn
  USE noncollin_module,   ONLY : noncolin
  USE control_flags,      ONLY : gamma_only, twfcollect
  USE environment,        ONLY : environment_start, environment_end
  USE symm_base,          ONLY : nrot, nsym, s, t_rev
  USE parameters,         ONLY : npk
  USE exx,                ONLY : nq1,nq2,nq3, ecutfock, igk_exx, xkq_collect, &
                                 nkqs, exxinit, exx_mp_init, use_ace, exxbuff,&
                                 index_xk, exx_fft, exx_grid_init 
  USE gvecw,              ONLY: ecutwfc, gcutw
  USE gvect,              ONLY : g, ngm
  USE gvecs,              ONLY : nls
  USE funct,              ONLY : dft_force_hybrid
  USE wvfct,              ONLY : nbnd, npwx, g2kin, et, wg
  USE wavefunctions_module, ONLY : evc
  USE buffers,            ONLY : save_buffer, open_buffer, close_buffer
  USE scf,                ONLY : rho
  USE lsda_mod,           ONLY : nspin, lsda
  USE io_rho_xml,         ONLY : write_scf
  USE input_parameters,   ONLY : nk1, nk2, nk3, k1, k2, k3, k_points, &
                              occupations, calculation !, nkstot,
  USE noncollin_module,   ONLY : nspin_mag, npol
  USE fft_interfaces,     ONLY : fwfft
  !
  USE qexsd_module,       ONLY : qexsd_input_obj
  USE qes_types_module,   ONLY : input_type
  USE fft_base,           ONLY : dffts
  !USE qexsd_input,        ONLY : qexsd_init_k_points_ibz
  USE control_flags,      ONLY : gamma_only, io_level
  USE start_k, ONLY : init_start_k
  ! 
  IMPLICIT NONE
  !
  INTERFACE
     SUBROUTINE   pw_init_qexsd_input(obj,obj_tagname)
     IMPORT                       :: input_type
     TYPE(input_type)             :: obj
     CHARACTER(LEN=*),INTENT(IN)  :: obj_tagname
     END SUBROUTINE
  END INTERFACE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, ik, ibnd, ik_idx, ik_idx_kpt, ik_idx_exx, is
  CHARACTER(len=4) :: spin_component
  CHARACTER(len=256) :: outdir
  !INTEGER :: nq(3)
  LOGICAL :: exst, opnd, exst_mem
  REAL(DP),ALLOCATABLE :: et0(:,:), wg0(:,:), yk(:,:), wk0(:)
  INTEGER, EXTERNAL  :: n_plane_waves
  COMPLEX(DP),ALLOCATABLE :: psic(:), evx(:,:)
  !
  LOGICAL           :: use_ace_back, exx_status_back
  REAL(DP)          :: ecutfock_back
  INTEGER           :: nq_back(3)

  ! these are in wannier module.....-> integer :: ispinw, ikstart, ikstop, iknum
  NAMELIST / inputpp / outdir, prefix !, nq
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
  CALL mp_bcast(outdir,ionode_id, world_comm)
  CALL mp_bcast(tmp_dir,ionode_id, world_comm)
  CALL mp_bcast(prefix,ionode_id, world_comm)
  !CALL mp_bcast(nq,ionode_id, world_comm)
  !
  WRITE(stdout,*)
  WRITE(stdout,*) ' Reading nscf_save data'
  CALL read_file
  !print*, "initial", nks, nkstot
  CALL open_buffer(iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst)
  !
  twfcollect = .false.
  !
  WRITE(stdout,*)
  IF ( npool > 1 .and. nspin_mag>1) CALL errore( 'open_grid', &
      'pools+spin not implemented, use wf_collect instead', npool )
  !
  IF(gamma_only) CALL errore("open_grid", &
      "not implemented, and pointless, for gamma-only",1)
  !
  ! Here we trap restarts from a different number of nodes.
  !
  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('open_grid', &
     'pw.x run on a different number of procs/pools. Use wf_collect=.true.',1)
  !
!  CALL openfil_pp()
  !
  ! Store some variables related to exx to be put back before writing
  use_ace_back = use_ace
  ecutfock_back = ecutfock
  nq_back = (/ nq1, nq2, nq3 /)
  exx_status_back = .true.
  CALL dft_force_hybrid(exx_status_back)
  !
  nq1 = -1
  nq2 = -1
  nq3 = -1
  ecutfock = 4*ecutwfc
  use_ace = .false.
  
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
  !print*, "spin:", nspin, nspin_mag
  nks = nkqs
  nkstot = nks
  ! Temporary order of xk points, they will be resorted later (if nspin>1) 
  xk(:,1:nks) = xkq_collect(:,1:nks)
  wk(1:nks) = 1._dp/DBLE(nks) !/DBLE(nspin_mag)
  npwx = n_plane_waves(gcutw, nks, xk, g, ngm)
  !
  DEALLOCATE(igk_k, ngk, et, wg)
  ALLOCATE(igk_k(npwx,nks), ngk(nks))
  ALLOCATE(et(nbnd,nks), wg(nbnd,nks))
  !
  DEALLOCATE(evc)
  ALLOCATE(evc(npwx,nbnd))
  !
  prefix = TRIM(prefix)//"_open"
  nwordwfc = nbnd * npwx * npol
  !WRITE(stdout,*) ' Nwordwfc:', nwordwfc, nbnd, npwx, npol
  CALL open_buffer(iunwfc, 'wfc', nwordwfc, +1, exst_mem, exst)
  !
  twfcollect = .false.
  CALL write_scf(rho, nspin)
  !
  ALLOCATE(psic(dffts%nnr), evx(npol*npwx, nbnd))
  DO ik = 1, nks !/nspin_mag
  !DO is = 1, nspin_mag
    ik_idx_kpt = ik !+ (is-1)*(nks/nspin_mag) !(ik-1)*nspin_mag + is
    ik_idx_exx = ik !+ (is-1)*(nks/nspin_mag)
!    print*, ik_idx_kpt, ik_idx_exx
    xk(:,ik_idx_kpt) = xkq_collect(:,ik_idx_exx)
!    wk(ik_idx_kpt) = 1._dp/nkqs*DBLE(nspin_mag)
!     print*, ik, is, xk(:,ik_idx_kpt)
    !
    CALL gk_sort (xk(:,ik_idx_kpt), ngm, g, ecutwfc / tpiba2, &
                  ngk(ik_idx_kpt), igk_k(:,ik_idx_kpt), g2kin)
!     print*, size(exxbuff,1), size(exxbuff,2), nwordwfc, npwx, &
!             dffts%nnr, exx_fft%dfftt%nnr
    DO ibnd = 1, nbnd
      psic(1:dffts%nnr) = exxbuff(1:dffts%nnr,ibnd,ik_idx_exx)
      CALL fwfft('Wave', psic, dffts)
      evx(1:ngk(ik_idx_kpt),ibnd) = psic(nls(igk_k(1:ngk(ik_idx_kpt),ik_idx_kpt)))
      !
      IF(noncolin)THEN
        psic(1:dffts%nnr) = exxbuff(dffts%nnr+1:2*dffts%nnr,ibnd,ik_idx_exx)
        CALL fwfft('Wave', psic, dffts)
        evx(npwx+1:npwx+ngk(ik_idx_kpt),ibnd) = &
                                  psic(nls(igk_k(1:ngk(ik_idx_kpt),ik_idx_kpt)))
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
  nk1 = nq1
  nk2 = nq2
  nk3 = nq3
  calculation = 'bands'
  k_points = "automatic"
  CALL init_start_k(nk1,nk2,nk3, k1, k2, k3, "automatic",nks/nspin_mag, xk, wk)
#if defined(__OLDXML)
#else
  !
  ! HACK: rebuild input structure, this uses unallocated stuff
  !print*, nk1, nk2, nk3, k1, k2, k3, k_points
  !CALL pw_init_qexsd_input(qexsd_input_obj, obj_tagname="input")
#endif
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
  WRITE(stdout,'(5x,a,3i4)') "Dimensions:", nk1, nk2, nk3
  WRITE(stdout,'(5x,a,3i4)') "Shift:     ", k1,k2,k3
  WRITE(stdout,'(5x,a)') "List to be put in the .win file of wannier90: &
                          &(already in crystal/fractionary coordinates):"
    
  DO ik = 1, nks/nspin_mag
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
