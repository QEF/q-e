!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!
! Original version by Andrea Ferretti
! Modified mainly by Layla Martin-Samos
! Modified by Joe Stenuit
!
!=----------------------------------------------------------------------------=!
  MODULE io_base_export
!=----------------------------------------------------------------------------=!

  USE kinds

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: file_version = 202
  INTEGER :: restart_module_verbosity = 0

  END MODULE

!-----------------------------------------------------------------------
program gwl_punch
  !-----------------------------------------------------------------------
  !
  ! read in PWSCF data in XML format using IOTK lib
  ! then prepare matrices for GWL calculation
  ! 
  ! input:  namelist "&inputpp", with variables
  !   prefix       prefix of input files saved by program pwscf
  !   outdir       temporary directory where files resides
  !   pp_file      output file. If it is omitted, a directory 
  !                "prefix.export/" is created in outdir and
  !                some output files are put there. Anyway all the data 
  !                are accessible through the "prefix.export/index.xml" file which
  !                contains implicit pointers to all the other files in the
  !                export directory. If reading is done by the IOTK library
  !                all data appear to be in index.xml even if physically it
  !                is not. 
  !   uspp_spsi    using US PP if set .TRUE. writes S | psi > 
  !                and | psi > separately in the output file 
  !   single_file  one-file output is produced
  !   ascii        ....
  !
  !   pseudo_dir   pseudopotential directory
  !   psfile(:)    name of the pp file for each species 
  !    

  USE kinds,     ONLY : i4b
  USE gvect,  ONLY : mill
  use pwcom
  USE constants,            ONLY : rytoev
  use io_global, ONLY : stdout, ionode, ionode_id
  use io_files,  ONLY : psfile, pseudo_dir
  use io_files,  ONLY : prefix, tmp_dir
  use ions_base, ONLY : ntype => nsp
  use iotk_module
  use mp_pools, ONLY : kunit
  use mp, ONLY: mp_bcast
  use mp_world, ONLY: world_comm, mpime
  use control_flags, ONLY : gamma_only
  use uspp, ONLY : okvan
  use ldaU, ONLY : lda_plus_u
  USE basis,                ONLY : swfcatom
  use scf, only : vrs, vltot, v, kedtau
  USE klist,                ONLY : xk, wk, nks, nkstot
  USE fft_base,             ONLY : dfftp
  USE wannier_gw,       ONLY :  lwannier, &
                                num_nbndv, &
                                num_nbnds, &
                                nset, &
                                l_truncated_coulomb, &
                                truncation_radius, &
                                remainder, &
                                restart_gww, &
                                numw_prod, &
                                l_gram,&
                                l_head,&
                                n_gauss,&
                                omega_gauss, &
                                l_exchange, &
                                l_zero, &
                                l_wing, &
                                grid_type, &
                                nset_overlap, &
                                nspace,&
                                ecutoff_global, &
                                maxiter2,&
                                diago_thr2, &
                                l_plot_mlwf,&
                                l_pmatrix,&
                                npcol,&
                                nprow,&
                                n_pola_lanczos,&
                                n_self_lanczos,&
                                nsteps_lanczos_pola,&
                                nsteps_lanczos_self,&
                                s_pola_lanczos,&
                                s_self_lanczos,&
                                s_g_lanczos,&
                                l_pmat_diago,&
                                pmat_ethr, &
                                pmat_cutoff,&
                                pmat_type,&
                                lanczos_restart,&
                                n_pola_lanczos_eff,&
                                n_self_lanczos_eff,&
                                n_pmat,&
                                s_pmat,&
                                n_fast_pmat,&
                                off_fast_pmat,&
                                l_fast_pola,&
                                l_v_basis,&
                                v_cutoff,&
                                l_iter_algorithm,&
                                dual_pb, &
                                l_t_wannier,&
                                dual_vt,&
                                dual_vs,&
                                wannier_thres,&
                                s_first_state,&
                                s_last_state,&
                                l_selfconsistent,&
                                l_whole_s,&
                                l_ts_eigen,&
                                l_frac_occ,&
                                num_nbndv_min,&
                                l_cond_pol_base,&
                                l_semicore,&
                                n_semicore,&
                                l_semicore_read,&
                                l_verbose,&
                                l_contour,&
                                l_real,&
                                l_bse,&
                                s_bse,&
                                dual_bse,&
                                l_big_system,&
                                extra_pw_cutoff,&
                                l_list,&
                                l_scissor,&
                                scissor,&
                                l_full,&
                                n_full,&
                                l_simple
                 
 
  USE exchange_custom, ONLY : exchange_fast_dual
  
  !
  implicit none
  integer :: i, kunittmp, ios

  character(len=200) :: pp_file
  character(len=iotk_attlenx) :: attr
  logical :: found, uspp_spsi, ascii, single_file, raw
!  INTEGER(i4b), EXTERNAL :: C_MKDIR
 CHARACTER(LEN=256), EXTERNAL :: trimcheck
 CHARACTER(LEN=256) :: outdir

  NAMELIST /inputpw4gww/ prefix, outdir, pp_file, uspp_spsi, ascii, single_file, raw, &
                     psfile, pseudo_dir, &
                      lwannier, num_nbndv, &
                               nset,num_nbnds, &
                               l_truncated_coulomb, &
                               truncation_radius, &
                               remainder, restart_gww, numw_prod, &
                               l_gram, l_head, n_gauss, omega_gauss, l_exchange, &
                               l_zero, l_wing, grid_type, &
                               nset_overlap, nspace, &
                               ecutoff_global,&
                               maxiter2,diago_thr2,l_plot_mlwf,&
                               l_pmatrix, npcol,nprow,&
                               n_pola_lanczos, nsteps_lanczos_pola,nsteps_lanczos_self,&
                               s_pola_lanczos,s_self_lanczos,n_self_lanczos,s_g_lanczos,&
                               l_pmat_diago,pmat_ethr,pmat_cutoff, pmat_type, lanczos_restart,&
                               n_pola_lanczos_eff,n_self_lanczos_eff,n_pmat,s_pmat,n_fast_pmat,&
                               off_fast_pmat,l_fast_pola,l_v_basis,v_cutoff,l_iter_algorithm,&
                               dual_pb, l_t_wannier, dual_vt, dual_vs,wannier_thres,s_first_state,&
                               s_last_state,l_selfconsistent,l_whole_s,l_ts_eigen,l_frac_occ,num_nbndv_min,&
                               l_cond_pol_base,l_semicore,n_semicore,l_semicore_read, l_verbose, l_contour,&
                               l_real,exchange_fast_dual,l_bse,s_bse,dual_bse,l_big_system,extra_pw_cutoff,&
                               l_list,l_scissor,scissor,l_full,n_full,l_simple
                    

  !
  call start_pw4gww( )
  !
  !   set default values for variables in namelist
  !
  prefix='export'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  pp_file= ' '
  uspp_spsi = .FALSE.
  ascii = .FALSE.
  single_file = .FALSE.
  raw = .FALSE.
!
!       nppstr   = 1
!
  lwannier = .false.
  wannier_thres=0.d0
  num_nbndv(1:2) = 1
  num_nbnds = 1
  nset = 250
  l_truncated_coulomb = .false.
  truncation_radius = 10.d0
  remainder=-1
  restart_gww=-1
  numw_prod=1
  l_gram=.false.
  l_head=.false.
  l_exchange=.false.
  n_gauss=79
  omega_gauss=20.d0
  l_zero=.true.
  l_wing=.false.
  grid_type=3
  nset_overlap=250
  nspace=1
  ecutoff_global = 400.d0
  maxiter2=0
  diago_thr2=0.d0
  l_plot_mlwf=.false.
  l_pmatrix=.false.
  npcol=1
  nprow=1
  n_pola_lanczos=400
  n_self_lanczos=600
  nsteps_lanczos_pola=20
  nsteps_lanczos_self=40
  s_pola_lanczos=0.5d0
  s_self_lanczos=1d-12
  s_g_lanczos=0.d0
  l_pmat_diago=.true.
  pmat_ethr=1d-5
  pmat_cutoff=3.d0
  pmat_type=4
  lanczos_restart=0
  n_pola_lanczos_eff=0
  n_self_lanczos_eff=0
  n_pmat=100
  s_pmat=0.01d0
  n_fast_pmat=0
  off_fast_pmat=0.d0
  l_fast_pola=.false.
  l_v_basis=.false.
  v_cutoff=0.01d0
  l_iter_algorithm=.true.
  dual_pb=1.d0
  dual_vt=1.d0
  dual_vs=1.d0
  l_t_wannier=.true.
  s_first_state=0
  s_last_state=0
  l_selfconsistent=.false.
  l_whole_s=.false.
  l_ts_eigen=.true.
  l_frac_occ=.false.
  num_nbndv_min(1:2)=1
  l_cond_pol_base=.false.
  l_semicore=.false.
  n_semicore=1
  l_semicore_read=.false.
  l_verbose=.false.
  l_contour=.false.
  l_real=.false.
  exchange_fast_dual=4.d0
  l_bse=.false.
  s_bse=0.d0
  dual_bse=1.d0
  l_big_system=.false.
  l_list=.false.
  extra_pw_cutoff=0.d0
  l_scissor=.false.
  scissor=0.d0
  l_full=.false.
  n_full=0
  l_simple=.false.
  !
  !    Reading input file
  !
  IF ( ionode ) THEN
      !
      CALL input_from_file ( )
      !
      READ(5,inputpw4gww,IOSTAT=ios)
      !
!      call read_namelists( 'PW4GWW' )
      !
      IF (ios /= 0) CALL errore ('pw4gww', 'reading inputpw4gww namelist', ABS(ios) )
      !
!-----------------------------------------------------------------------
!      IF( pp_file == ' ' ) THEN
          !
!          pp_file = TRIM(prefix)//".export/index.xml"
          !
!          if(ionode) ios = C_MKDIR( TRIM(outdir)//"/"//TRIM(prefix)// &
!                           ".export" , LEN(TRIM(outdir)//"/"//TRIM(prefix)//".export") )
!      ENDIF
      !
  ENDIF
#if !defined(__MPI)
  dual_pb=4.d0
  dual_vs=4.d0
  dual_vt=4.d0
#endif
!-------------------------------------------------------------------------
  ! ... Broadcasting variables
!------------------------------------------------------------------------
  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( pp_file, ionode_id, world_comm )
  CALL mp_bcast( uspp_spsi, ionode_id, world_comm )
  CALL mp_bcast( ascii, ionode_id, world_comm )
  CALL mp_bcast( single_file, ionode_id, world_comm )
  CALL mp_bcast( raw, ionode_id, world_comm )
  CALL mp_bcast( pseudo_dir, ionode_id, world_comm )
  CALL mp_bcast( psfile, ionode_id, world_comm )
  CALL mp_bcast( lwannier,      ionode_id, world_comm )
  CALL mp_bcast( wannier_thres, ionode_id, world_comm)
  CALL mp_bcast( num_nbndv,     ionode_id, world_comm )
  CALL mp_bcast( num_nbnds,     ionode_id, world_comm )
  CALL mp_bcast( nset,          ionode_id, world_comm )
  CALL mp_bcast(l_truncated_coulomb, ionode_id, world_comm)
  CALL mp_bcast(truncation_radius, ionode_id, world_comm)
  CALL mp_bcast(remainder, ionode_id, world_comm)
  CALL mp_bcast(restart_gww, ionode_id, world_comm)
  call mp_bcast(numw_prod, ionode_id, world_comm)
  CALL mp_bcast(l_gram, ionode_id, world_comm)
  CALL mp_bcast(l_head, ionode_id, world_comm)
  CALL mp_bcast(n_gauss, ionode_id, world_comm)
  CALL mp_bcast(omega_gauss, ionode_id, world_comm)
  CALL mp_bcast(l_exchange, ionode_id, world_comm)
  CALL mp_bcast(l_zero, ionode_id, world_comm)
  CALL mp_bcast(l_wing, ionode_id, world_comm)
  CALL mp_bcast(grid_type, ionode_id, world_comm)
  CALL mp_bcast(nset_overlap, ionode_id, world_comm)
  CALL mp_bcast(nspace, ionode_id, world_comm)
  CALL mp_bcast(ecutoff_global, ionode_id, world_comm)
  CALL mp_bcast(maxiter2, ionode_id, world_comm)
  CALL mp_bcast(diago_thr2, ionode_id, world_comm)
  CALL mp_bcast(l_plot_mlwf, ionode_id, world_comm)
  CALL mp_bcast(l_pmatrix, ionode_id, world_comm)
  CALL mp_bcast(npcol, ionode_id, world_comm)
  CALL mp_bcast(nprow, ionode_id, world_comm)
  CALL mp_bcast(n_pola_lanczos, ionode_id, world_comm)
  CALL mp_bcast(n_self_lanczos, ionode_id, world_comm)
  CALL mp_bcast(nsteps_lanczos_pola, ionode_id, world_comm)
  CALL mp_bcast(nsteps_lanczos_self, ionode_id, world_comm)
  CALL mp_bcast(s_pola_lanczos, ionode_id, world_comm)
  CALL mp_bcast(s_self_lanczos, ionode_id, world_comm)
  CALL mp_bcast(s_g_lanczos, ionode_id, world_comm)
  CALL mp_bcast(l_pmat_diago, ionode_id, world_comm)
  CALL mp_bcast(pmat_ethr, ionode_id, world_comm)
  CALL mp_bcast(pmat_cutoff, ionode_id, world_comm)
  CALL mp_bcast(pmat_type, ionode_id, world_comm)
  CALL mp_bcast(lanczos_restart, ionode_id, world_comm)
  CALL mp_bcast(n_pola_lanczos_eff, ionode_id, world_comm)
  CALL mp_bcast(n_self_lanczos_eff, ionode_id, world_comm)
  CALL mp_bcast(n_pmat, ionode_id, world_comm)
  CALL mp_bcast(s_pmat, ionode_id, world_comm)
  CALL mp_bcast(n_fast_pmat, ionode_id, world_comm)
  CALL mp_bcast(off_fast_pmat, ionode_id, world_comm)
  CALL mp_bcast(l_fast_pola, ionode_id, world_comm)
  CALL mp_bcast(l_v_basis, ionode_id, world_comm)
  CALL mp_bcast(v_cutoff, ionode_id, world_comm)
  CALL mp_bcast(l_iter_algorithm, ionode_id, world_comm)
  CALL mp_bcast(dual_pb, ionode_id, world_comm)
  CALL mp_bcast(dual_vt, ionode_id, world_comm)
  CALL mp_bcast(dual_vs, ionode_id, world_comm)
  CALL mp_bcast(l_t_wannier, ionode_id, world_comm)
  CALL mp_bcast(s_first_state, ionode_id, world_comm)
  CALL mp_bcast(s_last_state, ionode_id, world_comm)
  CALL mp_bcast(l_selfconsistent, ionode_id, world_comm)
  CALL mp_bcast(l_whole_s, ionode_id, world_comm)
  CALL mp_bcast(l_ts_eigen, ionode_id, world_comm)
  CALL mp_bcast(l_frac_occ, ionode_id, world_comm)
  CALL mp_bcast(num_nbndv_min, ionode_id, world_comm)
  CALL mp_bcast(l_cond_pol_base, ionode_id, world_comm)
  CALL mp_bcast(l_semicore, ionode_id, world_comm)
  CALL mp_bcast(n_semicore, ionode_id, world_comm)
  CALL mp_bcast(l_semicore_read, ionode_id, world_comm)
  CALL mp_bcast(l_verbose, ionode_id, world_comm)
  CALL mp_bcast(l_contour, ionode_id, world_comm)
  CALL mp_bcast(l_real, ionode_id, world_comm)
  CALL mp_bcast(exchange_fast_dual, ionode_id, world_comm)
  CALL mp_bcast(l_bse, ionode_id, world_comm)
  CALL mp_bcast(s_bse, ionode_id, world_comm)
  CALL mp_bcast(dual_bse, ionode_id, world_comm)
  CALL mp_bcast(l_big_system, ionode_id, world_comm)
  CALL mp_bcast(extra_pw_cutoff, ionode_id, world_comm)
  CALL mp_bcast(l_list, ionode_id, world_comm)
  CALL mp_bcast(l_scissor, ionode_id, world_comm)
  CALL mp_bcast(scissor, ionode_id, world_comm)
  CALL mp_bcast(l_full, ionode_id, world_comm)
  CALL mp_bcast(n_full, ionode_id, world_comm)
  CALL mp_bcast(l_simple, ionode_id, world_comm)

  call read_file 


#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  !
 

  call openfil_pw4gww


! read wave functions (direct access)
  call read_export(pp_file,kunittmp,uspp_spsi, ascii, single_file, raw)
  !

! after read_file everything is known
! realy?

  call summary()

!
! init some quantities ...
!
  CALL hinit0()
!
  if(lda_plus_u) then 
    CALL init_ns()
  endif

  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )

!-------------------------------------------------
! allocating wannier stuff (from init_run.f90)
!-----------------------------------------------------
  CALL allocate_wannier()

! This is something from hinit0.f90, qpointlist ????

!
  
!
! -----------------------------------------------------
! now calculating the first wannier stuff (first in non_scf.f90)
! -----------------------------------------------------

  if(l_verbose) write(stdout,*) 'To check, we print the KS eigenvalues:'
  FLUSH( stdout )
  !
  CALL print_ks_energies()
  !

!  IF(l_head .and. .not.gamma_only) THEN
!     write(stdout,*) 'BEFORE calculate_head'
!     FLUSH( stdout )
!     CALL calculate_head
!     write(stdout,*) 'AFTER calculate_head'
!     FLUSH( stdout )
!  ENDIF
  !

  IF(l_exchange) THEN
    IF(gamma_only) THEN
      call dft_exchange(num_nbndv(1),num_nbnds,nset) 
    ELSE
      !!! add this, since wk are used in dft_exchange_k
      !
      CALL weights  ( )
      !
      if(l_verbose) write(stdout,*) 'BEFORE dft_exchange_k'
      FLUSH( stdout )
      !call dft_exchange_k(num_nbndv,num_nbnds,ecutoff_global)
      if(l_verbose) write(stdout,*) 'AFTER dft_exchange_k'
      FLUSH( stdout )
    ENDIF
  ENDIF


 

  if(l_verbose) write(stdout,*) 'BEFORE produce_wannier_gamma'
  FLUSH( stdout )
  CALL produce_wannier_gamma
  if(l_verbose) write(stdout,*) 'AFTER produce_wannier_gamma'
  FLUSH( stdout )
!     ENDIF
 
!
!
!deallocate wannier stuff (clean_pw.f90)
!
  CALL deallocate_wannier()


  call stop_pp
  stop
end program gwl_punch
!
!-----------------------------------------------------------------------
subroutine read_export (pp_file,kunit,uspp_spsi, ascii, single_file, raw)
  !-----------------------------------------------------------------------
  !
  use iotk_module


  use kinds,          ONLY : DP 
  use pwcom  
  use gvecw,          ONLY : gcutw
  use control_flags,  ONLY : gamma_only  
  use becmod,         ONLY : bec_type, becp, calbec, &
                             allocate_bec_type, deallocate_bec_type
!  use symme,          ONLY : nsym, s, invsym, sname, irt, ftau
!  use symme,          ONLY : nsym, s, invsym, irt, ftau
!  use char,           ONLY : sname
! occhio sname is in symme which is now outside pwcom
  use  uspp,          ONLY : nkb, vkb
  use wavefunctions_module,  ONLY : evc
  use io_files,       ONLY : nd_nmbr, prefix, iunwfc, nwordwfc, iunsat, nwordatwfc
  use io_files,       ONLY : pseudo_dir, psfile
  use io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
  use mp_pools,       ONLY : nproc_pool, my_pool_id, intra_pool_comm, &
                             inter_pool_comm
  use mp,             ONLY : mp_sum, mp_max
  use mp_world,       ONLY : world_comm, nproc, mpime
  use ldaU,           ONLY : lda_plus_u
  USE basis,          ONLY : swfcatom

  implicit none

  CHARACTER(5), PARAMETER :: fmt_name="QEXPT"
  CHARACTER(5), PARAMETER :: fmt_version="1.1.0"

  integer, intent(in) :: kunit
  character(80), intent(in) :: pp_file
  logical, intent(in) :: uspp_spsi, ascii, single_file, raw

  integer :: i, j, k, ig, ik, ibnd, na, ngg,ig_, ierr
  integer, allocatable :: kisort(:)
  real(DP) :: xyz(3), tmp(3)
  integer :: ike, iks, npw_g, npwx_g, ispin, local_pw
  integer, allocatable :: ngk_g( : )
  integer, allocatable :: itmp_g( :, : )
  real(DP),allocatable :: rtmp_g( :, : )
  real(DP),allocatable :: rtmp_gg( : )
  integer, allocatable :: itmp1( : )
  integer, allocatable :: igwk( :, : )
  integer, allocatable :: l2g_new( : )
  integer, allocatable :: igk_l2g( :, : )
  integer, external :: global_kpoint_index

  real(DP) :: wfc_scal 
  logical :: twf0, twfm
  character(iotk_attlenx) :: attr
  complex(DP), allocatable :: sevc (:,:)

  write(stdout,*) "nkstot=", nkstot

  IF( nkstot > 0 ) THEN

     IF( ( kunit < 1 ) .OR. ( MOD( nkstot, kunit ) /= 0 ) ) &
       CALL errore( ' write_export ',' wrong kunit ', 1 )

     IF( ( nproc_pool > nproc ) .OR. ( MOD( nproc, nproc_pool ) /= 0 ) ) &
       CALL errore( ' write_export ',' nproc_pool ', 1 )

     iks = global_kpoint_index (nkstot, 1)
     ike = iks + nks - 1
     
  END IF

  write(stdout,*) "after first init"

  ! find out the global number of G vectors: ngm_g  
  ngm_g = ngm
  call mp_sum( ngm_g , intra_pool_comm )

  ! collect all G vectors across processors within the pools
  ! and compute their modules
  !
  allocate( itmp_g( 3, ngm_g ) )
  allocate( rtmp_g( 3, ngm_g ) )
  allocate( rtmp_gg( ngm_g ) )

  itmp_g = 0
  do  ig = 1, ngm
    itmp_g( 1, ig_l2g( ig ) ) = mill(1, ig )
    itmp_g( 2, ig_l2g( ig ) ) = mill(2, ig )
    itmp_g( 3, ig_l2g( ig ) ) = mill(3, ig )
  end do
  call mp_sum( itmp_g , intra_pool_comm )
  !
  ! here we are in crystal units
  rtmp_g(1:3,1:ngm_g) = REAL( itmp_g(1:3,1:ngm_g) )
  !
  ! go to cartesian units (tpiba)
  call cryst_to_cart( ngm_g, rtmp_g, bg , 1 )
  !
  ! compute squared moduli
  do  ig = 1, ngm_g 
     rtmp_gg(ig) = rtmp_g(1,ig)**2 + rtmp_g(2,ig)**2 + rtmp_g(3,ig)**2 
  enddo
  deallocate( rtmp_g )

  ! build the G+k array indexes
  allocate ( igk_l2g ( npwx, nks ) )
  allocate ( kisort( npwx ) )
  do ik = 1, nks
     kisort = 0
     npw = npwx
     call gk_sort (xk (1, ik+iks-1), ngm, g, gcutw, npw, kisort(1), g2kin)
     !
     ! mapping between local and global G vector index, for this kpoint
     !
     DO ig = 1, npw
        !
        igk_l2g(ig,ik) = ig_l2g( kisort(ig) )
        !
     END DO
     !
     igk_l2g( npw+1 : npwx, ik ) = 0
     !
     ngk (ik) = npw
  end do
  deallocate (kisort)

  ! compute the global number of G+k vectors for each k point
  allocate( ngk_g( nkstot ) )
  ngk_g = 0
  ngk_g( iks:ike ) = ngk( 1:nks )
  CALL mp_sum( ngk_g, world_comm )

  ! compute the Maximum G vector index among all G+k and processors
  npw_g = MAXVAL( igk_l2g(:,:) )
  CALL mp_max( npw_g, world_comm )

  ! compute the Maximum number of G vector among all k points
  npwx_g = MAXVAL( ngk_g( 1:nkstot ) )

  deallocate(rtmp_gg)

  allocate( igwk( npwx_g,nkstot ) )

  write(stdout,*) "after g stuff"

! wfc grids

  DO ik = 1, nkstot
    igwk(:,ik) = 0
    !
    ALLOCATE( itmp1( npw_g ), STAT= ierr )
    IF ( ierr/=0 ) CALL errore('pw_export','allocating itmp1', ABS(ierr) )
    itmp1 = 0
    ! 
    IF( ik >= iks .AND. ik <= ike ) THEN 
      DO  ig = 1, ngk( ik-iks+1 )
        itmp1( igk_l2g( ig, ik-iks+1 ) ) = igk_l2g( ig, ik-iks+1 ) 
      END DO
    END IF
    !
    CALL mp_sum( itmp1, world_comm )
    !
    ngg = 0
    DO  ig = 1, npw_g
      IF( itmp1( ig ) == ig ) THEN
        ngg = ngg + 1
        igwk( ngg , ik) = ig
      END IF
    END DO
    IF( ngg /= ngk_g( ik ) ) THEN
      WRITE( stdout,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
    END IF
    !
    DEALLOCATE( itmp1 )
    !
  ENDDO
  !
  deallocate( itmp_g )
  
  write(stdout,*)"after wfc waves"

#if defined(__MPI)
  call poolrecover (et, nbnd, nkstot, nks)
#endif
 
  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.

  do ik = 1, nkstot
     local_pw = 0
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN

       call davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), - 1)
       IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunsat, (ik-iks+1), -1 )
       local_pw = ngk(ik-iks+1)

     ENDIF


     allocate(l2g_new(local_pw))

     l2g_new = 0
     do ig = 1, local_pw
       ngg = igk_l2g(ig,ik-iks+1)
       do ig_ = 1, ngk_g(ik)
         if(ngg == igwk(ig_,ik)) then
           l2g_new(ig) = ig_
           exit
         end if
       end do
     end do


     ispin = isk( ik )
     !  WRITE(0,*) ' ### ', ik,nkstot,iks,ike,kunit,nproc,nproc_pool 
     deallocate(l2g_new)
  end do
  !  

  write(stdout,*) "after davcio"

  ! If specified and if USPP are used the wfcs S_psi are written  
  ! | spsi_nk > = \hat S | psi_nk >  
  ! where S is the overlap operator of US PP 
  !  
  IF ( uspp_spsi .AND. nkb > 0 ) THEN

       ALLOCATE( sevc(npwx,nbnd), STAT=ierr )
       IF (ierr/=0) CALL errore( ' read_export ',' Unable to allocate SEVC ', ABS(ierr) )

       CALL init_us_1
       CALL init_at_1

       CALL allocate_bec_type (nkb,nbnd,becp)

       do ik = 1, nkstot
 
           local_pw = 0
           IF( (ik >= iks) .AND. (ik <= ike) ) THEN
               
               CALL davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), - 1)
               CALL init_us_2(npw, igk_k(1,ik-iks+1), xk(1, ik), vkb)
               local_pw = ngk(ik-iks+1)
                            
               IF ( gamma_only ) THEN
                  if(nkb>0) CALL calbec ( ngk_g(ik), vkb, evc, becp )
               ELSE
                  CALL calbec ( npw, vkb, evc, becp )
               ENDIF
               CALL s_psi(npwx, npw, nbnd, evc, sevc)
           ENDIF

           ALLOCATE(l2g_new(local_pw))

           l2g_new = 0
           DO ig = 1, local_pw
             ngg = igk_l2g(ig,ik-iks+1)
             DO ig_ = 1, ngk_g(ik)
               IF(ngg == igwk(ig_,ik)) THEN
                 l2g_new(ig) = ig_
                 EXIT
               ENDIF
             ENDDO
           ENDDO

           ispin = isk( ik )
           DEALLOCATE(l2g_new)
       ENDDO
      
       DEALLOCATE( sevc, STAT=ierr )
       IF ( ierr/= 0 ) CALL errore('read_export','Unable to deallocate SEVC',ABS(ierr))
       CALL deallocate_bec_type ( becp )
  ENDIF

  DEALLOCATE( igk_l2g )
  DEALLOCATE( igwk )
  DEALLOCATE ( ngk_g )

end subroutine read_export
