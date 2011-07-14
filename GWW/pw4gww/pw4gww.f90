! Copyright (C) 2003-2006 Andrea Ferretti and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
program pp_punch
  !-----------------------------------------------------------------------
  !
  ! writes PWSCF data for postprocessing purposes in XML format using IOTK lib
  ! Wave-functions are collected and written using IO_BASE module.
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
  use pwcom
  use fft_base,  ONLY : dfftp
  USE constants,       ONLY : rytoev
  use io_global, ONLY : stdout, ionode, ionode_id
  use io_files,  ONLY : psfile, pseudo_dir
  use io_files,  ONLY : prefix, tmp_dir, outdir
  use ions_base, ONLY : ntype => nsp
  use iotk_module
  use mp_global, ONLY : mpime, kunit
  use mp, ONLY: mp_bcast
  use control_flags, ONLY : gamma_only
  use realus, ONLY : qpointlist
  use uspp, ONLY : okvan
  use ldaU, ONLY : lda_plus_u
  use scf, only : vrs, vltot, v, kedtau
  USE klist,                ONLY : xk, wk, nks, nkstot
  USE wannier_gw,       ONLY :  lwannier, &
                                cutoff_wsq, &
                                cutoff_wsq_c, &
                                cutoff_wpr, &
                                cutoff_wpr_wpr, &
                                cutoff_overlap,&
                                num_nbndv, &
                                num_nbnds, &
                                lsmallgrid, &
                                lnonorthogonal, &
                                no_radius, &
                                lggrid, &
                                nset, &
                                ultra_alpha_v, &
                                ultra_alpha_c, &
                                ultra_alpha_c2, &
                                ultra_alpha_c_prim, &
                                num_nbndc_set, &
                                l_truncated_coulomb, &
                                truncation_radius, &
                                r_cutoff_products, &
                                remainder, &
                                restart_gww, &
                                numw_prod, &
                                cutoff_wpr_vc,&
                                cutoff_wpr_vc2,&
                                num_nbnd_first,&
                                cutoff_wpr_prim,&
                                cutoff_wpr_prim2,&
                                l_gram,&
                                l_head,&
                                n_gauss,&
                                omega_gauss, &
                                l_exchange, &
                                tau_gauss, &
                                l_zero, &
                                l_wing, &
                                grid_type, &
                                cprim_type, &
                                cprim_first, &
                                cprim_last, &
                                l_vcw_overlap, &
                                nset_overlap, &
                                nspace,&
                                lambda_ene,&
                                e_min_cutoff, &
                                e_max_cutoff, &
                                v_min_cutoff, &
                                v_max_cutoff, &
                                l_orthonorm_products, &
                                cutoff_products, &
                                ecutoff_global, &
                                l_wpwt_terms,&
                                l_polarization_analysis,&
                                cutoff_polarization,&
                                nc_polarization_analysis,&
                                l_only_val_cond,&
                                l_no_val_cond_sec,&
                                maxiter2,&
                                diago_thr2, &
                                l_plot_mlwf,&
                                l_plot_ulwf,&
                                l_ultra_external,&
                                nbnd_normal,&
                                num_nbnd_delta,&
                                num_nbnd_upper,&
                                l_pmatrix,&
                                npcol,&
                                nprow,&
                                l_assume_ortho,&
                                l_coulomb_analysis,&
                                cutoff_coulomb_analysis,&
                                mem_per_core

  !
  implicit none
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  integer :: i, kunittmp, ios

  character(len=200) :: pp_file
  character(len=iotk_attlenx) :: attr
  logical :: found, uspp_spsi, ascii, single_file, raw
!  INTEGER(i4b), EXTERNAL :: C_MKDIR

  NAMELIST /inputpw4gww/ prefix, outdir, pp_file, uspp_spsi, ascii, single_file, raw, &
                     psfile, pseudo_dir, &
                     lwannier, cutoff_wsq, cutoff_wpr, num_nbndv, &
                               lsmallgrid, lnonorthogonal, no_radius, cutoff_overlap,&
                               cutoff_wsq_c,lggrid, nset,num_nbnds, ultra_alpha_v, &
                               ultra_alpha_c, num_nbndc_set, l_truncated_coulomb, &
                               truncation_radius, cutoff_wpr_wpr, r_cutoff_products, &
                               remainder, ultra_alpha_c_prim, restart_gww, numw_prod, &
                               cutoff_wpr_vc,cutoff_wpr_vc2,num_nbnd_first, cutoff_wpr_prim, &
                               l_gram, ultra_alpha_c2, l_head, n_gauss, omega_gauss, l_exchange, &
                               tau_gauss,l_zero,cutoff_wpr_prim2, l_wing, grid_type, cprim_type, &
                               cprim_first,cprim_last, l_vcw_overlap, nset_overlap, nspace, &
                               lambda_ene, e_min_cutoff, e_max_cutoff, v_min_cutoff, v_max_cutoff, &
                               l_orthonorm_products,cutoff_products,ecutoff_global,&
                               l_wpwt_terms, l_polarization_analysis,&
                               cutoff_polarization, nc_polarization_analysis,l_only_val_cond,&
                               l_no_val_cond_sec,maxiter2,diago_thr2,l_plot_mlwf,l_plot_ulwf,&
                               l_ultra_external,nbnd_normal,num_nbnd_delta,num_nbnd_upper,&
                               l_pmatrix, npcol,nprow,l_assume_ortho,l_coulomb_analysis,&
                               cutoff_coulomb_analysis, mem_per_core

  !
  call start_pw4gww( )
  !
  !   set default values for variables in namelist
  !
  prefix='export'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  pp_file= ' '
  uspp_spsi = .FALSE.
  ascii = .FALSE.
  single_file = .FALSE.
  raw = .FALSE.
!
!       nppstr   = 1
!
       lwannier = .true.
       cutoff_wsq = 0.9999d0
       cutoff_wsq_c = 0.9999d0
       cutoff_wpr = 1000d0
       cutoff_overlap = 0.d0
       cutoff_wpr_wpr = 0.02
       num_nbndv = 1
       num_nbnds = 1
       lsmallgrid = .false.
       lnonorthogonal = .true.
       lggrid=.true.
       no_radius = 8.d0
       nset = 200
       ultra_alpha_v = 1000d0
       ultra_alpha_c = 1000d0
       ultra_alpha_c2 = 1000d0
       ultra_alpha_c_prim = 1000d0
       num_nbndc_set = 0
       l_truncated_coulomb = .false.
       truncation_radius = 10.d0
       r_cutoff_products = 10.d0
       remainder=-1
       restart_gww=0
       numw_prod=1
       cutoff_wpr_vc=0.05d0
       cutoff_wpr_vc2=1000d0
       cutoff_wpr_prim=0.0d0
       cutoff_wpr_prim2=0.0d0
       num_nbnd_first=0
       l_gram=.false.
       l_head=.false.
       l_exchange=.true.
       n_gauss=79
       omega_gauss=20.d0
       tau_gauss=10.d0
       l_zero=.false.
       l_wing=.false.
       grid_type=2
       cprim_type=2
       cprim_first=1
       cprim_last=1
       l_vcw_overlap=.true.
       nset_overlap=1000
       nspace=1
       lambda_ene=0.d0
       e_min_cutoff = 0.d0
       e_max_cutoff = 50.d0
       v_min_cutoff = 0.d0
       v_max_cutoff = 10.d0
       l_orthonorm_products = .true.
       cutoff_products = 0.01d0
       ecutoff_global = 80.d0
       l_wpwt_terms=.true.
       l_polarization_analysis = .true.
       cutoff_polarization = 0.1d0
       nc_polarization_analysis = 10000
       l_only_val_cond = .true.
       l_no_val_cond_sec = .true.
       maxiter2=0
       diago_thr2=0.d0
       l_plot_mlwf=.false.
       l_plot_ulwf=.false.
       l_ultra_external=.false.
       nbnd_normal=0
       num_nbnd_delta=100
       num_nbnd_upper=10
       l_pmatrix=.false.
       npcol=1
       nprow=1
       l_assume_ortho=.true.
       l_coulomb_analysis=.false.
       cutoff_coulomb_analysis=0.d0
       mem_per_core=1600000000
!

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
!-------------------------------------------------------------------------
  ! ... Broadcasting variables
!------------------------------------------------------------------------
  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( outdir, ionode_id )
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( pp_file, ionode_id )
  CALL mp_bcast( uspp_spsi, ionode_id )
  CALL mp_bcast( ascii, ionode_id )
  CALL mp_bcast( single_file, ionode_id )
  CALL mp_bcast( raw, ionode_id )
  CALL mp_bcast( pseudo_dir, ionode_id )
  CALL mp_bcast( psfile, ionode_id )


  !
!      CALL mp_bcast( nppstr,        ionode_id )
!
       CALL mp_bcast( lwannier,      ionode_id )
       CALL mp_bcast( cutoff_wsq,    ionode_id )
       CALL mp_bcast( cutoff_wsq_c,    ionode_id )
       CALL mp_bcast( cutoff_wpr,    ionode_id )
       CALL mp_bcast( cutoff_wpr_wpr,    ionode_id )
       CALL mp_bcast( cutoff_overlap,ionode_id )
       CALL mp_bcast( num_nbndv,     ionode_id )
       CALL mp_bcast( num_nbnds,     ionode_id )
       CALL mp_bcast( lsmallgrid,    ionode_id )
       CALL mp_bcast( lggrid,        ionode_id )
       CALL mp_bcast( lnonorthogonal,ionode_id )
       CALL mp_bcast( no_radius,     ionode_id )
       CALL mp_bcast( nset,          ionode_id )
       CALL mp_bcast(ultra_alpha_c,  ionode_id )
       CALL mp_bcast(ultra_alpha_c2,  ionode_id )
       CALL mp_bcast(ultra_alpha_v,  ionode_id )
       CALL mp_bcast(ultra_alpha_c_prim, ionode_id)
       CALL mp_bcast(num_nbndc_set,  ionode_id )
       CALL mp_bcast(l_truncated_coulomb, ionode_id)
       CALL mp_bcast(truncation_radius, ionode_id)
       CALL mp_bcast(r_cutoff_products, ionode_id)
       CALL mp_bcast(remainder, ionode_id)
       CALL mp_bcast(restart_gww, ionode_id)
       call mp_bcast(numw_prod, ionode_id)
       call mp_bcast(cutoff_wpr_vc,  ionode_id)
       call mp_bcast(cutoff_wpr_vc2,  ionode_id)
       call mp_bcast(num_nbnd_first,  ionode_id)
       call mp_bcast(cutoff_wpr_prim,  ionode_id)
       call mp_bcast(cutoff_wpr_prim2,  ionode_id)
       CALL mp_bcast(l_gram, ionode_id)
       CALL mp_bcast(l_head, ionode_id)
       CALL mp_bcast(n_gauss, ionode_id)
       CALL mp_bcast(omega_gauss, ionode_id)
       CALL mp_bcast(l_exchange, ionode_id)
       CALL mp_bcast(tau_gauss, ionode_id)
       CALL mp_bcast(l_zero, ionode_id)
       CALL mp_bcast(l_wing, ionode_id)
       CALL mp_bcast(grid_type, ionode_id)
       CALL mp_bcast(cprim_type, ionode_id)
       CALL mp_bcast(cprim_first, ionode_id)
       CALL mp_bcast(cprim_last, ionode_id)
       CALL mp_bcast(l_vcw_overlap, ionode_id)
       CALL mp_bcast(nset_overlap, ionode_id)
       CALL mp_bcast(nspace, ionode_id)
       CALL mp_bcast(lambda_ene, ionode_id)
       CALL mp_bcast(e_min_cutoff, ionode_id)
       CALL mp_bcast(e_max_cutoff, ionode_id)
       CALL mp_bcast(v_min_cutoff, ionode_id)
       CALL mp_bcast(v_max_cutoff, ionode_id)
       CALL mp_bcast(l_orthonorm_products, ionode_id)
       CALL mp_bcast(cutoff_products, ionode_id)
       CALL mp_bcast(ecutoff_global, ionode_id)
       CALL mp_bcast(l_wpwt_terms, ionode_id)
       CALL mp_bcast(l_polarization_analysis, ionode_id)
       CALL mp_bcast(cutoff_polarization, ionode_id)
       CALL mp_bcast(nc_polarization_analysis, ionode_id)
       CALL mp_bcast(l_only_val_cond, ionode_id)
       CALL mp_bcast(l_no_val_cond_sec, ionode_id)
       CALL mp_bcast(maxiter2, ionode_id)
       CALL mp_bcast(diago_thr2, ionode_id)
       CALL mp_bcast(l_plot_mlwf, ionode_id)
       CALL mp_bcast(l_plot_ulwf, ionode_id)
       CALL mp_bcast(l_ultra_external, ionode_id)
       CALL mp_bcast(nbnd_normal, ionode_id)
       CALL mp_bcast(num_nbnd_delta, ionode_id)
       CALL mp_bcast(num_nbnd_upper, ionode_id)
       CALL mp_bcast(l_pmatrix, ionode_id)
       CALL mp_bcast(npcol, ionode_id)
       CALL mp_bcast(nprow, ionode_id)
       CALL mp_bcast(l_assume_ortho, ionode_id)
       CALL mp_bcast(l_coulomb_analysis, ionode_id)
       CALL mp_bcast(cutoff_coulomb_analysis, ionode_id)
       CALL mp_bcast(mem_per_core, ionode_id)
  !

  call read_file


#if defined __PARA
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
! init some quantities igk,....
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
  write(stdout,*) 'PRIMA QPOINT',l_exchange, okvan!ATTENZIONE
  IF ( (lwannier .or. l_head .or. l_exchange).and. okvan) CALL qpointlist()
!
! -----------------------------------------------------
! now calculating the first wannier stuff (first in non_scf.f90)
! -----------------------------------------------------

  write(stdout,*) 'To check, we print the KS eigenvalues:'
  CALL flush_unit( stdout )
  !
  CALL print_ks_energies()
  !

!  IF(l_head .and. .not.gamma_only) THEN
!     write(stdout,*) 'BEFORE calculate_head'
!     CALL flush_unit( stdout )
!     CALL calculate_head
!     write(stdout,*) 'AFTER calculate_head'
!     CALL flush_unit( stdout )
!  ENDIF
  !

  IF(l_exchange) THEN
    IF(gamma_only) THEN
      call dft_exchange(num_nbndv,num_nbnds,nset)
    ELSE
      !!! add this, since wk are used in dft_exchange_k
      !
      CALL weights  ( )
      !
      write(stdout,*) 'BEFORE dft_exchange_k'
      CALL flush_unit( stdout )
      call dft_exchange_k(num_nbndv,num_nbnds,ecutoff_global)
      write(stdout,*) 'AFTER dft_exchange_k'
      CALL flush_unit( stdout )
    ENDIF
  ENDIF


  IF(lwannier) THEN
!     IF(.not.gamma_only) THEN  ! not yet implemented
!       CALL produce_wannier
!     ELSE
        write(stdout,*) 'BEFORE produce_wannier_gamma'
        CALL flush_unit( stdout )
        CALL produce_wannier_gamma
        write(stdout,*) 'AFTER produce_wannier_gamma'
        CALL flush_unit( stdout )
!     ENDIF
  ENDIF
!
!
!deallocate wannier stuff (clean_pw.f90)
!
  CALL deallocate_wannier()
!  call write_export (pp_file, kunittmp, uspp_spsi, ascii, single_file, raw)

  call stop_pp
  stop
end program pp_punch
!
!-----------------------------------------------------------------------
subroutine read_export (pp_file,kunit,uspp_spsi, ascii, single_file, raw)
  !-----------------------------------------------------------------------
  !
  use iotk_module


  use kinds,          ONLY : DP
  use pwcom
  use control_flags,  ONLY : gamma_only
  use becmod,         ONLY : bec_type, becp, calbec, &
                             allocate_bec_type, deallocate_bec_type
  use symm_base,       ONLY : nsym, s, invsym, irt, ftau
  use  uspp,          ONLY : nkb, vkb
  use wavefunctions_module,  ONLY : evc
  use io_files,       ONLY : nd_nmbr, outdir, prefix, iunwfc, nwordwfc, iunsat, nwordatwfc
  use io_files,       ONLY : pseudo_dir, psfile
  use io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
  use mp_global,      ONLY : nproc, nproc_pool, mpime
  use mp_global,      ONLY : my_pool_id, intra_pool_comm, inter_pool_comm
  use mp,             ONLY : mp_sum, mp_max
  use ldaU,           ONLY : swfcatom, lda_plus_u

  implicit none

  CHARACTER(5), PARAMETER :: fmt_name="QEXPT"
  CHARACTER(5), PARAMETER :: fmt_version="1.1.0"

  integer, intent(in) :: kunit
  character(80), intent(in) :: pp_file
  logical, intent(in) :: uspp_spsi, ascii, single_file, raw

  integer :: i, j, k, ig, ik, ibnd, na, ngg,ig_, ierr
  integer, allocatable :: kisort(:)
  real(DP) :: xyz(3), tmp(3)
  integer :: npool, nkbl, nkl, nkr, npwx_g
  integer :: ike, iks, npw_g, ispin, local_pw
  integer, allocatable :: ngk_g( : )
  integer, allocatable :: itmp_g( :, : )
  real(DP),allocatable :: rtmp_g( :, : )
  real(DP),allocatable :: rtmp_gg( : )
  integer, allocatable :: itmp1( : )
  integer, allocatable :: igwk( :, : )
  integer, allocatable :: l2g_new( : )
  integer, allocatable :: igk_l2g( :, : )


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

     !  find out the number of pools
     npool = nproc / nproc_pool

     !  find out number of k points blocks
     nkbl = nkstot / kunit

     !  k points per pool
     nkl = kunit * ( nkbl / npool )

     !  find out the reminder
     nkr = ( nkstot - nkl * npool ) / kunit

     !  Assign the reminder to the first nkr pools
     IF( my_pool_id < nkr ) nkl = nkl + kunit

     !  find out the index of the first k point in this pool
     iks = nkl * my_pool_id + 1
     IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

     !  find out the index of the last k point in this pool
     ike = iks + nkl - 1

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
    itmp_g( 1, ig_l2g( ig ) ) = mill(1,ig )
    itmp_g( 2, ig_l2g( ig ) ) = mill(2,ig )
    itmp_g( 3, ig_l2g( ig ) ) = mill(3,ig )
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
     call gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, kisort(1), g2kin)
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
  CALL mp_sum( ngk_g )

  ! compute the Maximum G vector index among all G+k and processors
  npw_g = MAXVAL( igk_l2g(:,:) )
  CALL mp_max( npw_g )

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
    CALL mp_sum( itmp1 )
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

#ifdef __PARA
  call poolrecover (et, nbnd, nkstot, nks)
#endif

  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.

  do ik = 1, nkstot
     local_pw = 0
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN

       call davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)
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

               CALL gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
               CALL davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)

               CALL init_us_2(npw, igk, xk(1, ik), vkb)
               local_pw = ngk(ik-iks+1)

               IF ( gamma_only ) THEN
                  CALL calbec ( ngk_g(ik), vkb, evc, becp )
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
