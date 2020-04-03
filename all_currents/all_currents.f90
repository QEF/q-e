!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !!!!!!!!!!!!    THERMOCURRENTS    !!!!!!!!!!!!
!    !!! calculation of energy currents with QE !!!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  written by:
!     Aris Marcolongo
!  partially updated and modified by:
!     Loris Ercole, Federico Grasselli, Pietro Delugas
!
!  at SISSA, Via Bonomea 265, 34136 Trieste, Italy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This program is not part of the official QE release yet.
!  Please do not redistribute.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
!  use grid_dimensions, ONLY : nrxx
!changed
  USE fft_base,        ONLY : dfftp
  USE constants,       ONLY : rytoev
  use io_global, ONLY : stdout, ionode, ionode_id
  use io_files,  ONLY : psfile, pseudo_dir, diropn
  use io_files,  ONLY : prefix, tmp_dir
  use input_parameters, ONLY : outdir 
  use ions_base, ONLY : ntype => nsp
  use iotk_module
  USE environment, ONLY : environment_start, environment_end
  use mp_global, ONLY : kunit, mp_startup
  use mp_world, ONLY : mpime
  use mp, ONLY: mp_bcast, mp_barrier
  use mp_pools, ONLY : intra_pool_comm
  use control_flags, ONLY : gamma_only
  use realus, ONLY : qpointlist
  use uspp, ONLY : okvan
  use ldaU, ONLY : lda_plus_u
  use scf, only : vrs, vltot, v, kedtau
  use klist,                ONLY : xk, wk, nks, nkstot
  use zero_mod
  use hartree_mod
  use wavefunctions_module,  ONLY : evc
  use io_files,       ONLY : nd_nmbr, tmp_dir, prefix, nwordwfc, iunwfc
  use pw_restart_new, ONLY : read_collected_to_evc
  USE wrappers,      ONLY : f_mkdir_safe
  !
  implicit none
  integer :: i, kunittmp, ios
  integer ::handle
  logical :: exst

  character(len=200) :: pp_file, dirname
  character(len=iotk_attlenx) :: attr
  logical :: found, uspp_spsi, ascii, single_file, raw
!  INTEGER(i4b), EXTERNAL :: C_MKDIR
   CHARACTER(LEN=256), EXTERNAL :: trimcheck
   INTEGER, EXTERNAL            :: find_free_unit
!
  INTERFACE
    subroutine diropn_due(pref, unit, extension, recl, exst, tmp_dir_)
       implicit none
       character(len=*),intent(in)   :: pref, extension
       character(len=*), optional :: tmp_dir_
       integer            :: unit,recl
       logical            :: exst 
    end subroutine 
  END INTERFACE 

  NAMELIST /inputhartree/ prefix_uno, prefix_due, delta_t, init_linear, &
file_output, file_dativel,outdir, thermodir

  NAMELIST /inputzero/ eta, n_max,status,l_zero
 

!  call start_all_currents()
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'HARTRIS' )
  !
  !   set default values for variables in namelist
  !
  prefix='export'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  thermodir='./'
  pp_file= ' '
  uspp_spsi = .FALSE.
  ascii = .FALSE.
  single_file = .FALSE.
  raw = .FALSE.
  prefix_due='pwscf'
  delta_t=1.d0
  n_max=5
  eta=1.0
  status="undefined"
  init_linear="scratch"
  file_output="corrente_def"
  file_dativel="velocita_def"

!
!       nppstr   = 1
!
    
  !    Reading input file
  !
  IF ( ionode ) THEN
      !
      CALL input_from_file ( )
      !
      READ(5,inputhartree,IOSTAT=ios)
      !
      READ(5,inputzero,IOSTAT=ios)
      !
!      call read_namelists( 'PW4GWW' )
      !
      IF (ios /= 0) CALL errore ('main', 'reading inputp namelists', ABS(ios) )
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
  call mp_barrier(intra_pool_comm)
!-------------------------------------------------------------------------
  ! ... Broadcasting variables
!------------------------------------------------------------------------
  tmp_dir = trimcheck( outdir )
  thermodir = trimcheck( thermodir )
  ! ... create thermodir
  IF ( ionode ) ios = f_mkdir_safe( TRIM(thermodir) )
  IF ( (ios==-1) > 0 ) CALL errore ('check_thermodir','thermodir cannot be opened',1)
  if (ionode) print *, "    thermodir = "// thermodir

  CALL mp_bcast( tmp_dir,     ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( thermodir,    ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( prefix_uno,  ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( prefix_due,  ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( delta_t,     ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( eta,         ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( n_max,       ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( status,      ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( init_linear, ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( file_output, ionode_id, MPI_COMM_WORLD )
  CALL mp_bcast( file_dativel,ionode_id, MPI_COMM_WORLD )
  !
  prefix = trim(prefix_uno)
  call read_file
  dirname = TRIM(tmp_dir) //TRIM(prefix_due) //'.save'//'/'
  iunwfc  = 1000+iunwfc
  call diropn_due(trim(prefix_due), iunwfc, 'wfc', 2 * nwordwfc, exst, tmp_dir ) 
  
  call read_collected_to_evc( dirname) 
  IF ( ionode ) WRITE( stdout, '(/,5x,A,/,5x,A)') &
     'Reading data from directory:'// TRIM( dirname )
  close(iunwfc)
  iunwfc = iunwfc -1000

#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  !
  call openfil_all_currents

! read wave functions (direct access)
  call read_export(pp_file,kunittmp,uspp_spsi, ascii, single_file, raw)
  !

! after read_file everything is known
! realy?

  call summary()
  call setup_nbnd_occ()
!
! init some quantities igk,....
!
  CALL hinit0()
!
  if(lda_plus_u) then
    CALL init_ns()
  endif

!changed
  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, .false. )


! This is something from hinit0.f90, qpointlist ????

!
!  write(stdout,*) 'PRIMA QPOINT',l_exchange, okvan!ATTENZIONE
!  IF ( (lwannier .or. l_head .or. l_exchange).and. okvan) CALL qpointlist()
!
! -----------------------------------------------------
! now calculating the first wannier stuff (first in non_scf.f90)
! -----------------------------------------------------

  write(stdout,*) 'To check, we print the KS eigenvalues:'
  !CALL flush_unit( stdout )
  !
  CALL print_ks_energies()
  !
!  handle=find_free_unit()
!  CALL diropn ( handle, 'evc', 2*nwordwfc, exst )
!  CALL davcio (evc, 2*nwordwfc, handle, 1, + 1)

!  do ipw=1,npw
!     do inbd=1,nbnd
!        CALL davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), + 1)
!     end do
!  end do
!  IF(l_head .and. .not.gamma_only) THEN
!     write(stdout,*) 'BEFORE calculate_head'
!     CALL flush_unit( stdout )
!     CALL calculate_head
!     write(stdout,*) 'AFTER calculate_head'
!     CALL flush_unit( stdout )
!  ENDIF
  !
  
  call start_clock('all_currents')
  if (status=='compute') then
      call routine_hartree()
  end if
! 
  call routine_zero()
!
  call stop_clock('all_currents')
  call print_clock('all_currents')
  !
  call environment_end ( 'HARTRIS' )
  !
  call stop_pp()
  !
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
  use gvect,          ONLY : ngm, ngm_g, mill, ig_l2g
  use control_flags,  ONLY : gamma_only
  use becmod,         ONLY : bec_type, becp, calbec, &
                             allocate_bec_type, deallocate_bec_type
  use symm_base,       ONLY : nsym, s, invsym, irt, ftau
  use  uspp,          ONLY : nkb, vkb
  use wavefunctions_module,  ONLY : evc
  use io_files,       ONLY : nd_nmbr, tmp_dir, prefix, iunwfc, nwordwfc, iunsat, nwordatwfc
  use io_files,       ONLY : pseudo_dir, psfile
  use io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
  use mp_world,       ONLY : nproc, mpime
  use mp_global,      ONLY : nproc_pool, my_pool_id, intra_pool_comm, inter_pool_comm
  use mp,             ONLY : mp_sum, mp_max
  use ldaU,           ONLY : lda_plus_u
  use basis,          ONLY : swfcatom 
  use gvect,          ONLY : ig_l2g
  use cell_base,      ONLY : bg, tpiba2
  use gvecw,          ONLY : ecutwfc
  use gvect,          ONLY : g
  use klist,          ONLY : igk_k

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
  CALL mp_sum( ngk_g, intra_pool_comm )

  ! compute the Maximum G vector index among all G+k and processors
  npw_g = MAXVAL( igk_l2g(:,:) )
  CALL mp_max( npw_g , intra_pool_comm)

  ! compute the Maximum number of G vector among all k points
  npwx_g = MAXVAL( ngk_g( 1:nkstot )  )

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
    CALL mp_sum( itmp1, intra_pool_comm )
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

#ifdef __MPI
  call poolrecover (et, nbnd, nkstot, nks)
#endif

  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.

  do ik = 1, nkstot
     local_pw = 0
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN
!*mio*
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

               CALL gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik), g2kin)
!
               CALL davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), - 1)

               CALL init_us_2(npw, igk_k(1,ik), xk(1, ik), vkb)
               local_pw = ngk(ik-iks+1)

               IF ( gamma_only ) THEN
                  CALL calbec ( npw, vkb, evc, becp )
               ELSE
                  CALL calbec ( ngk_g(ik), vkb, evc, becp )
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
