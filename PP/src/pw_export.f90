!
! Copyright (C) 2003-2009 Andrea Ferretti and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE io_base_export
!=----------------------------------------------------------------------------=!

! do i = 1, nk             !                                                   !
!   WAVEFUNCTIONS( i )     !  write_restart_wfc         read_restart_wfc       !
! end do                   !                                                   !

  USE io_global,  ONLY : stdout
  USE kinds
  USE parameters, ONLY: nsx

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: file_version = 202
  INTEGER :: restart_module_verbosity = 0

  INTERFACE write_restart_wfc
    MODULE PROCEDURE write_restart_wfc1, write_restart_wfc2
  END INTERFACE


!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write wavefunctions to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_wfc1(iuni, &
      ik, nk, kunit, ispin, nspin, scal, wf0, t0, wfm, tm, ngw, gamma_only, nbnd, igl, ngwl )
!
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_get, mp_max
      USE mp_pools, ONLY: me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool, npool
      USE mp_world,  ONLY: mpime, nproc, root, world_comm
      USE io_global, ONLY: ionode, ionode_id
      USE iotk_module
!
      IMPLICIT NONE
!
      INTEGER, INTENT(in) :: iuni
      INTEGER, INTENT(in) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP), INTENT(in) :: wf0(:,:)
      COMPLEX(DP), INTENT(in) :: wfm(:,:)
      INTEGER, INTENT(in) :: ngw   !
      LOGICAL, INTENT(in) :: gamma_only
      INTEGER, INTENT(in) :: nbnd
      INTEGER, INTENT(in) :: ngwl
      INTEGER, INTENT(in) :: igl(:)
      REAL(DP), INTENT(in) :: scal
      LOGICAL, INTENT(in) :: t0, tm

      INTEGER :: i, j, ierr, idum = 0
      INTEGER :: iks, ike, nkt, ikt, igwx
      INTEGER :: ipmask( nproc ), ipsour
      INTEGER, EXTERNAL :: global_kpoint_index
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)

      CHARACTER(len=20) :: section_name = 'wfc'

      LOGICAL :: twrite = .true.

      INTEGER :: ierr_iotk
      CHARACTER(len=iotk_attlenx) :: attr

!
! ... Subroutine Body
!

        ! set working variables for k point index (ikt) and k points number (nkt)
        ikt = ik
        nkt = nk

        iks = global_kpoint_index (nkt, 1)
        ike = iks + nk - 1

        ipmask = 0
        ipsour = ionode_id

        !  find out the index of the processor which collect the data in the pool of ik
        IF( npool > 1 ) THEN
          IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
            IF( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
          ENDIF
          CALL mp_sum( ipmask, world_comm )
          DO i = 1, nproc
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
          ENDDO
        ENDIF

        igwx = 0
        ierr = 0
        IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
          IF( ngwl > size( igl ) ) THEN
            ierr = 1
          ELSE
            igwx = maxval( igl(1:ngwl) )
          ENDIF
        ENDIF

        ! get the maximum index within the pool
        !
        CALL mp_max( igwx, intra_pool_comm )

        ! now notify all procs if an error has been found
        !
        CALL mp_max( ierr, world_comm )

        IF( ierr > 0 ) &
          CALL errore(' write_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipsour /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipsour, 1, world_comm )
        ENDIF

        IF(ionode) THEN
          CALL iotk_write_begin(iuni,"Kpoint"//iotk_index(ik))
          CALL iotk_write_attr (attr,"ngw",ngw,first=.true.)
          CALL iotk_write_attr (attr,"nbnd",nbnd)
          CALL iotk_write_attr (attr,"gamma_only",gamma_only)
          CALL iotk_write_attr (attr,"ik",ik)
          CALL iotk_write_attr (attr,"nk",nk)
          CALL iotk_write_attr (attr,"kunit",kunit)
          CALL iotk_write_attr (attr,"ispin",ispin)
          CALL iotk_write_attr (attr,"nspin",nspin)
          CALL iotk_write_attr (attr,"scal",scal)
          CALL iotk_write_attr (attr,"igwx",igwx)
          CALL iotk_write_empty(iuni,"Info",attr)
        ENDIF


        ALLOCATE( wtmp( max(igwx,1) ) )
        wtmp = 0.0d0

        DO j = 1, nbnd
          IF( t0 ) THEN
            IF( npool > 1 ) THEN
              IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
                CALL mergewf(wf0(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm)
              ENDIF
              IF( ipsour /= ionode_id ) THEN
                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j, world_comm )
              ENDIF
            ELSE
              CALL mergewf(wf0(:,j), wtmp, ngwl, igl, mpime, nproc, &
                           ionode_id, world_comm )
            ENDIF

            IF( ionode ) THEN
              CALL iotk_write_dat(iuni,"Wfc"//iotk_index(j),wtmp(1:igwx))
            ENDIF
          ELSE
          ENDIF
        ENDDO

        DO j = 1, nbnd
          IF( tm ) THEN
            IF( npool > 1 ) THEN
              IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
                CALL mergewf(wfm(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm)
              ENDIF
              IF( ipsour /= ionode_id ) THEN
                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j, world_comm )
              ENDIF
            ELSE
              CALL mergewf(wfm(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id, world_comm )
            ENDIF
            IF( ionode ) THEN
              CALL iotk_write_dat(iuni,"Wfcm"//iotk_index(j),wtmp(1:igwx))
            ENDIF
          ELSE
          ENDIF
        ENDDO
        IF(ionode) CALL iotk_write_end  (iuni,"Kpoint"//iotk_index(ik))

        DEALLOCATE( wtmp )

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!
!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE write_restart_wfc2(iuni, nbnd)
      USE io_global, ONLY: ionode, ionode_id
      IMPLICIT NONE
      INTEGER, INTENT(in) :: iuni, nbnd
      LOGICAL :: twrite = .false.
      INTEGER :: idum, i
      CHARACTER(len=20) :: section_name = 'wfc'
      idum = nbnd
      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!
!
!
!=----------------------------------------------------------------------------=!
  END MODULE
!=----------------------------------------------------------------------------=!



!-----------------------------------------------------------------------
PROGRAM pw_export
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
  USE wrappers,  ONLY : f_mkdir_safe
  USE pwcom
  USE fft_base,  ONLY : dfftp
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE io_files,  ONLY : psfile, pseudo_dir, prefix, tmp_dir
  USE ions_base, ONLY : ntype => nsp
  USE iotk_module
  USE mp_global, ONLY : mp_startup
  USE mp_pools,  ONLY : kunit
  USE mp_world,  ONLY: world_comm
  USE mp,        ONLY: mp_bcast
  USE environment,   ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=256) :: outdir
  !
  INTEGER :: ik, i, kunittmp, ios

  CHARACTER(len=200) :: pp_file
  CHARACTER(len=iotk_attlenx) :: attr
  LOGICAL :: found, uspp_spsi, ascii, single_file, raw

  NAMELIST /inputpp/ prefix, outdir, pp_file, uspp_spsi, ascii, single_file, &
                     raw, psfile, pseudo_dir
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PW_EXPORT' )
  !
  !   set default values for variables in namelist
  !
  prefix='export'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  pp_file= ' '
  uspp_spsi = .false.
  ascii = .false.
  single_file = .false.
  raw = .false.

  !
  !    Reading input file
  !
  IF ( ionode ) THEN
      !
      CALL input_from_file ( )
      !
      READ(5,inputpp,IOSTAT=ios)
      !
      IF (ios /= 0) CALL errore ('pw_export', 'reading inputpp namelist', abs(ios) )
      !
      IF( pp_file == ' ' ) THEN
          !
          pp_file = trim(prefix)//".export/index.xml"
          !
          IF(ionode) ios = f_mkdir_safe( trim(outdir)//"/"//trim(prefix)//".export" )
      ENDIF
      !
  ENDIF
  !
  ! ... Broadcasting variables
  !
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


  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  !

#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  CALL write_export (pp_file, kunittmp, uspp_spsi, ascii, single_file, raw)

  CALL environment_end ( 'PW_EXPORT' )

  CALL stop_pp
  STOP

CONTAINS

!
!-----------------------------------------------------------------------
SUBROUTINE write_export (pp_file,kunit,uspp_spsi, ascii, single_file, raw)
  !-----------------------------------------------------------------------
  !
  USE iotk_module


  USE kinds,          ONLY : DP
  USE pwcom
  USE gvecw,          ONLY : ecutwfc, gcutw
  USE start_k,        ONLY : nk1, nk2, nk3, k1, k2, k3
  USE control_flags,  ONLY : gamma_only
  USE global_version, ONLY : version_number
  USE becmod,         ONLY : bec_type, becp, calbec, &
                             allocate_bec_type, deallocate_bec_type
  USE symm_base,      ONLY : nsym, s, invsym, sname, irt, ftau
  USE  uspp,          ONLY : nkb, vkb
  USE wavefunctions_module,  ONLY : evc
  USE io_files,       ONLY : nd_nmbr, tmp_dir, prefix, iunwfc, nwordwfc
  USE io_files,       ONLY : pseudo_dir, psfile
  USE io_base_export, ONLY : write_restart_wfc
  USE io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
  USE mp_pools,       ONLY : my_pool_id, intra_pool_comm, inter_pool_comm, &
                             nproc_pool
  USE mp,             ONLY : mp_sum, mp_max
  USE mp_world,       ONLY : world_comm, nproc

  IMPLICIT NONE

  CHARACTER(5), PARAMETER :: fmt_name="QEXPT"
  CHARACTER(5), PARAMETER :: fmt_version="1.1.0"

  INTEGER, INTENT(in) :: kunit
  CHARACTER(80), INTENT(in) :: pp_file
  LOGICAL, INTENT(in) :: uspp_spsi, ascii, single_file, raw

  INTEGER :: i, j, k, ig, ik, ibnd, na, ngg,ig_, ierr
  real(DP) :: xyz(3), tmp(3)
  INTEGER :: ike, iks, npw_g, npwx_g, ispin, local_pw
  INTEGER, EXTERNAL :: global_kpoint_index
  INTEGER, ALLOCATABLE :: ngk_g( : )
  INTEGER, ALLOCATABLE :: itmp_g( :, : )
  real(DP),ALLOCATABLE :: rtmp_g( :, : )
  real(DP),ALLOCATABLE :: rtmp_gg( : )
  INTEGER, ALLOCATABLE :: itmp1( : )
  INTEGER, ALLOCATABLE :: igwk( :, : )
  INTEGER, ALLOCATABLE :: l2g_new( : )
  INTEGER, ALLOCATABLE :: igk_l2g( :, : )


  real(DP) :: wfc_scal
  LOGICAL :: twf0, twfm
  CHARACTER(iotk_attlenx) :: attr
  COMPLEX(DP), ALLOCATABLE :: sevc (:,:)

  REAL(DP), ALLOCATABLE :: raux(:)

  IF( nkstot > 0 ) THEN

     IF( ( kunit < 1 ) .or. ( mod( nkstot, kunit ) /= 0 ) ) &
       CALL errore( ' write_export ',' wrong kunit ', 1 )

     IF( ( nproc_pool > nproc ) .or. ( mod( nproc, nproc_pool ) /= 0 ) ) &
       CALL errore( ' write_export ',' nproc_pool ', 1 )

     iks = global_kpoint_index (nkstot, 1)
     ike = iks + nks - 1

  ENDIF

  ! find out the global number of G vectors: ngm_g
  ngm_g = ngm
  CALL mp_sum( ngm_g , intra_pool_comm )


  !  Open file PP_FILE

  IF( ionode ) THEN
    WRITE(0,*) "Opening file "//trim(pp_file)
    CALL iotk_open_write(50,file=trim(tmp_dir) // trim(pp_file))
    WRITE(0,*) "Reconstructing the main grid"
  ENDIF

  ! collect all G vectors across processors within the pools
  ! and compute their modules
  !
  ALLOCATE( itmp_g( 3, ngm_g ) )
  ALLOCATE( rtmp_g( 3, ngm_g ) )
  ALLOCATE( rtmp_gg( ngm_g ) )

  itmp_g = 0
  DO  ig = 1, ngm
    itmp_g( 1, ig_l2g( ig ) ) = mill(1,ig )
    itmp_g( 2, ig_l2g( ig ) ) = mill(2,ig )
    itmp_g( 3, ig_l2g( ig ) ) = mill(3,ig )
  ENDDO
  CALL mp_sum( itmp_g , intra_pool_comm )
  !
  ! here we are in crystal units
  rtmp_g(1:3,1:ngm_g) = REAL( itmp_g(1:3,1:ngm_g) )
  !
  ! go to cartesian units (tpiba)
  CALL cryst_to_cart( ngm_g, rtmp_g, bg , 1 )
  !
  ! compute squared moduli
  DO  ig = 1, ngm_g
     rtmp_gg(ig) = rtmp_g(1,ig)**2 + rtmp_g(2,ig)**2 + rtmp_g(3,ig)**2
  ENDDO
  DEALLOCATE( rtmp_g )

  ! build the G+k array indexes
  ALLOCATE ( igk_l2g ( npwx, nks ) )
  DO ik = 1, nks
     !
     ! mapping between local and global G vector index, for this kpoint
     !
     npw = ngk(ik)
     DO ig = 1, npw
        !
        igk_l2g(ig,ik) = ig_l2g( igk_k(ig,ik) )
        !
     ENDDO
     !
     igk_l2g( npw+1 : npwx, ik ) = 0
     !
  ENDDO

  ! compute the global number of G+k vectors for each k point
  ALLOCATE( ngk_g( nkstot ) )
  ngk_g = 0
  ngk_g( iks:ike ) = ngk( 1:nks )
  CALL mp_sum( ngk_g, world_comm )

  ! compute the Maximum G vector index among all G+k and processors
  npw_g = maxval( igk_l2g(:,:) )
  CALL mp_max( npw_g, world_comm )

  ! compute the Maximum number of G vector among all k points
  npwx_g = maxval( ngk_g( 1:nkstot ) )


  IF( ionode ) THEN
    !
    WRITE(0,*) "Writing header"
    CALL iotk_write_begin(50,"Header")
    CALL iotk_write_attr (attr,"name",trim(fmt_name),FIRST=.true.)
    CALL iotk_write_attr (attr,"version",trim(fmt_version))
    CALL iotk_write_empty(50,"format", ATTR=attr)
    !
    CALL iotk_write_attr (attr,"name","Quantum ESPRESSO",FIRST=.true.)
    CALL iotk_write_attr (attr,"version",trim(version_number))
    CALL iotk_write_empty(50,"creator", ATTR=attr)
    CALL iotk_write_end(50,"Header")
    !
    WRITE(0,*) "Writing dimensions"
    CALL iotk_write_begin(50,"Dimensions")
    CALL iotk_write_attr (attr,"nktot",nkstot,first=.true.)
    CALL iotk_write_attr (attr,"nspin",nspin)
    CALL iotk_write_attr (attr,"nk1",nk1)
    CALL iotk_write_attr (attr,"nk2",nk2)
    CALL iotk_write_attr (attr,"nk3",nk3)
    CALL iotk_write_attr (attr,"s1",k1)
    CALL iotk_write_attr (attr,"s2",k2)
    CALL iotk_write_attr (attr,"s3",k3)
    CALL iotk_write_empty(50,"Kpoints",attr)
    CALL iotk_write_attr (attr,"nbnd",nbnd,first=.true.)
    CALL iotk_write_empty(50,"Bands",attr)
    CALL iotk_write_attr (attr,"gamma_only",gamma_only,first=.true.)
    CALL iotk_write_empty(50,"Gamma_tricks",attr)
    CALL iotk_write_attr (attr,"npw",ngm_g,first=.true.)
    CALL iotk_write_empty(50,"Main_grid",attr)
    CALL iotk_write_attr (attr,"npwx",npwx_g,first=.true.)
    CALL iotk_write_empty(50,"Wfc_grid",attr)
    CALL iotk_write_attr (attr,"natoms",nat,first=.true.)
    CALL iotk_write_empty(50,"Atoms",attr=attr)
    CALL iotk_write_attr (attr,"nsym",nsym,first=.true.)
    CALL iotk_write_empty(50,"Symmops",attr=attr)
    CALL iotk_write_end  (50,"Dimensions")

    WRITE(0,*) "Writing cell"
    CALL iotk_write_attr (attr,"units","a.u.",first=.true.)
    CALL iotk_write_begin(50,"Cell",attr=attr)
    CALL iotk_write_attr (attr,"alat",alat,first=.true.)
    CALL iotk_write_attr (attr,"omega",omega)
    CALL iotk_write_attr (attr,"tpiba",tpiba)
    CALL iotk_write_attr (attr,"tpiba2",tpiba2)
    CALL iotk_write_empty(50,"Data",attr=attr)
    CALL iotk_write_attr (attr,"xyz",at(:,1)*alat,first=.true.)
    CALL iotk_write_empty(50,"a1",attr=attr)
    CALL iotk_write_attr (attr,"xyz",at(:,2)*alat,first=.true.)
    CALL iotk_write_empty(50,"a2",attr=attr)
    CALL iotk_write_attr (attr,"xyz",at(:,3)*alat,first=.true.)
    CALL iotk_write_empty(50,"a3",attr=attr)
    CALL iotk_write_attr (attr,"xyz",bg(:,1)*tpiba,first=.true.)
    CALL iotk_write_empty(50,"b1",attr=attr)
    CALL iotk_write_attr (attr,"xyz",bg(:,2)*tpiba,first=.true.)
    CALL iotk_write_empty(50,"b2",attr=attr)
    CALL iotk_write_attr (attr,"xyz",bg(:,3)*tpiba,first=.true.)
    CALL iotk_write_empty(50,"b3",attr=attr)
    CALL iotk_write_end(50,"Cell")

    WRITE(0,*) "Writing atoms"
    CALL iotk_write_begin(50,"Atoms")
    CALL iotk_write_attr (attr,"natoms",nat,FIRST=.true.)
    CALL iotk_write_attr (attr,"nspecies",nsp)
    CALL iotk_write_empty(50,"Data",attr=attr)
    CALL iotk_write_attr (attr,"units","alat",FIRST=.true.)
    CALL iotk_write_begin(50,"Positions",attr=attr)
    DO i = 1, nat
      xyz = tau(:,i)
!
! this line convert to crystal representation
!      call cryst_to_cart(1,xyz,bg,-1)
!
      CALL iotk_write_attr (attr,"type",trim(atm(ityp(i))),first=.true.)
      CALL iotk_write_attr (attr,"xyz",xyz)
      CALL iotk_write_empty(50,"atom"//trim(iotk_index(i)),attr=attr)
    ENDDO
    CALL iotk_write_end(50,"Positions")
    CALL iotk_write_begin(50,"Types")
    CALL iotk_write_attr (attr,"pseudo_dir",trim(pseudo_dir),FIRST=.true.)
    CALL iotk_write_empty(50,"Data",attr=attr)
    DO i=1, nsp
        CALL iotk_write_attr (attr,"type",trim(atm(i)),FIRST=.true.)
        CALL iotk_write_attr (attr,"pseudo_file",trim(psfile(i)) )
        CALL iotk_write_empty(50,"specie"//trim(iotk_index(i)), ATTR=attr )
    ENDDO
    CALL iotk_write_end  (50,"Types")
    CALL iotk_write_end  (50,"Atoms")

    WRITE(0,*) "Writing symmetry operations"
    CALL iotk_write_begin(50,"Symmetry")
    CALL iotk_write_attr(attr,"nsym",nsym,first=.true.)
    CALL iotk_write_attr(attr,"invsym",invsym)
    CALL iotk_write_empty(50,"symmops",attr)
    !
    ! The matrix s is the transpose of the symmetry matrix in direct space,
    ! in units of a_i.
    !
    DO i=1,nsym
       !
       CALL iotk_write_attr ( attr,"name", trim(sname(i)), FIRST=.true. )
       CALL iotk_write_empty(50,"info"//trim(iotk_index(i)), ATTR=attr )
       !
       tmp(1) = ftau(1,i) / dble( dfftp%nr1 )
       tmp(2) = ftau(2,i) / dble( dfftp%nr2 )
       tmp(3) = ftau(3,i) / dble( dfftp%nr3 )
       !
       CALL iotk_write_attr(attr,"units","crystal",first=.true.)
       !
       CALL iotk_write_dat (50,"sym"//trim(iotk_index(i)), &
                                s(1:3,1:3,i), ATTR=attr, COLUMNS=3)
       CALL iotk_write_dat (50,"trasl"//trim(iotk_index(i)), tmp(:), ATTR=attr )
       !
    ENDDO
    !
    CALL iotk_write_end  (50,"Symmetry")

    WRITE(0,*) "Writing k-mesh"
    CALL iotk_write_attr (attr,"nk",nkstot,first=.true.)
    CALL iotk_write_begin(50,"Kmesh",attr=attr)
    CALL iotk_write_dat  (50,"weights",wk(1:nkstot))
    CALL iotk_write_dat  (50,"k",xk(1:3,1:nkstot),fmt="(3f15.9)")
    CALL iotk_write_end  (50,"Kmesh")

    WRITE(0,*) "Writing other parameters"
    CALL iotk_write_begin(50,"Other_parameters")
    CALL iotk_write_attr(attr,"wfc",ecutwfc,first=.true.)
    CALL iotk_write_attr(attr,"rho",dual*ecutwfc)
    CALL iotk_write_attr(attr,"units","Rydberg")
    CALL iotk_write_empty(50,"Cutoff",attr)
    CALL iotk_write_attr(attr,"nr1",dfftp%nr1,first=.true.)
    CALL iotk_write_attr(attr,"nr2",dfftp%nr2)
    CALL iotk_write_attr(attr,"nr3",dfftp%nr3)
    CALL iotk_write_empty(50,"Space_grid",attr)
    CALL iotk_write_attr(attr,"nelec",nelec,first=.true.)
    CALL iotk_write_empty(50,"Charge",attr)
    CALL iotk_write_end   (50,"Other_parameters")

    WRITE(0,*) "Writing main grid"
    CALL iotk_write_attr(attr,"npw",   ngm_g,first=.true.)
    CALL iotk_write_attr(attr,"gamma_only", gamma_only )
    CALL iotk_write_attr(attr,"cutoff","NOT AVAILABLE")
    IF(.not.single_file) &
      CALL iotk_link(50,"Main_grid","mgrid",create=.true.,binary=.not.ascii,raw=raw)
    CALL iotk_write_begin(50,"Main_grid",attr=attr)
    CALL iotk_write_attr(attr,"units", "crystal",first=.true.)
    CALL iotk_write_dat(50,"g",itmp_g(1:3,1:ngm_g),fmt="(3i5)", attr=attr)
    CALL iotk_write_attr(attr,"units", "tpiba^2",first=.true.)
    CALL iotk_write_dat(50,"gg",rtmp_gg(1:ngm_g),attr=attr)
    CALL iotk_write_end(50,"Main_grid")
  ENDIF
  DEALLOCATE( rtmp_gg )

  ! for each k point build and write the global G+k indexes array
  ALLOCATE( igwk( npwx_g,nkstot ) )
  WRITE(0,*) "Writing grids for wfc"
  CALL iotk_write_attr (attr,"npwx",npwx_g,first=.true.)
  IF(ionode) CALL iotk_write_begin(50,"Wfc_grids",ATTR=attr)


  DO ik = 1, nkstot
    igwk(:,ik) = 0
    !
    ALLOCATE( itmp1( npw_g ), STAT= ierr )
    IF ( ierr/=0 ) CALL errore('pw_export','allocating itmp1', abs(ierr) )
    itmp1 = 0
    !
    IF( ik >= iks .and. ik <= ike ) THEN
      DO  ig = 1, ngk( ik-iks+1 )
        itmp1( igk_l2g( ig, ik-iks+1 ) ) = igk_l2g( ig, ik-iks+1 )
      ENDDO
    ENDIF
    !
    CALL mp_sum( itmp1, world_comm )
    !
    ngg = 0
    DO  ig = 1, npw_g
      IF( itmp1( ig ) == ig ) THEN
        ngg = ngg + 1
        igwk( ngg , ik) = ig
      ENDIF
    ENDDO
    IF( ngg /= ngk_g( ik ) ) THEN
      WRITE( stdout,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
    ENDIF
    !
    DEALLOCATE( itmp1 )
    !
    IF( ionode ) THEN
      CALL iotk_write_attr (attr,"npw",ngk_g(ik),first=.true.)
      CALL iotk_write_attr(attr,"gamma_only", gamma_only )
      CALL iotk_write_attr (attr,"kcry",xk(1:3,ik))
      IF(.not.single_file) &
          CALL iotk_link(50,"Kpoint"//iotk_index(ik),"grid"//iotk_index(ik), &
                         create=.true.,binary=.not.ascii,raw=raw)
      CALL iotk_write_begin(50,"Kpoint"//iotk_index(ik),attr)
      CALL iotk_write_dat  (50,"index",igwk(1:ngk_g(ik),ik))
      CALL iotk_write_dat  (50,"grid",itmp_g(1:3,igwk(1:ngk_g(ik),ik)),fmt="(3i5)")
      CALL iotk_write_end  (50,"Kpoint"//iotk_index(ik))
    ENDIF
  ENDDO

  IF(ionode) CALL iotk_write_end(50,"Wfc_grids")
  DEALLOCATE( itmp_g )


#if defined(__MPI)
  CALL poolrecover (et, nbnd, nkstot, nks)
#endif
!
  ALLOCATE(raux(1:nbnd))
!

  WRITE(0,*) "Writing band structure"

  IF( ionode ) THEN
    CALL iotk_write_attr (attr,"nspin",nspin,first=.true.)
    CALL iotk_write_attr (attr,"nk",nkstot)
    CALL iotk_write_attr (attr,"nbnd",nbnd)
    CALL iotk_write_attr (attr,"efermi",ef)
    CALL iotk_write_attr (attr,"units","Rydberg")
    CALL iotk_write_begin(50,"Eigenvalues",attr=attr)
    DO ik=1,nkstot
      CALL iotk_write_dat(50,"e"//iotk_index(ik),et(1:nbnd,ik))
    ENDDO
    CALL iotk_write_end  (50,"Eigenvalues")
  ENDIF

  IF( ionode ) THEN
    CALL iotk_write_attr (attr,"nspin",nspin,first=.true.)
    CALL iotk_write_attr (attr,"nk",nkstot)
    CALL iotk_write_attr (attr,"nbnd",nbnd)
    CALL iotk_write_begin(50,"OCCUPATIONS",attr=attr)
    DO ik=1,nkstot
      IF ( wk(ik) == 0.D0 ) THEN
        !
        raux = wg(:,ik)
        !
      ELSE
        !
        raux = wg(:,ik) / wk(ik)
        !
      END IF
      CALL iotk_write_dat(50,"wg"//iotk_index(ik),raux(1:nbnd))
    ENDDO
    CALL iotk_write_end  (50,"OCCUPATIONS")
  ENDIF
  !
  DEALLOCATE(raux)
  !
  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.


  WRITE(0,*) "Writing Eigenvectors"
  IF( ionode ) CALL iotk_write_begin(50, "Eigenvectors")

  DO ik = 1, nkstot
    IF(.not.single_file .and. ionode) &
       CALL iotk_link(50,"Kpoint"//iotk_index(ik),"wfc"//iotk_index(ik), &
                         create=.true.,binary=.not.ascii,raw=raw)

     local_pw = 0
     IF( (ik >= iks) .and. (ik <= ike) ) THEN

       CALL davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), - 1)
       local_pw = ngk(ik-iks+1)

     ENDIF


     ALLOCATE(l2g_new(local_pw))

     l2g_new = 0
     DO ig = 1, local_pw
       ngg = igk_l2g(ig,ik-iks+1)
       DO ig_ = 1, ngk_g(ik)
         IF(ngg == igwk(ig_,ik)) THEN
           l2g_new(ig) = ig_
           exit
         ENDIF
       ENDDO
     ENDDO


     ispin = isk( ik )
     !  WRITE(0,*) ' ### ', ik,nkstot,iks,ike,kunit,nproc,nproc_pool
     CALL write_restart_wfc(50, ik, nkstot, kunit, ispin, nspin, &
         wfc_scal, evc, twf0, evc, twfm, npw_g, gamma_only, nbnd, &
         l2g_new(:),local_pw )
     DEALLOCATE(l2g_new)
  ENDDO
  IF( ionode ) CALL iotk_write_end  (50, "Eigenvectors")


  !
  ! If specified and if USPP are used the wfcs S_psi are written
  ! | spsi_nk > = \hat S | psi_nk >
  ! where S is the overlap operator of US PP
  !
  IF ( uspp_spsi .and. nkb > 0 ) THEN

       ALLOCATE( sevc(npwx,nbnd), STAT=ierr )
       IF (ierr/=0) CALL errore( ' write_export ',' Unable to allocate SEVC ', abs(ierr) )

       WRITE(0,*) "Writing Eigenvectors_Spsi"
       IF( ionode ) CALL iotk_write_begin(50, "Eigenvectors_Spsi")

       CALL init_us_1
       CALL init_at_1

       CALL allocate_bec_type (nkb,nbnd, becp)

       DO ik = 1, nkstot

           IF(.not.single_file .and. ionode) &
                 CALL iotk_link(50,"Kpoint"//iotk_index(ik),"swfc"//iotk_index(ik), &
                       create=.true.,binary=.not.ascii,raw=raw)

           local_pw = 0
           IF( (ik >= iks) .and. (ik <= ike) ) THEN
 
               CALL davcio (evc, 2*nwordwfc, iunwfc, (ik-iks+1), - 1)

               npw = ngk(ik-iks+1)
               local_pw = npw
	       CALL init_us_2(npw, igk_k(1,ik-iks+1), xk(1, ik), vkb)

               IF ( gamma_only ) THEN
                  CALL calbec ( npw, vkb, evc, becp )
                  WRITE(0,*) 'Gamma only PW_EXPORT not yet tested'
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
                 exit
               ENDIF
             ENDDO
           ENDDO

           ispin = isk( ik )
           CALL write_restart_wfc(50, ik, nkstot, kunit, ispin, nspin, &
               wfc_scal, sevc, twf0, sevc, twfm, npw_g, gamma_only, nbnd, &
               l2g_new(:),local_pw )
           DEALLOCATE(l2g_new)
       ENDDO
       IF( ionode ) CALL iotk_write_end  (50, "Eigenvectors_Spsi")

       DEALLOCATE( sevc, STAT=ierr )
       IF ( ierr/= 0 ) CALL errore('pw_export','Unable to deallocate SEVC',abs(ierr))
       CALL deallocate_bec_type ( becp )
  ENDIF

  DEALLOCATE( igk_l2g )
  DEALLOCATE( igwk )
  DEALLOCATE ( ngk_g )

  IF( ionode ) THEN
    CALL iotk_close_write(50)
  ENDIF

END SUBROUTINE write_export

END PROGRAM pw_export

