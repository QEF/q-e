!
! Copyright (C) 2003-2006 Andrea Ferretti and Quantum-Espresso group
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
      USE mp, ONLY: mp_sum, mp_get, mp_bcast, mp_max
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool, world_comm
      USE io_global, ONLY: ionode, ionode_id
      USE iotk_module
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP), INTENT(IN) :: wf0(:,:)
      COMPLEX(DP), INTENT(IN) :: wfm(:,:)
      INTEGER, INTENT(IN) :: ngw   ! 
      LOGICAL, INTENT(IN) :: gamma_only
      INTEGER, INTENT(IN) :: nbnd
      INTEGER, INTENT(IN) :: ngwl
      INTEGER, INTENT(IN) :: igl(:)
      REAL(DP), INTENT(IN) :: scal
      LOGICAL, INTENT(IN) :: t0, tm

      INTEGER :: i, j, ierr, idum = 0
      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx
      INTEGER :: npool, ipmask( nproc ), ipsour
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)

      CHARACTER(LEN=20) :: section_name = 'wfc'

      LOGICAL :: twrite = .TRUE.

      INTEGER :: ierr_iotk
      CHARACTER(LEN=iotk_attlenx) :: attr

!
! ... Subroutine Body
!

        ! set working variables for k point index (ikt) and k points number (nkt)
        ikt = ik
        nkt = nk

        !  find out the number of pools
        npool = nproc / nproc_pool 

        !  find out number of k points blocks
        nkbl = nkt / kunit  

        !  k points per pool
        nkl = kunit * ( nkbl / npool )

        !  find out the reminder
        nkr = ( nkt - nkl * npool ) / kunit

        !  Assign the reminder to the first nkr pools
        IF( my_pool_id < nkr ) nkl = nkl + kunit

        !  find out the index of the first k point in this pool
        iks = nkl * my_pool_id + 1
        IF( my_pool_id >= nkr ) iks = iks + nkr * kunit
      
        !  find out the index of the last k point in this pool
        ike = iks + nkl - 1

        ipmask = 0
        ipsour = ionode_id

        !  find out the index of the processor which collect the data in the pool of ik
        IF( npool > 1 ) THEN
          IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            IF( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
          END IF
          CALL mp_sum( ipmask )
          DO i = 1, nproc
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
          END DO
        END IF

        igwx = 0
        ierr = 0
        IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
          IF( ngwl > SIZE( igl ) ) THEN
            ierr = 1
          ELSE
            igwx = MAXVAL( igl(1:ngwl) )
          END IF
        END IF

        ! get the maximum index within the pool
        !
        CALL mp_max( igwx, intra_pool_comm ) 

        ! now notify all procs if an error has been found 
        !
        CALL mp_max( ierr ) 

        IF( ierr > 0 ) &
          CALL errore(' write_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipsour /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipsour, 1 )
        END IF

        if(ionode) then
          call iotk_write_begin(iuni,"Kpoint"//iotk_index(ik))
          call iotk_write_attr (attr,"ngw",ngw,first=.true.)
          call iotk_write_attr (attr,"nbnd",nbnd)
          call iotk_write_attr (attr,"gamma_only",gamma_only)
          call iotk_write_attr (attr,"ik",ik)
          call iotk_write_attr (attr,"nk",nk)
          call iotk_write_attr (attr,"kunit",kunit)
          call iotk_write_attr (attr,"ispin",ispin)
          call iotk_write_attr (attr,"nspin",nspin)
          call iotk_write_attr (attr,"scal",scal)
          call iotk_write_attr (attr,"igwx",igwx)
          call iotk_write_empty(iuni,"Info",attr)
        end if


        ALLOCATE( wtmp( MAX(igwx,1) ) )
        wtmp = 0.0d0

        DO j = 1, nbnd
          IF( t0 ) THEN
            IF( npool > 1 ) THEN
              IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
                CALL mergewf(wf0(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm)
              END IF
              IF( ipsour /= ionode_id ) THEN
                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j, world_comm )
              END IF
            ELSE
              CALL mergewf(wf0(:,j), wtmp, ngwl, igl, mpime, nproc, &
                           ionode_id, world_comm )
            END IF

            if( ionode ) then
              call iotk_write_dat(iuni,"Wfc"//iotk_index(j),wtmp(1:igwx))
            end if
          ELSE
          END IF
        END DO

        DO j = 1, nbnd
          IF( tm ) THEN
            IF( npool > 1 ) THEN
              IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
                CALL mergewf(wfm(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm)
              END IF
              IF( ipsour /= ionode_id ) THEN
                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j, world_comm )
              END IF
            ELSE
              CALL mergewf(wfm(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id, world_comm )
            END IF
            if( ionode ) then
              call iotk_write_dat(iuni,"Wfcm"//iotk_index(j),wtmp(1:igwx))
            end if
          ELSE
          END IF
        END DO
        if(ionode) call iotk_write_end  (iuni,"Kpoint"//iotk_index(ik))

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
      USE mp, ONLY: mp_bcast
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iuni, nbnd
      LOGICAL :: twrite = .FALSE.
      INTEGER :: idum, i
      CHARACTER(LEN=20) :: section_name = 'wfc'
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
  

  use pwcom
  use io_global, ONLY : stdout, ionode, ionode_id
  use io_files,  ONLY : psfile, pseudo_dir, trimcheck
  use io_files,  ONLY : nd_nmbr, prefix, tmp_dir, outdir
  use ions_base, ONLY : ntype => nsp
  use iotk_module
  use mp_global, ONLY : mpime, kunit
  use mp, ONLY: mp_bcast

  !
  implicit none
  integer :: ik, i, kunittmp, ios

  character(len=200) :: pp_file
  character(len=iotk_attlenx) :: attr
  logical :: found, uspp_spsi, ascii, single_file, raw
  INTEGER, EXTERNAL :: C_MKDIR

  NAMELIST /inputpp/ prefix, outdir, pp_file, uspp_spsi, ascii, single_file, raw, &
                     psfile, pseudo_dir

  !
  call start_postproc (nd_nmbr)
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
  !    Reading input file
  !
  IF ( ionode ) THEN
      !
      CALL input_from_file ( )
      !
      READ(5,inputpp,IOSTAT=ios)
      !
      IF (ios /= 0) CALL errore ('pw_export', 'reading inputpp namelist', ABS(ios) )
      !
      IF( pp_file == ' ' ) THEN
          !
          pp_file = TRIM(prefix)//".export/index.xml"
          !
          if(ionode) ios = C_MKDIR( TRIM(outdir)//"/"//TRIM(prefix)// &
                           ".export" , LEN(TRIM(outdir)//"/"//TRIM(prefix)//".export") )
      ENDIF
      !
  ENDIF
  !
  ! ... Broadcasting variables
  !
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
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file  
  call openfil_pp
  !

#if defined __PARA
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  call write_export (pp_file, kunittmp, uspp_spsi, ascii, single_file, raw)

  call stop_pp
  stop
end program pp_punch
!
!-----------------------------------------------------------------------
subroutine write_export (pp_file,kunit,uspp_spsi, ascii, single_file, raw)
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  use iotk_module


  use kinds,          ONLY : DP 
  use pwcom  
  use control_flags,  ONLY : gamma_only  
  use global_version, ONLY : version_number
  use becmod,         ONLY : becp, rbecp, calbec, allocate_bec, deallocate_bec
  use wavefunctions_module,  ONLY : evc
  use io_files,       ONLY : nd_nmbr, outdir, prefix, iunwfc, nwordwfc
  use io_files,       ONLY : pseudo_dir, psfile
  use io_base_export, ONLY : write_restart_wfc
  use io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
  use mp_global,      ONLY : nproc, nproc_pool, mpime
  use mp_global,      ONLY : my_pool_id, intra_pool_comm, inter_pool_comm
  use mp,             ONLY : mp_sum, mp_max

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

  ! find out the global number of G vectors: ngm_g  
  ngm_g = ngm
  call mp_sum( ngm_g , intra_pool_comm )


  !  Open file PP_FILE

  if( ionode ) then
    write(0,*) "Opening file "//trim(pp_file)
    call iotk_open_write(50,file=TRIM(outdir)//'/'//TRIM(pp_file))
    write(0,*) "Reconstructing the main grid"
  end if

  ! collect all G vectors across processors within the pools
  ! and compute their modules
  !
  allocate( itmp_g( 3, ngm_g ) )
  allocate( rtmp_g( 3, ngm_g ) )
  allocate( rtmp_gg( ngm_g ) )

  itmp_g = 0
  do  ig = 1, ngm
    itmp_g( 1, ig_l2g( ig ) ) = ig1( ig )
    itmp_g( 2, ig_l2g( ig ) ) = ig2( ig )
    itmp_g( 3, ig_l2g( ig ) ) = ig3( ig )
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


  if( ionode ) then
    !
    write(0,*) "Writing header"
    CALL iotk_write_begin(50,"Header")
    CALL iotk_write_attr (attr,"name",TRIM(fmt_name),FIRST=.true.)
    CALL iotk_write_attr (attr,"version",TRIM(fmt_version))
    CALL iotk_write_empty(50,"format", ATTR=attr)
    !
    CALL iotk_write_attr (attr,"name","Quantum-ESPRESSO",FIRST=.true.)
    CALL iotk_write_attr (attr,"version",TRIM(version_number))
    CALL iotk_write_empty(50,"creator", ATTR=attr)
    CALL iotk_write_end(50,"Header")
    !
    write(0,*) "Writing dimensions"
    call iotk_write_begin(50,"Dimensions")
    call iotk_write_attr (attr,"nktot",nkstot,first=.true.)
    call iotk_write_attr (attr,"nspin",nspin)
    call iotk_write_attr (attr,"nk1",nk1)
    call iotk_write_attr (attr,"nk2",nk2)
    call iotk_write_attr (attr,"nk3",nk3)
    call iotk_write_attr (attr,"s1",k1)
    call iotk_write_attr (attr,"s2",k2)
    call iotk_write_attr (attr,"s3",k3)
    call iotk_write_empty(50,"Kpoints",attr)
    call iotk_write_attr (attr,"nbnd",nbnd,first=.true.)
    call iotk_write_empty(50,"Bands",attr)
    call iotk_write_attr (attr,"gamma_only",gamma_only,first=.true.)
    call iotk_write_empty(50,"Gamma_tricks",attr)
    call iotk_write_attr (attr,"npw",ngm_g,first=.true.)
    call iotk_write_empty(50,"Main_grid",attr)
    call iotk_write_attr (attr,"npwx",npwx_g,first=.true.)
    call iotk_write_empty(50,"Wfc_grid",attr)
    call iotk_write_attr (attr,"natoms",nat,first=.true.)
    call iotk_write_empty(50,"Atoms",attr=attr)
    call iotk_write_attr (attr,"nsym",nsym,first=.true.)
    call iotk_write_empty(50,"Symmops",attr=attr)
    call iotk_write_end  (50,"Dimensions")

    write(0,*) "Writing cell"
    call iotk_write_attr (attr,"units","a.u.",first=.true.)
    call iotk_write_begin(50,"Cell",attr=attr)
    call iotk_write_attr (attr,"alat",alat,first=.true.)
    call iotk_write_attr (attr,"omega",omega)
    call iotk_write_attr (attr,"tpiba",tpiba)
    call iotk_write_attr (attr,"tpiba2",tpiba2)
    call iotk_write_empty(50,"Data",attr=attr)
    call iotk_write_attr (attr,"xyz",at(:,1)*alat,first=.true.)
    call iotk_write_empty(50,"a1",attr=attr)
    call iotk_write_attr (attr,"xyz",at(:,2)*alat,first=.true.)
    call iotk_write_empty(50,"a2",attr=attr)
    call iotk_write_attr (attr,"xyz",at(:,3)*alat,first=.true.)
    call iotk_write_empty(50,"a3",attr=attr)
    call iotk_write_attr (attr,"xyz",bg(:,1)*tpiba,first=.true.)
    call iotk_write_empty(50,"b1",attr=attr)
    call iotk_write_attr (attr,"xyz",bg(:,2)*tpiba,first=.true.)
    call iotk_write_empty(50,"b2",attr=attr)
    call iotk_write_attr (attr,"xyz",bg(:,3)*tpiba,first=.true.)
    call iotk_write_empty(50,"b3",attr=attr)
    call iotk_write_end(50,"Cell")

    write(0,*) "Writing atoms"
    call iotk_write_begin(50,"Atoms")
    call iotk_write_attr (attr,"natoms",nat,FIRST=.TRUE.)
    call iotk_write_attr (attr,"nspecies",nsp)
    call iotk_write_empty(50,"Data",attr=attr)
    call iotk_write_attr (attr,"units","alat",FIRST=.TRUE.)
    call iotk_write_begin(50,"Positions",attr=attr)
    DO i = 1, nat
      xyz = tau(:,i) 
!
! this line convert to crystal representation
!      call cryst_to_cart(1,xyz,bg,-1)
!
      call iotk_write_attr (attr,"type",TRIM(atm(ityp(i))),first=.true.)
      call iotk_write_attr (attr,"xyz",xyz)
      call iotk_write_empty(50,"atom"//trim(iotk_index(i)),attr=attr)
    ENDDO
    call iotk_write_end(50,"Positions")
    call iotk_write_begin(50,"Types")
    call iotk_write_attr (attr,"pseudo_dir",TRIM(pseudo_dir),FIRST=.TRUE.)
    call iotk_write_empty(50,"Data",attr=attr)
    DO i=1, nsp 
        call iotk_write_attr (attr,"type",TRIM(atm(i)),FIRST=.TRUE.)
        call iotk_write_attr (attr,"pseudo_file",TRIM(psfile(i)) )
        call iotk_write_empty(50,"specie"//TRIM(iotk_index(i)), ATTR=attr )
    ENDDO
    call iotk_write_end  (50,"Types")
    call iotk_write_end  (50,"Atoms")

    write(0,*) "Writing symmetry operations"
    call iotk_write_begin(50,"Symmetry")
    call iotk_write_attr(attr,"nsym",nsym,first=.true.)
    call iotk_write_attr(attr,"invsym",invsym)
    call iotk_write_empty(50,"symmops",attr)
    !
    ! The matrix s is the transpose of the symmetry matrix in direct space,
    ! in units of a_i.
    !
    DO i=1,nsym
       !
       CALL iotk_write_attr ( attr,"name", TRIM(sname(i)), FIRST=.TRUE. )
       CALL iotk_write_empty(50,"info"//TRIM(iotk_index(i)), ATTR=attr )
       !
       tmp(1) = ftau(1,i) / DBLE( nr1 )
       tmp(2) = ftau(2,i) / DBLE( nr2 )
       tmp(3) = ftau(3,i) / DBLE( nr3 )
       !
       CALL iotk_write_attr(attr,"units","crystal",first=.TRUE.)
       !
       CALL iotk_write_dat (50,"sym"//TRIM(iotk_index(i)), &
                                s(1:3,1:3,i), ATTR=attr, COLUMNS=3)
       CALL iotk_write_dat (50,"trasl"//TRIM(iotk_index(i)), tmp(:), ATTR=attr )
       !
    ENDDO
    !
    CALL iotk_write_end  (50,"Symmetry")

    write(0,*) "Writing k-mesh"
    call iotk_write_attr (attr,"nk",nkstot,first=.true.)
    call iotk_write_begin(50,"Kmesh",attr=attr)
    call iotk_write_dat  (50,"weights",wk(1:nkstot))
    call iotk_write_dat  (50,"k",xk(1:3,1:nkstot),fmt="(3f15.9)")
    call iotk_write_end  (50,"Kmesh")

    write(0,*) "Writing other parameters"
    call iotk_write_begin(50,"Other_parameters")
    call iotk_write_attr(attr,"wfc",ecutwfc,first=.true.)
    call iotk_write_attr(attr,"rho",dual*ecutwfc)
    call iotk_write_attr(attr,"units","Rydberg")
    call iotk_write_empty(50,"Cutoff",attr)
    call iotk_write_attr(attr,"nr1",nr1,first=.true.)
    call iotk_write_attr(attr,"nr2",nr2)
    call iotk_write_attr(attr,"nr3",nr3)
    call iotk_write_empty(50,"Space_grid",attr)
    call iotk_write_attr(attr,"nelec",nelec,first=.true.)
    call iotk_write_empty(50,"Charge",attr)
    call iotk_write_end   (50,"Other_parameters")

    write(0,*) "Writing main grid"
    call iotk_write_attr(attr,"npw",   ngm_g,first=.true.)
    call iotk_write_attr(attr,"gamma_only", gamma_only )
    call iotk_write_attr(attr,"cutoff","NOT AVAILABLE")
    if(.not.single_file) &
      call iotk_link(50,"Main_grid","mgrid",create=.true.,binary=.not.ascii,raw=raw)
    call iotk_write_begin(50,"Main_grid",attr=attr)
    call iotk_write_attr(attr,"units", "crystal",first=.true.)
    call iotk_write_dat(50,"g",itmp_g(1:3,1:ngm_g),fmt="(3i5)", attr=attr)
    call iotk_write_attr(attr,"units", "tpiba^2",first=.true.)
    call iotk_write_dat(50,"gg",rtmp_gg(1:ngm_g),attr=attr)
    call iotk_write_end(50,"Main_grid")
  end if
  deallocate( rtmp_gg )

  ! for each k point build and write the global G+k indexes array
  allocate( igwk( npwx_g,nkstot ) )
  write(0,*) "Writing grids for wfc"
  call iotk_write_attr (attr,"npwx",npwx_g,first=.true.)
  if(ionode) call iotk_write_begin(50,"Wfc_grids",ATTR=attr)


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
    IF( ionode ) THEN
      CALL iotk_write_attr (attr,"npw",ngk_g(ik),first=.true.)
      call iotk_write_attr(attr,"gamma_only", gamma_only )
      CALL iotk_write_attr (attr,"kcry",xk(1:3,ik))
      IF(.NOT.single_file) &
          CALL iotk_link(50,"Kpoint"//iotk_index(ik),"grid"//iotk_index(ik), &
                         create=.true.,binary=.not.ascii,raw=raw)
      CALl iotk_write_begin(50,"Kpoint"//iotk_index(ik),attr)
      CALL iotk_write_dat  (50,"index",igwk(1:ngk_g(ik),ik))
      CALL iotk_write_dat  (50,"grid",itmp_g(1:3,igwk(1:ngk_g(ik),ik)),fmt="(3i5)")
      CALL iotk_write_end  (50,"Kpoint"//iotk_index(ik))
    END IF
  END DO

  if(ionode) call iotk_write_end(50,"Wfc_grids")
  deallocate( itmp_g )


#ifdef __PARA
  call poolrecover (et, nbnd, nkstot, nks)
#endif
 

  Write(0,*) "Writing band structure"

  if( ionode ) then
    call iotk_write_attr (attr,"nspin",nspin,first=.true.)
    call iotk_write_attr (attr,"nk",nkstot)
    call iotk_write_attr (attr,"nbnd",nbnd)
    call iotk_write_attr (attr,"efermi",ef)
    call iotk_write_attr (attr,"units","Rydberg")
    call iotk_write_begin(50,"Eigenvalues",attr=attr)
    do ik=1,nkstot
      call iotk_write_dat(50,"e"//iotk_index(ik),et(1:nbnd,ik))
    end do
    call iotk_write_end  (50,"Eigenvalues")
  end if


  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.


  Write(0,*) "Writing Eigenvectors"
  if( ionode ) call iotk_write_begin(50, "Eigenvectors")

  do ik = 1, nkstot
    if(.not.single_file .and. ionode) &
       call iotk_link(50,"Kpoint"//iotk_index(ik),"wfc"//iotk_index(ik), &
                         create=.true.,binary=.not.ascii,raw=raw)

     local_pw = 0
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN

       call davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)
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
     CALL write_restart_wfc(50, ik, nkstot, kunit, ispin, nspin, &
         wfc_scal, evc, twf0, evc, twfm, npw_g, gamma_only, nbnd, &
         l2g_new(:),local_pw )
     deallocate(l2g_new)
  end do
  if( ionode ) call iotk_write_end  (50, "Eigenvectors")


  !  
  ! If specified and if USPP are used the wfcs S_psi are written  
  ! | spsi_nk > = \hat S | psi_nk >  
  ! where S is the overlap operator of US PP 
  !  
  IF ( uspp_spsi .AND. nkb > 0 ) THEN

       ALLOCATE( sevc(npwx,nbnd), STAT=ierr )
       IF (ierr/=0) CALL errore( ' write_export ',' Unable to allocate SEVC ', ABS(ierr) )

       Write(0,*) "Writing Eigenvectors_Spsi"
       if( ionode ) call iotk_write_begin(50, "Eigenvectors_Spsi")

       CALL init_us_1
       CALL init_at_1

       CALL allocate_bec (nkb,nbnd)

       do ik = 1, nkstot
 
           IF(.NOT.single_file .AND. ionode) &
                 call iotk_link(50,"Kpoint"//iotk_index(ik),"swfc"//iotk_index(ik), &
                       create=.true.,binary=.not.ascii,raw=raw)

           local_pw = 0
           IF( (ik >= iks) .AND. (ik <= ike) ) THEN
               
               CALL gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
               CALL davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)

               CALL init_us_2(npw, igk, xk(1, ik), vkb)
               local_pw = ngk(ik-iks+1)
                            
               IF ( gamma_only ) THEN
                  CALL calbec ( ngk_g(ik), vkb, evc, rbecp )
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
                 EXIT
               ENDIF
             ENDDO
           ENDDO

           ispin = isk( ik )
           CALL write_restart_wfc(50, ik, nkstot, kunit, ispin, nspin, &
               wfc_scal, sevc, twf0, sevc, twfm, npw_g, gamma_only, nbnd, &
               l2g_new(:),local_pw )
           DEALLOCATE(l2g_new)
       ENDDO
       if( ionode ) call iotk_write_end  (50, "Eigenvectors_Spsi")
      
       DEALLOCATE( sevc, STAT=ierr )
       IF ( ierr/= 0 ) CALL errore('pw_export','Unable to deallocate SEVC',ABS(ierr))
       CALL deallocate_bec ()
  ENDIF

  DEALLOCATE( igk_l2g )
  DEALLOCATE( igwk )
  DEALLOCATE ( ngk_g )

  IF( ionode ) THEN
    call iotk_close_write(50)
  END IF

end subroutine write_export

