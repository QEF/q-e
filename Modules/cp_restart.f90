!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE cp_restart
!------------------------------------------------------------------------------!

    !  This module contains subroutines to write and read data required to
    !  restart a calculation from the disk  
!
      IMPLICIT NONE
      SAVE

      PRIVATE :: write_wfc, read_wfc

!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    SUBROUTINE create_directory( dirname )
      USE io_global, ONLY: ionode, ionode_id, stdout
      USE mp, ONLY: mp_bcast
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=256) :: dirname
      INTEGER :: ios, ik
      !
      INTEGER, EXTERNAL :: C_MKDIR

      IF ( ionode ) THEN
          WRITE( stdout, * ) 'Creating dir : ',TRIM( dirname )
          !
          ios = C_MKDIR( TRIM( dirname ), LEN_TRIM( dirname ) )
      END IF
      CALL mp_bcast( ios, ionode_id )
      IF( ios /= 0 ) THEN
        CALL errore(' cp_writefile ', ' unable to create directory '//dirname , ios )
      END IF  

      IF ( ionode ) THEN
        OPEN( UNIT = 4, FILE = TRIM( dirname ) // '/' // 'tst' , &
            STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
      END IF
      !
      CALL mp_bcast( ios, ionode_id )
      !
      IF ( ios /= 0 ) &
           CALL errore( ' cp_writefile : ', TRIM( dirname ) // &
                      & ' non existent or non writable', 1 )

      RETURN
    END SUBROUTINE

!------------------------------------------------------------------------------!

    SUBROUTINE cp_writefile( scradir, ascii, nfi, simtime, acc, nk, xk, wk, &
        ht, htm, htm2, htvel, taui, cdmi, stau0, svel0, staum, svelm, force, &
        occ0, occm, lambda0, lambdam, b1, b2, b3, c0, cm  )
      !
      USE iotk_module
      USE kinds, ONLY: dbl
      USE io_global, ONLY: ionode, ionode_id, stdout
      USE parser, ONLY: int_to_char
      USE control_flags, ONLY: ndw, gamma_only, force_pairing
      USE io_files, ONLY: iunpun, xmlpun
      USE printout_base, ONLY: title
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b
      USE gvecp, ONLY: ngm, ngmt, ecutp
      USE gvecs, ONLY: ngs, ngst
      USE gvecw, ONLY: ngw, ngwt, ecutw
      USE electrons_base, ONLY: nspin, nbnd, nelt, nel, nupdwn, iupdwn
      USE electrons_nose, ONLY: xnhe0, xnhem, xnhep, xnhem2, vnhe
      USE cell_base, ONLY: ibrav, alat, celldm
      USE cell_nose, ONLY: xnhh0, xnhhm, xnhhp, vnhh
      USE ions_base, ONLY: nsp, nat, na, atm, zv, pmass
      USE ions_nose, ONLY: vnhp, xnhp0, xnhpm, xnhpp, xnhpm2 
      USE reciprocal_vectors, ONLY: ngwt, ngw, ig_l2g, mill_l
      USE mp, ONLY: mp_sum

      
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: scradir
      LOGICAL, INTENT(IN) :: ascii
      INTEGER, INTENT(IN) :: nfi    !  index of the current step
      REAL(dbl), INTENT(IN) :: simtime   !  simulated time
      REAL(dbl), INTENT(IN) :: acc(:)   !  
      INTEGER, INTENT(IN) :: nk    !  number of kpoints
      REAL(dbl), INTENT(IN) :: xk(:,:) ! k points coordinates 
      REAL(dbl), INTENT(IN) :: wk(:)   ! k points weights
      REAL(dbl), INTENT(IN) :: ht(3,3) ! 
      REAL(dbl), INTENT(IN) :: htm(3,3) ! 
      REAL(dbl), INTENT(IN) :: htm2(3,3) ! 
      REAL(dbl), INTENT(IN) :: htvel(3,3) ! 
      REAL(dbl), INTENT(IN) :: taui(:,:) ! 
      REAL(dbl), INTENT(IN) :: cdmi(:) ! 
      REAL(dbl), INTENT(IN) :: stau0(:,:) ! 
      REAL(dbl), INTENT(IN) :: svel0(:,:) ! 
      REAL(dbl), INTENT(IN) :: staum(:,:) ! 
      REAL(dbl), INTENT(IN) :: svelm(:,:) ! 
      REAL(dbl), INTENT(IN) :: force(:,:) ! 
      REAL(dbl), INTENT(IN) :: occ0(:,:,:) ! 
      REAL(dbl), INTENT(IN) :: occm(:,:,:) ! 
      REAL(dbl), INTENT(IN) :: lambda0(:,:) ! 
      REAL(dbl), INTENT(IN) :: lambdam(:,:) ! 
      REAL(dbl), INTENT(IN) :: b1(3) ! 
      REAL(dbl), INTENT(IN) :: b2(3) ! 
      REAL(dbl), INTENT(IN) :: b3(3) ! 
      COMPLEX(dbl), INTENT(IN) :: c0(:,:,:,:) ! 
      COMPLEX(dbl), INTENT(IN) :: cm(:,:,:,:) ! 
      !
      CHARACTER(LEN=256) :: dirname, kdirname, filename
      CHARACTER(LEN=5) :: kindex
      INTEGER :: strlen
      CHARACTER(iotk_attlenx) :: attr
      INTEGER :: kunit
      INTEGER :: k1, k2, k3
      INTEGER :: nk1, nk2, nk3
      INTEGER :: i, ispin, ig, nspin_wfc
      INTEGER, ALLOCATABLE :: mill(:,:)
      REAL(dbl) :: wfscal, omega, htm1(3,3)

      dirname = 'RESTART' // int_to_char( ndw )
      IF ( LEN( scradir ) > 1 ) THEN
         strlen  = index(scradir,' ') - 1
         dirname = scradir(1:strlen) // '/' // dirname
      END IF 

      CALL create_directory( dirname )

      DO i = 1, nk
        WRITE( kindex, fmt='( I5.5 )' ) i
        kdirname = TRIM( dirname ) // '/K' // kindex
        CALL create_directory( kdirname )
      END DO

      !  Some default values
      !
      kunit = 1
      k1 = 0
      k2 = 0
      k3 = 0
      nk1 = 1
      nk2 = 1
      nk3 = 1

      CALL invmat( 3, ht, htm1, omega )

      !  Open XML descriptor

      IF ( ionode ) THEN
        write(stdout,*) "Opening file "//trim(xmlpun)
        IF( ascii ) THEN
          call iotk_open_write(iunpun,file=TRIM(dirname)//'/'//TRIM(xmlpun),binary=.false.)
        ELSE
          call iotk_open_write(iunpun,file=TRIM(dirname)//'/'//TRIM(xmlpun),binary=.true.)
        ENDIF
      END IF

      IF( ionode ) THEN
        !
        call iotk_write_begin(iunpun,"STATUS")
        call iotk_write_attr (attr,"nfi",nfi,first=.true.)
        call iotk_write_empty(iunpun,"STEP",attr)
        call iotk_write_dat (iunpun, "TIME", simtime)
        call iotk_write_dat (iunpun, "TITLE", TRIM(title) )
        call iotk_write_end(iunpun,"STATUS")
        !
        call iotk_write_begin(iunpun,"DIMENSIONS")
        call iotk_write_attr (attr,"nx",nr1,first=.true.)
        call iotk_write_attr (attr,"ny",nr2)
        call iotk_write_attr (attr,"nz",nr3)
        call iotk_write_empty(iunpun,"FFT_GRID",attr)
        call iotk_write_attr (attr,"nx",nr1s,first=.true.)
        call iotk_write_attr (attr,"ny",nr2s)
        call iotk_write_attr (attr,"nz",nr3s)
        call iotk_write_empty(iunpun,"SMOOTH_FFT_GRID",attr)
        call iotk_write_attr (attr,"nx",nr1b,first=.true.)
        call iotk_write_attr (attr,"ny",nr2b)
        call iotk_write_attr (attr,"nz",nr3b)
        call iotk_write_empty(iunpun,"SMALLBOX_FFT_GRID",attr)
        call iotk_write_attr (attr,"ng",ngmt,first=.true.)
        call iotk_write_empty(iunpun,"GVECTORS",attr)
        call iotk_write_attr (attr,"ng",ngst,first=.true.)
        call iotk_write_empty(iunpun,"SMOOTH_GVECTORS",attr)
        call iotk_write_end(iunpun,"DIMENSIONS")
        !
        call iotk_write_begin(iunpun,"ELECTONS")
        call iotk_write_attr (attr,"nspin",nspin,first=.true.)
        call iotk_write_attr (attr,"nbnd",nbnd)
        call iotk_write_attr (attr,"nel",nelt)
        call iotk_write_empty(iunpun,"BANDS",attr)
        call iotk_write_attr (attr,"nel",nel(1),first=.true.)
        call iotk_write_attr (attr,"nbnd",nupdwn(1))
        call iotk_write_empty(iunpun,"SPINUP",attr)
        call iotk_write_attr (attr,"nel",nel(2),first=.true.)
        call iotk_write_attr (attr,"nbnd",nupdwn(2))
        call iotk_write_empty(iunpun,"SPINDW",attr)
        call iotk_write_end(iunpun,"ELECTONS")
        ! 
        call iotk_write_begin(iunpun,"AVERAGES")
        call iotk_write_dat (iunpun, "ACCUMULATORS", acc)
        call iotk_write_end(iunpun,"AVERAGES")
        !
        call iotk_write_attr (attr,"unit","Hartree",first=.true.)
        call iotk_write_begin(iunpun,"CUTOFFS",attr)
        call iotk_write_dat (iunpun, "ECUTRHO", ecutp/2.0d0)
        call iotk_write_dat (iunpun, "ECUTWFC", ecutw/2.0d0)
        call iotk_write_end(iunpun,"CUTOFFS")
        !
        call iotk_write_begin(iunpun,"LATTICE")
        call iotk_write_attr (attr, "ibrav",ibrav,first=.true.)
        call iotk_write_attr (attr, "alat", alat )
        call iotk_write_empty(iunpun,"TYPE",attr)
        call iotk_write_dat (iunpun, "CELLDM", celldm(1:6))
        call iotk_write_dat (iunpun, "ht", ht)
        call iotk_write_dat (iunpun, "htm", htm)
        call iotk_write_dat (iunpun, "htm2", htm2)
        call iotk_write_dat (iunpun, "htvel", htvel)
        call iotk_write_begin(iunpun,"NOSE")
          call iotk_write_dat (iunpun, "xnhh0", xnhh0)
          call iotk_write_dat (iunpun, "xnhhm", xnhhm)
          call iotk_write_dat (iunpun, "xnhhp", xnhhp)
          call iotk_write_dat (iunpun, "vnhh", vnhh)
        call iotk_write_end(iunpun,"NOSE")
        call iotk_write_begin(iunpun,"RECIPROCAL_BASIS")
          call iotk_write_dat (iunpun, "b1", b1)
          call iotk_write_dat (iunpun, "b2", b2)
          call iotk_write_dat (iunpun, "b3", b3)
        call iotk_write_end(iunpun,"RECIPROCAL_BASIS")
        call iotk_write_end(iunpun,"LATTICE")
        ! 
        call iotk_write_attr (attr, "nsp", nsp, first=.true.)
        call iotk_write_attr (attr, "nat", nat )
        call iotk_write_begin(iunpun,"IONS", attr)
        DO i = 1, nsp
          call iotk_write_attr (attr, "label",atm(i),first=.true.)
          call iotk_write_attr (attr, "na",na(i) )
          call iotk_write_attr (attr, "zv", zv(i) )
          call iotk_write_attr (attr, "mass", pmass(i) )
          call iotk_write_empty(iunpun,"SPECIE",attr)
        END DO
        call iotk_write_dat(iunpun, "taui", taui(1:3,1:nat) )
        call iotk_write_dat(iunpun, "cdmi", cdmi(1:3) )
        call iotk_write_dat(iunpun, "stau0", stau0(1:3,1:nat) )
        call iotk_write_dat(iunpun, "svel0", svel0(1:3,1:nat) )
        call iotk_write_dat(iunpun, "staum", staum(1:3,1:nat) )
        call iotk_write_dat(iunpun, "svelm", svelm(1:3,1:nat) )
        call iotk_write_dat(iunpun, "force", force(1:3,1:nat) )
        call iotk_write_begin(iunpun,"NOSE")
          call iotk_write_dat (iunpun, "xnhp0", xnhp0)
          call iotk_write_dat (iunpun, "xnhpm", xnhpm)
          call iotk_write_dat (iunpun, "xnhpm2", xnhpm2)
          call iotk_write_dat (iunpun, "xnhpp", xnhpp)
          call iotk_write_dat (iunpun, "vnhp", vnhp)
        call iotk_write_end(iunpun,"NOSE")
        call iotk_write_end(iunpun,"IONS")
        ! 
        call iotk_write_begin(iunpun,"ELECTRONS")
        call iotk_write_begin(iunpun,"NOSE")
          call iotk_write_dat (iunpun, "xnhe0", xnhe0)
          call iotk_write_dat (iunpun, "xnhem", xnhem)
          call iotk_write_dat (iunpun, "xnhem2", xnhem2)
          call iotk_write_dat (iunpun, "xnhep", xnhep)
          call iotk_write_dat (iunpun, "vnhe", vnhe)
        call iotk_write_end(iunpun,"NOSE")
        call iotk_write_end(iunpun,"ELECTRONS")
        !
        call iotk_write_begin(iunpun,"K_POINTS")
        call iotk_write_attr (attr, "nk",nk,first=.true.)
        call iotk_write_attr (attr, "gamma_only", gamma_only )
        call iotk_write_attr (attr, "kunit", kunit )
        call iotk_write_attr (attr, "k1", k1 )
        call iotk_write_attr (attr, "k2", k2 )
        call iotk_write_attr (attr, "k3", k3 )
        call iotk_write_attr (attr, "nk1", nk1 )
        call iotk_write_attr (attr, "nk2", nk2 )
        call iotk_write_attr (attr, "nk3", nk3 )
        call iotk_write_empty(iunpun,"GRID",attr)
        DO i = 1, nk
          call iotk_write_attr (attr, "ik",i,first=.true.)
          call iotk_write_attr (attr,"xyz",xk(:,i))
          call iotk_write_attr (attr,"w",wk(i))
          call iotk_write_empty(iunpun,"K",attr)
        END DO
        call iotk_write_end(iunpun,"K_POINTS")
        !
      END IF

      !
      !  K points dependent quantities stored in directory K00000
      !
      DO i = 1, nk
        ! 
        !  Occupation number
        !
        WRITE( kindex, fmt='( I5.5 )' ) i
        kdirname = TRIM( dirname ) // '/K' // kindex
        !
        filename = TRIM( kdirname ) // '/occupations.dat'
        IF( ionode ) THEN
          OPEN( unit = 10, file = TRIM(filename), status = 'UNKNOWN', form = 'UNFORMATTED' )
        END IF
        !
        IF( ionode ) THEN
          DO ispin = 1, nspin
            WRITE( 10 ) occ0( :, i, ispin )
            WRITE( 10 ) occm( :, i, ispin )
          END DO
        END IF
        !
        IF( ionode ) CLOSE( unit = 10 )
        ! 
        !  Wave functions
        !
        nspin_wfc = nspin
        IF( force_pairing ) nspin_wfc = 1
        wfscal = SQRT( omega )
        !
        DO ispin = 1, nspin_wfc
          !
          WRITE( filename, fmt='( I1 )' ) ispin
          filename = TRIM( kdirname ) // '/wf' // TRIM( filename ) // '_t0.dat'
          IF( ionode ) THEN
            OPEN( unit = 10, file = TRIM(filename), status = 'UNKNOWN', form = 'UNFORMATTED' )
          END IF
          CALL write_wfc( 10, i, nk, kunit, ispin, nspin_wfc, wfscal, c0(:,:,i,ispin), ngwt, nbnd, ig_l2g, ngw )
          IF( ionode ) CLOSE( unit = 10 )
          !
          WRITE( filename, fmt='( I1 )' ) ispin
          filename = TRIM( kdirname ) // '/wf' // TRIM( filename ) // '_tm.dat'
          IF( ionode ) THEN
            OPEN( unit = 10, file = TRIM(filename), status = 'UNKNOWN', form = 'UNFORMATTED' )
          END IF
          CALL write_wfc( 10, i, nk, kunit, ispin, nspin_wfc, wfscal, cm(:,:,i,ispin), ngwt, nbnd, ig_l2g, ngw )
          IF( ionode ) CLOSE( unit = 10 )
          !
        END DO
        !
        !
      END DO

      !   
      !  Write G vectors (up to ecutrho) to file
      !   
      ALLOCATE( mill(3,ngmt) )
      mill = 0
      DO ig = 1, ngm
        mill(:,ig_l2g(ig)) = mill_l(:,ig)
      END DO
      CALL mp_sum( mill )
      !
      filename = TRIM( dirname ) // '/gvectors.dat'
      IF( ionode ) THEN
        OPEN( unit = 10, file = TRIM(filename), status = 'UNKNOWN', form = 'UNFORMATTED' )
        WRITE( 10 ) mill
        CLOSE( unit = 10 )
      END IF
      DEALLOCATE( mill )

      !
      !  Write matrix lambda to file
      !
      filename = TRIM( dirname ) // '/lambda.dat'
      IF( ionode ) THEN
        OPEN( unit = 10, file = TRIM(filename), status = 'UNKNOWN', form = 'UNFORMATTED' )
        WRITE( 10 ) lambda0
        WRITE( 10 ) lambdam
        CLOSE( unit = 10 )
      END IF

      if( ionode ) then
        call iotk_close_write( iunpun )
      end if 

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!

    SUBROUTINE cp_readfile( scradir, ascii, nfi, simtime, acc, nk, xk, wk, &
        ht, htm, htm2, htvel, taui, cdmi, stau0, svel0, staum, svelm, force, &
        occ0, occm, lambda0, lambdam, b1, b2, b3, c0, cm  )
      !
      USE iotk_module
      USE kinds, ONLY: dbl
      USE io_global, ONLY: ionode, ionode_id, stdout
      USE parser, ONLY: int_to_char
      USE control_flags, ONLY: ndr, gamma_only, force_pairing
      USE io_files, ONLY: iunpun, xmlpun
      USE printout_base, ONLY: title
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b
      USE gvecp, ONLY: ngm, ngmt, ecutp
      USE gvecs, ONLY: ngs, ngst
      USE gvecw, ONLY: ngw, ngwt, ecutw
      USE electrons_base, ONLY: nspin, nbnd, nelt, nel, nupdwn, iupdwn
      USE electrons_nose, ONLY: xnhe0, xnhem, xnhep, xnhem2, vnhe
      USE cell_base, ONLY: ibrav, alat, celldm
      USE cell_nose, ONLY: xnhh0, xnhhm, xnhhp, vnhh
      USE ions_base, ONLY: nsp, nat, na, atm, zv, pmass
      USE ions_nose, ONLY: vnhp, xnhp0, xnhpm, xnhpp, xnhpm2
      USE reciprocal_vectors, ONLY: ngwt, ngw, ig_l2g, mill_l
      USE mp, ONLY: mp_sum, mp_bcast

      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: scradir
      LOGICAL, INTENT(IN) :: ascii
      INTEGER, INTENT(INOUT) :: nfi    !  index of the current step
      REAL(dbl), INTENT(INOUT) :: simtime   !  simulated time
      REAL(dbl), INTENT(INOUT) :: acc(:)   !
      INTEGER, INTENT(IN) :: nk    !  number of kpoints
      REAL(dbl), INTENT(INOUT) :: xk(:,:) ! k points coordinates
      REAL(dbl), INTENT(INOUT) :: wk(:)   ! k points weights
      REAL(dbl), INTENT(INOUT) :: ht(3,3) !
      REAL(dbl), INTENT(INOUT) :: htm(3,3) !
      REAL(dbl), INTENT(INOUT) :: htm2(3,3) !
      REAL(dbl), INTENT(INOUT) :: htvel(3,3) !
      REAL(dbl), INTENT(INOUT) :: taui(:,:) !
      REAL(dbl), INTENT(INOUT) :: cdmi(:) !
      REAL(dbl), INTENT(INOUT) :: stau0(:,:) !
      REAL(dbl), INTENT(INOUT) :: svel0(:,:) !
      REAL(dbl), INTENT(INOUT) :: staum(:,:) !
      REAL(dbl), INTENT(INOUT) :: svelm(:,:) !
      REAL(dbl), INTENT(INOUT) :: force(:,:) ! 
      REAL(dbl), INTENT(INOUT) :: occ0(:,:,:) !
      REAL(dbl), INTENT(INOUT) :: occm(:,:,:) !
      REAL(dbl), INTENT(INOUT) :: lambda0(:,:) !
      REAL(dbl), INTENT(INOUT) :: lambdam(:,:) !
      REAL(dbl), INTENT(INOUT) :: b1(3) !
      REAL(dbl), INTENT(INOUT) :: b2(3) !
      REAL(dbl), INTENT(INOUT) :: b3(3) !
      COMPLEX(dbl), INTENT(INOUT) :: c0(:,:,:,:) !
      COMPLEX(dbl), INTENT(INOUT) :: cm(:,:,:,:) !
      !
      CHARACTER(LEN=256) :: dirname, kdirname, filename
      CHARACTER(LEN=5) :: kindex
      INTEGER :: strlen
      CHARACTER(iotk_attlenx) :: attr
      INTEGER :: kunit
      INTEGER :: k1, k2, k3
      INTEGER :: nk1, nk2, nk3
      INTEGER :: i, ispin, ig, nspin_wfc, ierr
      INTEGER, ALLOCATABLE :: mill(:,:)
      REAL(dbl) :: wfscal, omega, htm1(3,3)
      !
      ! Variables read for testing pourposes
      !
      INTEGER :: ibrav_
      INTEGER :: nat_ , nsp_, na_
      INTEGER :: nk_ , ik_
      LOGICAL :: gamma_only_
      REAL(dbl) :: alat_
      REAL(dbl) :: pmass_ , zv_
      REAL(dbl) :: celldm_ ( 6 )
      REAL(dbl) :: wfscal_
      INTEGER :: ispin_ , nspin_ , ngwt_ , nbnd_ , nelt_

      dirname = 'RESTART' // int_to_char( ndr )
      IF ( LEN( scradir ) > 1 ) THEN
         strlen  = index(scradir,' ') - 1
         dirname = scradir(1:strlen) // '/' // dirname
      END IF
      filename = TRIM( dirname ) // '/' // 'restart.xml'
     
      IF( ionode ) THEN
        CALL iotk_open_read( iunpun, file = TRIM( filename ), binary = .FALSE., root = attr )
      END IF

      ierr = 0
      IF( ionode ) THEN
        !
        call iotk_scan_begin(iunpun,"STATUS")
        call iotk_scan_empty(iunpun,"STEP",attr)
        call iotk_scan_attr (attr,"nfi",nfi)
        call iotk_scan_dat (iunpun, "TIME", simtime )
        call iotk_scan_dat (iunpun, "TITLE", title )
        call iotk_scan_end(iunpun,"STATUS")
        !
        call iotk_scan_begin(iunpun,"ELECTONS")
        call iotk_scan_empty(iunpun,"BANDS",attr)
        call iotk_scan_attr (attr,"nspin",nspin_ )
        call iotk_scan_attr (attr,"nbnd",nbnd_ )
        call iotk_scan_attr (attr,"nel",nelt_ )
        !call iotk_scan_empty(iunpun,"SPINUP",attr)
        !call iotk_scan_attr (attr,"nel",nel(1))
        !call iotk_scan_attr (attr,"nbnd",nupdwn(1))
        !call iotk_scan_empty(iunpun,"SPINDW",attr)
        !call iotk_scan_attr (attr,"nel",nel(2))
        !call iotk_scan_attr (attr,"nbnd",nupdwn(2))
        call iotk_scan_end(iunpun,"ELECTONS")
        !
        IF( ( nspin_ /= nspin ) .OR. ( nbnd_ /= nbnd ) .OR. ( nelt_ /= nelt ) ) THEN 
          ierr = 30
          GOTO 100
        END IF
        ! 
        !
        call iotk_scan_begin(iunpun,"AVERAGES")
        call iotk_scan_dat (iunpun, "ACCUMULATORS", acc)
        call iotk_scan_end(iunpun,"AVERAGES")
        !
        call iotk_scan_begin(iunpun,"LATTICE")
        call iotk_scan_empty(iunpun,"TYPE",attr)
        call iotk_scan_attr (attr, "ibrav",ibrav_ )
        call iotk_scan_attr (attr, "alat", alat_ )
        call iotk_scan_dat (iunpun, "CELLDM", celldm_ (1:6) )
        call iotk_scan_dat (iunpun, "ht", ht)
        call iotk_scan_dat (iunpun, "htm", htm)
        call iotk_scan_dat (iunpun, "htm2", htm2)
        call iotk_scan_dat (iunpun, "htvel", htvel)
        call iotk_scan_begin(iunpun,"NOSE")
          call iotk_scan_dat (iunpun, "xnhh0", xnhh0)
          call iotk_scan_dat (iunpun, "xnhhm", xnhhm)
          call iotk_scan_dat (iunpun, "xnhhp", xnhhp)
          call iotk_scan_dat (iunpun, "vnhh", vnhh)
        call iotk_scan_end(iunpun,"NOSE")
        call iotk_scan_begin(iunpun,"RECIPROCAL_BASIS")
          call iotk_scan_dat (iunpun, "b1", b1)
          call iotk_scan_dat (iunpun, "b2", b2)
          call iotk_scan_dat (iunpun, "b3", b3)
        call iotk_scan_end(iunpun,"RECIPROCAL_BASIS")
        call iotk_scan_end(iunpun,"LATTICE")
        !
        call iotk_scan_begin(iunpun,"IONS", attr)
        call iotk_scan_attr (attr, "nsp", nsp_ )
        call iotk_scan_attr (attr, "nat", nat_ )
        IF( nsp_ /= nsp .OR. nat_ /= nat ) THEN 
          ierr = 10
          GOTO 100
        END IF
        !
        DO i = 1, nsp
          call iotk_scan_empty(iunpun, "SPECIE", attr)
          ! call iotk_scan_attr (attr, "label", atm(i) )
          call iotk_scan_attr (attr, "na", na_ )
          IF( na_ /= na( i ) ) THEN 
            ierr = 11
            GOTO 100
          END IF
          call iotk_scan_attr (attr, "zv", zv_ )
          call iotk_scan_attr (attr, "mass", pmass_ )
        END DO
        call iotk_scan_dat(iunpun, "taui", taui(1:3,1:nat) )
        call iotk_scan_dat(iunpun, "cdmi", cdmi(1:3) )
        call iotk_scan_dat(iunpun, "stau0", stau0(1:3,1:nat) )
        call iotk_scan_dat(iunpun, "svel0", svel0(1:3,1:nat) )
        call iotk_scan_dat(iunpun, "staum", staum(1:3,1:nat) )
        call iotk_scan_dat(iunpun, "svelm", svelm(1:3,1:nat) )
        call iotk_scan_dat(iunpun, "force", force(1:3,1:nat) )
        call iotk_scan_begin(iunpun,"NOSE")
          call iotk_scan_dat (iunpun, "xnhp0", xnhp0)
          call iotk_scan_dat (iunpun, "xnhpm", xnhpm)
          call iotk_scan_dat (iunpun, "xnhpm2", xnhpm2)
          call iotk_scan_dat (iunpun, "xnhpp", xnhpp)
          call iotk_scan_dat (iunpun, "vnhp", vnhp)
        call iotk_scan_end(iunpun,"NOSE")
        call iotk_scan_end(iunpun,"IONS")
        ! 
        call iotk_scan_begin(iunpun,"ELECTRONS")
        call iotk_scan_begin(iunpun,"NOSE")
          call iotk_scan_dat (iunpun, "xnhe0", xnhe0)
          call iotk_scan_dat (iunpun, "xnhem", xnhem)
          call iotk_scan_dat (iunpun, "xnhem2", xnhem2)
          call iotk_scan_dat (iunpun, "xnhep", xnhep)
          call iotk_scan_dat (iunpun, "vnhe", vnhe)
        call iotk_scan_end(iunpun,"NOSE")
        call iotk_scan_end(iunpun,"ELECTRONS")
        !
        call iotk_scan_begin(iunpun,"K_POINTS")
        call iotk_scan_empty(iunpun,"GRID",attr)
        call iotk_scan_attr (attr, "nk",nk_ )
        call iotk_scan_attr (attr, "gamma_only", gamma_only_ )
        !call iotk_scan_attr (attr, "kunit", kunit )
        !call iotk_scan_attr (attr, "k1", k1 )
        !call iotk_scan_attr (attr, "k2", k2 )
        !call iotk_scan_attr (attr, "k3", k3 )
        !call iotk_scan_attr (attr, "nk1", nk1 )
        !call iotk_scan_attr (attr, "nk2", nk2 )
        !call iotk_scan_attr (attr, "nk3", nk3 )
        IF( ( nk_ /= nk ) .OR.  &
            ( gamma_only .AND. .NOT. gamma_only_ ) .OR. &
            ( .NOT. gamma_only .AND. gamma_only_ ) ) THEN 
          ierr = 20
          GOTO 100
        END IF
        DO i = 1, nk
          call iotk_scan_empty(iunpun,"K",attr)
          call iotk_scan_attr (attr, "ik",ik_ )
          ! call iotk_scan_attr (attr,"xyz",xk(:,i))
          ! call iotk_scan_attr (attr,"w",wk(i))
        END DO
        call iotk_scan_end(iunpun,"K_POINTS")
        !
      END IF
      !
 100  CONTINUE
      !
      CALL mp_bcast( ierr, ionode_id )
      CALL mp_bcast( attr, ionode_id )
      IF( ierr /= 0 ) THEN
        CALL errore( " cp_readfile ", attr, ierr )
      END IF
      !
      CALL mp_bcast( nfi, ionode_id )
      CALL mp_bcast( simtime, ionode_id )
      CALL mp_bcast( title, ionode_id )
      CALL mp_bcast( acc, ionode_id )
      !
      CALL mp_bcast( ht, ionode_id )
      CALL mp_bcast( htm, ionode_id )
      CALL mp_bcast( htm2, ionode_id )
      CALL mp_bcast( htvel, ionode_id )
      CALL mp_bcast( xnhh0, ionode_id )
      CALL mp_bcast( xnhhm, ionode_id )
      CALL mp_bcast( xnhhp, ionode_id )
      CALL mp_bcast( vnhh, ionode_id )
      CALL mp_bcast( b1, ionode_id )
      CALL mp_bcast( b2, ionode_id )
      CALL mp_bcast( b3, ionode_id )
      !
      CALL mp_bcast(stau0, ionode_id)
      CALL mp_bcast(svel0, ionode_id)
      CALL mp_bcast(staum, ionode_id)
      CALL mp_bcast(svelm, ionode_id)
      CALL mp_bcast(taui, ionode_id)
      CALL mp_bcast(force, ionode_id)
      CALL mp_bcast(cdmi, ionode_id)
      CALL mp_bcast(xnhp0, ionode_id)
      CALL mp_bcast(xnhpm, ionode_id)
      CALL mp_bcast(xnhpm2, ionode_id)
      CALL mp_bcast(xnhpp, ionode_id)
      CALL mp_bcast(vnhp, ionode_id)
      !
      CALL mp_bcast(xnhe0, ionode_id)
      CALL mp_bcast(xnhem, ionode_id)
      CALL mp_bcast(xnhem2, ionode_id)
      CALL mp_bcast(xnhep, ionode_id)
      CALL mp_bcast(vnhe, ionode_id)

      IF( ionode ) THEN
        CALL iotk_close_read( iunpun )
      END IF

      !
      !  K points dependent quantities stored in directory K00000
      !
      DO i = 1, nk
        ! 
        !  Occupation number
        !
        WRITE( kindex, fmt='( I5.5 )' ) i
        kdirname = TRIM( dirname ) // '/K' // kindex
        !
        filename = TRIM( kdirname ) // '/occupations.dat'
        IF( ionode ) THEN
          OPEN( unit = 10, file = TRIM(filename), status = 'OLD', form = 'UNFORMATTED' )
        END IF
        !
        IF( ionode ) THEN
          DO ispin = 1, nspin
            READ( 10 ) occ0( :, i, ispin )
            READ( 10 ) occm( :, i, ispin )
          END DO
        END IF
        !
        IF( ionode ) CLOSE( unit = 10 )
        ! 
        !  Wave functions
        !
        nspin_wfc = nspin
        IF( force_pairing ) nspin_wfc = 1
        !
        DO ispin = 1, nspin_wfc
          !
          WRITE( filename, fmt='( I1 )' ) ispin
          filename = TRIM( kdirname ) // '/wf' // TRIM( filename ) // '_t0.dat'
          IF( ionode ) THEN
            OPEN( unit = 10, file = TRIM(filename), status = 'OLD', form = 'UNFORMATTED' )
          END IF
          CALL read_wfc( 10, i, nk, kunit, ispin_ , nspin_ , wfscal_ , c0(:,:,i,ispin), ngwt_ , nbnd_ , ig_l2g, ngw )
          IF( ionode ) CLOSE( unit = 10 )
          !
          WRITE( filename, fmt='( I1 )' ) ispin
          filename = TRIM( kdirname ) // '/wf' // TRIM( filename ) // '_tm.dat'
          IF( ionode ) THEN
            OPEN( unit = 10, file = TRIM(filename), status = 'OLD', form = 'UNFORMATTED' )
          END IF
          CALL read_wfc( 10, i, nk, kunit, ispin_ , nspin_ , wfscal_ , cm(:,:,i,ispin), ngwt_ , nbnd_ , ig_l2g, ngw )
          IF( ionode ) CLOSE( unit = 10 )
          !
        END DO
        !
        !
      END DO


      RETURN
    END SUBROUTINE



!
!
!=----------------------------------------------------------------------------=!

!
! ..  This subroutine write wavefunctions to the disk
!

    SUBROUTINE write_wfc(iuni, ik, nk, kunit, ispin, nspin, scal, wf, ngw, nbnd, igl, ngwl )
!
      USE kinds, ONLY: dbl
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_get, mp_bcast, mp_max
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool, my_image_id
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      INTEGER, INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(dbl), INTENT(IN) :: wf(:,:)
      INTEGER, INTENT(IN) :: ngw   ! 
      INTEGER, INTENT(IN) :: nbnd
      INTEGER, INTENT(IN) :: ngwl
      INTEGER, INTENT(IN) :: igl(:)
      REAL(dbl), INTENT(IN) :: scal

      INTEGER :: i, j, ierr
      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx
      INTEGER :: npool, ipmask( nproc ), ipsour
      COMPLEX(dbl), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)

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
        !
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

        ! get the maximum G vector index within the pool
        !
        CALL mp_max( igwx, intra_pool_comm ) 

        ! now notify all procs if an error has been found 
        !
        CALL mp_max( ierr ) 

        IF( ierr > 0 ) &
          CALL errore(' write_wfc ',' wrong size ngl ', ierr )

        IF( ipsour /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipsour, 1 )
        END IF

        IF( ionode ) WRITE(iuni) ngw, nbnd, ik, nk, kunit, ispin, nspin, scal
        IF( ionode ) WRITE(iuni) igwx

        ALLOCATE( wtmp( MAX(igwx,1) ) )
        wtmp = 0.0d0

        DO j = 1, nbnd
          IF( npool > 1 ) THEN
            IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
              CALL mergewf( wf(:,j), wtmp, ngwl, igl, me_pool, nproc_pool, root_pool, intra_pool_comm )
            END IF
            IF( ipsour /= ionode_id ) THEN
              CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j )
            END IF
          ELSE
            CALL mergewf( wf(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id)
          END IF
          IF( ionode ) WRITE(iuni) ( wtmp( i ), i=1,igwx )
        END DO

        DEALLOCATE( wtmp )

      RETURN
    END SUBROUTINE

!
!
!
!=----------------------------------------------------------------------------=!

    SUBROUTINE read_wfc(iuni, ik, nk, kunit, ispin, nspin, scal, wf, ngw, nbnd, igl, ngwl )
!
      USE kinds, ONLY: dbl
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_put, mp_bcast, mp_max, mp_get
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool, my_image_id
      USE io_global, ONLY: ionode, ionode_id
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: iuni
      COMPLEX(dbl), INTENT(INOUT) :: wf(:,:)
      INTEGER, INTENT(IN) :: ik, nk, kunit
      INTEGER, INTENT(OUT) :: ngw, nbnd, ispin, nspin
      REAL(dbl), INTENT(OUT) :: scal
      INTEGER, INTENT(IN) :: ngwl
      INTEGER, INTENT(IN) :: igl(:)

      INTEGER :: i, j
      COMPLEX(dbl), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)
      INTEGER :: ierr

      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx, igwx_
      INTEGER :: ik_, nk_, kunit_
      INTEGER :: npool, ipmask( nproc ), ipdest
!
! ... Subroutine Body
!

        IF( ionode ) READ(iuni) ngw, nbnd, ik_, nk_, kunit_, ispin, nspin, scal
        IF( ionode ) READ(iuni) igwx_

        CALL mp_bcast( ngw, ionode_id )
        CALL mp_bcast( nbnd, ionode_id )
        CALL mp_bcast( ik_, ionode_id )
        CALL mp_bcast( nk_, ionode_id )
        CALL mp_bcast( kunit_, ionode_id )
        CALL mp_bcast( ispin, ionode_id )
        CALL mp_bcast( nspin, ionode_id )
        CALL mp_bcast( scal, ionode_id )
        CALL mp_bcast( igwx_, ionode_id )

        ! set working variables for k point index (ikt) and k points number (nkt)
        ikt = ik
        nkt = nk

        !  find out the number of pools
        npool = nproc / nproc_pool

        !  find out number of k points blocks (each block contains kunit k points)
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
        ipdest = ionode_id

        !  find out the index of the processor which collect the data in the pool of ik
        IF( npool > 1 ) THEN
          IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            IF( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
          END IF
          CALL mp_sum( ipmask )
          DO i = 1, nproc
            IF( ipmask(i) == 1 ) ipdest = ( i - 1 )
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
          CALL errore(' read_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipdest /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipdest, 1 )
        END IF

        ALLOCATE( wtmp( MAX(igwx_, igwx) ) )

        DO j = 1, nbnd

            IF( j <= SIZE( wf, 2 ) ) THEN

              IF( ionode ) READ(iuni) ( wtmp(i), i=1,igwx_ )
              IF( igwx > igwx_ ) wtmp( (igwx_ + 1) : igwx ) = 0.0d0
 
              IF( npool > 1 ) THEN
                IF( ipdest /= ionode_id ) THEN
                  CALL mp_put( wtmp, wtmp, mpime, ionode_id, ipdest, j )
                END IF
                IF( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
                  CALL splitwf(wf(:,j), wtmp, ngwl, igl, me_pool, nproc_pool, root_pool, intra_pool_comm)
                END IF
              ELSE
                CALL splitwf(wf(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id)
              END IF

            END IF

        END DO

        DEALLOCATE( wtmp )

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------------=!
!

!------------------------------------------------------------------------------!
  END MODULE cp_restart
!------------------------------------------------------------------------------!
