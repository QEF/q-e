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

!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    SUBROUTINE create_directories( dirname )
      USE io_global, ONLY: ionode, ionode_id, stdout
      USE mp, ONLY: mp_bcast
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=256) :: dirname
      INTEGER :: ios
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

    SUBROUTINE cp_writefile( scradir, ascii, nfi, simtime, acc, nk, xk, wk )
      USE iotk_module
      USE io_global, ONLY: ionode, ionode_id, stdout
      USE parser, ONLY: int_to_char
      USE control_flags, ONLY: ndw, gamma_only
      USE io_files, ONLY: iunpun, xmlpun
      USE printout_base, ONLY: title
      USE kinds, ONLY: dbl
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b
      USE gvecp, ONLY: ngm, ecutp
      USE gvecs, ONLY: ngs
      USE gvecw, ONLY: ecutw
      USE electrons_base, ONLY: nspin, nbnd, nelt, nel, nupdwn, iupdwn
      USE cell_base, ONLY: ibrav, alat, celldm
      USE ions_base, ONLY: nsp, nat, na, atm, zv, pmass
      
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: scradir
      LOGICAL, INTENT(IN) :: ascii
      INTEGER, INTENT(IN) :: nfi    !  index of the current step
      REAL(dbl), INTENT(IN) :: simtime   !  simulated time
      REAL(dbl), INTENT(IN) :: acc(:)   !  
      INTEGER, INTENT(IN) :: nk    !  number of kpoints
      REAL(dbl), INTENT(IN) :: xk(:,:) ! k points coordinates 
      REAL(dbl), INTENT(IN) :: wk(:)   ! k points weights
      !
      CHARACTER(LEN=256) :: dirname
      INTEGER :: strlen
      CHARACTER(iotk_attlenx) :: attr
      INTEGER :: kunit
      INTEGER :: k1, k2, k3
      INTEGER :: nk1, nk2, nk3
      INTEGER :: i

      dirname = 'RESTART' // int_to_char( ndw )
      IF ( LEN( scradir ) > 1 ) THEN
         strlen  = index(scradir,' ') - 1
         dirname = scradir(1:strlen) // '/' // dirname
      END IF 

      CALL create_directories( dirname )

      !  Some default values
      !
      kunit = 1
      k1 = 0
      k2 = 0
      k3 = 0
      nk1 = 1
      nk2 = 1
      nk3 = 1


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
        call iotk_write_attr (attr,"ng",ngm,first=.true.)
        call iotk_write_empty(iunpun,"GVECTORS",attr)
        call iotk_write_attr (attr,"ng",ngs,first=.true.)
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
        call iotk_write_begin(iunpun,"LATTICE")
        call iotk_write_attr (attr, "ibrav",ibrav,first=.true.)
        call iotk_write_attr (attr, "alat", alat )
        call iotk_write_dat (iunpun, "CELLDM", celldm, attr)
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
        call iotk_write_end(iunpun,"IONS")
        

      
      END IF

      if( ionode ) then
        call iotk_close_write( iunpun )
      end if 

      RETURN
    END SUBROUTINE


!------------------------------------------------------------------------------!
  END MODULE cp_restart
!------------------------------------------------------------------------------!
