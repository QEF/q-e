!
! Copyright (C) 2002-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE trpmd_io_routines
  !----------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines used for I/O in path
  ! ... optimisations
  !
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : pi, autoev, bohr_radius_angs, eV_to_kelvin, rytoev
  USE io_global,  ONLY : meta_ionode, meta_ionode_id
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: write_output
  PUBLIC :: new_image_init, get_new_image, stop_other_images
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_output()
       !-----------------------------------------------------------------------
       !
       USE ring_io_units_module,  ONLY : iunpath
       USE ring_variables, ONLY :  path_length, & !num_of_images,
                                   pos !pes
       USE ring_formats,   ONLY : run_info, run_output
       USE ions_base,             ONLY : zv, ityp, nat
       USE fcp_variables,         ONLY : lfcpopt, fcp_mu
       USE fcp_opt_routines, ONLY : fcp_neb_ef, fcp_neb_nelec
       use pimd_variables, ONLY: pes, nbeadMD
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER   :: image
       REAL (DP) :: inter_image_distance
       !
       !
       IF ( .NOT. meta_ionode ) RETURN
       !
     !!!  WRITE( UNIT = iunpath, &    !!! <----my mod.
     !!!         FMT = '(/,5X,"activation energy (->) = ",F10.6," eV")' ) &  !!! <----my mod.
     !!!         activation_energy   !!! <----my mod.
     !!!  WRITE( UNIT = iunpath, &   !!! <----my mod.
     !!!         FMT = '(5X,"activation energy (<-) = ",F10.6," eV",/)' ) &   !!! <----my mod.
     !!!         activation_energy + ( pes(1) - pes(num_of_images) ) * autoev    !!! <----my mod.
       !
       WRITE( UNIT = iunpath, FMT = run_info )
       !
       path_length = 0.0_DP
       !
       DO image = 1, nbeadMD
          !
          IF ( image > 1 ) &
             path_length = path_length + &
                           norm( pos(:,image) - pos(:,image-1) )
          !
          WRITE( UNIT = iunpath, FMT = run_output ) &
              image, pes(image) * autoev!!!, error(image), frozen(image)   !!! <----my mod.
          !
       END DO
       !
       IF ( lfcpopt ) THEN
          WRITE(iunpath,'(/,5X,"image",2X,"Fermi energy (eV)",11X, &
                        & "error (V)",4X,"tot_charge",/)')
          DO image = 1, nbeadMD
          !
             WRITE(iunpath,'(5X,I5,9X,F10.6,10X,F10.6,4X,F10.6)') &
                 image, fcp_neb_ef(image)*rytoev, &
                 (fcp_mu-fcp_neb_ef(image))*rytoev, &
                 SUM( zv(ityp(1:nat)) ) - fcp_neb_nelec(image)
          !
          END DO
       END IF
       !
       inter_image_distance = path_length / DBLE( nbeadMD - 1 )
       !
     END SUBROUTINE write_output
     !
     !-----------------------------------------------------------------------
     SUBROUTINE new_image_init( nimage, fii, outdir )
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine initializes the file needed for the
       ! ... parallelization among images
       !
       USE ring_io_units_module, ONLY : iunnewimage
       USE io_files, ONLY : prefix
       USE ring_variables, ONLY : tune_load_balance
       !
       IMPLICIT NONE
       !
       INTEGER,          INTENT(IN) :: nimage, fii
       CHARACTER(LEN=*), INTENT(IN) :: outdir
       !
       !
       IF ( nimage == 1 .OR. .NOT.tune_load_balance ) RETURN
       !
       OPEN( UNIT = iunnewimage, FILE = TRIM( outdir ) // &
           & TRIM( prefix ) // '.newimage' , STATUS = 'UNKNOWN' )
       !
       WRITE( iunnewimage, * ) fii + nimage
       !
       CLOSE( UNIT = iunnewimage, STATUS = 'KEEP' )
       !
       RETURN
       !
     END SUBROUTINE new_image_init
     !
     !-----------------------------------------------------------------------
     SUBROUTINE get_new_image( nimage, image, outdir )
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine is used to get the new image to work on
       ! ... the "prefix.LOCK" file is needed to avoid (when present) that
       ! ... other jobs try to read/write on file "prefix.newimage"
       !
       USE io_files,       ONLY : iunnewimage, iunlock, prefix
       USE io_global,      ONLY : ionode
       USE ring_variables, ONLY : tune_load_balance
       !
       IMPLICIT NONE
       !
       INTEGER,          INTENT(IN)    :: nimage
       INTEGER,          INTENT(INOUT) :: image
       CHARACTER(LEN=*), INTENT(IN)    :: outdir
       !
       INTEGER            :: ioerr
       CHARACTER(LEN=256) :: filename
       LOGICAL            :: opened
       !
       !
       IF ( .NOT.ionode ) RETURN
       !
       IF ( nimage > 1 ) THEN
          !
          IF ( tune_load_balance ) THEN
             !
             filename = TRIM( outdir ) // TRIM( prefix ) // '.LOCK'
             !
             open_loop: DO
                !
                OPEN( UNIT = iunlock, FILE = TRIM( filename ), &
                     & IOSTAT = ioerr, STATUS = 'NEW' )
                !
                IF ( ioerr > 0 ) CYCLE open_loop
                !
                INQUIRE( UNIT = iunnewimage, OPENED = opened )
                !
                IF ( .NOT. opened ) THEN
                   !
                   OPEN( UNIT = iunnewimage, FILE = TRIM( outdir ) // &
                       & TRIM( prefix ) // '.newimage' , STATUS = 'OLD' )
                   !
                   READ( iunnewimage, * ) image
                   !
                   CLOSE( UNIT = iunnewimage, STATUS = 'DELETE' )
                   !
                   OPEN( UNIT = iunnewimage, FILE = TRIM( outdir ) // &
                       & TRIM( prefix ) // '.newimage' , STATUS = 'NEW' )
                   !
                   WRITE( iunnewimage, * ) image + 1
                   !
                   CLOSE( UNIT = iunnewimage, STATUS = 'KEEP' )
                   !
                   EXIT open_loop
                   !
                END IF
                !
             END DO open_loop
             !
             CLOSE( UNIT = iunlock, STATUS = 'DELETE' )
             !
          ELSE
             !
             image = image + nimage
             !
          END IF
          !
       ELSE
          !
          image = image + 1
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE get_new_image
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stop_other_images()
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine is used to send a stop signal to other images
       ! ... this is done by creating the exit_file on the working directory
       !
       USE io_files,  ONLY : iunexit, exit_file
       USE io_global, ONLY : ionode
       !
       IMPLICIT NONE
       !
       !
       IF ( .NOT. ionode ) RETURN
       !
       OPEN( UNIT = iunexit, FILE = TRIM( exit_file ) )
       CLOSE( UNIT = iunexit, STATUS = 'KEEP' )
       !
       RETURN
       !
     END SUBROUTINE stop_other_images
     !
END MODULE trpmd_io_routines

FUNCTION input_images_getarg( ) RESULT(input_images)
  !-----------------------------------------------------------------------------
  !
  ! check for command-line option "-input_images N" or "--input_images N",
  ! return N (0 if not found)
  !
  USE kinds,         ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER :: input_images
  CHARACTER(len=256) ::  myname
  INTEGER :: iiarg, nargs, i, i0
  !
  nargs = command_argument_count()
  input_images = 0
  !
  DO iiarg = 1, nargs
     !
     CALL get_command_argument( iiarg, myname)
     !
     IF ( TRIM( myname ) == '-input_images' .OR. &
          TRIM( myname ) == '--input_images' ) THEN
        !
        CALL get_command_argument( ( iiarg + 1 ) , myname )
        !
        READ(myname,*) input_images
        RETURN
        !
     END IF
     !
  ENDDO
  !
  RETURN
  !
END FUNCTION input_images_getarg

!----------------------------------------------------------------------------
SUBROUTINE close_io_units(myunit)
  !-----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: myunit
  !
  LOGICAL :: opnd
  !
  INQUIRE( UNIT = myunit, OPENED = opnd )
  IF ( opnd ) CLOSE( UNIT = myunit )
  !
END SUBROUTINE close_io_units
!
!----------------------------------------------------------------------------
SUBROUTINE open_io_units(myunit,file_name,lappend)
  !-----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: myunit
  CHARACTER(LEN=256), intent(in) :: file_name
  LOGICAL, intent(in) :: lappend
  !
  LOGICAL :: opnd
  !
  INQUIRE( UNIT = myunit, OPENED = opnd )
  IF ( opnd ) CLOSE( UNIT = myunit )
  OPEN( UNIT = myunit, FILE = TRIM(file_name), &
  STATUS = 'UNKNOWN', POSITION = 'APPEND' )
  !
END SUBROUTINE open_io_units
