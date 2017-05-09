!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM importexportbinary
  !-----------------------------------------------------------------------
  !
  ! Program for exporting the charge density to a non-binary, text XML format
  ! and to import it back to binary for restarting.
  ! This code converts the charge density file, the spin polarization file, 
  ! and manually copies in the new output directory the data-file.xml and the
  ! UPF files stored in the prefix.save directory.
  ! All the rest (in particular: K??????/eigvals.xml files) is not copied
  ! and must be copied by hand, if needed.
  !
  ! Routine written by G. Pizzi, EPFL (2015) 
  !
  ! The input format is described both below and in INPUT_IMPORTEXPORT_BINARY
  ! in the PP/doc folder.
  !
  ! It accepts in input a &INPUTPP namelist, with the following variables:
  ! * outdir: the outdir directory that contains the charge density you want to
  !           convert 
  ! * prefix: the prefix used in the pw.x jobs
  ! * newoutdir: the new outdir to which you want to copy the files in the new
  !           (binary or non-binary) format
  ! * direction: a string, must be either 'export' or 'import', any other value
  !           is not valid. 'export' means that in outdir there is the (standard)
  !           binary charge density, and in newoutdir the text charge density
  !           will be written; 'import' instead is used to read from outdir
  !           a previously exported charge density and write back in newoutdir
  !           the binary format to be used to restart a QE calculation.
  !           Default: 'export'
  !
  USE io_global,  ONLY : ionode
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start, environment_end

  !
  IMPLICIT NONE
  !
  ! Initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'IMPORTEXPORT' )
  !
  ! Read the input file
  ! 
  IF ( ionode )  CALL input_from_file ( )
  !
  ! Do the actual job
  !
  CALL impexp ()
  !
  ! Exit
  ! 
  CALL environment_end ( 'IMPORTEXPORT' )
  !
  CALL stop_pp()
  !
END PROGRAM importexportbinary
!
!-----------------------------------------------------------------------
SUBROUTINE impexp ()
  !-----------------------------------------------------------------------

  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, prefix, psfile, pseudo_dir, xmlpun
  USE ions_base, ONLY : nsp
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE io_rho_xml,    ONLY : write_scf
  USE scf,           ONLY : rho
  USE lsda_mod,      ONLY : nspin
  USE xml_io_base,   ONLY : rho_binary, create_directory
  USE wrappers, ONLY: f_copy


  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  character (len=1), external :: capital

  INTEGER :: ios, len, l

  CHARACTER(len=256)  :: outdir, newoutdir, new_tmp_dir, direction
  CHARACTER(len=256)  :: directionupper, filename, old_tmp_dir
  CHARACTER(len=1024) :: sourcef, destf

  LOGICAL :: export ! Internal value to track if we are importing or exporting

  NAMELIST / inputpp / outdir, prefix, newoutdir, direction

  !
  !   set default values for variables in namelist
  !
  direction = 'export'
  prefix = 'pwscf'
  newoutdir = ' '
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     new_tmp_dir = trimcheck (newoutdir)
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, world_comm)
  !
  IF ( ios /= 0) CALL errore ('importexport', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( new_tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( direction, ionode_id, world_comm )

  ! Store original tmp_dir to preserve it, we will need it later
  old_tmp_dir = tmp_dir

  ! check direction, set default newoutdir if not specified
  len = len_trim(direction)
  directionupper = ' '
  do l = 1, len
     directionupper (l:l) = capital (direction(l:l) )
  enddo

  if (trim(directionupper) .eq. "IMPORT") then
     export = .FALSE.
  else if (trim(directionupper) .eq. "EXPORT") then
     export = .TRUE.
  else
     CALL errore ('importexport', 'wrong direction (must be either IMPORT or EXPORT)', 1)
  end if
 
  ! Set a sensible newoutdir if none was specified
  if ( trim( newoutdir ) == ' ' ) then
     if (export .eqv. .TRUE. ) then
        newoutdir = './export/'
     else
        newoutdir = './import/'
     end if
  end if

  ! Set the binary flag to the proper value for READING
  ! NOTE! requires to change the definition of rho_binary
  ! from LOGICAL, PARAMETER to LOGICAL, SAVE!
  if ( export .eqv. .TRUE. ) then
     rho_binary = .TRUE.
  else     
     rho_binary = .FALSE. 
  end if

  ! Read XML and the charge densities; don't read the eigenval.xml
  ! files that are not needed for import/export purposes
  CALL read_xml_file_nobs
  
  ! A trick here:
  ! I now change the out dir to the new one, to store there the outputs
  ! (this is done on all processors, I already broadcasted the values)
  tmp_dir = new_tmp_dir

  ! If the new output directory does not exist, I create it
  CALL create_directory( tmp_dir )
  
  ! Set the binary flag to the proper value for WRITING
  ! (i.e., I invert the value of rho_binary)
  if ( export .eqv. .TRUE. ) then
     rho_binary = .FALSE.
  else     
     rho_binary = .TRUE.
  end if

  ! Now I can store the new charge density in the proper binary/non-binary format
  CALL write_scf(rho, nspin)
  
  ! I need to copy XML file
  filename =  TRIM( xmlpun )
  sourcef = TRIM( old_tmp_dir ) // TRIM( prefix ) // '.save/' // TRIM( filename )
  destf   = TRIM( new_tmp_dir ) // TRIM( prefix ) // '.save/' // TRIM( filename )
  ios = f_copy( TRIM( sourcef ), TRIM( destf ))
  IF ( ios /= 0) CALL errore ('importexport', 'copying the '//TRIM(filename)//' file', abs(ios))

  ! I also need to copy the UPF files
  do l=1, nsp
     sourcef = TRIM( old_tmp_dir ) // TRIM( prefix ) // '.save/' // TRIM(psfile(l))
     destf = TRIM( new_tmp_dir ) // TRIM( prefix ) // '.save/' // TRIM(psfile(l))
     ios = f_copy( TRIM( sourcef ), TRIM( destf ))
     IF ( ios /= 0) CALL errore ('importexport', 'copying the ' // TRIM(psfile(l)) // ' pseudo', abs(ios))
  end do

  ! I am not copying the K?????/eigenval.xml files, this must be done by hand if needed
  if (ionode) then
     write(*,'(3X,A)') "NOTE! The band structure (files eigenval.xml) are not copied!"
  end if

  ! Rho and the various files have been written, exit

END SUBROUTINE impexp
