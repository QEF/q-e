!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE metadyn_io
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the I/O methods used by meta-dynamics
  !
  ! ... code written by Carlo Sbraccia (2005)
  !
  USE kinds, ONLY : DP
  !
  USE iotk_module
  !
  USE xml_io_base, ONLY : create_directory
  !
  IMPLICIT NONE
  !
  CHARACTER(iotk_attlenx) :: attr
  !
  PRIVATE
  !
  PUBLIC :: write_metadyn_restart, &
            read_metadyn_restart, &
            write_axsf_file
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_metadyn_restart( dirname, iter, tau, energy, pos_unit )
      !------------------------------------------------------------------------
      !
      USE metadyn_vars,       ONLY : max_metadyn_iter, g_amplitude, &
                                     gaussian_pos, fe_grad, fe_step
      USE constraints_module, ONLY : nconstr, target
      USE io_global,          ONLY : ionode, ionode_id
      USE mp_global,          ONLY : MPIME
      USE mp,                 ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      INTEGER,          INTENT(IN) :: iter
      REAL(DP),         INTENT(IN) :: tau(:,:)
      REAL(DP),         INTENT(IN) :: energy
      REAL(DP),         INTENT(IN) :: pos_unit
      !
      INTEGER            :: i, ia
      CHARACTER(LEN=256) :: filename, metadyn_dir
      INTEGER            :: iunit, ierr
      !
      !
      IF ( ionode ) THEN
         !
         ! ... look for an empty unit (only ionode needs it)
         !
         CALL iotk_free_unit( iunit, ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      CALL errore( 'write_metadyn_restart', &
                   'no free units to write the restart file', ierr )
      !
      ! ... the restart information is written in a sub-directory of
      ! .. the 'save' directory
      !
      CALL create_directory( dirname )
      !
      metadyn_dir = TRIM( dirname ) // '/meta-dynamics'
      !
      CALL create_directory( metadyn_dir )
      !
      filename = TRIM( metadyn_dir ) // '/' // "metadyn-descriptor.xml"
      !
      ! ... only ionode writes the file
      !
      IF ( .NOT. ionode ) RETURN
      !
      ! ... descriptor file
      !
      CALL iotk_open_write( iunit, FILE = filename, &
                            ROOT = "METADYNAMICS", BINARY = .FALSE. )
      !
      CALL iotk_write_dat( iunit, &
                           "NUM_OF_COLLECTIVE_VARIABLES", nconstr )
      !
      CALL iotk_write_dat( iunit, &
                           "NUM_OF_STEPS", max_metadyn_iter )
      !
      CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
      CALL iotk_write_dat( iunit, &
                           "GAUSSIAN_AMPLITUDE", g_amplitude, ATTR = attr )
      !
      CALL iotk_write_attr( attr, "UNITS", "depend on the" // &
                          & "type of collective variables", FIRST = .TRUE. )
      CALL iotk_write_dat( iunit, "GAUSSIAN_SPREAD", fe_step(:), ATTR = attr )
      !
      CALL iotk_write_dat( iunit, "STEP", iter )
      !
      DO i = 1, iter
         !
         filename = TRIM( metadyn_dir ) // '/' // &
                  & 'iteration' // TRIM( iotk_index( i ) ) // '.xml'
         !
         CALL iotk_link( iunit, "ITERATION" // TRIM( iotk_index( i ) ), &
                         filename, CREATE = .FALSE., BINARY = .FALSE. )
         !
      END DO
      !
      CALL iotk_close_write( iunit )
      !
      ! ... information about the last step
      !
      CALL iotk_open_write( iunit, FILE = filename, ROOT = 'iteration' // &
                          & TRIM( iotk_index( iter ) ), BINARY = .FALSE. )
      !
      CALL iotk_write_begin( iunit, "IONS" )
      !
      CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_empty( iunit, "UNITS_FOR_IONIC_POS", attr )
      !
      DO i = 1, SIZE( tau, DIM = 2 )
         !
         CALL iotk_write_attr( attr, "tau", tau(:,i)*pos_unit, FIRST = .TRUE. )
         CALL iotk_write_empty( iunit, &
                              & "ATOM" // TRIM( iotk_index( i ) ), attr )
         !
      END DO
      !
      CALL iotk_write_end( iunit, "IONS" )
      !
      CALL iotk_write_dat( iunit, "COLLECTIVE_VARIABLES", target(:) )
      !
      CALL iotk_write_dat( iunit, "GAUSSIAN_CENTERS", gaussian_pos(:) )
      !
      CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
      CALL iotk_write_dat( iunit, "POTENTIAL_ENERGY", energy, ATTR = attr )
      !
      CALL iotk_write_attr( attr, "UNITS", "Hartree / Bohr", FIRST = .TRUE. )
      CALL iotk_write_dat( iunit, &
                           "POTENTIAL_OF_MEAN_FORCE", fe_grad(:), ATTR = attr )
      !
      CALL iotk_close_write( iunit )
      !
      RETURN
      !
    END SUBROUTINE write_metadyn_restart
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_metadyn_restart( dirname, tau, pos_unit )
      !------------------------------------------------------------------------
      !
      USE metadyn_vars,       ONLY : g_amplitude, gaussian_pos, fe_grad, &
                                     metadyn_history, first_metadyn_iter
      USE constraints_module, ONLY : nconstr, target
      USE io_global,          ONLY : ionode, ionode_id
      USE mp,                 ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      REAL(DP),         INTENT(OUT) :: tau(:,:)
      REAL(DP),         INTENT(IN)  :: pos_unit
      !
      INTEGER            :: nconstr_in
      INTEGER            :: i, ia
      CHARACTER(LEN=256) :: filename, tag
      INTEGER            :: iunit, ierr
      !
      !
      ! ... look for an empty unit
      !
      IF ( ionode ) THEN
         !
         CALL iotk_free_unit( iunit, ierr )
         !
         CALL errore( 'read_metadyn_restart', &
                      'no free units to read the restart file', ierr )
         !
         filename = TRIM( dirname ) // &
                  & '/meta-dynamics/' // "metadyn-descriptor.xml"
         !
         ! ... descriptor file
         ! 
         CALL iotk_open_read( iunit, FILE = filename, IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      CALL errore( 'read_metadyn_restart', &
                   'restart file ' // TRIM( filename ) // ' not found', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_dat( iunit, &
                             "NUM_OF_COLLECTIVE_VARIABLES", nconstr_in )
         !
         IF ( nconstr_in == nconstr ) THEN
            !
            nconstr = nconstr_in
            !
         ELSE
            !
            CALL errore( 'read_metadyn_restart', &
                         'wrong number of collective variables', 1 )
            !
         END IF
         !
         CALL iotk_scan_dat( iunit, "STEP", first_metadyn_iter )
         !
         DO i = 1, first_metadyn_iter
            !
            tag = "ITERATION" // TRIM( iotk_index( i ) )
            !
            CALL iotk_scan_begin( iunit, TRIM( tag ) )
            !
            CALL iotk_scan_dat( iunit, &
                                "GAUSSIAN_CENTERS", metadyn_history(:,i) )
            !
            CALL iotk_scan_end( iunit, TRIM( tag ) )
            !
         END DO
         !
         CALL iotk_close_read( iunit )
         !
         ! ... information about the last step
         !
         CALL iotk_open_read( iunit, FILE = filename )
         !
         tag = "ITERATION" // TRIM( iotk_index( first_metadyn_iter ) )
         !
         CALL iotk_scan_begin( iunit, TRIM( tag ) )
         !
         CALL iotk_scan_begin( iunit, "IONS" )
         !
         DO i = 1, SIZE( tau, DIM = 2 )
            !
            CALL iotk_scan_empty( iunit, &
                                  "ATOM" // TRIM( iotk_index( i ) ), attr )
            CALL iotk_scan_attr( attr, "tau", tau(:,i) )
            !
         END DO
         !
         CALL iotk_scan_end( iunit, "IONS" )
         !
         CALL iotk_scan_dat( iunit, "COLLECTIVE_VARIABLES", target(:) )
         CALL iotk_scan_dat( iunit, "GAUSSIAN_CENTERS", gaussian_pos(:) )
         CALL iotk_scan_dat( iunit, "POTENTIAL_OF_MEAN_FORCE", fe_grad(:) )
         !
         CALL iotk_scan_end( iunit, TRIM( tag ) )
         !
         CALL iotk_close_read( iunit )
         !
         ! ... positions are converted to internal units
         !
         tau(:,:) = tau(:,:) / pos_unit
         !
      END IF
      !
      CALL mp_bcast( nconstr,            ionode_id )
      CALL mp_bcast( first_metadyn_iter, ionode_id )
      CALL mp_bcast( metadyn_history,    ionode_id )
      CALL mp_bcast( tau,                ionode_id )
      CALL mp_bcast( target,             ionode_id )
      CALL mp_bcast( gaussian_pos,       ionode_id )
      CALL mp_bcast( fe_grad,            ionode_id )
      !
      RETURN
      !
    END SUBROUTINE read_metadyn_restart
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_axsf_file( image, tau, tau_units )
      !------------------------------------------------------------------------
      !
      USE input_parameters, ONLY : atom_label
      USE io_files,         ONLY : iunaxsf
      USE constants,        ONLY : bohr_radius_angs
      USE ions_base,        ONLY : nat, ityp
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN) :: image
      REAL(DP), INTENT(IN) :: tau(:,:)
      REAL(DP), INTENT(IN) :: tau_units
      !
      INTEGER :: atom
      !
      !
      WRITE( UNIT = iunaxsf, FMT = '(" PRIMCOORD ",I5)' ) image
      WRITE( UNIT = iunaxsf, FMT = '(I5,"  1")' ) nat
      !
      DO atom = 1, nat
         !
         WRITE( UNIT = iunaxsf, FMT = '(A2,3(2X,F18.10))' ) &
                TRIM( atom_label(ityp(atom)) ), &
             tau(1,atom) * tau_units * bohr_radius_angs, &
             tau(2,atom) * tau_units * bohr_radius_angs, &
             tau(3,atom) * tau_units * bohr_radius_angs
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE write_axsf_file
    !
END MODULE metadyn_io
