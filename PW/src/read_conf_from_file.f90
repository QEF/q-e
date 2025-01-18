!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE read_conf_from_file( stop_on_error, nat, nsp, tau, alat, at, &
                                is_tau_read )
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE constants,       ONLY : eps8
  USE io_global,       ONLY : stdout, ionode, ionode_id
  USE io_files,        ONLY : restart_dir, xmlfile
  USE mp,              ONLY : mp_bcast
  USE mp_images,       ONLY : intra_image_comm
  USE qexsd_module,    ONLY : qexsd_readschema
  USE qexsd_copy,      ONLY : qexsd_copy_atomic_structure
  USE qes_types_module,ONLY : output_type
  USE qes_libs_module, ONLY : qes_reset
  USE qes_bcast_module,ONLY : qes_bcast
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in)    :: stop_on_error
  INTEGER, INTENT(in)    :: nat
  INTEGER, INTENT(in)    :: nsp
  REAL(DP),INTENT(out)   :: alat
  REAL(DP),INTENT(out)   :: at(3,3)
  REAL(DP),INTENT(inout) :: tau(3,nat)
  LOGICAL,INTENT(out)    :: is_tau_read
  !
  ! ... local variables
  !
  TYPE ( output_type) :: output_obj
  !
  INTEGER :: ierr, nat_, ibrav_, natomwfc_
  INTEGER, ALLOCATABLE :: ityp_(:)
  REAL(dp), ALLOCATABLE :: tau_(:,:)
  CHARACTER (LEN=6) :: atm_(nsp)
  !
  is_tau_read = .FALSE.
  !
  WRITE( stdout, '(/5X,"Atomic positions and unit cell read from directory:", &
                &  /,5X,A)') TRIM(restart_dir())
  !
  ! ... check if restart file is present, if so read config parameters
  !
  IF (ionode) CALL qexsd_readschema ( xmlfile(), ierr, output_obj )
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF ( ierr /= 0 .AND. stop_on_error ) CALL errore ( 'read_conf_from_file', &
       'fatal error reading xml file', ABS(ierr) )
  IF (ierr /= 0 ) THEN
     !
     WRITE( stdout, '(5X,"Nothing found: ", &
                       & "using input atomic positions and unit cell",/)' )
     !
  ELSE
     !
     CALL qes_bcast(output_obj, ionode_id, intra_image_comm)
     CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
          atm_, nat_, tau_, ityp_, alat, at(:,1), at(:,2), at(:,3), ibrav_ , natomwfc_)
     CALL qes_reset (output_obj)
     IF ( nat_ /= nat ) CALL errore('read_conf_from_file','bad number of atoms',1)
     at(:,:) = at(:,:) / alat
     tau_(:,1:nat) = tau_(:,1:nat)/alat
     IF ( SUM ( (tau_(:,1:nat)-tau(:,1:nat))**2 ) > eps8 ) THEN
        WRITE( stdout, '(5X,"Atomic positions from file used, from input discarded")' )
        tau(:,1:nat) = tau_(:,1:nat)
        !
        is_tau_read = .TRUE.
        !
     END IF
     DEALLOCATE ( tau_, ityp_ )
     WRITE( stdout, * )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE read_conf_from_file
