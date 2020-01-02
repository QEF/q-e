!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE read_conf_from_file( stop_on_error, nat, nsp, tau, at )
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE constants,       ONLY : eps8
  USE io_global,       ONLY : stdout, ionode, ionode_id
  USE io_files,        ONLY : restart_dir, xmlfile, &
                              psfile, pseudo_dir, pseudo_dir_cur
  USE mp,              ONLY : mp_bcast
  USE mp_images,       ONLY : intra_image_comm
  USE qexsd_module,    ONLY : qexsd_readschema
  USE qexsd_copy,      ONLY : qexsd_copy_atomic_species, &
                              qexsd_copy_atomic_structure
  USE qes_types_module,ONLY : output_type
  USE qes_libs_module, ONLY : qes_reset
  USE qes_bcast_module,ONLY : qes_bcast
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in)    :: stop_on_error
  INTEGER, INTENT(in)    :: nat
  INTEGER, INTENT(in)    :: nsp
  REAL(DP),INTENT(out)   :: at(3,3)
  REAL(DP),INTENT(inout) :: tau(3,nat)
  !
  TYPE ( output_type) :: output_obj
  !
  INTEGER :: ierr  
  REAL(dp) :: alat_
  !
  CHARACTER (LEN=3) :: atm_(nsp)
  INTEGER, ALLOCATABLE :: ityp_(:)
  REAL(dp), ALLOCATABLE :: tau_(:,:)
  REAL(dp):: amass_(nsp)
  INTEGER :: nat_, ibrav_
  !
  pseudo_dir_cur = restart_dir () 
  WRITE( stdout, '(/5X,"Atomic positions and unit cell read from directory:", &
                &  /,5X,A)') pseudo_dir_cur
  !
  ! ... check if restart file is present, if so read config parameters
  !
  IF (ionode) CALL qexsd_readschema ( xmlfile(), ierr, output_obj )
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF ( ierr > 0 .OR. (ierr < 0 .AND. stop_on_error) ) &
       CALL errore ( 'read_conf_from_file', &
       'fatal error reading xml file', ABS(ierr) ) 
  CALL qes_bcast(output_obj, ionode_id, intra_image_comm)
  !
  IF (ierr == 0 ) THEN
     ! FIXME: what to do with the following data?
     ! CALL qexsd_copy_atomic_species ( output_obj%atomic_species, &
     !      nsp, atm, amass, PSFILE=psfile, PSEUDO_DIR=pseudo_dir )
     !IF ( pseudo_dir == ' ' ) pseudo_dir=pseudo_dir_cur
     !
     CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
          atm_, nat_, tau_, ityp_, alat_, at(:,1), at(:,2), at(:,3), ibrav_ )
     CALL qes_reset (output_obj)
     IF ( nat_ /= nat ) CALL errore('read_conf_from_file','bad number of atoms',1)
     at(:,:) = at(:,:) / alat_
     tau_(:,1:nat) = tau_(:,1:nat)/alat_
     IF ( SUM ( (tau_(:,1:nat)-tau(:,1:nat))**2 ) > eps8 ) THEN
        WRITE( stdout, '(5X,"Atomic positions from file used, from input discarded")' )
        tau(:,1:nat) = tau_(:,1:nat)
     END IF
     DEALLOCATE ( tau_, ityp_ )
     !
  ELSE
     !
     WRITE( stdout, '(5X,"Nothing found: ", &
                       & "using input atomic positions and unit cell",/)' )
     RETURN
     !
  END IF
  !
  WRITE( stdout, * )
  !
  RETURN
  !
END SUBROUTINE read_conf_from_file
