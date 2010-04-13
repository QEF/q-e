! 
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE wfcinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes an estimate of the starting wavefunctions
  ! ... from superposition of atomic wavefunctions and/or random wavefunctions.
  !
  USE io_global,            ONLY : stdout
  USE basis,                ONLY : natomwfc, starting_wfc
  USE bp,                   ONLY : lelfield
  USE klist,                ONLY : xk, nks, ngk
  USE control_flags,        ONLY : io_level, lscf
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE io_files,             ONLY : nwordwfc, nwordatwfc, iunwfc, iunigk, iunsat
  USE buffers,              ONLY : get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, current_k, igk
  USE wannier_new,              ONLY : use_wannier
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !
  !
  CALL start_clock( 'wfcinit' )
  !
  ! ... Needed for LDA+U
  !
  IF ( lda_plus_u .OR. use_wannier .OR. one_atom_occupations ) CALL orthoatwfc()
  !
  ! ... state what is going to happen
  !
  IF ( TRIM(starting_wfc) == 'file' ) THEN
     !
     WRITE( stdout, '(5X,"Starting wfc from file")' )
     !
  ELSE IF ( starting_wfc == 'atomic' ) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        !
        WRITE( stdout, '(5X,"Starting wfc are ",I4," atomic wfcs")' ) natomwfc
        !
     ELSE
        !
        WRITE( stdout, '(5X,"Starting wfc are ",I4," atomic + ", &
             &           I4," random wfc")' ) natomwfc, nbnd-natomwfc
        !
     END IF
     !
  ELSE IF ( TRIM(starting_wfc) == 'atomic+random' ) THEN
     !
     WRITE( stdout, '(5X,"Starting wfc are ",I4," randomized atomic wfcs")' ) &
             natomwfc
     !
  ELSE
     !
     WRITE( stdout, '(5X,"Starting wfc are random")' )
     !
  END IF
  !
  ! ... for non-scf calculations, the starting wavefunctions are not 
  ! ... calculated here but immediately before diagonalization
  !
  IF ( .NOT. lscf .AND. .NOT. lelfield ) THEN
     !
     CALL stop_clock( 'wfcinit' )
     !
     RETURN
     !
  END IF
  !
  IF ( TRIM(starting_wfc) == 'file' ) THEN
     !
     ! ... wavefunctions are to be read from file: store wavefunction into
     ! ... memory if c_bands will not do it (for a single k-point);
     ! ... return and do nothing otherwise (c_bands will read wavefunctions)
     !
     IF ( nks == 1 .AND. (io_level < 2) ) &
        CALL get_buffer ( evc, nwordwfc, iunwfc, 1 )
     !
     CALL stop_clock( 'wfcinit' )
     !
     RETURN
     !
  END IF
  !
  IF ( nks > 1 ) REWIND( iunigk )
  !
  ! ... calculate and write all starting wavefunctions to file
  !
  DO ik = 1, nks
     !
     ! ... various initializations: k, spin, number of PW, indices
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk (ik)
     IF ( nks > 1 ) READ( iunigk ) igk
     !
     call g2_kin (ik)
     !
     ! ... Calculate nonlocal pseudopotential projectors |beta>
     !
     IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
     !
     ! ... LDA+U: read atomic wavefunctions for U term in Hamiltonian
     !
     IF ( lda_plus_u ) &
         CALL davcio( swfcatom, nwordatwfc, iunsat, ik, - 1 )
     !
     ! ... calculate starting wavefunctions
     !
     CALL init_wfc ( ik )
     !
     ! ... write  starting wavefunctions to file
     !
     IF ( nks > 1 .OR. (io_level > 1) ) &
         CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
  END DO
  !
  CALL stop_clock( 'wfcinit' )
  !
  RETURN
  !
END SUBROUTINE wfcinit
!
!----------------------------------------------------------------------------
SUBROUTINE init_wfc ( ik )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes starting wavefunctions for k-point ik
  !
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : lelfield
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp
  USE constants,            ONLY : tpi
  USE cell_base,            ONLY : tpiba2
  USE basis,                ONLY : natomwfc, starting_wfc
  USE gvect,                ONLY : g, gstart
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, et
  USE uspp,                 ONLY : nkb, okvan
  USE noncollin_module,     ONLY : npol
  USE wavefunctions_module, ONLY : evc
  USE random_numbers,       ONLY : randy
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !
  INTEGER :: ibnd, ig, ipol, n_starting_wfc, n_starting_atomic_wfc
  LOGICAL :: lelfield_save
  !
  REAL(DP) :: rr, arg
  REAL(DP), ALLOCATABLE :: etatom(:) ! atomic eigenvalues
  !
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:,:) ! atomic wfcs for initialization
  !
  !
  IF ( starting_wfc(1:6) == 'atomic' ) THEN
     !
     n_starting_wfc = MAX( natomwfc, nbnd )
     n_starting_atomic_wfc = natomwfc
     !
  ELSE IF ( starting_wfc == 'random' ) THEN
     !
     n_starting_wfc = nbnd
     n_starting_atomic_wfc = 0
     !
  ELSE
     !
     ! ...case 'file' should not be done here
     !
     CALL errore ( 'init_wfc', &
          'invalid value for startingwfc: ' // TRIM ( starting_wfc ) , 1 )
     !
  END IF
  !
  ALLOCATE( wfcatom( npwx, npol, n_starting_wfc ) )
  !
  IF ( starting_wfc(1:6) == 'atomic' ) THEN
     !
     CALL atomic_wfc( ik, wfcatom )
     !
     IF ( starting_wfc == 'atomic+random' .AND. &
         n_starting_wfc == n_starting_atomic_wfc ) THEN
         !
         ! ... in this case, introduce a small randomization of wavefunctions
         ! ... to prevent possible "loss of states"
         !
         DO ibnd = 1, n_starting_atomic_wfc
            !
            DO ipol = 1, npol
               !
               DO ig = 1, npw
                  !
                  rr  = randy()
                  arg = tpi * randy()
                  !
                  wfcatom(ig,ipol,ibnd) = wfcatom(ig,ipol,ibnd) * &
                     ( 1.0_DP + 0.05_DP * CMPLX( rr*COS(arg), rr*SIN(arg) ,kind=DP) ) 
                  !
               END DO
               !
            END DO
            !
         END DO
         !
     END IF
     !
  END IF
  !
  ! ... if not enough atomic wfc are available,
  ! ... fill missing wfcs with random numbers
  !
  DO ibnd = n_starting_atomic_wfc + 1, n_starting_wfc
     !
     DO ipol = 1, npol
        ! 
        wfcatom(:,ipol,ibnd) = (0.0_dp, 0.0_dp)
        !
        DO ig = 1, npw
           !
           rr  = randy()
           arg = tpi * randy()
           !
           wfcatom(ig,ipol,ibnd) = &
                CMPLX( rr*COS( arg ), rr*SIN( arg ) ,kind=DP) / &
                       ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                         ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                         ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.0_DP )
        END DO
        !
     END DO
     !
  END DO
  !
  ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
  !
  ALLOCATE( etatom( n_starting_wfc ) )
  !
  ! ... Allocate space for <beta|psi>
  !
  CALL allocate_bec_type ( nkb, n_starting_wfc, becp )
  !
  ! ... the following trick is for electric fields with Berry's phase:
  ! ... by setting lelfield = .false. one prevents the calculation of
  ! ... electric enthalpy in the Hamiltonian (cannot be calculated
  ! ... at this stage: wavefunctions at previous step are missing)
  ! 
  lelfield_save = lelfield
  lelfield = .FALSE.
  !
  CALL rotate_wfc ( npwx, npw, n_starting_wfc, gstart, &
                    nbnd, wfcatom, npol, okvan, evc, etatom )
  !
  lelfield = lelfield_save
  !
  ! ... copy the first nbnd eigenvalues
  ! ... eigenvectors are already copied inside routine rotate_wfc
  !
  et(1:nbnd,ik) = etatom(1:nbnd)
  !
  CALL deallocate_bec_type ( becp )
  DEALLOCATE( etatom )
  DEALLOCATE( wfcatom )
  !
  RETURN
  !
END SUBROUTINE init_wfc
