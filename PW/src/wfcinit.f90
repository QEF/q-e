! 
! Copyright (C) 2001-2013 Quantum ESPRESSO group
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
  ! ... It also open needed files or memory buffers
  !
  USE io_global,            ONLY : stdout
  USE basis,                ONLY : natomwfc, starting_wfc
  USE bp,                   ONLY : lelfield
  USE klist,                ONLY : xk, nks, ngk, igk_k
  USE control_flags,        ONLY : io_level, lscf, twfcollect 
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE io_files,             ONLY : nwordwfc, nwordwfcU, iunhub, iunwfc,&
                                   diropn, tmp_dir, prefix
  USE buffers,              ONLY : open_buffer, get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, current_k
  USE wannier_new,          ONLY : use_wannier
#if defined (__XSD)
  USE pw_restart_new,       ONLY : pw_readschema_file, read_collected_to_evc 
#else
  USE pw_restart,           ONLY : pw_readfile
#endif
  USE mp_bands,             ONLY : nbgrp, root_bgrp,inter_bgrp_comm
  USE mp,                   ONLY : mp_bcast
  USE qes_types_module,            ONLY : output_type
  USE qes_libs_module,             ONLY : qes_reset_output
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ierr
  LOGICAL :: exst, exst_mem, exst_file, opnd_file, twfcollect_file = .FALSE.
  CHARACTER (256 )                        :: dirname
#if defined (__XSD) 
  TYPE ( output_type )                    :: output_obj
#endif 
  !
  !
  !
  CALL start_clock( 'wfcinit' )
  !
  ! ... Orthogonalized atomic functions needed for LDA+U and other cases
  !
  IF ( use_wannier .OR. one_atom_occupations ) CALL orthoatwfc ( use_wannier )
  IF ( lda_plus_u ) CALL orthoUwfc()
  !
  ! ... open files/buffer for wavefunctions (nwordwfc set in openfil)
  ! ... io_level > 1 : open file, otherwise: open buffer
  !
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
  !
  ! ... now the various possible wavefunction initializations
  ! ... first a check: is "tmp_dir"/"prefix".wfc found on disk?
  !
  IF ( TRIM(starting_wfc) == 'file' .AND. .NOT. exst_file) THEN
     !
     ! ... "tmp_dir"/"prefix".wfc not found on disk: try to read
     ! ... wavefunctions in "collected" format from "prefix".save/, 
     ! ... rewrite them (in pw_readfile) using the internal format
     !
     ierr = 1
#if defined(__XSD)
     !ALLOCATE ( output_obj) 
     CALL pw_readschema_file(IERR = ierr, RESTART_OUTPUT = output_obj )
     IF ( ierr == 0 ) THEN 
        twfcollect_file = output_obj%band_structure%wf_collected   
        dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save' 
        IF ( twfcollect_file ) CALL read_collected_to_evc(dirname )
     END IF 
     CALL qes_reset_output ( output_obj ) 
     !DEALLOCATE ( output_obj )    
#else
     CALL pw_readfile( 'wave', ierr )
#endif
     IF ( ierr > 0 ) THEN
        WRITE( stdout, '(5X,"Cannot read wfc : file not found")' )
        starting_wfc = 'atomic+random'
     END IF
     !
     ! ... workaround: with k-point parallelization and 1 k-point per pool,
     ! ... pw_readfile does not leave evc properly initialized on all pools
     !
     IF ( nks == 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
     !
  ELSE IF ( TRIM(starting_wfc) == 'file' .AND. exst_file) THEN
     !
     ! ... wavefunctions are read from file (or buffer) in routine 
     ! ... c_bands, but not if there is a single k-point. In such
     ! ... a case, we read wavefunctions (directly from file in 
     ! ... order to avoid a useless buffer allocation) here
     !
     IF ( nks == 1 ) THEN
         inquire (unit = iunwfc, opened = opnd_file)
         if (.not.opnd_file) CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
         CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, -1 )
         if(.not.opnd_file) CLOSE ( UNIT=iunwfc, STATUS='keep' )
     END IF
     !
  END IF
  !
  ! ... state what will happen
  !
  IF ( TRIM(starting_wfc) == 'file' ) THEN
     !
     WRITE( stdout, '(5X,"Starting wfc from file")' )
     !
  ELSE IF ( starting_wfc == 'atomic' ) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfc are ",I4," atomic wfcs")' ) natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfc are ",I4," atomic + ", &
             &           I4," random wfc")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE IF ( TRIM(starting_wfc) == 'atomic+random' .AND. natomwfc > 0) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfc are ",I4," randomized atomic wfcs")')&
             natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfc are ",I4," randomized atomic wfcs + "&
             &          ,I4," random wfc")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE
     !
     WRITE( stdout, '(5X,"Starting wfc are random")' )
     !
  END IF
  !
  ! ... exit here if starting from file or for non-scf calculations.
  ! ... In the latter case the starting wavefunctions are not 
  ! ... calculated here but just before diagonalization (to reduce I/O)
  !
  IF (  ( .NOT. lscf .AND. .NOT. lelfield ) .OR. TRIM(starting_wfc) == 'file' ) THEN
     !
     CALL stop_clock( 'wfcinit' )
     RETURN
     !
  END IF
  !
  ! ... calculate and write all starting wavefunctions to file
  !
  DO ik = 1, nks
     !
     ! ... Hpsi initialization: k-point index, spin, kinetic energy
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     call g2_kin (ik)
     !
     ! ... More Hpsi initialization: nonlocal pseudopotential projectors |beta>
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... Needed for LDA+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
        CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )
     !
     ! ... calculate starting wavefunctions (calls Hpsi)
     !
     CALL init_wfc ( ik )
     !
     ! ... write  starting wavefunctions to file
     !
     IF ( nks > 1 .OR. (io_level > 1) .OR. lelfield ) &
         CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
  END DO
  !
  CALL stop_clock( 'wfcinit' )
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
  USE klist,                ONLY : xk, ngk, igk_k
  USE wvfct,                ONLY : nbnd, npwx, et
  USE uspp,                 ONLY : nkb, okvan
  USE noncollin_module,     ONLY : npol
  USE wavefunctions_module, ONLY : evc
  USE random_numbers,       ONLY : randy
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm, my_bgrp_id
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ik
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
     CALL start_clock( 'wfcinit:atomic' )
     CALL atomic_wfc( ik, wfcatom )
     CALL stop_clock( 'wfcinit:atomic' )
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
               DO ig = 1, ngk(ik)
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
        DO ig = 1, ngk(ik)
           !
           rr  = randy()
           arg = tpi * randy()
           !
           wfcatom(ig,ipol,ibnd) = &
                CMPLX( rr*COS( arg ), rr*SIN( arg ) ,kind=DP) / &
                       ( ( xk(1,ik) + g(1,igk_k(ig,ik)) )**2 + &
                         ( xk(2,ik) + g(2,igk_k(ig,ik)) )**2 + &
                         ( xk(3,ik) + g(3,igk_k(ig,ik)) )**2 + 1.0_DP )
        END DO
        !
     END DO
     !
  END DO
  
  ! when band parallelization is active, the first band group distributes
  ! the wfcs to the others making sure all bgrp have the same starting wfc

  if (my_bgrp_id > 0) wfcatom(:,:,:) = (0.d0,0.d0)
  call mp_sum(wfcatom,inter_bgrp_comm)
  !
  ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
  !
  ALLOCATE( etatom( n_starting_wfc ) )
  !
  ! ... Allocate space for <beta|psi>
  !
  CALL allocate_bec_type ( nkb, n_starting_wfc, becp, intra_bgrp_comm )
  !
  ! ... the following trick is for electric fields with Berry's phase:
  ! ... by setting lelfield = .false. one prevents the calculation of
  ! ... electric enthalpy in the Hamiltonian (cannot be calculated
  ! ... at this stage: wavefunctions at previous step are missing)
  ! 
  lelfield_save = lelfield
  lelfield = .FALSE.
  !
  ! ... subspace diagonalization (calls Hpsi)
  !
  CALL start_clock( 'wfcinit:wfcrot' )
  CALL rotate_wfc ( npwx, ngk(ik), n_starting_wfc, gstart, &
                    nbnd, wfcatom, npol, okvan, evc, etatom )
  CALL stop_clock( 'wfcinit:wfcrot' )
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
