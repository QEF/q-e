!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
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
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE basis,                ONLY : natomwfc, starting_wfc
  USE bp,                   ONLY : lelfield
  USE klist,                ONLY : xk, nks, ngk, igk_k
  USE control_flags,        ONLY : io_level, lscf
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU, lda_plus_u_kind
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE io_files,             ONLY : nwordwfc, nwordwfcU, iunhub, iunwfc,&
                                   diropn, xmlfile, restart_dir
  USE buffers,              ONLY : open_buffer, close_buffer, get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb
  USE wavefunctions,        ONLY : evc
  USE wvfct,                ONLY : nbnd, current_k
  USE wannier_new,          ONLY : use_wannier
  USE pw_restart_new,       ONLY : read_collected_wfc
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_images,            ONLY : intra_image_comm
  USE qexsd_module,         ONLY : qexsd_readschema
  USE qes_types_module,     ONLY : output_type
  USE qes_libs_module,      ONLY : qes_reset
  USE wavefunctions_gpum,   ONLY : using_evc
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ierr, exst_sum
  LOGICAL :: exst, exst_mem, exst_file, opnd_file, twfcollect_file
  CHARACTER (LEN=256)  :: dirname
  TYPE ( output_type ) :: output_obj
  !
  CALL start_clock( 'wfcinit' )
  CALL using_evc(0) ! this may be removed
  !
  ! ... Orthogonalized atomic functions needed for DFT+U and other cases
  !
  IF ( (use_wannier .OR. one_atom_occupations ) .AND. lda_plus_u ) &
       CALL errore ( 'wfcinit', 'currently incompatible options', 1 )
  IF ( use_wannier .OR. one_atom_occupations ) CALL orthoatwfc ( use_wannier )
  IF ( lda_plus_u ) CALL orthoUwfc(.FALSE.)
  !
  ! ... open files/buffer for wavefunctions (nwordwfc set in openfil)
  ! ... io_level > 1 : open file, otherwise: open buffer
  !
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
  !
  IF ( TRIM(starting_wfc) == 'file') THEN
     ! Check whether all processors have found a file when opening a buffer
     IF (exst_file) THEN
        exst_sum = 0
     ELSE
        exst_sum = 1
     END IF
     CALL mp_sum (exst_sum, intra_image_comm)
     !
     ! Check whether wavefunctions are collected (info in xml file)
     dirname = restart_dir ( )
     IF (ionode) CALL qexsd_readschema ( xmlfile(), ierr, output_obj )
     CALL mp_bcast(ierr, ionode_id, intra_image_comm)
     IF ( ierr <= 0 ) THEN
        ! xml file is valid
        IF (ionode) twfcollect_file = output_obj%band_structure%wf_collected
        CALL mp_bcast(twfcollect_file, ionode_id, intra_image_comm)
        CALL qes_reset  ( output_obj )
     ELSE
        ! xml file not found or not valid
        twfcollect_file = .FALSE.
     END IF
     !
     IF ( twfcollect_file ) THEN
        !
        DO ik = 1, nks
           CALL read_collected_wfc ( dirname, ik, evc )
           CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
        END DO
        !
     ELSE IF ( exst_sum /= 0 ) THEN
        !
        WRITE( stdout, '(5X,"Cannot read wfcs: file not found")' )
        IF (exst_file) THEN
           CALL close_buffer(iunwfc, 'delete')
           CALL open_buffer(iunwfc,'wfc', nwordwfc, io_level, exst_mem, exst_file)
        END IF
        starting_wfc = 'atomic+random'
        !
     ELSE
        !
        ! ... wavefunctions are read from file (or buffer) not here but
        !  ...in routine c_bands. If however there is a single k-point,
        ! ... c_bands doesn't read wavefunctions, so we read them here
        ! ... (directly from file to avoid a useless buffer allocation)
        !
        IF ( nks == 1 ) THEN
           INQUIRE (unit = iunwfc, opened = opnd_file)
           IF ( .NOT.opnd_file ) CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
           CALL using_evc(2)
           CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, -1 )
           IF ( .NOT.opnd_file ) CLOSE ( UNIT=iunwfc, STATUS='keep' )
        END IF
     END IF
  END IF
  !
  ! ... state what will happen
  !
  IF ( TRIM(starting_wfc) == 'file' ) THEN
     !
     WRITE( stdout, '(5X,"Starting wfcs from file")' )
     !
  ELSE IF ( starting_wfc == 'atomic' ) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," atomic wfcs")' ) natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," atomic + ", &
             &           I4," random wfcs")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE IF ( TRIM(starting_wfc) == 'atomic+random' .AND. natomwfc > 0) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," randomized atomic wfcs")')&
             natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," randomized atomic wfcs + "&
             &          ,I4," random wfcs")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE
     !
     WRITE( stdout, '(5X,"Starting wfcs are random")' )
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
  ! ... calculate and write all starting wavefunctions to buffer
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
     ! ... Needed for DFT+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
        CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )
     !
     ! DFT+U+V: calculate the phase factor at a given k point
     !
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) CALL phase_factor(ik)
     !
     ! ... calculate starting wavefunctions (calls Hpsi)
     !
     CALL init_wfc ( ik )
     !
     ! ... write  starting wavefunctions to file
     !
     IF ( nks > 1 .OR. (io_level > 1) .OR. lelfield ) CALL using_evc(0)
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
  USE basis,                ONLY : natomwfc, starting_wfc
  USE gvect,                ONLY : g, gstart
  USE klist,                ONLY : xk, ngk, igk_k
  USE wvfct,                ONLY : nbnd, npwx, et
  USE uspp,                 ONLY : nkb, okvan
  USE noncollin_module,     ONLY : npol
  USE wavefunctions,        ONLY : evc
  USE random_numbers,       ONLY : randy
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm, &
                                   nbgrp, root_bgrp_id
  USE mp,                   ONLY : mp_bcast
  USE xc_lib,               ONLY : xclib_dft_is, stop_exx
  !
  USE wavefunctions_gpum,   ONLY : using_evc
  USE wvfct_gpum,           ONLY : using_et
  USE becmod_subs_gpum,     ONLY : using_becp_auto
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
  IF ( n_starting_atomic_wfc > 0 ) THEN
     !
     CALL start_clock( 'wfcinit:atomic' ); !write(*,*) 'start wfcinit:atomic' ; FLUSH(6)
     CALL atomic_wfc( ik, wfcatom )
     CALL stop_clock( 'wfcinit:atomic' ); !write(*,*) 'stop wfcinit:atomic' ; FLUSH(6)
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
  ! FIXME: maybe this should be done once evc are computed, not here?
  !
  IF( nbgrp > 1 ) CALL mp_bcast( wfcatom, root_bgrp_id, inter_bgrp_comm )
  !
  ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
  !
  ALLOCATE( etatom( n_starting_wfc ) )
  !
  ! ... Allocate space for <beta|psi>
  !
  CALL allocate_bec_type ( nkb, n_starting_wfc, becp, intra_bgrp_comm )
  CALL using_becp_auto (2)
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
  IF ( xclib_dft_is('hybrid')  ) CALL stop_exx()
  CALL start_clock( 'wfcinit:wfcrot' ); !write(*,*) 'start wfcinit:wfcrot' ; FLUSH(6)
  CALL rotate_wfc ( npwx, ngk(ik), n_starting_wfc, gstart, nbnd, wfcatom, npol, okvan, evc, etatom )
  CALL using_evc(1)  ! rotate_wfc (..., evc, etatom) -> evc : out (not specified)
  CALL stop_clock( 'wfcinit:wfcrot' ); !write(*,*) 'stop wfcinit:wfcrot' ; FLUSH(6)
  !
  lelfield = lelfield_save
  !
  ! ... copy the first nbnd eigenvalues
  ! ... eigenvectors are already copied inside routine rotate_wfc
  !
  CALL using_et(1)
  et(1:nbnd,ik) = etatom(1:nbnd)
  !
  CALL deallocate_bec_type ( becp )
  CALL using_becp_auto (2)
  DEALLOCATE( etatom )
  DEALLOCATE( wfcatom )
  !
  RETURN
  !
END SUBROUTINE init_wfc
