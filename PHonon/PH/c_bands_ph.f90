!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
! Questa e' una copia di c_bands_nscf intesa per un confronto con 
! thermo_pw. 
!
!
SUBROUTINE c_bands_nscf_ph( )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine for Hamiltonian diagonalization routines
  ! ... specialized to non-self-consistent calculations (no electric field)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer, open_buffer
  USE basis,                ONLY : starting_wfc
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_lr,           ONLY : lgamma
  USE control_flags,        ONLY : ethr, restart, isolve, io_level, iverbosity
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions,        ONLY : evc
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE check_stop,           ONLY : check_stop_now
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : domag
  USE save_ph,              ONLY : tmp_dir_save
  USE io_files,             ONLY : tmp_dir, prefix
  !
  IMPLICIT NONE
  !
  REAL(DP) :: avg_iter, ethr_
  ! average number of H*psi products
  INTEGER :: ik_, ik, nkdum, ios, iuawfc, lrawfc
  ! ik_: k-point already done in a previous run
  ! ik : counter on k points
  LOGICAL :: exst, exst_mem
  !
  REAL(DP), EXTERNAL :: get_clock
  !
  CALL start_clock( 'c_bands' )
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands(ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  !
  DO ik = 1, ik_
     CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
  END DO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSE IF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")')
  ELSE IF ( isolve == 2 ) THEN
     WRITE( stdout, '(5X,"PPCG style diagonalization")')
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve)
  END IF
  IF (tmp_dir /= tmp_dir_save) THEN
     iuawfc = 20
     lrawfc = nbnd * npwx * npol
     CALL open_buffer (iuawfc, 'wfc', lrawfc, io_level, exst_mem, exst, &
                                                         tmp_dir_save)
     IF (.NOT.exst.AND..NOT.exst_mem) THEN
        CALL errore ('c_bands_ph', 'file '//trim(prefix)//'.wfc not found', 1)
     END IF
  ENDIF
  !
  ! ... For each k point (except those already calculated if restarting)
  ! ... diagonalizes the hamiltonian
  !
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     call g2_kin( ik )
     ! 
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... Needed for LDA+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
          CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )
     !
     ! ... calculate starting  wavefunctions
     !
     IF ( iverbosity > 0 ) WRITE( stdout, 9001 ) ik
     !
     IF ( TRIM(starting_wfc) == 'file' ) THEN
        !
        CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
        !
     ELSE
        !
        CALL init_wfc ( ik )
        !
     END IF
     !
     ! ... diagonalization of bands for k-point ik
     !
     call diag_bands ( 1, ik, avg_iter )
     !
     !  In the noncolinear magnetic case we have k, k+q, -k -k-q and
     !  to the last two wavefunctions we must apply t_rev.
     !  When lgamma is true we have only k and -k
     !
     IF (noncolin.AND.domag) THEN
        IF (lgamma.AND. MOD(ik,2)==0) THEN
           CALL apply_trev(evc, ik, ik-1)
        ELSEIF (.NOT.lgamma.AND.(MOD(ik,4)==3.OR.MOD(ik,4)==0)) THEN
           CALL apply_trev(evc, ik, ik-2)
        ENDIF
     ENDIF
     !
     ! ... save wave-functions (unless disabled in input)
     !
     IF ( io_level > -1 ) CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
     ! ... beware: with pools, if the number of k-points on different
     ! ... pools differs, make sure that all processors are still in
     ! ... the loop on k-points before checking for stop condition
     !
     nkdum  = kunit * ( nkstot / kunit / npool )
     IF (ik .le. nkdum) THEN
        !
        ! ... stop requested by user: save restart information,
        ! ... save wavefunctions to file
        !
        IF (check_stop_now()) THEN
           CALL save_in_cbands(ik, ethr, avg_iter, et )
           RETURN
        END IF
     ENDIF
     !
     ! report about timing
     !
     IF ( iverbosity > 0 ) THEN
        WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
        FLUSH( stdout )
     ENDIF
     !
  END DO k_loop
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  avg_iter = avg_iter / nkstot
  !
  WRITE( stdout, '(/,5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)' ) &
       ethr, avg_iter
  IF (tmp_dir /= tmp_dir_save) CALL close_buffer(iuawfc,'keep')
  !
  CALL stop_clock( 'c_bands' )
  !
  RETURN
  !
  ! formats
  !
9001 FORMAT(/'     Computing kpt #: ',I5 )
9000 FORMAT( '     total cpu time spent up to now is ',F10.1,' secs' )
  !
END SUBROUTINE c_bands_nscf_ph
