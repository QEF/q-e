!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE read_file()
  !----------------------------------------------------------------------------
  !
  ! ... This routine allocates space for all quantities already computed
  ! ... in the pwscf program and reads them from the data file.
  !
  !
  USE kinds,            ONLY : DP
  USE parameters,       ONLY : natx
  USE ions_base,        ONLY : nat, nsp, ityp, tau, if_pos
  USE basis,            ONLY : natomwfc
  USE cell_base,        ONLY : tpiba2, at, bg
  USE force_mod,        ONLY : force
  USE klist,            ONLY : nkstot, nks, xk, wk
  USE lsda_mod,         ONLY : lsda, nspin, current_spin, isk
  USE wvfct,            ONLY : nbnd, nbndx, et, wg
  USE symme,            ONLY : irt, nsym, ftau, s
  USE ktetra,           ONLY : tetra, ntetra 
  USE extfield,         ONLY : forcefield, tefield
  USE cellmd,           ONLY : cell_factor, lmovecell
  USE gvect,            ONLY : gg, ecutwfc, ngm, g, nr1, nr2, nr3, &
                               eigts1, eigts2, eigts3
  USE gsmooth,          ONLY : ngms, nls, nrx1s, nr1s, nr2s, nr3s
  USE scf,              ONLY : rho, vr
  USE vlocal,           ONLY : strf
  USE io_files,         ONLY : tmp_dir, prefix, iunpun
  USE restart_module,   ONLY : readfile_new
  USE noncollin_module, ONLY : noncolin, npol
  USE mp_global,        ONLY : kunit
  USE pw_restart,       ONLY : pw_readfile
  
  !
  IMPLICIT NONE
  !
  INTEGER               :: i, ik, ibnd, ios, ierr
  REAL(DP), ALLOCATABLE :: et_g(:,:), wg_g(:,:)
  REAL(DP)              :: rdum(1,1)
  !
  !
  ! ... a value of zero cause the parameter to be read from the ".save" file
  !
  kunit = 0
  !
  ! ... here we read the variables that dimension the system
  ! ... in parallel execution, only root proc reads the file
  ! ... and then broadcasts the values to all other procs
  !
#if defined(__NEWPUNCH)
  !
  ! ... a reset of the internal flgas is necessary because some codes call
  ! ... read_file() more than once
  !
  CALL pw_readfile( 'reset', ierr )
  CALL pw_readfile( 'dim',   ierr )
  !
  CALL errore( 'read_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.new-save', ierr )
  !
#else
  !
  CALL readfile_new( 'dim', iunpun, rdum, rdum, kunit, 0, 0, ierr )
  !
#endif
  !
  CALL errore( 'read_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  ! ... allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  IF ( nat <= 0 .OR. nat > natx ) &
     CALL errore( 'read_file', 'wrong number of atoms', 1 )
  !
  ! ... allocation
  !
  ALLOCATE( ityp( nat ) )
  !
  ALLOCATE( tau(    3, nat ) )
  ALLOCATE( if_pos( 3, nat ) )
  ALLOCATE( force(  3, nat ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  !
  ALLOCATE( irt( 48, nat ) )    
  ALLOCATE( tetra( 4, MAX( ntetra, 1 ) ) )    
  !
  ! ... here we read all the variables defining the system
  ! ... in parallel execution, only root proc read the file
  ! ... and then broadcast the values to all other procs
  !
#if defined(__NEWPUNCH)
  !
!-------------------------------------------------------------------------------
! ... XML punch-file
!-------------------------------------------------------------------------------
  !
  CALL set_dimensions()
  !
  ! ... parallel execution: distribute across pools k-points and
  ! ... related variables (not a smart implementation):
  ! ... nks and nkstot are redefined by the following routine
  !
  CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !
  current_spin = 1
  !
  cell_factor = 1.D0
  lmovecell = .FALSE.
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  CALL ggen()
  !
  ! ... allocate wavefunctions and related quantities (including et and wg)
  !
  nbndx = nbnd
  !
  CALL allocate_wfc()
  !
  CALL pw_readfile( 'nowave', ierr )
  !
  CALL poolscatter( nbnd , nkstot, et, nks, et )
  CALL poolscatter( nbnd , nkstot, wg, nks, wg )
  !
  CALL checkallsym( nsym, s, nat, tau, &
                    ityp, at, bg, nr1, nr2, nr3, irt, ftau )
  !
  ! ... read pseudopotentials
  !
  CALL pw_readfile( 'pseudo', ierr )
  !
  CALL readpp()
  !
  ! ... allocate the potential
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  !
#else
  !
!-------------------------------------------------------------------------------
! ... standard punch-file
!-------------------------------------------------------------------------------
  !
  ALLOCATE( et_g( nbnd, nkstot ), wg_g( nbnd, nkstot ) )
  !
  CALL readfile_new( 'nowave', iunpun, et_g, wg_g, kunit, 0, 0, ierr )
  !
  CALL errore( 'read_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  ! ... parallel execution: distribute across pools k-points and
  ! ... related variables (not a smart implementation)
  !
  nks = nkstot
  !
  ! ... nks and nkstot are redefined by the following routine
  !
  CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  cell_factor = 1.D0
  lmovecell = .FALSE.
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  CALL ggen()
  !
  ! ... allocate the potential
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  !
  ! ... allocate wavefunctions and related quantities (including et and wg)
  !
  nbndx = nbnd
  !
  CALL allocate_wfc()
  !
  et = et_g
  wg = wg_g
  !
  DEALLOCATE( et_g, wg_g )
  !
  CALL poolscatter( nbnd , nkstot, et, nks, et )
  CALL poolscatter( nbnd , nkstot, wg, nks, wg )
  !
#endif
  !
  ! ... read the charge density
  !
  CALL io_pot( - 1, 'rho', rho, nspin )
  !
  ! read the potential
  !
  CALL io_pot( - 1, 'pot', vr, nspin )
  !
  ! ... re-calculate the local part of the pseudopotential vltot
  ! ... and the core correction charge (if any) - This is done here
  ! ... for compatibility with the previous version of read_file
  !
  CALL init_vloc()
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                   nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
  !
  CALL setlocal()
  !
  CALL set_rhoc()
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_dimensions()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE gvect,     ONLY : ecutwfc, dual, gcutm
      USE gsmooth,   ONLY : gcutms, doublegrid
      USE klist,     ONLY : nks, nkstot
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutm = dual * ecutwfc / tpiba2
      !
      doublegrid = ( dual > 4.D0 )
      !
      IF ( doublegrid ) THEN
         !
         gcutms = 4.D0 * ecutwfc / tpiba2
         !
      ELSE
         !
         gcutms = gcutm
         !
      END IF
      !
      nks = nkstot
      !
    END SUBROUTINE set_dimensions
    !
END SUBROUTINE read_file
