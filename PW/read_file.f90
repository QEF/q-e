!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE read_file
  !-----------------------------------------------------------------------
  !
  !     This routine allocates space for all quantities already computed
  !     in the pwscf program and reads them from the data file.
  !
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,            ONLY : natomwfc
  USE cell_base,        ONLY : tpiba2, bg
  USE force_mod,        ONLY : force
  USE klist,            ONLY : nkstot, nks, xk, wk
  USE lsda_mod,         ONLY : lsda, nspin, current_spin, isk
  USE wvfct,            ONLY : nbnd, nbndx, et, wg
  USE symme,            ONLY : irt 
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
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: nax =1000 ! an unlikely large number of atoms
  INTEGER :: i, ik, ibnd, ios, ierr
  !
  REAL(kind=DP), ALLOCATABLE :: et_g(:,:), wg_g(:,:)
  REAL(kind=DP) :: rdum(1,1)
  INTEGER :: kunittmp
  !
  ! choose the fortran unit to attach to the file
  !
  iunpun = 4
  !
  !  a value of zero cause the parameter to be read from the ".save" file
  !
  kunittmp = 0

  !  here we read the variables that dimension the system
  !  in parallel execution, only root proc reads the file
  !  and then broadcasts the values to all other procs
  !
  CALL readfile_new( 'dim', iunpun, rdum, rdum, kunittmp, 0, 0, ierr )
  IF( ierr /= 0 ) THEN
    CALL errore ('read_file', 'problem reading file '// &
      &      TRIM(tmp_dir)//TRIM(prefix)//'.save', ierr)
  END IF
  !
#ifdef __PARA
  kunit = kunittmp
#endif
  !
  !  allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  IF ( nat <= 0 .OR. nat > nax ) &
       CALL errore ('read_file', 'wrong number of atoms', 1)
  !
  ALLOCATE( et_g(nbnd,  nkstot), wg_g(nbnd,  nkstot) )

  ALLOCATE(tau (3, nat) )
  ALLOCATE(ityp (nat) )
  ALLOCATE(force (3, nat) )
  IF (tefield) ALLOCATE(forcefield (3, nat) )
  ALLOCATE (irt( 48, nat))    
  ALLOCATE (tetra(4, MAX(ntetra,1)))    
  !
  !     here we read all the variables defining the system
  !     in parallel execution, only root proc read the file
  !     and then broadcast the values to all ather procs
  !
  CALL readfile_new( 'nowave', iunpun, et_g, wg_g, kunittmp, 0, 0, ierr )
  IF( ierr /= 0 ) THEN
    CALL errore ('read_file', 'problem reading file '// &
      &      TRIM(tmp_dir)//TRIM(prefix)//'.save', ierr)
  END IF
  !
  !
#ifdef __PARA
  kunit = kunittmp
  ! parallel execution: distribute across pools k-points and
  ! related variables (not a smart implementation)
  nks = nkstot
  ! nks and nkstot are redefined by the following routine
  CALL divide_et_impera (xk, wk, isk, lsda, nkstot, nks)
#endif
  !
  !  check whether LSDA
  !
  IF (lsda) THEN
     nspin = 2
     npol=1
  ELSEIF (noncolin) THEN
     nspin=4
     npol = 2
     current_spin=1
  ELSE
     nspin = 1
     npol=1
     current_spin = 1
  ENDIF
  cell_factor = 1.d0
  lmovecell = .FALSE.
  !
  !   allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft
  CALL ggen
  !
  !    allocate the potential
  !
  CALL allocate_locpot
  CALL allocate_nlpot
  !
  !    allocate wavefunctions and related quantities (including et and wg)
  !
  nbndx = nbnd
  CALL allocate_wfc
  !
  et = et_g
  wg = wg_g
  !
  DEALLOCATE( et_g, wg_g )
  !
#ifdef __PARA
  CALL poolscatter (nbnd , nkstot, et, nks, et)
  CALL poolscatter (nbnd , nkstot, wg, nks, wg)
#endif
  !
  ! read the charge density
  !
  CALL io_pot ( - 1, TRIM(prefix)//'.rho', rho, nspin)
  !
  ! read the potential
  !
  CALL io_pot ( - 1, TRIM(prefix)//'.pot', vr, nspin)
  !
  ! re-calculate the local part of the pseudopotential vltot
  ! and the core correction charge (if any) - This is done here
  ! for compatibility with the previous version of read_file
  !
  CALL init_vloc
  CALL struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, &
       nr3, strf, eigts1, eigts2, eigts3)
  CALL setlocal
  CALL set_rhoc
  !
  RETURN
END SUBROUTINE read_file
