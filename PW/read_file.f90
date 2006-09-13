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
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, nsp, ntyp => nsp, ityp, tau, if_pos
  USE basis,                ONLY : natomwfc
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg, npwx
  USE symme,                ONLY : irt, nsym, ftau, s
  USE ktetra,               ONLY : tetra, ntetra 
  USE extfield,             ONLY : forcefield, tefield
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE gvect,                ONLY : gg, ecutwfc, ngm, g, nr1, nr2, nr3, nrxx,&
                                   nrx1, nrx2, nrx3, eigts1, eigts2, eigts3, &
                                   nl, gstart
  USE gsmooth,              ONLY : ngms, nls, nrx1s, nr1s, nr2s, nr3s
  USE spin_orb,             ONLY : so, lspinorb
  USE scf,                  ONLY : rho, rhog, rho_core, rhog_core, vr
  USE wavefunctions_module, ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc
  USE uspp_param,           ONLY : nbeta, jjj, tvanp
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_global,            ONLY : kunit
  USE pw_restart,           ONLY : pw_readfile
  USE xml_io_base,          ONLY : pp_check_file
  USE uspp,                 ONLY : okvan
  !
  IMPLICIT NONE
  !
  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, charge
  REAL(DP) :: sr(3,3,48)
  LOGICAL  :: exst
  !
  !
  ! ... first we check if the file can be used for post-processing
  !
  IF ( .NOT. pp_check_file() ) &
     CALL infomsg( 'read_file', 'file ' // TRIM( tmp_dir ) // TRIM( prefix ) &
               & // '.save not guaranteed to be safe for post-processing', -1 )
  !
  ! ... here we read the variables that dimension the system
  ! ... in parallel execution, only root proc reads the file
  ! ... and then broadcasts the values to all other procs
  !
  ! ... a reset of the internal flags is necessary because some codes call
  ! ... read_file() more than once
  !
  CALL pw_readfile( 'reset', ierr )
  CALL pw_readfile( 'dim',   ierr )
  !
  CALL errore( 'read_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  ! ... allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  IF ( nat <= 0 ) &
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
  ! ... allocate memory for eigenvalues and weights (read from file)
  !
  nbndx = nbnd
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  CALL pw_readfile( 'nowave', ierr )
  !
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
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
  okvan = ANY ( tvanp(1:nsp) )
  !
  ! ... check for spin-orbit pseudopotentials
  !
  DO nt = 1, nsp
     !
     so(nt) = ( nbeta(nt) > 0 )
     !
     DO nb = 1, nbeta(nt)
        !
        so(nt) = so(nt) .AND. ( ABS( jjj(nb,nt) ) > 1.D-7 )
        !
     END DO
     !
  END DO
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( ntyp )
  !
  ! ... allocate the potential and wavefunctions
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  CALL allocate_wfc()
  !
  ! ... read the charge density
  !
  CALL pw_readfile( 'rho', ierr )
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
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho(:,is)
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     !
     rhog(:,is) = psic(nl(:))
     !
  END DO
  !
  ! ... recalculate the potential
  !
  CALL v_of_rho( rho, rhog, rho_core, rhog_core, &
                 ehart, etxc, vtxc, etotefield, charge, vr )
  !
  ! ... reads the wavefunctions and writes them in 'distributed' form 
  ! ... to unit iunwfc (for compatibility)
  !
  nwordwfc = 2*nbnd*npwx*npol
  !
  CALL diropn( iunwfc, 'wfc', nwordwfc, exst )
  !
  CALL pw_readfile( 'wave', ierr )
  !
  CLOSE( UNIT = iunwfc, STATUS = 'KEEP' )
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
