!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE potinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes the self consistent potential in the array
  ! ... vr. There are three possible cases:
  !
  ! ... a) In this run the code is restarting from a broken run
  ! ... b) The potential (or rho) is read from file
  ! ... c) if a and b are both false, the total charge is computed
  ! ...    as a sum of atomic charges, and the corresponding potential
  ! ...    is saved in vr
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE brilz,         ONLY : alat, omega
  USE basis,         ONLY : nat, startingpot
  USE klist,         ONLY : nelec
  USE lsda_mod,      ONLY : lsda, nspin
  USE gvect,         ONLY:  ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                            nrxx, nl, g, gg
  USE gsmooth,       ONLY : doublegrid
  USE control_flags, ONLY : imix, lscf
  USE scf,           ONLY : rho, rho_core, vltot, vr, vrs
  USE ener,          ONLY : ehart, etxc, vtxc
  USE ldaU,          ONLY : niter_with_fixed_ns
  USE ldaU,          ONLY : lda_plus_u, Hubbard_lmax, ns, nsnew
  USE io_files,      ONLY : prefix, iunocc, input_drho
  USE mp,            ONLY : mp_bcast
  USE mp_global,     ONLY : intra_image_comm, me_image, root_image
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL (KIND=DP) :: charge               ! the starting charge
  INTEGER        :: ios
  INTEGER        :: ldim                 ! integer variable for I/O control
  LOGICAL        :: exst 
  !
  ! ... end of local variables
  !
  !
  IF ( me_image == root_image ) THEN
     !
     IF ( imix >= 0 .AND. lscf ) THEN
        CALL seqopn( 4, TRIM( prefix )//'.rho', 'UNFORMATTED', exst )
     ELSE
        CALL seqopn( 4, TRIM( prefix )//'.pot', 'UNFORMATTED', exst )
     END IF
     IF ( exst ) THEN
        CLOSE( UNIT = 4, STATUS = 'KEEP' )
     ELSE
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
     END IF
     !
  END IF
  !
  CALL mp_bcast( exst, root_image, intra_image_comm )
  !
  IF ( startingpot == 'file' .AND. exst ) THEN
     ! 
     ! ... First case, the potential is read from file
     ! ... NB: this case applies also for a restarting run, in which case
     ! ...     potential and rho files have been read from the restart file
     !
     IF ( imix >= 0 .AND. lscf ) THEN
        !      
        CALL io_pot( -1, TRIM( prefix )//'.rho', rho, nspin )
        !       
        WRITE( stdout, '(/5X,"The initial density is read from file ", A14)' ) &
            TRIM( prefix ) // '.rho'
        !
        ! ... here we compute the potential which correspond to the 
        ! ... initial charge
        !
        CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                       nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                       ehart, etxc, vtxc, charge, vr )
        !       
        IF ( ABS( charge - nelec ) / charge > 1.0D-4 ) &
           WRITE( stdout, '(/5X,"starting charge =",F10.5)') charge
        !    
     ELSE
        !
        CALL io_pot( -1, TRIM( prefix )//'.pot', vr, nspin )
        !
        WRITE( stdout, '(/5X,"The initial potential is read from file ", A14)' ) &
            TRIM( prefix ) // '.pot'
        !    
     END IF
     !
     ! ... The occupations ns also need to be read in order to build up 
     ! ... the poten
     !
     IF ( lda_plus_u ) THEN  
        !
        ldim = 2 * Hubbard_lmax + 1
        !
        IF ( me_image == root_image ) THEN
           !
           CALL seqopn( iunocc, TRIM( prefix )//'.occup', 'FORMATTED', exst )
           READ( UNIT = iunocc, FMT = * ) ns
           CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
           !
        ELSE
           !  
           ns(:,:,:,:) = 0.D0
           !
        END IF
        !
        CALL reduce( ( ldim * ldim * nspin * nat ), ns )  
        CALL poolreduce( ( ldim * ldim * nspin * nat ), ns )  
        !
        CALL DCOPY( ( ldim * ldim * nspin * nat ), ns, 1, nsnew, 1 )
        !
     END IF
     !
  ELSE
     !   
     ! ... Second case, the potential is built from a superposition 
     ! ... of atomic charges contained in the array rho_at and already 
     ! ... set in readin-readva
     !     
     IF ( startingpot == 'file' .AND. .NOT. exst ) &
        WRITE( stdout, '(5X,"Cannot read pot/rho file: not found")' )
     !
     WRITE( stdout, '(/5X,"Initial potential from superposition of free atoms")' )
     !
     ! ... in the lda+U case set the initial value of ns
     !
     IF ( lda_plus_u ) THEN
        !
        ldim = 2 * Hubbard_lmax + 1
        CALL init_ns  
        CALL DCOPY( ( ldim * ldim * nspin * nat ), ns, 1, nsnew, 1 )
        !
     END IF
     !
     CALL atomic_rho( rho, nspin )
     !
     IF ( input_drho /= ' ' ) THEN
        !
        IF ( lsda ) CALL errore( 'potinit', ' lsda not allowed in drho', 1 )
        CALL io_pot( -1, input_drho, vr, nspin )
        WRITE( stdout, '(/5X,"a scf correction to at. rho is read from", a14)' ) &
            input_drho
        CALL DAXPY( nrxx, 1.D0, vr, 1, rho, 1 )
        !
     END IF
     !
     ! ... here we compute the potential which correspond to the 
     ! ... initial charge
     !
     CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                    nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                    ehart, etxc, vtxc, charge, vr )
     !   
     IF ( ABS( charge - nelec ) / charge > 1.0D-4 ) &
          WRITE( stdout, '(/5X,"starting charge =",F10.5)') charge
     !
  END IF
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, vr, nrxx, nspin, doublegrid )
  !
  ! ... write on output the parameters used in the lda+U calculation
  !
  IF ( lda_plus_u ) THEN
     !
     WRITE( stdout, '(/5X,"Parameters of the lda+U calculation:")')
     WRITE( stdout, '(5X,"Number of iteration with fixed ns =",I3)') &
         niter_with_fixed_ns
     WRITE( stdout, '(5X,"Starting ns and Hubbard U :")')
     CALL write_ns
     !
  END IF
  !
  IF ( imix >= 0 ) CALL io_pot( +1,  TRIM( prefix )//'.rho', rho, nspin )
  !
  CALL io_pot( +1,  TRIM( prefix )//'.pot', vr, nspin )
  !
  RETURN
  !
END SUBROUTINE potinit
