!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE weights()
  !----------------------------------------------------------------------------
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho,
  ! ... Fermi energies, HOMO and LUMO, "-TS" term (gaussian)
  !
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : demet, ef, ef_up, ef_dw
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec, nelup, neldw, &
                                   two_fermi_energies, ngk
  USE ktetra,               ONLY : ltetra, ntetra, tetra
  USE lsda_mod,             ONLY : nspin, current_spin, isk
  USE noncollin_module,     ONLY : bfield
  USE wvfct,                ONLY : nbnd, wg, et
  USE mp_global,            ONLY : intra_image_comm, me_image, &
                                   root_image, npool, my_pool_id, inter_pool_comm
  USE mp,                   ONLY : mp_bcast, mp_sum
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: is, ibnd, ik
    ! counter on spin polarizations
    ! counter on bands
    ! counter on k points  
  real (DP) demet_up, demet_dw
  !
  demet         = 0.D0
  !
  IF ( .NOT. lgauss .AND. .NOT. ltetra .AND. .NOT. tfixed_occ ) THEN
     !
     ! ... calculate weights for the insulator case
     !
     IF ( two_fermi_energies ) THEN
        !
        CALL iweights( nks, wk, nbnd, nelup, et, ef_up, wg, 1, isk )
        CALL iweights( nks, wk, nbnd, neldw, et, ef_dw, wg, 2, isk )
        !
     ELSE
        !
        CALL iweights( nks, wk, nbnd, nelec, et, ef,    wg, 0, isk )
        !
     END IF
     !
  ELSE IF ( ltetra ) THEN
     !
     ! ... calculate weights for the metallic case using tetrahedra
     ! ... Important notice: all eigenvalues (et) must be present on 
     ! ... the first pool (poolreduce must have been called for et)
     !
     IF ( me_image == root_image ) THEN
        !
        IF (two_fermi_energies) THEN
           !
           CALL tweights( nkstot, nspin, nbnd, nelup, &
                          ntetra, tetra, et, ef_up, wg, 1, isk )
           CALL tweights( nkstot, nspin, nbnd, neldw, &
                          ntetra, tetra, et, ef_dw, wg, 2, isk )
           !
        ELSE
           !
           CALL tweights( nkstot, nspin, nbnd, nelec, &
                          ntetra, tetra, et, ef, wg, 0, isk )
           !
        END IF
       !
     END IF
     !
     CALL poolscatter( nbnd, nkstot, wg, nks, wg )
     !
     CALL mp_bcast( ef, root_image, intra_image_comm )
     !
  ELSE IF ( lgauss ) THEN
     !
     IF ( two_fermi_energies ) THEN
        !
        CALL gweights( nks, wk, nbnd, nelup, degauss, &
                       ngauss, et, ef_up, demet_up, wg, 1, isk )
        CALL gweights( nks, wk, nbnd, neldw, degauss, &
                       ngauss, et, ef_dw, demet_dw, wg, 2, isk )
        !
        demet = demet_up + demet_dw
        !
        bfield(3) = 0.5D0*( ef_up - ef_dw )
        !
     ELSE
        !
        CALL gweights( nks, wk, nbnd, nelec, degauss, &
                       ngauss, et, ef, demet, wg, 0, isk)
     END IF
     !
     CALL mp_sum( demet, inter_pool_comm )
     !
  ELSE IF ( tfixed_occ ) THEN
     !
     IF ( npool == 1 ) THEN
        !
        wg = f_inp
        !
     ELSE
        !
        wg(:,1) = f_inp(:,my_pool_id+1)
        wg(:,2) = f_inp(:,my_pool_id+1)
        !
     END IF
     !
     ef = - 1.D+20
     !
     DO is = 1, nspin
        !
        DO ibnd = 1, nbnd
           !
           IF ( wg(ibnd,is) > 0.D0 ) ef = MAX( ef, et(ibnd,is) )
           !
        END DO
        !
     END DO
     !
  END IF
  !
  ! ... collect all weights on the first pool;
  ! ... not needed for calculation but useful for printout 
  !
  CALL poolrecover( wg, nbnd, nkstot, nks )
  !
  RETURN
  !
END SUBROUTINE weights
