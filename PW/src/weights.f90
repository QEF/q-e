!
! Copyright (C) 2001-2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE weights()
  !----------------------------------------------------------------------------
  !! Calculates weights of Kohn-Sham orbitals used in calculation of rho,
  !! Fermi energies, HOMO and LUMO, "-TS" term (gaussian).
  !
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : demet, ef, ef_up, ef_dw
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE klist,                ONLY : ltetra, lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec, nelup, neldw, &
                                   two_fermi_energies
  USE ktetra,               ONLY : ntetra, tetra, tetra_type, tetra_weights, &
                                   opt_tetra_weights
  USE lsda_mod,             ONLY : nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, wg, et
  USE mp_images,            ONLY : intra_image_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE io_global,            ONLY : ionode, ionode_id
  USE gcscf_module,         ONLY : lgcscf, gcscf_mu, gcscf_beta
  USE wvfct_gpum,           ONLY : using_et, using_wg, using_wg_d
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ibnd, ik ! counters: bands, k-points
  REAL(DP) :: demet_up, demet_dw
  !
  CALL using_et(0)
  CALL using_wg(2)
  !
  demet = 0.D0
  !
  IF ( tfixed_occ .OR. ltetra ) THEN
     !
     ! ... For these two cases, the weights are computed on one processor,
     ! ... broadcast to the other. All eigenvalues (et) must be present on
     ! ... the first pool: poolreduce must have been called for et
     !
     IF ( ionode ) THEN
        !
        IF ( tfixed_occ ) THEN 
           ! 
           ! ... occupancies are fixed to the values read from input
           !
           DO ik = 1, nkstot
              wg(:,ik) = f_inp(:,isk(ik)) * wk(ik)
              IF ( nspin == 1 ) wg(:,ik) = wg(:,ik)/2.0_dp
           ENDDO
           !
           ef = -1.0d10
           DO ik = 1, nkstot
              DO ibnd = 1, nbnd
                 IF ( wg(ibnd,ik) > 0.D0 ) ef = MAX( ef, et(ibnd,ik) )
              ENDDO
           ENDDO
           !
        ELSE
           !
           ! ... calculate weights for the metallic case using tetrahedra
           !
           IF (tetra_type == 0) THEN
              !
              ! Bloechl's tetrahedra
              !
              IF (two_fermi_energies) THEN
                 !
                 CALL tetra_weights( nkstot, nspin, nbnd, nelup, et, ef_up, wg, 1, isk )
                 CALL tetra_weights( nkstot, nspin, nbnd, neldw, et, ef_dw, wg, 2, isk )
                 !
              ELSE
                 !
                 CALL tetra_weights( nkstot, nspin, nbnd, nelec, et, ef, wg, 0, isk )
                 !
              ENDIF
              !
           ELSE
              !
              ! Linear or Optimized tetrahedra
              !
              IF (two_fermi_energies) THEN
                 !
                 CALL opt_tetra_weights( nkstot, nspin, nbnd, nelup, et, ef_up, wg, 1, isk )
                 CALL opt_tetra_weights( nkstot, nspin, nbnd, neldw, et, ef_dw, wg, 2, isk )
                 !
              ELSE
                 !             
                 CALL opt_tetra_weights ( nkstot, nspin, nbnd, nelec, et, ef, wg, 0, isk )
                 !
              ENDIF
              !
           ENDIF ! tetra_type
           !
        ENDIF
        !
     ENDIF
     !
     CALL poolscatter( nbnd, nkstot, wg, nks, wg )
     CALL mp_bcast( ef, ionode_id, intra_image_comm )
     !
  ELSE
     !
      IF ( lgauss ) THEN
        !
        ! ... calculate weights for the metallic case using smearing
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
        ELSE
           !
           IF ( lgcscf ) THEN
              !
              ef = gcscf_mu
              !
              CALL gweights_mix( nks, wk, nbnd, nelec, degauss, &
                                 ngauss, et, ef, demet, wg, 0, isk, gcscf_beta )
              !
           ELSE
              !
              CALL gweights( nks, wk, nbnd, nelec, degauss, &
                             ngauss, et, ef, demet, wg, 0, isk )
              !
           END IF
           !
        ENDIF
        !
        CALL mp_sum( demet, inter_pool_comm )
        !
     ELSE
        !
        ! ... calculate weights for the insulator case
        !
        IF ( two_fermi_energies ) THEN
           !
           CALL iweights( nks, wk, nbnd, nelup, et, ef_up, wg, 1, isk )
           CALL iweights( nks, wk, nbnd, neldw, et, ef_dw, wg, 2, isk )
           !
           ! the following line to prevent NaN in Ef
           !
           ef = ( ef_up + ef_dw ) / 2.0_dp
           !
        ELSE
           !
           CALL iweights( nks, wk, nbnd, nelec, et, ef, wg, 0, isk )
           !
        ENDIF
        !
     ENDIF
     !
     ! ... collect all weights on the first pool;
     ! ... not needed for calculation but useful for printout 
     !
     CALL poolrecover( wg, nbnd, nkstot, nks )
     !
  ENDIF
#if defined(__CUDA)
  ! Sync here. Shouldn't be done and will be removed ASAP.
  CALL using_wg_d(0)
#endif
  !
  RETURN
  !
END SUBROUTINE weights
!
!
!----------------------------------------------------------------------------
SUBROUTINE weights_only()
  !----------------------------------------------------------------------------
  !! Calculates only weights of Kohn-Sham orbitals, with Fermi energy 
  !! given in input.
  !
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : demet, ef, ef_up, ef_dw
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE klist,                ONLY : ltetra, lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec, nelup, neldw,  &
                                   two_fermi_energies
  USE ktetra,               ONLY : ntetra, tetra, tetra_type, &
                                   tetra_weights_only, opt_tetra_weights_only
  USE lsda_mod,             ONLY : nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, wg, et
  USE mp_images,            ONLY : intra_image_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE io_global,            ONLY : ionode, ionode_id
  !
  USE wvfct_gpum,           ONLY : using_et, using_wg, using_wg_d
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ibnd, ik ! counters: bands, k-points
  REAL(DP) :: demet_up, demet_dw
  !
  CALL using_et(0)
  CALL using_wg(2)
  !
  demet = 0.D0
  !
  IF ( tfixed_occ .OR. ltetra ) THEN
     !
     ! ... For these two cases, the weights are computed on one processor,
     ! ... broadcast to the other. All eigenvalues (et) must be present on
     ! ... the first pool: poolreduce must have been called for et
     !
     IF ( ionode ) THEN
        !
        IF ( tfixed_occ ) THEN 
           ! 
           ! ... occupancies are fixed to the values read from input
           !
           DO ik = 1, nkstot
              wg(:,ik) = f_inp(:,isk(ik)) * wk(ik)
              IF ( nspin == 1 ) wg(:,ik) = wg(:,ik)/2.0_dp
           ENDDO
           !
        ELSE
           !
           ! ... calculate weights for the metallic case using tetrahedra
           !
           IF (tetra_type == 0) then
              !
              ! Bloechl's tetrahedra
              !
              IF (two_fermi_energies) THEN
                 !
                 CALL tetra_weights_only( nkstot, nspin, 1, isk, nbnd, nelup, et, ef_up, wg )
                 CALL tetra_weights_only( nkstot, nspin, 2, isk, nbnd, neldw, et, ef_dw, wg )
                 !
              ELSE
                 !
                 CALL tetra_weights_only( nkstot, nspin, 0, isk, nbnd, nelec, et, ef, wg )
                 !
              ENDIF
              !
           ELSE ! tetra_type == 1 .or. 2
              !
              ! Linear or Optimized tetrahedra
              !
              IF (two_fermi_energies) THEN
                 !
                 CALL opt_tetra_weights_only( nkstot, nspin, nbnd, et, ef_up, wg, 1, isk )
                 CALL opt_tetra_weights_only( nkstot, nspin, nbnd, et, ef_dw, wg, 2, isk )
                 !
              ELSE
                 !
                 CALL opt_tetra_weights_only( nkstot, nspin, nbnd, et, ef, wg, 0, isk )
                 !
              ENDIF
              !
           ENDIF ! tetra_type
           !
        ENDIF
        !
     ENDIF
     !
     CALL poolscatter( nbnd, nkstot, wg, nks, wg )
     !
  ELSE
     !
     IF ( lgauss ) THEN
        !
        ! ... calculate weights for the metallic case using smearing
        !
        IF ( two_fermi_energies ) THEN
           !
           CALL gweights_only( nks, wk, 1, isk, nbnd, nelup, degauss, &
                               ngauss, et, ef_up, demet_up, wg )
           CALL gweights_only( nks, wk, 2, isk, nbnd, neldw, degauss, &
                               ngauss, et, ef_dw, demet_dw, wg )
           !
           demet = demet_up + demet_dw
           !
        ELSE
           !
           CALL gweights_only( nks, wk, 0, isk, nbnd, nelec, degauss, &
                               ngauss, et, ef, demet, wg )
        ENDIF
        !
        CALL mp_sum( demet, inter_pool_comm )
        !
     ELSE
        !
        ! ... calculate weights for the insulator case
        !
        IF ( two_fermi_energies ) THEN
           !
           CALL iweights_only( nks, wk, 1, isk, nbnd, nelup, wg )
           CALL iweights_only( nks, wk, 2, isk, nbnd, neldw, wg )
           !
        ELSE
           !
           CALL iweights_only( nks, wk, 0, isk, nbnd, nelec, wg )
           !
        ENDIF
        !
     ENDIF
     !
     ! ... collect all weights on the first pool;
     ! ... not needed for calculation but useful for printout 
     !
     CALL poolrecover( wg, nbnd, nkstot, nks )
     !
  ENDIF
#if defined(__CUDA)
  ! Sync here. Shouldn't be done and will be removed ASAP.
  CALL using_wg_d(0)
#endif
  !
  RETURN
  !
END SUBROUTINE weights_only
