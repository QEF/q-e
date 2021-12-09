!
! Copyright (C) 2001-2018 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE phqscf
  !-----------------------------------------------------------------------
  !! This subroutine is the main driver of the self consistent cycle
  !! which gives as output the change of the wavefunctions and the
  !! change of the self-consistent potential due to a phonon of
  !! fixed q.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat
  USE lsda_mod,         ONLY : nspin
  USE io_global,        ONLY : stdout, ionode
  USE fft_base,         ONLY : dfftp
  USE uspp,             ONLY : okvan
  USE efield_mod,       ONLY : zstarue0, zstarue0_rec
  USE control_ph,       ONLY : zue, convt, rec_code
  USE partial,          ONLY : done_irr, comp_irr
  USE modes,            ONLY : nirr, npert
  USE uspp_param,       ONLY : nhm
  USE paw_variables,    ONLY : okpaw
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE recover_mod,      ONLY : write_rec
  USE mp_pools,         ONLY : inter_pool_comm
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE lrus,             ONLY : int3, int3_nc, int3_paw
  USE eqv,              ONLY : drhoscfs
  USE dynmat,           ONLY : dyn_hub_scf
  USE ldaU,             ONLY : lda_plus_u, Hubbard_lmax
  USE ldaU_lr,          ONLY : dnsscf
  USE ldaU_ph,          ONLY : dnsscf_all_modes
  USE units_ph,         ONLY : iundnsscf
  USE control_flags,    ONLY : iverbosity
  USE write_hub

  IMPLICIT NONE

  INTEGER :: irr, irr1, imode0, npe
  ! counter on the representations
  ! counter on the representations
  ! counter on the modes
  ! npert(irr)

  REAL(DP) :: tcpu, get_clock
  ! timing variables

  EXTERNAL get_clock
  ! the change of density due to perturbations

  CALL start_clock ('phqscf')
  !
  ! DFPT+U
  !
  IF (lda_plus_u) THEN
     ALLOCATE (dnsscf_all_modes(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,3*nat))
     ALLOCATE (dyn_hub_scf(3*nat,3*nat))
     dnsscf_all_modes = (0.d0, 0.d0)
     dyn_hub_scf      = (0.d0, 0.d0)
  ENDIF
  !
  ! For each irreducible representation we compute the change
  ! of the wavefunctions
  !
  DO irr = 1, nirr
     IF ( (comp_irr (irr)) .AND. (.NOT.done_irr (irr)) ) THEN
        npe=npert(irr)
        ALLOCATE (drhoscfs( dfftp%nnr , nspin_mag, npe))
        imode0 = 0
        DO irr1 = 1, irr - 1
           imode0 = imode0 + npert (irr1)
        ENDDO
        IF (npe == 1) THEN
           WRITE( stdout, '(//,5x,"Representation #",i4," mode #",i4)') &
                              irr, imode0 + 1
        ELSE
           WRITE( stdout, '(//,5x,"Representation #",i4," modes #",8i4)') &
                              irr, (imode0+irr1, irr1=1,npe)
        ENDIF
        !
        !    then for this irreducible representation we solve the linear system
        !
        IF (okvan) THEN
           ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, npe))
           IF (okpaw) ALLOCATE (int3_paw (nhm, nhm, nat, nspin_mag, npe))
           IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, nat, nspin, npe))
        ENDIF
        !
        ! DFPT+U: dnsscf in the phonon calculation
        ! is the scf change of atomic occupations ns 
        ! induced by atomic displacements.
        !
        IF (lda_plus_u) THEN
           ALLOCATE (dnsscf(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,npe))
           dnsscf = (0.d0, 0.d0)
        ENDIF
        !
        WRITE( stdout, '(/,5x,"Self-consistent Calculation")')
        CALL solve_linter (irr, imode0, npe, drhoscfs)
        WRITE( stdout, '(/,5x,"End of self-consistent calculation")')
        !
        !   Add the contribution of this mode to the dynamical matrix
        !
        IF (convt) THEN
           CALL drhodv (imode0, npe, drhoscfs)
           !
           !   add the contribution of the modes imode0+1 -> imode+npe
           !   to the effective charges Z(Us,E) (Us=scf,E=bare)
           !
           IF (zue) CALL add_zstar_ue (imode0, npe )
           IF (zue.AND. okvan) CALL add_zstar_ue_us(imode0, npe )
           !
           ! DFPT+U: calculate the scf part of the Hubbard dynamical matrix
           IF (lda_plus_u) THEN
              CALL dynmat_hub_scf (irr, imode0, npe)
              DEALLOCATE (dnsscf)
           ENDIF
           !
           IF (zue) THEN
              CALL mp_sum ( zstarue0_rec, intra_bgrp_comm )
              CALL mp_sum ( zstarue0_rec, inter_pool_comm )
              zstarue0(:,:)=zstarue0(:,:)+zstarue0_rec(:,:)
           ENDIF
           !
           WRITE( stdout, '(/,5x,"Convergence has been achieved ")')
           done_irr (irr) = .TRUE.
        ELSE
           WRITE( stdout, '(/,5x,"No convergence has been achieved ")')
           CALL stop_smoothly_ph (.FALSE.)
        ENDIF
        rec_code=20
        CALL write_rec('done_drhod',irr,0.0_DP,-1000,.false.,npe,&
                        drhoscfs)
        !
        IF (okvan) THEN
           DEALLOCATE (int3)
           IF (okpaw) DEALLOCATE (int3_paw)
           IF (noncolin) DEALLOCATE(int3_nc)
        ENDIF
        tcpu = get_clock ('PHONON')
        !
        DEALLOCATE (drhoscfs)
     ENDIF

  ENDDO
  !
  ! DFPT+U
  !
  IF (lda_plus_u) THEN
     !
     ! Write dnsscf_all_modes in the pattern basis u to file,
     ! because it is needed for el-ph calculations
     !
     IF (ionode) WRITE(iundnsscf,*) dnsscf_all_modes
     !
     ! Write dnsscf_all_modes in the cartesian coordinates
     ! to the standard output
     !
     IF (iverbosity==1) CALL write_dnsscf_ph()
     !
     ! Write the SCF and total Hubbard dynamical matrix
     ! to the standard output
     !
     IF (iverbosity==1) CALL write_dynmat_hub()
     !
     DEALLOCATE (dnsscf_all_modes)
     DEALLOCATE (dyn_hub_scf)
     !
  ENDIF
  !
  CALL stop_clock ('phqscf')
  !
  RETURN
  !
END SUBROUTINE phqscf
