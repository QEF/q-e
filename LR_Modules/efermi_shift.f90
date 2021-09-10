!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

MODULE efermi_shift
  USE kinds, ONLY : DP
  COMPLEX(DP),SAVE,PUBLIC :: def(3)
  ! the change of the Fermi energy for each pert.
  ! NB: def(3) should be def (npertx), but it is used only at Gamma
  !     where the dimension of irreps never exceeds 3
  !
  ! Define an abstract interface to use a callback
  ABSTRACT INTERFACE
     SUBROUTINE def_symmetrization(def, irr)
        USE kinds, ONLY : DP
        INTEGER :: irr
        COMPLEX(DP) :: def(3)
     END SUBROUTINE
  END INTERFACE
  !
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE ef_shift (npert, dos_ef, ldos, drhoscf, dbecsum, becsum1, irr, sym_def)
  !-----------------------------------------------------------------------
  !! This routine takes care of the effects of a shift of Ef, due to the
  !! perturbation, that can take place in a metal at q=0
  !! Optionally, update dbecsum using becsum1.
  !
  USE kinds,                ONLY : DP
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : gg
  USE buffers,              ONLY : get_buffer, save_buffer
  USE uspp_param,           ONLY : nhm
  USE noncollin_module,     ONLY : nspin_mag, nspin_lsda
  !
  IMPLICIT NONE
  !
  ! input/output variables
  !
  INTEGER, INTENT(IN) :: npert
  !! the number of perturbation
  REAL(DP), INTENT(IN) :: dos_ef
  !! density of states at Ef
  COMPLEX(DP), INTENT(IN) :: ldos(dfftp%nnr, nspin_mag)
  !! local DOS at Ef (with augmentation)
  COMPLEX(DP), INTENT(INOUT) :: drhoscf(dfftp%nnr, nspin_mag, npert)
  !! the change of the charge (with augmentation)
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum((nhm*(nhm+1))/2, nat, nspin_mag, npert)
  !! input:  dbecsum = 2 <psi|beta> <beta|dpsi>
  !! output: dbecsum = 2 <psi|beta> <beta|dpsi> + def * becsum1
  REAL(DP), INTENT(IN), OPTIONAL :: becsum1((nhm*(nhm+1))/2, nat, nspin_mag)
  !! becsum1 = wdelta * <psi|beta> <beta|psi>
  !! (where wdelta is a Dirac-delta-like function)
  INTEGER, INTENT(IN), OPTIONAL :: irr
  !! index of the current irr. rep. Used only in sym_def.
  PROCEDURE(def_symmetrization), OPTIONAL :: sym_def
  !! Symmetrization routine for the fermi energy change
  !
  ! local variables
  !
  INTEGER :: is
  !! counter on spin polarizations
  INTEGER :: ipert
  !! counter on perturbations
  COMPLEX(DP) :: delta_n
  !! the change in electron number
  !! This may be complex since perturbation may be complex
  !
  REAL(DP), external :: w0gauss
  !! the smeared delta function
  !
  call start_clock ('ef_shift')
  !
  ! This routine is used only at q=Gamma where the dimension of irrep never exceeds 3
  IF (npert > 3) CALL errore("ef_shift", "npert exceeds 3", 1)
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  WRITE( stdout, * )
  do ipert = 1, npert
     delta_n = (0.d0, 0.d0)
     do is = 1, nspin_lsda
        CALL fwfft ('Rho', drhoscf(:,is,ipert), dfftp)
        if (gg(1) < 1.0d-8) delta_n = delta_n + omega*drhoscf(dfftp%nl(1),is,ipert)
        CALL invfft ('Rho', drhoscf(:,is,ipert), dfftp)
     enddo
     call mp_sum ( delta_n, intra_bgrp_comm )
     IF ( ABS(dos_ef) > 1.d-18 ) THEN
        def (ipert) = - delta_n / dos_ef
     ELSE
        def (ipert) = 0.0_dp
     ENDIF
  enddo
  !
  ! symmetrizes the Fermi energy shift
  !
  IF (present(sym_def)) CALL sym_def(def, irr)
  WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift (Ry) =",2es15.4)')&
       (ipert, def (ipert) , ipert = 1, npert )
  !
  ! corrects the density response accordingly...
  !
  do ipert = 1, npert
     call zaxpy (dfftp%nnr*nspin_mag, def(ipert), ldos, 1, drhoscf(1,1,ipert), 1)
  enddo
  !
  ! In the PAW case there is also a metallic term
  !
  IF (PRESENT(dbecsum) .AND. PRESENT(becsum1)) THEN
     DO ipert = 1, npert
        dbecsum(:,:,:,ipert) = dbecsum(:,:,:,ipert) &
           + def(ipert) * CMPLX(becsum1(:,:,:), 0.0_DP, KIND=DP)
     ENDDO
  ENDIF
  !
  CALL stop_clock ('ef_shift')
  !
END SUBROUTINE ef_shift
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
SUBROUTINE ef_shift_wfc(npert, ldoss, drhoscf)
  !-----------------------------------------------------------------------
  !! This routine takes care of the effects of a shift of Ef, due to the
  !! perturbation, that can take place in a metal at q=0, on the wavefunction
  !
  USE kinds,                ONLY : DP
  USE mp,                   ONLY : mp_sum
  USE wavefunctions,        ONLY : evc
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE buffers,              ONLY : get_buffer, save_buffer
  USE wvfct,                ONLY : npwx, et
  USE klist,                ONLY : degauss, ngauss, ngk, ltetra
  USE ener,                 ONLY : ef
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE qpoint,               ONLY : nksq
  USE control_lr,           ONLY : nbnd_occ
  USE units_lr,             ONLY : iuwfc, lrwfc, lrdwf, iudwf
  USE eqv,                  ONLY : dpsi
  USE dfpt_tetra_mod,       ONLY : dfpt_tetra_delta
  !
  IMPLICIT NONE
  !
  ! input/output variables
  !
  INTEGER, INTENT(IN) :: npert
  !! the number of perturbation
  COMPLEX(DP), INTENT(IN) :: ldoss(dffts%nnr, nspin_mag)
  !! local DOS at Ef without augmentation
  COMPLEX(DP), INTENT(INOUT) :: drhoscf(dfftp%nnr, nspin_mag, npert)
  !! the change of the charge (with augmentation)
  !
  ! local variables
  !
  INTEGER :: npw, ibnd, ik, is, ipert, nrec, ikrec
  ! counter on occupied bands
  ! counter on k-point
  ! counter on spin polarizations
  ! counter on perturbations
  ! record number
  ! record position of wfc at k
  ! auxiliary for spin
  COMPLEX(DP) :: wfshift
  !! the shift coefficient for the wavefunction
  !! This may be complex since perturbation may be complex
  !
  REAL(DP), external :: w0gauss
  ! the smeared delta function
  !
  call start_clock ('ef_shift_wfc')
  !
  ! This routine is used only at q=Gamma where the dimension of irrep never exceeds 3
  IF (npert > 3) CALL errore("ef_shift_wfc", "npert exceeds 3", 1)
  !
  ! Update the perturbed wavefunctions according to the Fermi energy shift
  !
  do ik = 1, nksq
     npw = ngk (ik)
     !
     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
     !
     ikrec = ik
     if (nksq > 1) call get_buffer (evc, lrwfc, iuwfc, ikrec)
     !
     ! reads delta_psi from iunit iudwf, k=kpoint
     !
     do ipert = 1, npert
        nrec = (ipert - 1) * nksq + ik
        IF (nksq > 1 .OR. npert > 1) CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
        do ibnd = 1, nbnd_occ (ik)
           !
           if(ltetra) then
              wfshift = 0.5d0 * def(ipert) * dfpt_tetra_delta(ibnd,ik)
           else
              wfshift = 0.5d0 * def(ipert) * w0gauss( (ef-et(ibnd,ik))/degauss, ngauss) / degauss
           end if
           !
           IF (noncolin) THEN
              call zaxpy (npwx*npol,wfshift,evc(1,ibnd),1,dpsi(1,ibnd),1)
           ELSE
              call zaxpy (npw, wfshift, evc(1,ibnd), 1, dpsi(1,ibnd), 1)
           ENDIF
        enddo
        !
        ! writes corrected delta_psi to iunit iudwf, k=kpoint,
        !
        IF (nksq > 1 .OR. npert > 1) CALL save_buffer(dpsi, lrdwf, iudwf, nrec)
     enddo
  enddo
  !
  do ipert = 1, npert
     do is = 1, nspin_mag
        call zaxpy (dffts%nnr, def(ipert), ldoss(1,is), 1, drhoscf(1,is,ipert), 1)
     enddo
  enddo
  !
  CALL stop_clock ('ef_shift_wfc')
  !
END SUBROUTINE ef_shift_wfc
!-------------------------------------------------------------------------
!
END MODULE efermi_shift
!-------------------------------------------------------------------------
