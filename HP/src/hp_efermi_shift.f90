!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

MODULE hp_efermi_shift

  USE kinds, ONLY : DP

  COMPLEX(DP), SAVE, PUBLIC :: def
  
CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE hp_ef_shift (drhoscf, ldos, ldoss, dos_ef, dbecsum, becsum1)
  !-----------------------------------------------------------------------
  !
  ! This routine takes care of the effects of a shift of Ef, due to the
  ! perturbation, that can take place in a metal at q=0.
  ! The pertubation is neutral.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat
  USE io_global,            ONLY : stdout
  USE wavefunctions,        ONLY : evc
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : gg
  USE buffers,              ONLY : get_buffer, save_buffer
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npw, npwx, et
  USE klist,                ONLY : degauss, ngauss, ngk
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, nspin_lsda
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE uspp_param,           ONLY : nhm

  IMPLICIT NONE
  !
  COMPLEX(DP) :: drhoscf(dfftp%nnr,nspin_mag), &
                 ldos(dfftp%nnr,nspin_mag),    &
                 ldoss(dffts%nnr,nspin_mag)
  ! input/output: the change of the charge (with augmentation)
  ! input: local DOS at Ef (with augmentation)
  ! input: local DOS at Ef 
  REAL(DP) :: dos_ef 
  ! input: density of states at Ef
  !
  ! The PAW case
  !
  COMPLEX(DP), OPTIONAL :: dbecsum ((nhm*(nhm+1))/2, nat, nspin_mag, 1)
  ! input:  dbecsum = 2 <psi|beta> <beta|dpsi>
  ! output: dbecsum = 2 <psi|beta> <beta|dpsi> + def * becsum1
  REAL(DP),    OPTIONAL :: becsum1 ((nhm*(nhm+1))/2, nat, nspin_mag)
  ! input: becsum1 = wdelta * <psi|beta> <beta|psi>
  !        (where wdelta is a Dirac-delta-like function)
  !
  COMPLEX(DP) :: delta_n  ! the change in electron number
  INTEGER :: is           ! counter on spin polarizations
  !
  CALL start_clock ('hp_ef_shift')
  !
  delta_n = (0.d0, 0.d0)
  !
  DO is = 1, nspin_lsda
     !
     ! FFT to G-space
     CALL fwfft ('Rho', drhoscf(:,is), dfftp)
     !
     IF (gg(1).lt.1.0d-8) delta_n = delta_n + omega * drhoscf(dfftp%nl(1),is)
     !
     ! FFT to R-space
     CALL invfft ('Rho', drhoscf(:,is), dfftp)
     !
  ENDDO
  !
  CALL mp_sum ( delta_n, intra_bgrp_comm )
  !
  def = - delta_n / dos_ef
  !
  WRITE( stdout, '(6x,"Fermi energy shift (Ry) =", 2e12.4)') def
  !
  ! Corrects the density response accordingly:
  ! drhoscfh = drhoscfh + ldos * def
  ! See Eq.(75) in Rev. Mod. Phys. 73, 515 (2001).
  !
  CALL zaxpy (dfftp%nnr*nspin_mag, def, ldos, 1, drhoscf, 1)
  !
  ! In the PAW case there is also a metallic term
  !
  IF (PRESENT(dbecsum) .AND. PRESENT(becsum1)) &
     dbecsum(:,:,:,1) = dbecsum(:,:,:,1) + &
                          def * CMPLX(becsum1(:,:,:), 0.0_DP, kind=DP)
  !
  CALL stop_clock ('hp_ef_shift')
  !
  RETURN
  !
END SUBROUTINE hp_ef_shift

END MODULE hp_efermi_shift
