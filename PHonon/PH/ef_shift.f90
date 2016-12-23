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
  CONTAINS

!-----------------------------------------------------------------------
subroutine ef_shift (drhoscf, ldos, ldoss, dos_ef, irr, npe, flag)
  !-----------------------------------------------------------------------
  !    This routine takes care of the effects of a shift of Ef, due to the
  !    perturbation, that can take place in a metal at q=0
  !

  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : gg, nl
  USE buffers,              ONLY : get_buffer, save_buffer
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npw, npwx, et
  USE klist,                ONLY : degauss, ngauss, ngk, ltetra
  USE ener,                 ONLY : ef
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, nspin_lsda
! modules from phcom
  USE qpoint,               ONLY : nksq
  USE control_lr,           ONLY : nbnd_occ
  USE control_ph,           ONLY : lgamma_gamma
  USE units_ph,             ONLY : lrwfc, iuwfc, lrdwf, iudwf
  USE eqv,                  ONLY : dpsi
  USE modes,                ONLY : npert
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE dfpt_tetra_mod,       ONLY : dfpt_tetra_delta

  implicit none
  !
  ! input/output variables
  !
  integer :: npe
  ! input: the number of perturbation

  complex(DP) :: drhoscf(dfftp%nnr,nspin_mag,npe), &
       ldos(dfftp%nnr,nspin_mag), ldoss(dffts%nnr,nspin_mag)
  ! inp/out:the change of the charge
  ! inp: local DOS at Ef
  ! inp: local DOS at Ef without augme
  real(DP) :: dos_ef
  ! inp: density of states at Ef
  integer :: irr
  ! inp: index of the current irr. rep.
  logical :: flag
  ! inp: if true the eigenfunctions are updated
  !
  ! local variables
  !
  !--> these quantities may be complex since perturbation may be

  complex(DP) :: delta_n, wfshift
  ! the change in electron number
  ! the shift coefficient for the wavefunction

  real(DP), external :: w0gauss
  ! the smeared delta function

  integer :: ibnd, ik, is, ipert, nrec, ikrec
  ! counter on occupied bands
  ! counter on k-point
  ! counter on spin polarizations
  ! counter on perturbations
  ! record number
  ! record position of wfc at k
  ! auxiliary for spin
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  call start_clock ('ef_shift')
  if (.not.flag) then
     WRITE( stdout, * )
     do ipert = 1, npert (irr)
        delta_n = (0.d0, 0.d0)
        do is = 1, nspin_lsda
           CALL fwfft ('Dense', drhoscf(:,is,ipert), dfftp)
           if (gg(1).lt.1.0d-8) delta_n = delta_n + omega*drhoscf(nl(1),is,ipert)
           CALL invfft ('Dense', drhoscf(:,is,ipert), dfftp)
        enddo
        call mp_sum ( delta_n, intra_bgrp_comm )
        def (ipert) = - delta_n / dos_ef
     enddo
     !
     ! symmetrizes the Fermi energy shift
     !
     if (.not.lgamma_gamma) call sym_def (def, irr)
     WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift (Ry) =",2es15.4)')&
          (ipert, def (ipert) , ipert = 1, npert (irr) )
     !
     ! corrects the density response accordingly...
     !
     do ipert = 1, npert (irr)
        call zaxpy (dfftp%nnr*nspin_mag, def(ipert), ldos, 1, drhoscf(1,1,ipert), 1)
     enddo
  else
     !
     ! does the same for perturbed wfc
     !
     do ik = 1, nksq
        npw = ngk (ik)
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        ikrec = ik
        if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ikrec)
        !
        ! reads delta_psi from iunit iudwf, k=kpoint
        !
        do ipert = 1, npert (irr)
           nrec = (ipert - 1) * nksq + ik
           if (nksq.gt.1.or.npert(irr).gt.1) &
                call get_buffer(dpsi, lrdwf, iudwf, nrec)
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
           if (nksq.gt.1.or.npert(irr).gt.1) &
                call save_buffer (dpsi, lrdwf, iudwf, nrec)
        enddo
     enddo
     do ipert = 1, npert (irr)
        do is = 1, nspin_mag
           call zaxpy (dffts%nnr, def(ipert), ldoss(1,is), 1, drhoscf(1,is,ipert), 1)
        enddo
     enddo
  endif
  call stop_clock ('ef_shift')
  return
end subroutine ef_shift

!-----------------------------------------------------------------------
subroutine ef_shift_paw (drhoscf, dbecsum, ldos, ldoss, becsum1, &
                         dos_ef, irr, npe, flag)
  !-----------------------------------------------------------------------
  !    This routine takes care of the effects of a shift of Ef, due to the
  !    perturbation, that can take place in a metal at q=0
  !    This routine updates also dbecsum
  !

  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE ions_base,            ONLY : nat
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : omega
  USE buffers,              ONLY : get_buffer, save_buffer
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : gg, nl
  USE lsda_mod,             ONLY : nspin
  USE uspp_param,           ONLY : nhm
  USE wvfct,                ONLY : npw, npwx, et
  USE klist,                ONLY : degauss, ngauss, ngk, ltetra
  USE ener,                 ONLY : ef
! modules from phcom
  USE qpoint,               ONLY : nksq
  USE control_lr,           ONLY : nbnd_occ
  USE control_ph,           ONLY : lgamma_gamma
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag
  USE units_ph,             ONLY : lrwfc, iuwfc, lrdwf, iudwf
  USE eqv,                  ONLY : dpsi
  USE modes,                ONLY : npert
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE dfpt_tetra_mod,       ONLY : dfpt_tetra_delta

  implicit none
  !
  ! input/output variables
  !
  integer :: npe
  ! input: the number of perturbation

  complex(DP) :: drhoscf(dfftp%nnr,nspin_mag,npe), &
       ldos(dfftp%nnr,nspin_mag), ldoss(dffts%nnr,nspin_mag), &
       dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin_mag, npe)
  ! inp/out:the change of the charge
  ! inp: local DOS at Ef
  ! inp: local DOS at Ef without augme
  real(DP) :: becsum1 ( (nhm * (nhm + 1))/2 , nat , nspin_mag)
  !
  real(DP) :: dos_ef
  ! inp: density of states at Ef
  integer :: irr
  ! inp: index of the current irr. rep.
  logical :: flag
  ! inp: if true the eigenfunctions are updated
  !
  ! local variables
  !
  !--> these quantities may be complex since perturbation may be

  complex(DP) :: delta_n, wfshift
  ! the change in electron number
  ! the shift coefficient for the wavefunction

  real(DP), external :: w0gauss
  ! the smeared delta function

  integer :: ibnd, ik, is, ipert, nrec, ikrec
  ! counter on occupied bands
  ! counter on k-point
  ! counter on spin polarizations
  ! counter on perturbations
  ! record number
  ! record position of wfc at k
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  call start_clock ('ef_shift')
  if (.not.flag) then
     WRITE( stdout, * )
     do ipert = 1, npert (irr)
        delta_n = (0.d0, 0.d0)
        do is = 1, nspin_lsda
           CALL fwfft ('Dense', drhoscf(:,is,ipert), dfftp)
           if (gg(1).lt.1.0d-8) delta_n = delta_n + omega*drhoscf(nl(1),is,ipert)
           CALL invfft ('Dense', drhoscf(:,is,ipert), dfftp)
        enddo
        call mp_sum ( delta_n, intra_bgrp_comm )
        def (ipert) = - delta_n / dos_ef
     enddo
     !
     ! symmetrizes the Fermi energy shift
     !
     if (.not.lgamma_gamma) call sym_def (def, irr)
     WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift (Ry) =", 2es15.4)')&
          (ipert, def (ipert) , ipert = 1, npert (irr) )
     !
     ! corrects the density response accordingly...
     !
     do ipert = 1, npert (irr)
         drhoscf(:,:,ipert)=drhoscf(:,:,ipert)+def(ipert)*ldos(:,:)
         dbecsum(:,:,:,ipert)=dbecsum(:,:,:,ipert)+def(ipert)*&
                                          CMPLX(becsum1(:,:,:)*0.5_DP,0.0_DP,kind=DP)
     enddo
  else
     !
     ! does the same for perturbed wfc
     !
     do ik = 1, nksq
        npw = ngk (ik)
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        ikrec = ik
        if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ikrec)
        !
        ! reads delta_psi from iunit iudwf, k=kpoint
        !
        do ipert = 1, npert (irr)
           nrec = (ipert - 1) * nksq + ik
           if (nksq.gt.1.or.npert(irr).gt.1) &
                call get_buffer(dpsi, lrdwf, iudwf, nrec)
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
           if (nksq.gt.1.or.npert(irr).gt.1) &
                call save_buffer(dpsi, lrdwf, iudwf, nrec)
        enddo
     enddo
     do ipert = 1, npert (irr)
        do is = 1, nspin_mag
           call zaxpy (dffts%nnr, def(ipert), ldoss(1,is), 1, drhoscf(1,is,ipert), 1)
        enddo
     enddo
  endif
  call stop_clock ('ef_shift')
  return
end subroutine ef_shift_paw

END MODULE efermi_shift
