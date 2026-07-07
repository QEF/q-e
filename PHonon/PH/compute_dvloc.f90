!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine compute_dvloc (uact, addnlcc, dvlocin)
  !----------------------------------------------------------------------
  !! This routine calculates \(dV_\text{bare}/\text{dtau}\) in real space
  !! for one perturbation with a given q. The displacements are described
  !! by a vector u. The result is stored in \(\text{dvlocin}\).
  !!
  !
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ityp
  USE cell_base,        ONLY : tpiba
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : invfft
  USE gvect,            ONLY : eigts1, eigts2, eigts3, mill, g, gg
  USE gvecs,            ONLY : ngms
  USE uspp,             ONLY : nlcc_any
  USE eqv,              ONLY : vlocq
  USE qpoint,           ONLY : xq, eigqts
  USE modes,            ONLY : nmodes
  USE control_ph,       ONLY : lmultipole
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: uact(nmodes)
  !! input: the pattern of displacements
  LOGICAL, INTENT(IN) :: addnlcc
  !! input: if .true. the nonlinear core correction contribution is added
  COMPLEX(DP), INTENT(INOUT) :: dvlocin(dffts%nnr)
  !! output: the change of the local potential
  !
  ! ... local variables
  !
  INTEGER :: na
  !! counter on atoms
  INTEGER :: mu
  !! counter on modes
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! the type of atom
  !!
  INTEGER :: nnr, nnp, itmp, itmpp
  !!
  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  COMPLEX(DP), EXTERNAL :: Vaeps_dvloc
  COMPLEX(DP) :: pot
  !
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nl_d(:)
  !
  nl_d  => dffts%nl_d
#else
  INTEGER, ALLOCATABLE :: nl_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
#endif
  !
  call start_clock ('com_dvloc')
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  nnr = dffts%nnr
  !
  !$acc data present_or_copy(dvlocin(1:nnr)) copyin(vlocq) deviceptr(nl_d)
  !$acc kernels present(dvlocin)
  dvlocin(:) = (0.d0, 0.d0)
  !$acc end kernels
  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     u1 = uact (mu + 1)
     u2 = uact (mu + 2)
     u3 = uact (mu + 3)
     if ( abs(u1) + abs(u2) + abs(u3) > 1.0d-12) then
        nt = ityp (na)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        !$acc parallel loop present(eigts1, eigts2, eigts3, mill, g, dvlocin)
        do ig = 1, ngms
           gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
                  eigts3 (mill(3,ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           itmp = nl_d (ig)
           dvlocin (itmp) = dvlocin (itmp) + vlocq (ig, nt) * gu * fact * gtau
        enddo
        !
     endif
  enddo
  !
  if (nlcc_any.and.addnlcc) CALL errore('compute_dvloc','volentieri!',1)
  !
  IF (lmultipole .AND. gg(1) < 1d-8) THEN
    ! FM: refer potential to all-electron calculation, see routine description
    pot = Vaeps_dvloc(uact, dffts%nl(1))
    !$acc kernels
    dvlocin(dffts%nl(1)) = dvlocin(dffts%nl(1)) + pot
    !$acc end kernels
  ENDIF
  !
  ! Now we compute dV_loc/dtau in real space
  !
  !$acc host_data use_device(dvlocin)
  CALL invfft ('Rho', dvlocin, dffts)
  !$acc end host_data
  !
  !$acc end data
  !
#if !defined(__CUDA)
  DEALLOCATE(nl_d)
#endif
  !
  call stop_clock ('com_dvloc')
  !
end subroutine compute_dvloc
!----------------------------------------------------------------------
subroutine compute_dvloc_cc (uact, dvlocin)
  !----------------------------------------------------------------------
  !! This routine adds to \(dV_\text{bare}/\text{dtau}\), in real space,
  !! the core correction contribution for one perturbation with a given q.
  !! Moved out from compute_dvloc
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ityp
  USE cell_base,        ONLY : tpiba
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvecs,            ONLY : ngms
  USE lsda_mod,         ONLY : lsda, current_spin
  USE noncollin_module, ONLY : nspin_mag
  USE uspp,             ONLY : nlcc_any
  USE modes,            ONLY : nmodes
  USE dv_of_drho_lr,    ONLY : dv_of_drho_xc
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: uact(nmodes)
  !! input: the pattern of displacements
  COMPLEX(DP), INTENT(INOUT) :: dvlocin(dffts%nnr)
  !! output: the change of the local potential
  !
  ! ... local variables
  !
  INTEGER :: mu
  !! counter on modes
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! the type of atom
  !!
  INTEGER :: nnr, nnp, itmp, itmpp
  !!
  complex(DP), allocatable :: aux (:,:)
  complex(DP), pointer :: auxs (:)
  COMPLEX(DP), ALLOCATABLE :: drhoc(:)
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nl_d(:), nlp_d(:)
  !
  nl_d  => dffts%nl_d
  nlp_d  => dfftp%nl_d
#else
  INTEGER, ALLOCATABLE :: nl_d(:)
  INTEGER, ALLOCATABLE :: nlp_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  ALLOCATE( nlp_d(dfftp%ngm) )
  nl_d  = dffts%nl
  nlp_d  = dfftp%nl
#endif
  !
  ! add NLCC when present
  !
  !$acc data present_or_copy(dvlocin) deviceptr(nl_d, nlp_d)
  if (nlcc_any) then
     allocate (aux( dfftp%nnr,nspin_mag))
     !$acc data create(aux)
     allocate (drhoc( dfftp%nnr))
     !$acc data create(drhoc)
     !
     CALL addcore (uact, drhoc)
     !
     !$acc kernels
     aux(:,:) = (0.0_dp, 0.0_dp)
     !$acc end kernels
     !
     CALL dv_of_drho_xc(aux, drhoc = drhoc)
     !
     !$acc end data
     deallocate (drhoc)
     !
     !$acc host_data use_device(aux)
     IF (lsda) THEN
        CALL fwfft ('Rho', aux(:,current_spin), dfftp)
     ELSE
        CALL fwfft ('Rho', aux(:,1), dfftp)
     ENDIF
     !$acc end host_data
     !
     !  This is needed also when the smooth and the thick grids coincide to
     !  cut the potential at the cut-off - PG: does not make sense to me
     !
     allocate (auxs(dffts%nnr))
     !$acc data create(auxs)
     !$acc kernels
     auxs(:) = (0.d0, 0.d0)
     !$acc end kernels
     IF (lsda) THEN
       !$acc parallel loop present(auxs,aux)
       do ig=1,ngms
          itmp = nl_d(ig)
          itmpp = nlp_d(ig)
          auxs(itmp) = aux(itmpp,current_spin)
       enddo
     ELSE
       !$acc parallel loop present(auxs,aux)
       do ig=1,ngms
          itmp = nl_d(ig)
          itmpp = nlp_d(ig)
          auxs(itmp) = aux(itmpp,1)
       enddo
     ENDIF
     !
     !$acc host_data use_device(auxs)
     CALL invfft ('Rho', auxs, dffts)
     !$acc end host_data
     !$acc kernels present(dvlocin,auxs)
     dvlocin(:) = dvlocin(:) + auxs(:)
     !$acc end kernels
     !
     !$acc end data
     deallocate (auxs)
     !$acc end data
     deallocate (aux)
  endif
  !
  !$acc end data
  !
#if !defined(__CUDA)
  DEALLOCATE(nl_d)
  DEALLOCATE(nlp_d)
#endif
  !
end subroutine compute_dvloc_cc
