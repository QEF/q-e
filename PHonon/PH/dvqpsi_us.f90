!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us (ik, uact, addnlcc, becp1, alphap)
  !----------------------------------------------------------------------
  !! This routine calculates \(dV_\text{bare}/d\tau \cdot \psi\) for one
  !! perturbation with a given q. The displacements are described by a
  !! vector u.
  !! The result is stored in \(\text{dvpsi}\). The routine is called for
  !! each k-point and for each pattern u. It computes simultaneously all
  !! the bands. It implements Eq. (B29) of PRB 64, 235118 (2001). The
  !! contribution of the local pseudopotential is calculated here, that
  !! of the nonlocal pseudopotential in \(\texttt{dvqpsi_us_only}\).
  !
  !
  USE kinds,            ONLY : DP
  USE funct,            ONLY : dft_is_nonlocc
  USE xc_lib,           ONLY : xclib_dft_is
  USE ions_base,        ONLY : nat
  USE fft_base,         ONLY : dffts
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE lsda_mod,         ONLY : lsda, isk, current_spin
  USE noncollin_module, ONLY : npol
  USE ldaU,             ONLY : lda_plus_u
  USE wvfct,            ONLY : nbnd, npwx
  USE wavefunctions,    ONLY : evc
  USE eqv,              ONLY : dvpsi
  USE qpoint,           ONLY : ikqs, ikks
  USE klist,            ONLY : ngk, igk_k
  USE qpoint,           ONLY : nksq
  USE becmod,           ONLY : bec_type
  USE gvect,            ONLY : gg
  USE control_lr,       ONLY : lmultipole
  !
  IMPLICIT NONE
  !
  !   The dummy variables
  !
  INTEGER, INTENT(in) :: ik
  !! input: the k point
  COMPLEX(DP) :: uact(3*nat)
  !! input: the pattern of displacements
  LOGICAL :: addnlcc
  TYPE(bec_type) :: becp1(nksq), alphap(3,nksq)
  !
  ! ... local variables
  !
  INTEGER :: npw
  !! Number of pw
  INTEGER :: ikk
  !! the point k
  INTEGER :: npwq
  !! Number of q
  INTEGER :: ikq
  !! k-q index
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: ibnd
  !! counter on bands
  INTEGER :: ir
  !! counter on real mesh
  INTEGER :: ip, nnr, itmp
  !!
  complex(DP) , allocatable :: dvlocin (:), aux2 (:)
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
  call start_clock_gpu ('dvqpsi_us')
  allocate (dvlocin(dffts%nnr))
  allocate (aux2(dffts%nnr))
  !
  nnr = dffts%nnr
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  IF (lsda) current_spin = isk(ikk)
  !
  !    We start by computing the contribution the local potential.
  !
  !$acc data create(dvlocin(1:nnr),aux2(1:nnr)) copyout(dvpsi) present( igk_k ) deviceptr(nl_d)
  !
  !$acc kernels present(dvlocin, dvpsi)
  dvlocin(:) = (0.d0, 0.d0)
  dvpsi(:,:) = (0.d0, 0.d0)
  !$acc end kernels
  !
  ! Compute dV_loc/dtau in real space
  !
  IF (.NOT. lmultipole) THEN
    CALL compute_dvloc (uact, addnlcc, dvlocin)
  ELSE
    ! Bring potential in reciprocal space
    !$acc host_data use_device(dvlocin)
    CALL fwfft ('Rho', dvlocin, dffts)
    !$acc end host_data
    IF (gg(1) < 1d-8) THEN
      !$acc kernels
      dvlocin(dffts%nl(1)) = 1d0
     !$acc end kernels
    ENDIF
    ! Bring potential in real space
   !$acc host_data use_device(dvlocin)
    CALL invfft ('Rho', dvlocin, dffts)
   !$acc end host_data
  ENDIF
  ! Now we compute dV_loc/dtau * psi in real space
  !
  !$acc update device(evc)
  !
  do ibnd = 1, nbnd
     do ip=1,npol
        !$acc kernels present(aux2)
        aux2(:) = (0.d0, 0.d0)
        !$acc end kernels
        if (ip==1) then
           !$acc parallel loop present(aux2, igk_k)
           do ig = 1, npw
              itmp = nl_d (igk_k (ig,ikk) )
              aux2 ( itmp ) = evc (ig, ibnd)
           enddo
        else
           !$acc parallel loop present(aux2, igk_k)
           do ig = 1, npw
              itmp = nl_d (igk_k (ig,ikk) )
              aux2 ( itmp ) = evc (ig+npwx, ibnd)
           enddo
        end if
        !
        !  This wavefunction is computed in real space
        !
        !$acc host_data use_device(aux2)
        CALL invfft ('Wave', aux2, dffts)
        !$acc end host_data
        !$acc parallel loop present(aux2, dvlocin)
        do ir = 1, nnr
           aux2 (ir) = aux2 (ir) * dvlocin (ir)
        enddo
        !
        !  and finally dV_loc/dtau * psi is transformed in reciprocal space
        !
        !$acc host_data use_device(aux2)
        CALL fwfft ('Wave', aux2, dffts)
        !$acc end host_data
        if (ip==1) then
           !$acc parallel loop present( aux2, igk_k, dvpsi )
           do ig = 1, npwq
              itmp = nl_d (igk_k (ig,ikq) )
              dvpsi (ig, ibnd) = aux2 ( itmp )
           enddo
        else
           !$acc parallel loop present( aux2, igk_k, dvpsi )
           do ig = 1, npwq
              itmp = nl_d (igk_k (ig,ikq) )
              dvpsi (ig+npwx, ibnd) = aux2 ( itmp )
           enddo
        end if
     enddo
  enddo
  !$acc end data
  !
  deallocate (aux2)
  deallocate (dvlocin)
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients.
  !
  IF (.NOT. lmultipole) call dvqpsi_us_only (ik, uact, becp1, alphap)
  !
  ! DFPT+U: add the bare variation of the Hubbard potential
  !
  IF (lda_plus_u) CALL dvqhub_barepsi_us(ik, uact)
  !
  call stop_clock_gpu ('dvqpsi_us')
  !
#if !defined(__CUDA)
  DEALLOCATE(nl_d)
#endif
  return
end subroutine dvqpsi_us
