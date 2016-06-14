!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine cch_psi_all (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !

  USE kinds, only : DP
  USE becmod, ONLY : becp, calbec
  USE uspp, ONLY: nkb, vkb
  USE wvfct, ONLY : npwx, nbnd, current_k
  USE klist, ONLY : igk_k
  USE noncollin_module, ONLY : noncolin, npol

  USE eqv,  ONLY : evq
  USE qpoint, ONLY : ikqs

  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  USE control_lr, ONLY : alpha_pv, nbnd_occ

  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  complex(kind=DP) :: e (m)
  ! input: the eigenvalue + iu

  complex(kind=DP) :: h (npwx*npol, m), ah (npwx*npol, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex(kind=DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  call start_clock ('ch_psi')
  allocate (ps  ( nbnd , m))    
  allocate (hpsi( npwx * npol, m))    
  allocate (spsi( npwx * npol, m))    
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  ikq = ikqs(ik)
  current_k = ikq  ! used in h_psi
  CALL h_psi (npwx, n, m, h, hpsi)
  CALL s_psi (npwx, n, m, h, spsi)

  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah=(0.0_DP, 0.0_DP)
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     enddo
     IF (noncolin) THEN
        do ig = 1, n
           ah (ig+npwx, ibnd) = hpsi (ig+npwx, ibnd) - e (ibnd) * &
                                spsi (ig+npwx, ibnd)
        enddo
     END IF
  enddo
  !
  !   Here we compute the projector in the valence band
  !
  ps (:,:) = (0.d0, 0.d0)

  IF (noncolin) THEN
     call zgemm ('C', 'N', nbnd_occ (ikq) , m, npwx*npol, (1.d0, 0.d0) , evq, &
          npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
  ELSE
     call zgemm ('C', 'N', nbnd_occ (ikq) , m, n, (1.d0, 0.d0) , evq, &
          npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
  ENDIF
  ps (:,:) = ps(:,:) * alpha_pv
  call mp_sum (ps, intra_bgrp_comm)

  hpsi (:,:) = (0.d0, 0.d0)
  IF (noncolin) THEN
     call zgemm ('N', 'N', npwx*npol, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
          npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
  ELSE
     call zgemm ('N', 'N', n, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
          npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
  END IF
  spsi(:,:) = hpsi(:,:)
  !
  !    And apply S again
  !
  call calbec (n, vkb, hpsi, becp, m)

  call s_psi (npwx, n, m, hpsi, spsi)
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     enddo
     IF (noncolin) THEN
        do ig = 1, n
           ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
        enddo
     END IF
  enddo

  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine cch_psi_all
