!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE ch_psi_all_vdw (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !

  USE pwcom
  USE becmod
  USE kinds, ONLY : DP
  USE phcom
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE wvfct,      ONLY : igk

  IMPLICIT NONE

  INTEGER :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  COMPLEX(DP) :: e (m)
  ! input: the eigenvalue plus imaginary freq.

  COMPLEX(DP) :: h (npwx, m), ah (npwx, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  INTEGER :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  COMPLEX(DP), ALLOCATABLE :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  CALL start_clock ('ch_psi')
  ALLOCATE (ps  ( nbnd , m))
  ALLOCATE (hpsi( npwx , m))
  ALLOCATE (spsi( npwx , m))
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  CALL h_psiq_vdw (npwx, n, m, h, hpsi, spsi, igk)
  !
  CALL start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  DO ibnd = 1, m
     DO ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd)  - e (ibnd) * spsi (ig, ibnd)
     ENDDO
  ENDDO
  !
  DEALLOCATE (spsi)
  DEALLOCATE (hpsi)
  DEALLOCATE (ps)
  CALL stop_clock ('last')
  CALL stop_clock ('ch_psi')
  RETURN
  !
  !   Here we compute the projector in the valence band
  !
  IF (lgamma) THEN
     ikq = ik
  ELSE
     ikq = 2 * ik
  ENDIF
  ps (:,:) = (0.d0, 0.d0)

  CALL zgemm ('C', 'N', nbnd_occ (ikq) , m, n, (1.d0, 0.d0) , evq, &
       npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
  ps (:,:) = ps(:,:) * alpha_pv
#ifdef __MPI
  CALL mp_sum ( ps, intra_pool_comm )
#endif

  hpsi (:,:) = (0.d0, 0.d0)
  CALL zgemm ('N', 'N', n, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
       npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
  spsi(:,:) = hpsi(:,:)
  !
  !    And apply S again
  !
  !call calbec (n, vkb, hpsi, becp, m)     ! not needed in TFvW
  !
  CALL s_psi (npwx, n, m, hpsi, spsi)

  DO ibnd = 1, m
     DO ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     ENDDO
  ENDDO
  !
  DEALLOCATE (spsi)
  DEALLOCATE (hpsi)
  DEALLOCATE (ps)
  CALL stop_clock ('last')
  CALL stop_clock ('ch_psi')
  RETURN
END SUBROUTINE ch_psi_all_vdw
