!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine ch_psi_all2 (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
#include "f_defs.h"

  use pwcom
  use becmod
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, only : DP
  use phcom
  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  real(DP) :: e (m)
  ! input: the eigenvalue

  complex(DP) :: h (npwx*npol, m), ah (npwx*npol, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex(DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  call start_clock ('ch_psi')
  allocate (ps  ( nbnd , m))    
  allocate (hpsi( npwx*npol , m))    
  allocate (spsi( npwx*npol , m))    
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  call h_psiq2 (npwx, n, m, h, hpsi, spsi)

  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah=(0.d0,0.d0)
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
     do ibnd = 1, m
        do ig = 1, n
           ah (ig+npwx,ibnd)=hpsi(ig+npwx,ibnd)-e(ibnd)*spsi(ig+npwx,ibnd)
        end do
     end do
  END IF
  !
  !   Here we compute the projector in the valence band
  !
  if (lgamma) then
     ikq = ik
  else
     ikq = 2 * ik
  endif
  ps (:,:) = (0.d0, 0.d0)

  IF (noncolin) THEN
     call ZGEMM ('C', 'N', nbnd_occ (ikq) , m, npwx*npol, (1.d0, 0.d0) , evq, &
       npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
  ELSE
     call ZGEMM ('C', 'N', nbnd_occ (ikq) , m, n, (1.d0, 0.d0) , evq, &
       npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
  ENDIF
  ps (:,:) = ps(:,:) * alpha_pv
#ifdef __PARA
  call reduce (2 * nbnd * m, ps)
#endif

  hpsi (:,:) = (0.d0, 0.d0)
  IF (noncolin) THEN
     call ZGEMM ('N', 'N', npwx*npol, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
          npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
  ELSE
     call ZGEMM ('N', 'N', n, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
          npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
  END IF
  spsi(:,:) = hpsi(:,:)
  !
  !    And apply S again
  !
  IF (noncolin) THEN
     call ccalbec_nc (nkb, npwx, n, npol, m, becp_nc, vkb, hpsi)
     call s_psi_nc (npwx, n, m, hpsi, spsi)
  ELSE
     call ccalbec (nkb, npwx, n, m, becp, vkb, hpsi)
     call s_psi (npwx, n, m, hpsi, spsi)
  END IF
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
     do ibnd = 1, m
        do ig = 1, n
           ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
        enddo
     enddo
  END IF

  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine ch_psi_all2
