!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine h_psiq (lda, n, m, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with a m  wavefunctions  contained
  !     in psi. It first computes the bec matrix of these
  !     wavefunctions and then with the routines hus_1psi and
  !     s_psi computes for each band the required products
  !

  use pwcom
  USE wavefunctions_module,  ONLY: psic, psic_nc
  USE becmod, ONLY: becp, becp_nc
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, only : DP
  use phcom
  implicit none
  !
  !     Here the local variables
  !
  integer :: ibnd
  ! counter on bands

  integer :: lda, n, m
  ! input: the leading dimension of the array psi
  ! input: the real dimension of psi
  ! input: the number of psi to compute
  integer :: j
  ! do loop index

  complex(DP) :: psi (lda*npol, m), hpsi (lda*npol, m), spsi (lda*npol, m)
  complex(DP) :: sup, sdwn
  ! input: the functions where to apply H and S
  ! output: H times psi
  ! output: S times psi (Us PP's only)


  call start_clock ('h_psiq')
  call start_clock ('init')

  IF (noncolin) THEN
     call ccalbec_nc (nkb, npwx, n, npol, m, becp_nc, vkb, psi)
  ELSE
     call ccalbec (nkb, npwx, n, m, becp, vkb, psi)
  END IF
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  hpsi=(0.d0,0.d0)
  do ibnd = 1, m
     do j = 1, n
        hpsi (j, ibnd) = g2kin (j) * psi (j, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
     DO ibnd = 1, m
        DO j = 1, n
           hpsi (j+lda, ibnd) = g2kin (j) * psi (j+lda, ibnd)
        ENDDO
     ENDDO
  ENDIF
  call stop_clock ('init')
  !
  ! the local potential V_Loc psi. First the psi in real space
  !

  do ibnd = 1, m
     call start_clock ('firstfft')
     IF (noncolin) THEN
        psic_nc = (0.d0, 0.d0)
        do j = 1, n
           psic_nc(nls(igkq(j)),1) = psi (j, ibnd)
           psic_nc(nls(igkq(j)),2) = psi (j+lda, ibnd)
        enddo
        call cft3s (psic_nc(1,1), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        call cft3s (psic_nc(1,2), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     ELSE
        psic(:) = (0.d0, 0.d0)
        do j = 1, n
           psic (nls(igkq(j))) = psi (j, ibnd)
        enddo
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     END IF
     call stop_clock ('firstfft')
     !
     !   and then the product with the potential vrs = (vltot+vr) on the smoo
     !
     if (noncolin) then
        if (domag) then
           do j=1, nrxxs
              sup = psic_nc(j,1) * (vrs(j,1)+vrs(j,4)) + &
                    psic_nc(j,2) * (vrs(j,2)-(0.d0,1.d0)*vrs(j,3))
              sdwn = psic_nc(j,2) * (vrs(j,1)-vrs(j,4)) + &
                    psic_nc(j,1) * (vrs(j,2)+(0.d0,1.d0)*vrs(j,3))
              psic_nc(j,1)=sup
              psic_nc(j,2)=sdwn
           end do
        else
           do j=1, nrxxs
              psic_nc(j,1)=psic_nc(j,1) * vrs(j,1)
              psic_nc(j,2)=psic_nc(j,2) * vrs(j,1)
           enddo
        endif
     else
        do j = 1, nrxxs
           psic (j) = psic (j) * vrs (j, current_spin)
        enddo
     endif
     !
     !   back to reciprocal space
     !
     call start_clock ('secondfft')
     IF (noncolin) THEN
        call cft3s(psic_nc(1,1),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
        call cft3s(psic_nc(1,2),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
     !
     !   addition to the total product
     !
        do j = 1, n
           hpsi (j, ibnd) = hpsi (j, ibnd) + psic_nc (nls(igkq(j)), 1)
           hpsi (j+lda, ibnd) = hpsi (j+lda, ibnd) + psic_nc (nls(igkq(j)), 2)
        enddo
     ELSE
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
     !
     !   addition to the total product
     !
        do j = 1, n
           hpsi (j, ibnd) = hpsi (j, ibnd) + psic (nls(igkq(j)))
        enddo
     END IF
     call stop_clock ('secondfft')
  enddo
  !
  !  Here the product with the non local potential V_NL psi
  !

  IF (noncolin) THEN
     call add_vuspsi_nc (lda, n, m, psi, hpsi)
     call s_psi_nc (lda, n, m, psi, spsi)
  ELSE
     call add_vuspsi (lda, n, m, psi, hpsi)
     call s_psi (lda, n, m, psi, spsi)
  END IF

  call stop_clock ('h_psiq')
  return
end subroutine h_psiq
