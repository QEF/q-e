!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine h_psi (lda, n, m, psi, hpsi, spsi)  
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with m wavefunctions contained in psi
  ! input:
  !     lda   leading dimension of arrays psi, spsi, hpsi
  !     n     true dimension of psi, spsi, hpsi
  !     m     number of states psi
  !     psi   (G=0 component of psi is forced to be real)
  ! output:
  !     hpsi  H*psi
  !     spsi  S*psi
  !
  use pwcom  
  use gamma
  use rbecmod
  implicit none
  !
  integer :: lda, n, m  
  complex(kind=DP) :: psi (lda, m), hpsi (lda, m), spsi (lda, m)  
  !
  integer :: ibnd, j

  call start_clock ('h_psi')  
  call start_clock ('init')  
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  do ibnd = 1, m  
     ! set to zero the imaginary part of psi at G=0
     ! absolutely needed for numerical stability
     if (gstart==2) psi(1,ibnd) = cmplx(real(psi(1,ibnd)),0.d0)
     do j = 1, n  
        hpsi (j, ibnd) = g2kin (j) * psi (j, ibnd)  
     enddo
  enddo
  call stop_clock ('init')  
  !
  ! Here we add the Hubbard potential times psi
  !
  if (lda_plus_u) call vhpsi (lda, n, m, psi, hpsi)  
  !
  ! the local potential V_Loc psi
  !
  call vloc_psi(lda, n, m, psi, vrs(1,current_spin), hpsi)
  !
  !  Here the product with the non local potential V_NL psi
  !
  call pw_gemm ('Y', nkb, m, n, vkb, npwx, psi, npwx, becp, nkb)
  if (nkb.gt.0) call add_vuspsi (lda, n, m, psi, hpsi)  
  call s_psi (lda, n, m, psi, spsi)  
  call stop_clock ('h_psi')  
  return  
end subroutine h_psi

