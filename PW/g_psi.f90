!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine g_psi (lda, n, m, psi, e)
  !-----------------------------------------------------------------------
  !
  !    This routine computes an estimate of the inverse Hamiltonian
  !    and applies it to m wavefunctions
  !
  use parameters
  use g_psi_mod
  implicit none
  integer :: lda, n, m
  ! input: the leading dimension of psi
  ! input: the real dimension of psi
  ! input: the number of bands
  real(kind=DP) :: e (m)
  ! input: the eigenvectors
  complex(kind=DP) :: psi (lda, m)
  ! inp/out: the psi vector
  !
  !    Local variables
  !
  real(kind=DP), parameter :: eps = 1.0d-4
  ! a small number
  real(kind=DP) :: x, scala, denm
  integer :: k, i
  ! counter on psi functions
  ! counter on G vectors
  !
  call start_clock ('g_psi')
  !
  if (test_new_preconditioning) then
     scala = 1.0
     do k = 1, m
        do i = 1, n
           x = (h_diag(i) - e(k)*s_diag(i))*scala
           denm = (1+x+sqrt(1+(x-1)*(x-1)))/scala
           !         denm = 1.d0 + 16*x*x*x*x/(27.d0+18*x+12*x*x+8*x*x*x)
           psi (i, k) = psi (i, k) / denm
        enddo
     enddo
  else

     do k = 1, m
        do i = 1, n
           denm = h_diag (i) - e (k) * s_diag (i)
           !
           ! denm = g2+v(g=0) - e(k)
           !
           if (abs (denm) .lt.eps) denm = sign (eps, denm)
           !
           ! denm = sign( max( abs(denm),eps ), denm )
           !
           psi (i, k) = psi (i, k) / denm
        enddo
     enddo

  end if
  call stop_clock ('g_psi')
  return
end subroutine g_psi
