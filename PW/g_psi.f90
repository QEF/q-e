!
! Copyright (C) 2001 PWSCF group
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
!    first the input variables
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
!   One parameter
!
real(kind=DP) :: eps  
                          ! a small number
parameter (eps = 1.0d-4)  
!
!    The local variables
!
real(kind=DP) :: denm     ! the denominator
                    
integer :: k, i  
                          ! counter on psi functions
                          ! counter on G vectors
call start_clock ('g_psi')  
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
call stop_clock ('g_psi')  
return  
end subroutine g_psi

