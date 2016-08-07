!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this subroutine applies (H-e(iv)) to an array of num_nbndv wavefunctions
!NO SPIN YET

  subroutine hpsi_pw4gww( ndim,psi,ppsi,et,ik,numv)   
! ch_psi_all (n, h, ah, e, ik, m)    
    USE kinds,    ONLY : DP
    USE wvfct,    ONLY : npwx, npw, nbnd
    USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
    USE mp_world, ONLY : mpime, nproc
   

    implicit none

    INTEGER, INTENT(in) :: ndim !leading dimension of psi and psip
    INTEGER, INTENT(in) :: numv!number of bands
    INTEGER, INTENT(in) ::ik!dumm integer
    COMPLEX(kind=DP), INTENT(in) :: psi(ndim,numv)
    COMPLEX(kind=DP), INTENT(out) :: ppsi(ndim,numv)
    REAL(kind=DP) ::  et(numv)

    INTEGER :: iv

!apply h_psi
    do iv=1,numv
       call pc_operator(psi(1,iv),1,.false.)
    enddo
    call h_psi( ndim, npw, numv, psi, ppsi )
    do iv=1,numv
       ppsi(1:npw,iv)=ppsi(1:npw,iv)-et(iv)*psi(1:npw,iv)
       !ppsi(1:npw,iv)=ppsi(1:npw,iv)-100.d0*psi(1:npw,iv) 
    enddo
    do iv=1,numv
       call pc_operator(ppsi(1,iv),1,.false.)
    enddo
    
    return

  end subroutine hpsi_pw4gww
!-----------------------------------------------------------------
subroutine cg_psi_pw4gww (lda, n, m, psi, h_diag)
  !-----------------------------------------------------------------
  !
  !    This routine gives a preconditioning to the linear system solver.
  !    The preconditioning is diagonal in reciprocal space
  !
  !
  USE kinds, only : DP
  USE wvfct,    ONLY : et

  implicit none

  integer :: lda, n, m
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors

  complex(DP) :: psi (lda, m)
  ! inp/out: the vector to be preconditioned

  real(DP) :: h_diag (lda, m)
  ! input: the preconditioning vector

  integer :: k, i
  ! counter on bands
  ! counter on the elements of the vector
  !
  do k = 1, m
     do i = 1, n
        psi (i, k) = psi (i, k) * 1.d0/(1.d0+h_diag (i, k))!-et(k,1))
        !psi (i, k) = psi (i, k) * 1.d0/(h_diag (i, k)-100.d0)
     enddo
  enddo

  return
end subroutine cg_psi_pw4gww
