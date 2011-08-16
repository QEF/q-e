!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine eigenchnl(nchanl, nchanr, tchan, vec, eigen)
!
! It performs the eigenchannel decomposition, diagonalizing
! the matrix amat=T^+T, where T is a transmission matrix
!
  USE kinds, only : DP
implicit none

  integer ::  info
  integer ::  nchanl, &  ! number of channels in the left tip
              nchanr     !     ------------         right tip
  real(DP) :: eigen(nchanl)       ! eigenvalues
  complex(DP) ::   &
              tchan(nchanr, nchanl), & ! T matrix
              vec(nchanl, nchanl),   & ! eigenvectors
              x1, x2
  complex(DP), allocatable :: amat(:,:)

  allocate( amat(  nchanl, nchanl ) )

! amat=T^+T
  x1=(1.d0, 0.d0)
  x2=(0.d0, 0.d0)
  call zgemm('c', 'n', nchanl, nchanl, nchanr, x1, tchan, nchanr,  &
              tchan, nchanr, x2, amat, nchanl)

! looking for eigenvalues of amat
  info=-1
  call hev_ab(nchanl, amat, nchanl, eigen, vec, 0.d0, 0.d0, info)

  deallocate(amat)

  return
end subroutine eigenchnl

