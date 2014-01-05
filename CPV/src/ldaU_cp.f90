!
! Copyright (C) 2011-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
MODULE ldaU_cp
!-------------------------------------------------------------------------
  USE parameters, ONLY: nsx
  USE kinds
  implicit none
  save
  real(DP) :: Hubbard_U(nsx)
  real(DP) :: e_hubbard = 0.d0
  real(DP), allocatable :: ns(:,:,:,:)
  integer :: Hubbard_l(nsx), Hubbard_lmax=0, ldmx=0, nwfcU
  logical :: lda_plus_u
  COMPLEX(DP), allocatable::  vupsi(:,:)
  !
contains
  !
  subroutine ldaU_init0 ( nsp, lda_plus_u_, Hubbard_U_ )
!-----------------------------------------------------------------------
!
      USE constants,        ONLY: autoev
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsp
      LOGICAL, INTENT(IN) :: lda_plus_u_
      REAL(DP),INTENT(IN) :: Hubbard_U_(nsp)

      lda_plus_u = lda_plus_u_
      Hubbard_U(1:nsp) = Hubbard_U_(1:nsp) / autoev
      !
  END SUBROUTINE ldaU_init0
  !
  subroutine deallocate_lda_plus_u()
     !
     IF( ALLOCATED( ns ) ) DEALLOCATE( ns )
     IF( ALLOCATED( vupsi ) ) DEALLOCATE( vupsi )
     !
     !
  end subroutine
  !
end module ldaU_cp
