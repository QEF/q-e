! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE biot_savart(jpol)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the induced magentic field via the Biot-Savart law
  ! ... B_ind(r) = (1/c) \int d^3r' j(r') (r-r')/|r-r'|
  ! ... which in reciprocal space reads:
  ! ... B_ind(G) = (4\pi/c) (i G \times j(G))/G^2
  ! ... the G=0 is not computed here and is given by chi_bare
  USE kinds,                ONLY : DP
  USE constants,            ONLY : fpi
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk  
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, nlm, g, gg, ecutwfc, gcutm
  USE pwcom
  USE nmr_module

  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: jpol

  !-- local variables ----------------------------------------------------
  COMPLEX(DP), allocatable :: aux(:), j_of_g(:,:)
  complex(dp) :: fact
  INTEGER :: ig, ipol, ispin

  call start_clock('biot_savart')

  ! allocate memory
  allocate(aux(nrxxs), j_of_g(1:ngm,3))  

  ! transform current to reciprocal space
  j_of_g(:,:) = 0.d0
  do ispin = 1, nspin
    do ipol = 1, 3
      aux(1:nrxxs) = j_bare(1:nrxxs,ipol,jpol,ispin)
      call cft3s(aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)
      j_of_g(1:ngm,ipol) = j_of_g(1:ngm,ipol) + aux(nl(1:ngm))
    enddo
  enddo

  ! compute induced field in reciprocal space
  do ig = gstart, ngm
    fact = (0.d0,1.d0) * (fpi/c) / (gg(ig) * tpiba)
    b_ind(ig,1,jpol) = fact * (g(2,ig)*j_of_g(ig,3) - g(3,ig)*j_of_g(ig,2))
    b_ind(ig,2,jpol) = fact * (g(3,ig)*j_of_g(ig,1) - g(1,ig)*j_of_g(ig,3))
    b_ind(ig,3,jpol) = fact * (g(1,ig)*j_of_g(ig,2) - g(2,ig)*j_of_g(ig,1))
  enddo

  ! transform induced field in real space
  do ipol = 1, 3
    aux = (0.d0,0.d0)
    aux(nl(1:ngm)) = b_ind(1:ngm,ipol,jpol)
    call cft3s(aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1)
    b_ind_r(1:nrxxs,ipol,jpol) = real(aux(1:nrxxs))
  enddo
  
  deallocate(aux, j_of_g)
  call stop_clock('biot_savart')

END SUBROUTINE biot_savart



SUBROUTINE field_to_reciprocal_space
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, nlm, g, gg, ecutwfc, gcutm
  USE pwcom
  USE nmr_module
 
  IMPLICIT NONE
  complex(dp), allocatable :: aux(:)
  integer :: ipol, jpol

  allocate(aux(nrxxs))
  b_ind(:,:,:) = 0.d0
  do ipol = 1, 3
    do jpol = 1, 3
      aux(1:nrxxs) = b_ind_r(1:nrxxs,ipol,jpol)
      call cft3s(aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)
      b_ind(1:ngm,ipol,jpol) = aux(nl(1:ngm))
    enddo
  enddo
  deallocate(aux)

END SUBROUTINE field_to_reciprocal_space
