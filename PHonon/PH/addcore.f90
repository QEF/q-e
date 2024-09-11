!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine addcore (uact, drhoc)
  !
  !! This routine computes the change of the core charge
  !! when the atoms moves along the given mode.
  !
  !
  USE kinds,      only : DP
  use uspp_param, only : upf
  use ions_base,  only : nat, ityp
  use cell_base,  only : tpiba
  use fft_base,   only : dfftp
  use fft_interfaces, only: invfft
  use gvect,   only : ngm, mill, eigts1, eigts2, eigts3, g
  use modes,   only : nmodes
  use qpoint,  only : eigqts, xq
  use nlcc_ph, only : drc
  use uspp,    only : nlcc_any
  implicit none

  !
  !   The dummy variables
  !
  COMPLEX(DP), INTENT(IN) :: uact(nmodes)
  !! input: the pattern of displacements
  COMPLEX(DP), INTENT(OUT) :: drhoc(dfftp%nnr)
  !! the change of the core charge
  !
  !   Local variables
  !
  integer :: nt, ig, mu, na, nnp, itmp
  complex(DP) :: fact, gu, gu0, u1, u2, u3, gtau
  !
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nlp_d(:)
  !
  nlp_d  => dfftp%nl_d
#else
  INTEGER, ALLOCATABLE :: nlp_d(:)
  !
  ALLOCATE( nlp_d(dfftp%ngm) )
  nlp_d  = dfftp%nl
#endif
  !
  CALL start_clock_gpu('addcore')
  !
  ! compute the derivative of the core charge  along the given mode
  !
  nnp = dfftp%nnr
  !
  !$acc data copyin(drc) deviceptr(nlp_d)
  !$acc data present_or_copyin(drhoc(1:nnp))
  !
  !$acc kernels
  drhoc(:) = (0.d0, 0.d0)
  !$acc end kernels
  !
  if (nlcc_any) then
     do na = 1, nat
        nt = ityp (na)
        if (upf(nt)%nlcc) then
           fact = tpiba * (0.d0, -1.d0) * eigqts(na)
           mu = 3 * (na - 1)
           u1 = uact(mu + 1)
           u2 = uact(mu + 2)
           u3 = uact(mu + 3)
           if ( abs(u1) + abs(u2) + abs(u3) > 1.0d-12) then
              gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
              !$acc parallel loop present(eigts1, eigts2, eigts3, g, mill, drhoc)
              do ig = 1, ngm
                 gtau = eigts1 (mill (1,ig), na) * eigts2 (mill (2,ig), na) &
                      * eigts3 (mill (3,ig), na)
                 gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
                 itmp = nlp_d(ig)
                 drhoc (itmp) = drhoc (itmp) + drc (ig, nt) * gu * fact * gtau
              enddo
           endif
        endif
     enddo
     !
     !   transform to real space
     !
     !$acc host_data use_device(drhoc)
     CALL invfft ('Rho', drhoc, dfftp)
     !$acc end host_data
     !
  endif ! nlcc_any
  !
  !$acc update host(drhoc)
  !$acc end data
  !$acc end data
  !
  CALL stop_clock_gpu('addcore')
  !
end subroutine addcore
