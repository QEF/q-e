!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusldos (ldos, becsum1)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the local DOS the part which is due to
  !  the US augmentation.
  !
  !
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use fft_base,  only: dfftp
  use fft_interfaces, only: invfft
  USE gvect, ONLY : nl, eigts1, eigts2, eigts3, mill, gg, g, ngm
  USE wavefunctions_module,  ONLY: psic
  USE uspp, ONLY: okvan
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm
  USE noncollin_module, ONLY : nspin_mag
  implicit none
  complex(DP) :: ldos (dfftp%nnr, nspin_mag)
  ! local density of states

  real(DP) :: becsum1 ( (nhm * (nhm + 1) ) / 2, nat, nspin_mag)
  ! input: the becsum1 ter
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ijh, is
  ! counters

  real(DP), allocatable :: ylmk0 (:,:), qmod (:)
  ! the spherical harmonics
  ! the modulus of G

  complex(DP), allocatable :: aux (:,:), qgm (:)
  ! work space

  allocate (aux ( ngm , nspin_mag))
  allocate (ylmk0(ngm , lmaxq * lmaxq))
  allocate (qgm ( ngm))
  allocate (qmod( ngm))

  aux (:,:) = (0.d0,0.d0)
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    !
                    !  Multiply becsum and qg with the correct structure factor
                    !
                    do is = 1, nspin_mag
                       do ig = 1, ngm
                          aux (ig, is) = aux (ig, is) + &
                                         qgm (ig) * becsum1 (ijh, na, is) * &
                                       ( eigts1 (mill(1,ig), na) * &
                                         eigts2 (mill(2,ig), na) * &
                                         eigts3 (mill(3,ig), na) )
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  !     convert aux to real space and adds to the charge density
  !
  if (okvan) then
     do is = 1, nspin_mag
        psic (:) = (0.d0,0.d0)
        do ig = 1, ngm
           psic (nl (ig) ) = aux (ig, is)
        enddo
        CALL invfft ('Dense', psic, dfftp)
        call daxpy (dfftp%nnr, 1.d0, psic, 2, ldos(1,is), 2 )
     enddo
  endif
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (aux)
  return
end subroutine addusldos
