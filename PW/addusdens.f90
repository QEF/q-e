!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusdens
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
#include "machine.h"
  use pwcom
  USE wavefunctions,    ONLY : psic
  implicit none
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ijh, is
  ! counters

  real(kind=DP), allocatable :: qmod (:), ylmk0 (:,:)
  ! the modulus of G
  ! the spherical harmonics

  complex(kind=DP) :: skk
  complex(kind=DP), allocatable ::  aux (:,:)
  ! work space for FFT
  ! work space for rho(G,nspin)

  if (.not.okvan) return

  call start_clock ('addusdens')

  allocate (aux ( ngm, nspin))    
  allocate (qmod( ngm))    
  allocate (ylmk0( ngm, lqx * lqx))    

  aux (:,:) = (0.d0, 0.d0)
  call ylmr2 (lqx * lqx, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
#ifdef DEBUG_ADDUSDENS
  call start_clock ('addus:qvan2')
#endif
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
#ifdef DEBUG_ADDUSDENS
  call stop_clock ('addus:qvan2')
#endif
              ijh = ijh + 1
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    !
                    !  Multiply becsum and qg with the correct structure factor
                    !
#ifdef DEBUG_ADDUSDENS
  call start_clock ('addus:aux')
#endif
                    do is = 1, nspin
                       do ig = 1, ngm
                          skk = eigts1 (ig1 (ig), na) * &
                                eigts2 (ig2 (ig), na) * &
                                eigts3 (ig3 (ig), na)
                          aux(ig,is)=aux(ig,is) + qgm(ig)*skk*becsum(ijh,na,is)
                       enddo
                    enddo
#ifdef DEBUG_ADDUSDENS
  call stop_clock ('addus:aux')
#endif
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  deallocate (ylmk0)
  deallocate (qmod)
  !
  !     convert aux to real space and add to the charge density
  !
  do is = 1, nspin
     psic(:) = (0.d0, 0.d0)
     psic( nl(:) ) = aux(:,is)
     if (gamma_only) psic( nlm(:) ) = conjg(aux(:,is))
     call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     rho (:, is) = rho (:, is) + DREAL (psic (:) )
  enddo
  deallocate (aux)

  call stop_clock ('addusdens')
  return
end subroutine addusdens

