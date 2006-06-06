!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine addusdens
  !----------------------------------------------------------------------
  USE realus, ONLY: tqr,  addusdens_r
  if (tqr) then
     call addusdens_r
  else
     call addusdens_g
  endif
  return
end subroutine addusdens

!----------------------------------------------------------------------
subroutine addusdens_g
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, nl, nlm, gg, g, eigts1, eigts2, &
                                   eigts3, ig1, ig2, ig3
  USE lsda_mod,             ONLY : nspin
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : noncolin
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : becsum, okvan
  USE uspp_param,           ONLY : lmaxq, tvanp, nh
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  implicit none
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ijh, is, nspin0
  ! counters

  real(DP), allocatable :: qmod (:), ylmk0 (:,:)
  ! the modulus of G
  ! the spherical harmonics

  complex(DP) :: skk
  complex(DP), allocatable ::  aux (:,:), qgm(:)
  ! work space for rho(G,nspin)
  ! Fourier transform of q

  if (.not.okvan) return

  call start_clock ('addusdens')

  nspin0=nspin
  if (noncolin.and..not.domag) nspin0=1

  allocate (aux ( ngm, nspin0))    
  allocate (qmod( ngm))    
  allocate (qgm( ngm))    
  allocate (ylmk0( ngm, lmaxq * lmaxq))    

  aux (:,:) = (0.d0, 0.d0)
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
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
                    do is = 1, nspin0
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
  deallocate (qgm)
  deallocate (qmod)
  !
  !     convert aux to real space and add to the charge density
  !
  do is = 1, nspin0
     psic(:) = (0.d0, 0.d0)
     psic( nl(:) ) = aux(:,is)
     if (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is))
     call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     rho (:, is) = rho (:, is) +  DBLE (psic (:) )
  enddo
  deallocate (aux)

  call stop_clock ('addusdens')
  return
end subroutine addusdens_g

