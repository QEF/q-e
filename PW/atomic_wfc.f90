!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine atomic_wfc (ik, wfcatom)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the  superposition of atomic wavefunctions for a
  ! given k-point.
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nchix
  USE atom,       ONLY : nchi, lchi, chi, oc
  USE constants,  ONLY : tpi, fpi
  USE cell_base,  ONLY : omega, tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : ig1, ig2, ig3, eigts1, eigts2, eigts3, g
  USE klist,      ONLY : xk
  USE wvfct,      ONLY : npwx, npw, nbnd, igk
  USE us,         ONLY : tab_at, dq
  !
  implicit none
  !
  integer :: ik
  ! input: k-point
  complex(DP) :: wfcatom (npwx, natomwfc) ! output: atomic wavefunctions
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3
  !
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
  complex(DP), allocatable :: sk (:)
  real(DP) :: arg, px, ux, vx, wx
  complex(DP) :: kphase  , lphase

  call start_clock ('atomic_wfc')

  allocate ( qg(npw), chiq(npw,nchix,ntyp), gk(3,npw), sk(npw))

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  do nt = 1, ntyp
     do nb = 1, nchi (nt)
        lmax_wfc = max (lmax_wfc, lchi (nb, nt) )
     enddo
  enddo
  !
  allocate(ylm (npw,(lmax_wfc+1)**2) )
  !
  do ig = 1, npw
     gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
     gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
     gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = 1, npw
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  do nt = 1, ntyp
     do nb = 1, nchi (nt)
        if ( oc (nb, nt) >= 0.d0) then
           do ig = 1, npw
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           enddo
        endif
     enddo
  enddo

  do na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX (cos (arg), - sin (arg) )
     !
     !     sk is the structure factor
     !
     do ig = 1, npw
        iig = igk (ig)
        sk (ig) = kphase * eigts1 (ig1 (iig), na) * eigts2 (ig2 (iig), na) * &
                           eigts3 (ig3 (iig), na)
     enddo
     !
     nt = ityp (na)
     do nb = 1, nchi (nt)
        if (oc (nb, nt) >= 0.d0) then
           l = lchi (nb, nt)
           lphase = (0.d0,1.d0)**l
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
            do m = 1, 2 * l + 1
              lm = l**2 + m
              n_starting_wfc = n_starting_wfc + 1
              if (n_starting_wfc.gt.natomwfc) &
                   call errore ('atomic_wfc', 'too many wfcs', 1)
              do ig = 1, npw
                 wfcatom (ig, n_starting_wfc) = lphase * &
                      sk (ig) * ylm (ig, lm) * chiq (ig, nb, nt)
              enddo
           enddo
        endif
     enddo
  enddo

  if (n_starting_wfc.ne.natomwfc) call errore ('atomic_wfc', &
       'something wrong', 1)

  deallocate(qg, chiq ,gk, sk ,ylm)

  call stop_clock ('atomic_wfc')
  return
end subroutine atomic_wfc

