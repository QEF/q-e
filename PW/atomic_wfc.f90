!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine atomic_wfc (ik, wfcatom)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the  superposition of atomic wavefunctions for a
  ! given k-point.
  !
#include "machine.h"

  use parameters, only: dp, ndm, nchix
  USE atom, ONLY: nchi, lchi, chi, oc, r, rab, msh
  USE constants, ONLY: tpi, fpi
  USE brilz, ONLY: omega, tpiba
  USE basis, ONLY: nat, ntyp, natomwfc, ityp, tau
  USE gvect, ONLY: ig1, ig2, ig3, eigts1, eigts2, eigts3, g
  USE klist, ONLY: xk
  USE wvfct, ONLY: npwx, npw, nbnd, igk
  USE varie, ONLY: newpseudo
  implicit none
  integer :: ik
  ! input: k-point
  complex(kind=DP) :: wfcatom (npwx, natomwfc) ! output: atomic wavefunctions
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, ir, na, m, lm, ig, iig, nw
  !
  real(kind=DP), allocatable :: q (:), ylm (:,:), chiq (:,:,:), aux (:), &
       gk (:,:), vchi (:)
  complex(kind=DP), allocatable :: sk (:)
  real(kind=DP) :: vqint, arg
  complex(kind=DP) :: kphase  , lphase

  call start_clock ('atomic_wfc')

  allocate (q(npw),chiq(npw,nchix,ntyp),gk(3,npw),aux(ndm),vchi(ndm),sk(npw))

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
     q (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, q, ylm)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = 1, npw
     q (ig) = sqrt(q(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  do nt = 1, ntyp
     do nb = 1, nchi (nt)
        if (.not.newpseudo (nt) .or.oc (nb, nt) .gt.0.d0) then
           l = lchi (nb, nt)
           !
           !     here the  first term
           !
           call sph_bes (msh (nt), r (1, nt), q (1), l, aux)
           do ir = 1, msh (nt)
              vchi (ir) = chi (ir, nb, nt) * aux (ir) * r (ir, nt)
           enddo
           call simpson (msh (nt), vchi, rab (1, nt), vqint)
           chiq (1, nb, nt) = vqint
           !
           !    here the other terms
           !
           do ig = 2, npw
              ! dirty trick to speed up calculation:
              !    do not recalculate if |k+G_i|=|k+G_{i-1}|
              if (abs (q (ig) - q (ig - 1) ) .gt.1.0d-8) then
                 call sph_bes (msh (nt), r (1, nt), q (ig), l, aux)
                 do ir = 1, msh (nt)
                    vchi (ir) = chi (ir, nb, nt) * aux (ir) * r (ir, nt)
                 enddo
                 call simpson (msh (nt), vchi, rab (1, nt), vqint)
              endif
              chiq (ig, nb, nt) = vqint
           enddo
        endif
     enddo
  enddo

  do na = 1, nat
     arg = (xk (1, ik) * tau (1, na) + xk (2, ik) * tau (2, na) &
          + xk (3, ik) * tau (3, na) ) * tpi
     kphase = DCMPLX (cos (arg), - sin (arg) )
     !
     !     sk is the structure factor
     !
     do ig = 1, npw
        iig = igk (ig)
        sk (ig) = kphase * eigts1 (ig1 (iig), na) * eigts2 (ig2 (iig), &
             na) * eigts3 (ig3 (iig), na)
     enddo
     !
     nt = ityp (na)
     do nb = 1, nchi (nt)
        if (.not.newpseudo (nt) .or.oc (nb, nt) .gt.0.d0) then
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

  ! normalize atomic wfcs (not a bad idea in general and necessary to
  ! compute correctly lda+U projections)

  do nw = 1,natomwfc
     call DSCAL(2*npw,fpi/sqrt(omega),wfcatom(1,nw),1)
  end do

  deallocate(q, chiq ,gk ,aux ,vchi ,sk ,ylm)

  call stop_clock ('atomic_wfc')
  return
end subroutine atomic_wfc

