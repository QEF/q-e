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
subroutine init_at_1()
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform 
  ! of the atomic wavefunctions.
  !
  USE parameters, ONLY : nchix
  USE kinds,      ONLY : dp
  USE atom,       ONLY : nchi, lchi, chi, oc, rgrid, msh
  USE constants,  ONLY : fpi
  USE cell_base,  ONLY : omega
  USE ions_base,  ONLY : ntyp => nsp
  USE us,         ONLY : tab_at, nqx, dq
  !
  implicit none
  !
  integer :: n_starting_wfc, nt, nb, iq, ir, l, startq, lastq, ndm
  !
  real(DP), allocatable :: aux (:), vchi (:)
  real(DP) :: vqint, pref, q

  call start_clock ('init_at_1')

  ndm = MAXVAL (msh(1:ntyp))
  allocate (aux(ndm),vchi(ndm))

  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  pref = fpi/sqrt(omega)
  !
  ! needed to normalize atomic wfcs (not a bad idea in general and 
  ! necessary to compute correctly lda+U projections)
  !
  call divide (nqx, startq, lastq)
  tab_at(:,:,:) = 0.d0
  do nt = 1, ntyp
     do nb = 1, nchi (nt)
        if (oc (nb, nt) >= 0.d0) then
           l = lchi (nb, nt)
           do iq = startq, lastq
              q = dq * (iq - 1)
              call sph_bes (msh(nt), rgrid(nt)%r, q, l, aux)
              do ir = 1, msh(nt)
                 vchi(ir) = chi(ir,nb,nt) * aux(ir) * rgrid(nt)%r(ir)
              enddo
              call simpson (msh(nt), vchi, rgrid(nt)%rab, vqint)
              tab_at (iq, nb, nt) = vqint * pref
           enddo
        endif
     enddo
 enddo
#ifdef __PARA
  call reduce (nqx * nchix * ntyp, tab_at)
#endif

  deallocate(aux ,vchi)

  call stop_clock ('init_at_1')
  return

end subroutine init_at_1

