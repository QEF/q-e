!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine wfcinit
  !-----------------------------------------------------------------------
  !
  ! This routine computes an estimate of the starting wavefunctions
  ! from superposition of atomic wavefunctions.
  !
#include "machine.h"
  use pwcom
  implicit none
  !
  integer :: ik, ibnd, ig, ipol, n_starting_wfc
  ! counter on k points
  !    "     "   bands
  !    "     "  plane waves
  !    "     "  polarization
  ! number of starting wavefunctions
  complex(kind=DP), allocatable  ::  wfcatom(:,:)
  ! atomic wfcs for initialization
  real(kind=DP) ::rr, arg
  real(kind=DP), external :: rndm
  ! random function generation
  !
  ! state what is going to happen
  !
  if (startingwfc == 'file') then
     write (6, '(5x,a)') 'Starting wfc from file'
     !
     ! read the wavefunction into memory (if it is not done in c_bands)
     !
     if (nks.eq.1.and.reduce_io) call davcio(evc,nwordwfc,iunwfc,1,-1)
     return
  endif

  call start_clock ('wfcinit')
  if (startingwfc == 'atomic') then
     if (natomwfc >= nbnd) then
        write (6, '(5x,a)') 'Starting wfc are atomic'
     else
        write (6, '(5x,a,i3,a)') 'Starting wfc are atomic + ',&
             nbnd-natomwfc, ' random wfc'
     endif
     n_starting_wfc = max (natomwfc, nbnd)
  else
     write (6, '(5x,a)') 'Starting wfc are random'
     n_starting_wfc = nbnd
  endif
  !
  ! Needed for LDA+U
  !
  !!! if (lda_plus_u) call orthoatwfc
  if (nks > 1) rewind (iunigk)
  !
  !    we start a loop on k points
  !
  allocate (wfcatom( npwx, n_starting_wfc))
  !
  do ik = 1, nks
     if (lsda) current_spin = isk (ik)
     if (nks > 1) read (iunigk) npw, igk
     !
     !     here we compute the kinetic energy
     !
     do ig = 1, npw
        g2kin (ig) = (xk (1, ik) + g (1, igk (ig) ) ) **2 + &
                     (xk (2, ik) + g (2, igk (ig) ) ) **2 + &
                     (xk (3, ik) + g (3, igk (ig) ) ) **2
     enddo
     !
     ! Put the correct units on the kinetic energy
     !

     g2kin(:) = g2kin(:)*tpiba2

     !!! if (lda_plus_u) call davcio (swfcatom, nwordatwfc, iunat, ik, - 1)

     if (startingwfc == 'atomic') then
        call atomic_wfc (ik, wfcatom)

        !
        !  if not enough atomic wfc are available, complete with random wfcs
        !
        do ibnd = natomwfc + 1, nbnd
           do ig = 1, npw
              rr = rndm ()
              arg = tpi * rndm ()
              wfcatom (ig, ibnd) = DCMPLX (rr * cos (arg), rr * sin (arg) )&
                   / ( (xk (1, ik) + g (1, igk (ig) ) ) **2 + &
                   (xk (2, ik) + g (2, igk (ig) ) ) **2 + &
                   (xk (3, ik) + g (3, igk (ig) ) ) **2 + 1.0d0)
           enddo
        enddo

     else
        do ibnd = 1, nbnd
           do ig = 1, npw
              rr = rndm ()
              arg = tpi * rndm ()
              wfcatom (ig, ibnd) = DCMPLX (rr * cos (arg), rr * sin (arg) ) &
                   / ( (xk (1, ik) + g (1, igk (ig) ) ) **2 + &
                       (xk (2, ik) + g (2, igk (ig) ) ) **2 + &
                       (xk (3, ik) + g (3, igk (ig) ) ) **2 + 1.0d0)
           enddo
        enddo

     endif
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     !
     !   Diagonalize the Hamiltonian on the basis of atomic wfcs
     !
     !!! if (isolve.eq.1) then
     !!!   call cinitcgg &
     !!!        (npwx, npw, n_starting_wfc, nbnd, wfcatom, evc, et (1, ik))
     !!! else
        call rotate_wfc &
             (npwx, npw, n_starting_wfc, gstart, nbnd, wfcatom, evc, et(1, ik))
     !!! endif

     do ibnd = 1, nbnd
        do ig = npw + 1, npwx
           evc (ig, ibnd) = (0.d0, 0.d0)
        enddo
     enddo
     if (nks.gt.1.or..not.reduce_io) call davcio (evc, nwordwfc, iunwfc, ik, 1)

  enddo
  deallocate (wfcatom)
  if (iprint.eq.1) then
#ifdef __PARA
     call poolrecover (et, nbnd, nkstot, nks)
#endif
     do ik = 1, nkstot
        write (6, 9010) (xk (ipol, ik), ipol = 1, 3)
        write (6, '(2x,8f9.4)') (et (ibnd, ik) * rytoev, ibnd = 1, nbnd)
     enddo
  endif
#ifdef FLUSH
  call flush (6)
#endif
  call stop_clock ('wfcinit')

  return
9010 format (/'          k =',3f7.4,'     band energies (ev):'/)
end subroutine wfcinit

