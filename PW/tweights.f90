!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------


subroutine tweights (nks, nspin, nbndx, nbnd, nelec, ntetra, &
     tetra, et, ef, wg)
  !--------------------------------------------------------------------
  ! calculates weights with the tetrahedron method (Bloechl version)
  use parameters
  implicit none
  !
  integer :: nks, nspin, nbndx, nbnd, ntetra, tetra (4, ntetra)
  real(kind=DP) :: et (nbndx, nks), nelec

  real(kind=DP) :: wg (nbnd, nks), ef
  real(kind=DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra (4), dosef



  integer :: ik, ibnd, nt, nk, ns, i, kp1, kp2, kp3, kp4, itetra (4)
  ! Calculate the Fermi energy ef

  call efermit (et, nbndx, nbnd, nks, nelec, nspin, ntetra, tetra, &
       ef)
  do ik = 1, nks
     do ibnd = 1, nbnd
        wg (ibnd, ik) = 0.d0
     enddo

  enddo
  do ns = 1, nspin
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then
        nk = 0
     else
        nk = nks / 2
     endif
     do nt = 1, ntetra
        do ibnd = 1, nbnd
           !
           ! etetra are the energies at the vertexes of the nt-th tetrahedron
           !
           do i = 1, 4
              etetra (i) = et (ibnd, tetra (i, nt) + nk)
           enddo
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           !
           ! ...sort in ascending order: e1 < e2 < e3 < e4
           !
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           !
           ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
           !
           kp1 = tetra (itetra (1), nt) + nk
           kp2 = tetra (itetra (2), nt) + nk
           kp3 = tetra (itetra (3), nt) + nk
           kp4 = tetra (itetra (4), nt) + nk
           !
           ! calculate weights wg
           !
           if (ef.ge.e4) then
              wg (ibnd, kp1) = wg (ibnd, kp1) + 0.25d0 / ntetra
              wg (ibnd, kp2) = wg (ibnd, kp2) + 0.25d0 / ntetra
              wg (ibnd, kp3) = wg (ibnd, kp3) + 0.25d0 / ntetra
              wg (ibnd, kp4) = wg (ibnd, kp4) + 0.25d0 / ntetra
           elseif (ef.lt.e4.and.ef.ge.e3) then
              c4 = 0.25d0 / ntetra * (e4 - ef) **3 / (e4 - e1) / (e4 - e2) &
                   / (e4 - e3)
              dosef = 3.d0 / ntetra * (e4 - ef) **2 / (e4 - e1) / (e4 - e2) &
                   / (e4 - e3)
              wg (ibnd, kp1) = wg (ibnd, kp1) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + 0.25d0 / ntetra - c4 * &
                   (4.d0 - (e4 - ef) * (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) &
                   + 1.d0 / (e4 - e3) ) ) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * &
                   et (ibnd, kp4) ) / 40.d0
           elseif (ef.lt.e3.and.ef.ge.e2) then
              c1 = 0.25d0 / ntetra * (ef - e1) **2 / (e4 - e1) / (e3 - e1)
              c2 = 0.25d0 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef) &
                   / (e4 - e1) / (e3 - e2) / (e3 - e1)
              c3 = 0.25d0 / ntetra * (ef - e2) **2 * (e4 - ef) / (e4 - e2) &
                   / (e3 - e2) / (e4 - e1)
              dosef = 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * (3.d0 * &
                   (e2 - e1) + 6.d0 * (ef - e2) - 3.d0 * (e3 - e1 + e4 - e2) &
                   * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )
              wg (ibnd, kp1) = wg (ibnd, kp1) + c1 + (c1 + c2) * (e3 - ef) &
                   / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + c1 + c2 + c3 + (c2 + c3) &
                   * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + (c1 + c2) * (ef - e1) &
                   / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + (c1 + c2 + c3) * (ef - e1) &
                   / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 + &
                   e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
           elseif (ef.lt.e2.and.ef.ge.e1) then
              c4 = 0.25d0 / ntetra * (ef - e1) **3 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              dosef = 3.d0 / ntetra * (ef - e1) **2 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              wg (ibnd, kp1) = wg (ibnd, kp1) + c4 * (4.d0 - (ef - e1) &
                   * (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) ) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + c4 * (ef - e1) / (e2 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + c4 * (ef - e1) / (e3 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + c4 * (ef - e1) / (e4 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
           endif
        enddo
     enddo


  enddo
  ! add correct spin normalization : 2 for LDA, 1 for LSDA calculations
  do ik = 1, nks
     do ibnd = 1, nbnd
        wg (ibnd, ik) = wg (ibnd, ik) * 2.d0 / nspin
     enddo

  enddo
  return
end subroutine tweights
