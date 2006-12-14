!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine dos_t (et, nspin, nbnd, nks, ntetra, tetra, e, dost)
  !------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  integer :: nspin, nbnd, nks, ntetra, tetra (4, ntetra)

  real(DP) :: et (nbnd, nks), e, dost (2)
  integer :: itetra (4), nk, ns, nt, ibnd, i

  real(DP) :: etetra (4), e1, e2, e3, e4
  integer :: nspin0

  if (nspin==4) then
     nspin0=1
  else 
     nspin0=nspin
  endif
  do ns = 1, nspin0
     dost (ns) = 0.d0
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
           ! these are the energies at the vertexes of the nt-th tetrahedron
           do i = 1, 4
              etetra (i) = et (ibnd, tetra (i, nt) + nk)
           enddo
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           if (e.lt.e4.and.e.ge.e3) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * (3.0d0 * (e4 - e) **2 / &
                   (e4 - e1) / (e4 - e2) / (e4 - e3) )
           elseif (e.lt.e3.and.e.ge.e2) then
              dost (ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                   * (3.0d0 * (e2 - e1) + 6.0d0 * (e-e2) - 3.0d0 * (e3 - e1 + e4 - e2) &
                   / (e3 - e2) / (e4 - e2) * (e-e2) **2)
           elseif (e.lt.e2.and.e.gt.e1) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * 3.0d0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
           endif
        enddo

     enddo

     ! add correct spin normalization : 2 for LDA, 1 for LSDA or
     ! noncollinear calculations 

     if ( nspin == 1 ) dost (ns) = dost (ns) * 2.d0

  enddo
  return
end subroutine dos_t
