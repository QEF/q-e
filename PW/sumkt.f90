!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
FUNCTION sumkt (et, nbnd, nks, nspin, ntetra, tetra, e, is, isk)
  !--------------------------------------------------------------------
  !
  ! ... Sum over all states with tetrahedron method
  ! ... At Fermi energy e=E_F, sumkt(e) == number of electrons
  ! ... Generalization to noncollinear case courtesy of Yurii Timrov
  !
  USE kinds
  implicit none
  ! output variable
  real(DP) :: sumkt
  ! input variable
  integer, intent(in) :: nbnd, nks, nspin, ntetra, tetra (4, ntetra)
  real(DP), intent(in) :: et (nbnd, nks), e
  integer, intent(in) :: is, isk
  ! local variables
  real(DP) :: etetra (4), e1, e2, e3, e4
  integer :: nt, nk, ns, ibnd, i, nspin_lsda

  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  sumkt = 0.0d0
  do ns = 1, nspin_lsda
     if (is /= 0) then
        if ( ns .ne. is) cycle
     end if
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
           call piksort (4, etetra)
           !
           ! ...sort in ascending order: e1 < e2 < e3 < e4
           !
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           !
           ! calculate sum over k of the integrated charge
           !
           if (e.ge.e4) then
              sumkt = sumkt + 1.d0 / ntetra
           elseif (e.lt.e4.and.e.ge.e3) then
              sumkt = sumkt + 1.d0 / ntetra * (1.0d0 - (e4 - e) **3 / (e4 - e1) &
                   / (e4 - e2) / (e4 - e3) )
           elseif (e.lt.e3.and.e.ge.e2) then
              sumkt = sumkt + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * &
                   ( (e2 - e1) **2 + 3.0d0 * (e2 - e1) * (e-e2) + 3.0d0 * (e-e2) **2 - &
                   (e3 - e1 + e4 - e2) / (e3 - e2) / (e4 - e2) * (e-e2) **3)
           elseif (e.lt.e2.and.e.ge.e1) then
              sumkt = sumkt + 1.d0 / ntetra * (e-e1) **3 / (e2 - e1) / &
                   (e3 - e1) / (e4 - e1)
           endif
        enddo
     enddo
  enddo

  ! add correct spin normalization (2 for LDA, 1 for other cases)

  IF ( nspin == 1 ) sumkt = sumkt * 2.d0

  return

end function sumkt

subroutine piksort (n, a)
  USE kinds
  implicit none
  integer :: n

  real(DP) :: a (n)
  integer :: i, j
  real(DP) :: temp
  !
  do j = 2, n
     temp = a (j)
     do i = j - 1, 1, - 1
        if (a (i) .le.temp) goto 10
        a (i + 1) = a (i)
     enddo
     i = 0
10   a (i + 1) = temp
  enddo
  !
  return
end subroutine piksort
