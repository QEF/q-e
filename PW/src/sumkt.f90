!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
FUNCTION sumkt( et, nbnd, nks, nspin, ntetra, tetra, e, is, isk )
  !------------------------------------------------------------------
  !! Sum over all states with tetrahedron method.
  !
  !! At Fermi energy e=E_F, sumkt(e) == number of electrons.
  !
  !! Generalization to noncollinear case courtesy of Iurii Timrov.
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sumkt
  !! Output: see routine comments
  INTEGER, INTENT(IN) :: nbnd
  !! the number of bands
  INTEGER, INTENT(IN) :: nks
  !! the number of k points
  INTEGER, INTENT(IN) :: nspin
  !! the number of spin components
  INTEGER, INTENT(IN) :: ntetra
  !! the number of tetrahedra
  INTEGER, INTENT(IN) :: tetra(4,ntetra)
  !! the tetrahedron vertices
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! the energy eigenvalues
  REAL(DP), INTENT(IN) :: e
  !! see routine comments
  INTEGER, INTENT(IN) :: is
  !! spin component
  INTEGER, INTENT(IN) :: isk
  !! spin component for each k-point
  !
  ! ... local variables
  !
  REAL(DP) :: etetra(4), e1, e2, e3, e4
  INTEGER :: nt, nk, ns, ibnd, i, nspin_lsda
  !
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  !
  sumkt = 0.0_DP
  !
  DO ns = 1, nspin_lsda
     !
     IF ( is /= 0 ) THEN
        IF ( ns /= is) CYCLE
     ENDIF
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns==1) THEN
        nk = 0
     ELSE
        nk = nks / 2
     ENDIF
     !
     DO nt = 1, ntetra
        DO ibnd = 1, nbnd
           !
           ! etetra are the energies at the vertexes of the nt-th tetrahedron
           !
           DO i = 1, 4
              etetra(i) = et( ibnd, tetra(i,nt) + nk )
           ENDDO
           !
           CALL piksort( 4, etetra )
           !
           ! ...sort in ascending order: e1 < e2 < e3 < e4
           !
           e1 = etetra(1)
           e2 = etetra(2)
           e3 = etetra(3)
           e4 = etetra(4)
           !
           ! calculate sum over k of the integrated charge
           !
           IF ( e >= e4 ) THEN
              !
              sumkt = sumkt + 1.0_DP / ntetra
              !
           ELSEIF ( e<e4 .AND. e>=e3 ) THEN
              !
              sumkt = sumkt + 1.0_DP / ntetra * (1.0_DP - (e4 - e)**3 / (e4 - e1) &
                      / (e4 - e2) / (e4 - e3) )
              !
           ELSEIF ( e<e3 .AND. e>=e2) THEN
              !
              sumkt = sumkt + 1.0_DP / ntetra / (e3 - e1) / (e4 - e1) * &
                      ( (e2 - e1)**2 + 3.0_DP * (e2 - e1) * (e-e2) + 3.0_DP * (e-e2)**2 &
                      - (e3 - e1 + e4 - e2) / (e3 - e2) / (e4 - e2) * (e-e2)**3 )
              !
           ELSEIF (e<e2 .AND. e>=e1) THEN
              !
              sumkt = sumkt + 1.0_DP / ntetra * (e-e1)**3 / (e2 - e1) &
                      / (e3 - e1) / (e4 - e1)
              !
           ENDIF
           !
        ENDDO
     ENDDO
     !
  ENDDO
  !
  ! add correct spin normalization (2 for LDA, 1 for other cases)
  !
  IF ( nspin == 1 ) sumkt = sumkt * 2.0_DP
  !
  !
  RETURN
  !
END FUNCTION sumkt
!
!
!---------------------------------------------------
SUBROUTINE piksort( n, a )
  !-------------------------------------------------
  !! Sorts in ascending order.
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER :: n
  !! number of values
  REAL(DP) :: a(n)
  !! array of values
  !
  ! ... local variables
  !
  INTEGER :: i, j
  REAL(DP) :: temp
  !
  DO j = 2, n
     temp = a(j)
     DO i = j - 1, 1, - 1
        IF (a(i) <= temp) GOTO 10
        a(i+1) = a(i)
     ENDDO
     i = 0
10   a(i+1) = temp
  ENDDO
  !
  !
  RETURN
  !
END SUBROUTINE piksort
