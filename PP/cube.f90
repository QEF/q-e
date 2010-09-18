!
! Copyright (C) 2004 Axel Kohlmeyer
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! This file holds gaussian cube generation subroutines.
! Adapted by Axel Kohlmeyer from xsf.f90.
! updated by Axel Kohlmeyer on Sep 27, 2004.
! updated by PG on Sep. 15, 2005 to account for the case in which
! nr1x,nr2x,nr3x (the physical dimensions of array rho) differ from
!  nr1, nr2, nr3 (the true dimensions)
!

! -------------------------------------------------------------------
! this routine writes a gaussian 98 like formatted cubefile.
! atoms outside the supercell are wrapped back according to PBC.
! plain dumping of the data. no re-gridding or transformation to an
! orthorhombic box (needed for most .cube aware programs :-/).
! -------------------------------------------------------------------
SUBROUTINE write_cubefile ( alat, at, bg, nat, tau, atm, ityp, rho, &
     nr1, nr2, nr3, nr1x, nr2x, nr3x, ounit )

  USE kinds,  ONLY : DP

  IMPLICIT NONE
  INTEGER          :: nat, ityp(nat), ounit,nr1, nr2, nr3, nr1x, nr2x, nr3x
  CHARACTER(len=3) :: atm(*)
  real(DP)    :: alat, tau(3,nat), at(3,3), bg(3,3), rho(nr1x,nr2x,nr3x)

  ! --
  INTEGER          :: i, nt, i1, i2, i3, at_num
  INTEGER, EXTERNAL:: atomic_number
  real(DP)    :: at_chrg, tpos(3), inpos(3)

!C     WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
!C     TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
!C     THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
!C
!C     LINE   FORMAT      CONTENTS
!C     ===============================================================
!C      1     A           TITLE
!C      2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
!C      3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
!C      4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
!C      #ATOMS LINES OF ATOM COORDINATES:
!C      ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
!C      REST: 6E13.5      CUBE DATA
!C
!C     ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.

  WRITE(ounit,*) 'Cubfile created from PWScf calculation'
  WRITE(ounit,*) ' Total SCF Density'
!                        origin is forced to (0.0,0.0,0.0)
  WRITE(ounit,'(I5,3F12.6)') nat, 0.0d0, 0.0d0, 0.0d0
  WRITE(ounit,'(I5,3F12.6)') nr1, (alat*at(i,1)/dble(nr1),i=1,3)
  WRITE(ounit,'(I5,3F12.6)') nr2, (alat*at(i,2)/dble(nr2),i=1,3)
  WRITE(ounit,'(I5,3F12.6)') nr3, (alat*at(i,3)/dble(nr3),i=1,3)

  DO i=1,nat
     nt = ityp(i)
     ! find atomic number for this atom.
     at_num = atomic_number(trim(atm(nt)))
     at_chrg= dble(at_num)
     ! at_chrg could be alternatively set to valence charge
     ! positions are in cartesian coordinates and a.u.
     !
     ! wrap coordinates back into cell.
     tpos = matmul( transpose(bg), tau(:,i) )
     tpos = tpos - nint(tpos - 0.5d0)
     inpos = alat * matmul( at, tpos )
     WRITE(ounit,'(I5,5F12.6)') at_num, at_chrg, inpos
  ENDDO

  DO i1=1,nr1
     DO i2=1,nr2
        WRITE(ounit,'(6E13.5)') (rho(i1,i2,i3),i3=1,nr3)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE write_cubefile

