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
! nrx1,nrx2,nrx3 (the physical dimensions of array rho) differ from
!  nr1, nr2, nr3 (the true dimensions)
!

! -------------------------------------------------------------------
! this routine writes a gaussian 98 like formatted cubefile.
! atoms outside the supercell are wrapped back according to PBC.
! plain dumping of the data. no re-gridding or transformation to an
! orthorhombic box (needed for most .cube aware programs :-/).
! -------------------------------------------------------------------
subroutine write_cubefile ( alat, at, bg, nat, tau, atm, ityp, rho, &
     nr1, nr2, nr3, nrx1, nrx2, nrx3, ounit )

  USE kinds,  only : DP

  implicit none
  integer          :: nat, ityp(nat), ounit,nr1, nr2, nr3, nrx1, nrx2, nrx3
  character(len=3) :: atm(*)
  real(DP)    :: alat, tau(3,nat), at(3,3), bg(3,3), rho(nrx1,nrx2,nrx3)

  ! --
  integer          :: i, nt, i1, i2, i3, at_num
  integer, external:: atomic_number
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

  write(ounit,*) 'Cubfile created from PWScf calculation'
  write(ounit,*) ' Total SCF Density'
!                        origin is forced to (0.0,0.0,0.0)
  write(ounit,'(I5,3F12.6)') nat, 0.0d0, 0.0d0, 0.0d0
  write(ounit,'(I5,3F12.6)') nr1, (alat*at(i,1)/DBLE(nr1),i=1,3)
  write(ounit,'(I5,3F12.6)') nr2, (alat*at(i,2)/DBLE(nr2),i=1,3)
  write(ounit,'(I5,3F12.6)') nr3, (alat*at(i,3)/DBLE(nr3),i=1,3)

  do i=1,nat
     nt = ityp(i)
     ! find atomic number for this atom. 
     at_num = atomic_number(TRIM(atm(nt)))
     at_chrg= DBLE(at_num)
     ! at_chrg could be alternatively set to valence charge
     ! positions are in cartesian coordinates and a.u.
     !
     ! wrap coordinates back into cell.
     tpos = MATMUL( TRANSPOSE(bg), tau(:,i) )
     tpos = tpos - NINT(tpos - 0.5d0)
     inpos = alat * MATMUL( at, tpos )
     write(ounit,'(I5,5F12.6)') at_num, at_chrg, inpos
  enddo
  
  do i1=1,nr1
     do i2=1,nr2
        write(ounit,'(6E13.5)') (rho(i1,i2,i3),i3=1,nr3)
     enddo
  enddo
  return
end subroutine write_cubefile

