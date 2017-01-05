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
! updated by PG on Feb, 2013 upon suggestion by Thomas Gruber:
! workaround for VESTA - reverted to previous behavior in Oct 2013
! since workaround is no longer needed
! -------------------------------------------------------------------
! this routine writes a gaussian 98 like formatted cubefile.
! atoms outside the supercell are wrapped back according to PBC.
! plain dumping of the data. no re-gridding or transformation to an
! orthorhombic box (needed for most .cube aware programs :-/).
! -------------------------------------------------------------------
SUBROUTINE write_cubefile ( alat, at, bg, nat, tau, atm, ityp, rho, &
     nr1, nr2, nr3, nr1x, nr2x, nr3x, ounit )

  USE kinds,  ONLY : DP
  USE run_info, ONLY: title

  IMPLICIT NONE
  INTEGER, INTENT(IN):: nat, ityp(nat), ounit, nr1,nr2,nr3, nr1x,nr2x,nr3x
  CHARACTER(len=3), INTENT(IN) :: atm(*)
  REAL(DP), INTENT(IN) :: alat, tau(3,nat), at(3,3),bg(3,3), rho(nr1x,nr2x,nr3x)
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

  WRITE(ounit,*) 'Cubefile created from PWScf calculation'
  IF ( LEN_TRIM(title) > 1 ) THEN
     WRITE(ounit,*) TRIM(title) ! perhaps there is a better option...
  ELSE
     WRITE(ounit,'("Contains the selected quantity on a FFT grid")')
  END IF
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


! -------------------------------------------------------------------
! this routine instead writes a re-gridded cubefile (i.e. by B-spline
! interpolation)
! -------------------------------------------------------------------
SUBROUTINE write_cubefile_new (alat, nat, tau, atm, ityp, x0, &
               m1, m2, m3, e1, e2, e3, nx, ny, nz, carica, ounit)
  USE kinds,     ONLY : dp
  USE io_global, ONLY : stdout
  USE cell_base, ONLY : at
  implicit none
  integer, intent(in)  :: nat, ityp(nat), ounit, nx, ny, nz
  real(dp), intent(in) :: alat, tau(3,nat)
  character(len=3)     :: atm(*)
  real(dp), intent(in) :: m1, m2, m3, x0(3), e1(3), e2(3), e3(3), carica(nx,ny,nz)
  integer              :: ia, i, j, k, at_num
  integer, external    :: atomic_number
  real(dp)             :: at_chrg, tpos(3), inpos(3)
  real(dp)             :: bbmin(3), bbmax(3)
  integer, parameter   :: natomsmax = 10000
  real(dp)             :: taupos(3,natomsmax), pos(3)
  integer              :: natoms, taupostyp(natomsmax)

  ! generate bounding box
  bbmin(:) = 1d30
  bbmax(:) = -1d30
  call bbox(x0, bbmin, bbmax)
  call bbox(x0+e1, bbmin, bbmax)
  call bbox(x0+e2, bbmin, bbmax)
  call bbox(x0+e3, bbmin, bbmax)
  call bbox(x0+e1+e2, bbmin, bbmax)
  call bbox(x0+e2+e3, bbmin, bbmax)
  call bbox(x0+e3+e1, bbmin, bbmax)
  call bbox(x0+e1+e2+e3, bbmin, bbmax)
  write(stdout,'(5X,''Bounding box= ['',F12.4,'','',F12.4,'']'')') bbmin(1)*alat, bbmax(1)*alat
  write(stdout,'(5X,''              ['',F12.4,'','',F12.4,'']'')') bbmin(2)*alat, bbmax(2)*alat
  write(stdout,'(5X,''              ['',F12.4,'','',F12.4,'']'')') bbmin(3)*alat, bbmax(3)*alat

  ! generate atoms in bounding box
  natoms = 0
  do i = -5, 5
     do j = -5, 5
        do k = -5, 5
           do ia = 1, nat
              pos = tau(:,ia) + i*at(:,1) + j*at(:,2) + k*at(:,3)
              if (all(pos >= bbmin) .and. all(pos <= bbmax)) then
                 natoms = natoms + 1
                 if (natoms > natomsmax) &
                    call errore('write_cubefile_new', 'increase natomsmax', natoms)
                 taupos(:,natoms) = pos(:)
                 taupostyp(natoms) = ityp(ia)
              endif
           enddo
        enddo
     enddo
  enddo
  write(stdout,'(5X,I6,'' atoms inside bounding box'')') natoms

  write(ounit,*) 'cubfile created from pwscf calculation'
  write(ounit,*) 'total scf density'
  write(ounit,'(i5,3f12.6)') natoms, x0(:)*alat
  write(ounit,'(i5,3f12.6)') nx, alat*m1*e1(:)/dble(nx)
  write(ounit,'(i5,3f12.6)') ny, alat*m2*e2(:)/dble(ny)
  write(ounit,'(i5,3f12.6)') nz, alat*m3*e3(:)/dble(nz)

  do ia = 1, natoms
     at_num = atomic_number(trim(atm(taupostyp(ia))))
     at_chrg = dble(at_num)
     write(ounit,'(i5,5f12.6)') at_num, at_chrg, alat*taupos(:,ia)
  enddo

  do i=1,nx
     do j=1,ny
        write(ounit,'(6e13.5)') (carica(i,j,k),k=1,nz)
     enddo
  enddo
  return

END SUBROUTINE write_cubefile_new


SUBROUTINE bbox(r, bbmin, bbmax)
   USE kinds, only: dp
   implicit none
   real(dp), intent(in) :: r(3)
   real(dp), intent(inout) :: bbmin(3), bbmax(3)
   integer :: i
   do i = 1, 3
      bbmin(i) = min(bbmin(i), r(i))
      bbmax(i) = max(bbmax(i), r(i))
   enddo
END SUBROUTINE bbox

