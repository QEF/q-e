!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE write_tensor_field(name, ispin, field)
  !-----------------------------------------------------------------------
  !
  ! ... write the tensor field in xcrysden format
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout
  USE mp_global,                   ONLY : me_pool  
  USE cell_base,                   ONLY : at, bg, alat
  USE ions_base,                   ONLY : nat, tau, atm, ityp
  USE pwcom
  USE gipaw_module
  !--------------------------------------------------------------------
  character*(*) name
  integer :: ispin
  real(dp) :: field(nrx1,nrx2,nrx3,3)
  !--------------------------------------------------------------------
  integer, parameter :: ounit = 48
  character*80 :: fname
  integer :: ios, ipol
  
  if (me_pool /= 0) return

  do ipol = 1, 3
    ! form the name of the output file
    if (ispin == 0) fname = trim(name)//'_'
    if (ispin == 1) fname = trim(name)//'_UP_'
    if (ispin == 2) fname = trim(name)//'_DW_'

    if (ipol == 1) fname = trim(fname)//'X.xsf'
    if (ipol == 2) fname = trim(fname)//'Y.xsf'
    if (ipol == 3) fname = trim(fname)//'Z.xsf'
    write(stdout, '(5X,''write_tensor_field: '',A40)') fname

    open(unit=ounit, file=fname, iostat=ios, form='formatted', &
         status='unknown')
    if (ios /= 0) &
      call errore('write_tensor_field', 'error opening '//fname, ounit)

    call xsf_struct (alat, at, nat, tau, atm, ityp, nr1*nr2*nr3, ounit)
    call xsf_vector_3d(field(1,1,1,ipol), &
                       nr1, nr2, nr3, nrx1, nrx2, nrx3, at, bg, alat, ounit)
    close(unit=48)
  enddo
end subroutine write_tensor_field



!
! Copyright (C) 2003 Tone Kokalj
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file holds XSF (=Xcrysden Structure File) utilities.
! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
!
! -------------------------------------------------------------------
!   this routine writes the crystal structure in XSF format
! -------------------------------------------------------------------
subroutine xsf_struct (alat, at, nat, tau, atm, ityp, nr, ounit)
  USE kinds, only : DP
  implicit none
  integer          :: nat, ityp (nat), nr, ounit
  character(len=3) :: atm(*)
  real(DP)    :: alat, tau (3, nat), at (3, 3)
  ! --
  integer          :: i, j, n
  real(DP)    :: at1 (3, 3)
  ! convert lattice vectors to ANGSTROM units ...
  do i=1,3
     do j=1,3
        at1(j,i) = at(j,i)*alat*0.529177d0
     enddo
  enddo

  write(ounit,*) 'CRYSTAL'
  write(ounit,*) 'PRIMVEC'
  write(ounit,'(2(3F15.9/),3f15.9)') at1
  write(ounit,*) 'PRIMCOORD'
  write(ounit,*) nat + nr, 1

  do n=1,nat
     ! positions are in Angstroms
     write(ounit,'(a3,3x,3f15.9)') atm(ityp(n)), &
          tau(1,n)*alat*0.529177d0, &
          tau(2,n)*alat*0.529177d0, &
          tau(3,n)*alat*0.529177d0
  enddo
  return
end subroutine xsf_struct



! -------------------------------------------------------------------
!   this routine writes a 3D vector field
! -------------------------------------------------------------------
subroutine xsf_vector_3d(v, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                         at, bg, alat, ounit)
  USE kinds, only : DP
  implicit none
  integer  :: nrx1, nrx2, nrx3, nr1, nr2, nr3, ounit
  real(DP) :: at(3,3), bg(3,3), x(3), alat, v(nrx1,nrx2,nrx3,3)
  integer  :: i1, i2, i3

  do i1 = 1, nr1
    do i2 = 1, nr2
      do i3 = 1, nr3
        ! coordinate in angstrom
        x(1) = dble(i1)/dble(nr1)     
        x(2) = dble(i2)/dble(nr2)     
        x(3) = dble(i3)/dble(nr3)
        ! crystal to cartesian
        call trnvect (x, at, bg, 1)
        x = x * alat * 0.529177d0
        write(ounit,'(''X  '',3x,3f15.9,2x,3e12.4)') x, v(i1,i2,i3,1:3)
      enddo
    enddo
  enddo
end subroutine xsf_vector_3d
