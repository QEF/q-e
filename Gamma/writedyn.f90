!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine writedyn ( )
  !
  USE ions_base, ONLY : nat, tau, ityp, ntyp => nsp, atm, amass
  use cgcom
  use pwcom
  implicit none
  integer :: iudyn, nt, na, nb, i, j
  !
  iudyn = 20
  open(unit=iudyn,file=fildyn,form='formatted',status='unknown')
  !
  !  write the dynamical matrix on on file
  !
  write(iudyn,'(a)') title
  write(iudyn,'(a)') title_ph
  write(iudyn,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
  do nt = 1,ntyp
     write(iudyn,*) nt," '",atm(nt),"' ",amass(nt)
  end do
  do na=1,nat
     write(iudyn,'(2i5,3f15.7)') na,ityp(na),(tau(j,na),j=1,3)
  end do
  write (iudyn, '(/,5x,"Dynamical  Matrix in cartesian axes", &
       &         //,5x,"q = ( ",3f14.9," ) ",/)') 0.0d0,0.0d0,0.0d0
  do na = 1, nat
     do nb = 1, nat
        write(iudyn, '(2i3)') na, nb
        write(iudyn,'(3e24.12)') &
             ( (dyn(3*(na-1)+i,3*(nb-1)+j),0.d0,j=1,3),i=1,3)
     end do
  end do
  !
  !   as above, for dielectric tensor and effective charges
  !
  if (epsil) then
     write (iudyn, '(/,5x,"Dielectric Tensor:",/)')
     write (iudyn, '(3e24.12)') ( (epsilon0(i,j) , j=1,3), i=1,3)
     write (iudyn, '(/5x, "Effective Charges E-U: Z_{alpha}{s,beta}",/)')
     do na = 1,nat
        write (iudyn, '(5x,"atom # ",i4)') na
        write (iudyn, '(3e24.12)') ( (zstar (i,j, na),j=1,3),i=1,3)
     end do
  end if
  close (unit=iudyn)
  return
end subroutine writedyn
