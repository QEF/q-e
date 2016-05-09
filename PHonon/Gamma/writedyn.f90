!
! Copyright (C) 2003-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE writedyn ( )
  !
  USE ions_base, ONLY : nat, tau, ityp, ntyp => nsp, atm, amass
  USE run_info,  ONLY : title
  USE cell_base, ONLY : ibrav, celldm, at

  USE constants, ONLY : amu_ry
  USE cgcom
  IMPLICIT NONE
  INTEGER :: iudyn, nt, na, nb, i, j
  !
  iudyn = 20
  OPEN(unit=iudyn,file=fildyn,form='formatted',status='unknown')
  !
  !  write the dynamical matrix to file
  !
  WRITE(iudyn,'(a)') title
  WRITE(iudyn,'(a)') title_ph
  WRITE(iudyn,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
  IF (ibrav==0) THEN
     WRITE (iudyn,'("Basis vectors")')
     WRITE (iudyn,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
  END IF
  DO nt = 1,ntyp
     WRITE(iudyn,*) nt," '",atm(nt),"' ",amu_ry*amass(nt)
  ENDDO
  DO na=1,nat
     WRITE(iudyn,'(2i5,3f15.7)') na,ityp(na),(tau(j,na),j=1,3)
  ENDDO
  WRITE (iudyn, '(/,5x,"Dynamical  Matrix in cartesian axes", &
       &         //,5x,"q = ( ",3f14.9," ) ",/)') 0.0d0,0.0d0,0.0d0
  DO na = 1, nat
     DO nb = 1, nat
        WRITE(iudyn, '(2i5)') na, nb
        WRITE(iudyn,'(3e24.12)') &
             ( (dyn(3*(na-1)+i,3*(nb-1)+j),0.d0,j=1,3),i=1,3)
     ENDDO
  ENDDO
  !
  !   as above, for dielectric tensor and effective charges
  !
  IF (epsil) THEN
     WRITE (iudyn, '(/,5x,"Dielectric Tensor:",/)')
     WRITE (iudyn, '(3e24.12)') ( (epsilon0(i,j) , j=1,3), i=1,3)
     WRITE (iudyn, '(/5x, "Effective Charges E-U: Z_{alpha}{s,beta}",/)')
     DO na = 1,nat
        WRITE (iudyn, '(5x,"atom # ",i4)') na
        WRITE (iudyn, '(3e24.12)') ( (zstar (i,j, na),j=1,3),i=1,3)
     ENDDO
  ENDIF
  CLOSE (unit=iudyn)
  RETURN
END SUBROUTINE writedyn
