!
! Copyright (C) 2006 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! GNU License
!
! Eyvaz Isaev

! Theoretical Physics Department,
! Moscow State Institute of Steel and Alloys
! (Technological University)
!
! Condensed Matter Theory Group,
! Uppsala University, Sweden
!
! Eyvaz.Isaev@fysik.uu.se, 
! eyvaz_isaev@yahoo.com
!
      SUBROUTINE K_BRILLOUIN
      INCLUDE 'parameters.h'
      OPEN(UNIT=7,FILE='ttrinp',ACCESS='SEQUENTIAL') 
      REWIND (7)
      WRITE(6,9)
      READ(7,*) NPNT0,NTET0,NDIV
      WRITE(6,'(3I4)') NPNT0,NTET0,NDIV
      DO 5 KI=1,NPNT0
      READ(7,*)(PNT0(I,KI),I=1,3)
      WRITE(6,'(3f10.6)') (PNT0(I,KI),I=1,3)
5     CONTINUE
	npt=npnt0
      DO 6 KI=1,NTET0
      READ(7,*) (TTR0(I,KI),I=1,4)
6     CONTINUE
      WRITE(6,9)
9     FORMAT(6('*'),' input tetrahedra for BZ-integration ',6('*'))
      RETURN
      END 
